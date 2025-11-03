use anyhow::{Context, Result};
use clap::Parser;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use log::{error, info, warn, LevelFilter};
use rayon::prelude::*;
use serde::Deserialize;
use std::fs::{self, File};
use std::io::{self, BufReader, Read};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
// --- 新增：为 enable_steady_tick 引入 Duration ---
use std::time::Duration;

// --- 数据结构 ---

/// 用于从 JSON 文件中反序列化数据的结构体
// --- 修改：适配新的 JSON 报告格式 (PE + SE) ---
#[derive(Deserialize, Debug)]
struct RenamingReportEntry {
    sample_name: String,
    // --- 修复：重命名字段以消除 dead_code 警告 ---
    #[serde(rename = "library_type")] // 告诉 serde 匹配 JSON 中的 "library_type"
    _library_type: String, // "PE" 或 "SE" (下划线前缀告诉 Rust 编译器该字段有意未使用)

    // PE Fields (必须是 Option，因为 SE 样本没有这些字段)
    new_r1_path_relative: Option<String>,
    md5_r1: Option<String>,
    new_r2_path_relative: Option<String>,
    md5_r2: Option<String>,

    // SE Fields (必须是 Option，因为 PE 样本没有这些字段)
    new_se_path_relative: Option<String>,
    md5_se: Option<String>,
}

/// 用于在程序内部处理的验证任务
#[derive(Debug, Clone)]
struct VerificationTask {
    file_to_check: PathBuf,
    expected_md5: String,
    sample_name: String,
}

/// 存储单个文件验证结果的结构体
#[derive(Debug)]
struct VerificationResult {
    timestamp: String,
    sample_name: String,
    file_path: String,
    expected_md5: String,
    actual_md5: String,
    status: &'static str,
    message: String,
}

// --- 命令行参数 ---

/// 基于 JSON 报告，使用多线程并发验证文件 MD5 校验和。
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// 由 seq_preprocessor 生成的 JSON 报告文件。
    #[arg(short, long)]
    input: PathBuf,

    /// 包含已整理数据的根目录 (JSON 报告中相对路径的基准路径)。
    #[arg(short, long, default_value = ".")]
    base_dir: PathBuf,

    /// 用于并发验证的线程数 (0 表示使用 Rayon 的默认值)。
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// 生成 TSV 格式验证报告的输出文件路径。
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// 指定日志文件的路径。
    #[arg(long, default_value = "verifier.log")]
    log_file: PathBuf,
}

// --- 核心功能函数 ---

/// 配置日志记录器
fn setup_logger(log_path: &Path) -> Result<()> {
    if let Some(parent_dir) = log_path.parent() {
        fs::create_dir_all(parent_dir)?;
    }
    fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "{} [{}] {}",
                chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                record.level(),
                message
            ))
        })
        .level(LevelFilter::Info)
        .chain(io::stdout())
        .chain(fern::log_file(log_path)?)
        .apply()?;
    Ok(())
}

/// 计算文件的 MD5 值
fn calculate_md5(filepath: &Path) -> io::Result<String> {
    let file = File::open(filepath)?;
    let mut reader = BufReader::new(file);
    let mut context = md5::Context::new();
    let mut buffer = [0; 65536]; // 64KB buffer

    loop {
        let bytes_read = reader.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        context.consume(&buffer[..bytes_read]);
    }
    Ok(format!("{:x}", context.compute()))
}

/// 执行单个文件的验证任务
fn verify_file_task(task: &VerificationTask) -> VerificationResult {
    let timestamp = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
    let file_path_str = task.file_to_check.to_string_lossy().to_string();

    if !task.file_to_check.exists() {
        return VerificationResult {
            timestamp,
            sample_name: task.sample_name.clone(),
            file_path: file_path_str,
            expected_md5: task.expected_md5.clone(),
            actual_md5: "N/A".to_string(),
            status: "FAIL",
            message: "File not found".to_string(),
        };
    }

    match calculate_md5(&task.file_to_check) {
        Ok(actual_md5) => {
            if actual_md5.eq_ignore_ascii_case(&task.expected_md5) {
                VerificationResult { timestamp, sample_name: task.sample_name.clone(), file_path: file_path_str, expected_md5: task.expected_md5.clone(), actual_md5, status: "PASS", message: "MD5 match".to_string() }
            } else {
                VerificationResult { timestamp, sample_name: task.sample_name.clone(), file_path: file_path_str, expected_md5: task.expected_md5.clone(), actual_md5, status: "FAIL", message: "MD5 mismatch".to_string() }
            }
        }
        Err(e) => VerificationResult { timestamp, sample_name: task.sample_name.clone(), file_path: file_path_str, expected_md5: task.expected_md5.clone(), actual_md5: "N/A".to_string(), status: "FAIL", message: format!("Read error: {}", e) },
    }
}

/// 将验证结果写入 TSV 报告
fn generate_report(results: &[VerificationResult], output_file: &Path) -> Result<()> {
    info!("--- 正在生成验证报告至: {} ---", output_file.display());
    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(output_file)?;
    writer.write_record(&["CheckTime", "SampleName", "FilePath", "ExpectedMD5", "ActualMD5", "Status", "Message"])?;
    for res in results {
        writer.write_record(&[&res.timestamp, &res.sample_name, &res.file_path, &res.expected_md5, &res.actual_md5, res.status, &res.message])?;
    }
    writer.flush()?;
    info!("报告生成完成。");
    Ok(())
}

// --- 主函数 ---

fn main() -> Result<()> {
    let cli = Cli::parse();
    setup_logger(&cli.log_file)?;

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(cli.threads).build_global()?;
    }

    info!("--- 开始 MD5 验证流程 ---");
    info!("正在读取并解析 JSON 报告: {}", cli.input.display());

    let json_content = fs::read_to_string(&cli.input)
        .context(format!("无法读取 JSON 文件: {}", cli.input.display()))?;
    let report_entries: Vec<RenamingReportEntry> = serde_json::from_str(&json_content)
        .context("解析 JSON 文件失败，请检查文件格式是否正确。")?;

    let mut tasks = Vec::new();
    
    // --- 修改：更新任务收集逻辑以支持 PE 和 SE ---
    info!("正在从报告中收集验证任务...");
    for entry in &report_entries {
        // 检查 R1 (PE)
        // 只有当 MD5 和 路径 *同时* 存在时，才添加任务
        if let (Some(md5_r1), Some(path_r1)) = (&entry.md5_r1, &entry.new_r1_path_relative) {
            tasks.push(VerificationTask {
                file_to_check: cli.base_dir.join(path_r1),
                expected_md5: md5_r1.clone(),
                sample_name: entry.sample_name.clone(),
            });
        }
        
        // 检查 R2 (PE)
        if let (Some(md5_r2), Some(path_r2)) = (&entry.md5_r2, &entry.new_r2_path_relative) {
            tasks.push(VerificationTask {
                file_to_check: cli.base_dir.join(path_r2),
                expected_md5: md5_r2.clone(),
                sample_name: entry.sample_name.clone(),
            });
        }
        
        // 检查 SE (Long-Read)
        if let (Some(md5_se), Some(path_se)) = (&entry.md5_se, &entry.new_se_path_relative) {
            tasks.push(VerificationTask {
                file_to_check: cli.base_dir.join(path_se),
                expected_md5: md5_se.clone(),
                sample_name: entry.sample_name.clone(),
            });
        }
    }

    let num_tasks = tasks.len();
    if num_tasks == 0 {
        warn!("JSON 报告中未包含任何有效的 MD5 记录，无需验证。");
        return Ok(());
    }

    info!(
        "找到 {} 个文件待验证，开始使用 {} 个线程进行处理。",
        num_tasks,
        rayon::current_num_threads()
    );

    let pb = ProgressBar::new(num_tasks as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) | {per_sec} | ETA: {eta}")?
        .progress_chars("##-"));

    // --- 修改：启用稳定计时器以确保旋转动画持续 ---
    // 这将强制进度条大约每 100ms 重绘一次，从而使 {spinner} 动画保持活动
    pb.enable_steady_tick(Duration::from_millis(100));

    let has_failures = AtomicBool::new(false);
    
    let results: Vec<VerificationResult> = tasks
        .par_iter()
        .progress_with(pb)
        .map(|task| {
            let result = verify_file_task(task);
            if result.status == "FAIL" {
                has_failures.store(true, Ordering::Relaxed);
                error!("[FAIL] 样本: {}, 文件: {}, 原因: {}", task.sample_name, task.file_to_check.display(), result.message);
                if result.message == "MD5 mismatch" {
                    error!("    - 预期: {}", result.expected_md5);
                    error!("    - 实际:   {}", result.actual_md5);
                }
            } 
            result
        })
        .collect();

    if let Some(output_path) = &cli.output {
        generate_report(&results, output_path)?;
    }

    info!("======================================================");
    if has_failures.load(Ordering::Relaxed) {
        error!("验证过程中发现错误。请检查日志和报告文件以获取详细信息。");
        // --- 【修复】: 将 . 改为 :: ---
        std::process::exit(1);
    } else {
        info!("所有 {} 个文件均成功通过验证！", results.len());
    }
    info!("======================================================");

    Ok(())
}

