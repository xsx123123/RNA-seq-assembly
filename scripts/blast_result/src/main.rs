use clap::Parser;
use serde::Serialize;
use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;

/// 定义一个结构体来存储单次比对命中的信息
#[derive(Serialize, Debug)]
struct Hit {
    subject_id: String,
    score: f64,
    e_value: String,
}

/// 使用 clap 定义命令行参数
#[derive(Parser, Debug)]
#[command(author, version, about="一个用 Rust 编写的高性能 BLASTp 结果解析器", long_about = None)]
struct Cli {
    /// 输入的 blastp 结果文件路径
    #[arg(short, long)]
    input: PathBuf,

    /// 输出的 JSON 文件路径 (可选, 默认输出到屏幕)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// 为每个查询提取的最佳比对数量
    #[arg(short, long, default_value_t = 3)]
    top_n: usize,
}

/// 解析 BLAST 文件的核心函数
fn parse_blast_file(
    path: &PathBuf,
    top_n: usize,
) -> Result<HashMap<String, Vec<Hit>>, Box<dyn Error>> {
    let file = File::open(path)?;
    // 使用 BufReader 提高文件读取效率
    let reader = BufReader::new(file);

    let mut results: HashMap<String, Vec<Hit>> = HashMap::new();
    let mut current_query_id: Option<String> = None;
    let mut in_summary_section = false;
    let mut lines = reader.lines();

    while let Some(Ok(line)) = lines.next() {
        if line.starts_with("Query=") {
            // 从 "Query= query_id ..." 中提取 query_id
            if let Some(id) = line.split_whitespace().nth(1) {
                let query_id = id.to_string();
                // 确保为新的 query 创建一个空的 vector
                results.entry(query_id.clone()).or_insert_with(Vec::new);
                current_query_id = Some(query_id);
                in_summary_section = false;
            }
            continue;
        }

        if line.contains("Sequences producing significant alignments:") {
            in_summary_section = true;
            // 跳过标题后的两行（标题自身和空行）
            lines.next();
            lines.next();
            continue;
        }

        // 遇到详细比对区域或空行，则退出摘要区域
        if line.starts_with('>') || line.trim().is_empty() {
            in_summary_section = false;
            continue;
        }

        if in_summary_section {
            if let Some(ref query_id) = current_query_id {
                // 获取当前 query 对应的 hits 列表
                if let Some(hits) = results.get_mut(query_id) {
                    // 如果还没满 top_n，则解析该行
                    if hits.len() < top_n {
                        // 使用 rsplitn 从右边分割，这非常高效
                        let parts: Vec<&str> = line.trim().rsplitn(3, char::is_whitespace).collect();
                        if parts.len() == 3 {
                            // parts 的顺序是 [描述和ID, score, e_value]，但是是反的
                            let e_value = parts[0].to_string();
                            let score_str = parts[1];
                            let subject_part = parts[2];

                            // 从 subject_part 中获取第一个词作为 ID
                            if let Some(subject_id) = subject_part.split_whitespace().next() {
                                if let Ok(score) = score_str.parse::<f64>() {
                                    hits.push(Hit {
                                        subject_id: subject_id.to_string(),
                                        score,
                                        e_value,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(results)
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    println!(
        "正在从 '{}' 文件中解析 Top {} 的比对结果...",
        cli.input.display(),
        cli.top_n
    );

    let blast_results = parse_blast_file(&cli.input, cli.top_n)?;

    if blast_results.is_empty() {
        println!("未找到任何结果或解析失败。");
        return Ok(());
    }

    // 将结果序列化为格式化的 JSON 字符串
    let json_output = serde_json::to_string_pretty(&blast_results)?;

    match cli.output {
        Some(output_path) => {
            // 写入文件
            fs::write(&output_path, json_output)?;
            println!("结果已成功保存到 '{}'。", output_path.display());
        }
        None => {
            // 打印到标准输出
            println!("\n--- 解析结果 ---");
            println!("{}", json_output);
            println!("---   结束   ---");
        }
    }

    Ok(())
}