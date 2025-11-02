


files <- c(file.path(root_dir, “short_summary.specific.embryophyta_odb10.length_filtered_busco.json”),
			file.path(root_dir, “short_summary.specific.embryophyta_odb10.length_filtered_cd-hit-est_busco.json”),
			file.path(root_dir, “short_summary.specific.embryophyta_odb10.length_filtered_cd-hit-est_E90_transcript_busco.json”),
			file.path(root_dir, “short_summary.specific.embryophyta_odb10.length_filtered_cd-hit-est_remove_transcrip_busco.json”),
			file.path(root_dir, “short_summary.specific.embryophyta_odb10.TD2_CDS_busco.json”)) 
names <- c(“Origin Transcript”,“CD-HIT”,“E90 Transcript”,“E90 Remove Transcript”,“TD2 CDS Transcript”)
df <- read_busco_jsons(files, names)
p <- plot_busco_alluvial(df, sample_order = names, palette = c(”#ffb703”,”#fb8500”,”#8ecae6”,”#219ebc”)
save_busco_plots(p, file.path(root_dir, “BUSCO_evalution”), width = 7, height = 6.5, dpi = 1000)