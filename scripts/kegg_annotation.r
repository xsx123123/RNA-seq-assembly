# Author : JZHANG
# Date   : 2025-10-30
# Version: 1.0v
# ----- loading packages ----- #
library(qs)
library(KEGGREST)
library(dplyr)
# ----- run kegg annotation ----- #
kegg_dir <- '/home/zj/reference/KEGG_Datebase'
# ----- KEGG KO ANN ----- #
ko2path_raw <- keggLink("pathway", "ko")
ko2path <- data.frame(
  KO_ID = gsub("ko:", "", names(ko2path_raw)),
  Pathway_ID = gsub("path:", "", ko2path_raw)
)
rownames(ko2path) <- NULL
# save ko dataset
qsave(ko2path,file.path(kegg_dir,'ko_kegg_path.qs'))
# ----- KEGG KOPATH ANN ----- #
path2name_raw <- keggList("pathway")
path2name <- data.frame(
  Pathway_ID = gsub("path:", "", names(path2name_raw)),
  Pathway_Name = path2name_raw
)
rownames(path2name) <- NULL
qsave(path2name,file.path(kegg_dir,'kegg_path_kegg_name.qs'))