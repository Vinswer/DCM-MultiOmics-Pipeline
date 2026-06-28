rm(list = ls()); gc()

ORIGINAL_DIR <- ""
work_dir <- file.path(ORIGINAL_DIR, "00_rawdata")

if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}

setwd(ORIGINAL_DIR)
dirs <- c(
  file.path(ORIGINAL_DIR, "training_GSE57338"),
  file.path(ORIGINAL_DIR, "validation_GSE5406"),
  file.path(ORIGINAL_DIR, "singlecell_GSE183852"),
  file.path(ORIGINAL_DIR, "gene_sets")
)
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  cat("目录已就绪:", d, "\n")
}

library(GEOquery)
library(Biobase)
library(msigdbr)
library(data.table)
library(dplyr)
library(stringr)



download_geo_matrix <- function(geo_id, dest_dir) {
  gz_file <- paste0(geo_id, "_series_matrix.txt.gz")
  dest_file <- file.path(dest_dir, gz_file)

  if (file.exists(dest_file)) {
    cat("  文件已存在，跳过下载:", dest_file, "\n")
    return(dest_file)
  }

  prefix <- substr(geo_id, 1, nchar(geo_id) - 3)
  prefix <- paste0(prefix, "nnn")

  url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", prefix, "/",
                geo_id, "/matrix/", gz_file)

  cat("  下载URL:", url, "\n")
  cat("  保存至:", dest_file, "\n")

  tryCatch({
    download.file(url, destfile = dest_file, mode = "wb", method = "auto", quiet = FALSE)
    cat("  下载成功!\n")
    return(dest_file)
  }, error = function(e1) {
    cat("  方法1失败，尝试 curl...\n")
    tryCatch({
      download.file(url, destfile = dest_file, mode = "wb", method = "curl", quiet = FALSE)
      cat("  下载成功!\n")
      return(dest_file)
    }, error = function(e2) {
      cat("  方法2失败，尝试 wget...\n")
      tryCatch({
        download.file(url, destfile = dest_file, mode = "wb", method = "wget", quiet = FALSE)
        cat("  下载成功!\n")
        return(dest_file)
      }, error = function(e3) {
        cat("  所有下载方法均失败！请手动下载:", url, "\n")
        return(NULL)
      })
    })
  })
}







































































probe_to_gene <- function(expr_data, feature_data, gene_col = "Gene Symbol") {
  if (!gene_col %in% colnames(feature_data)) {
    possible_cols <- c("Gene Symbol", "gene_assignment", "GENE_SYMBOL",
                       "Symbol", "gene_symbol", "ILMN_Gene", "Gene symbol")
    gene_col <- intersect(possible_cols, colnames(feature_data))
    if (length(gene_col) == 0) {
      cat("  警告：无法找到基因符号列！可用列：\n")
      cat("  ", paste(colnames(feature_data), collapse = ", "), "\n")
      return(NULL)
    }
    gene_col <- gene_col[1]
  }
  cat("  使用基因符号列:", gene_col, "\n")

  probe_gene_map <- data.frame(
    probe_id    = rownames(feature_data),
    gene_symbol = as.character(feature_data[[gene_col]]),
    stringsAsFactors = FALSE
  )
  probe_gene_map <- probe_gene_map %>%
    filter(!is.na(gene_symbol) &
             gene_symbol != "" &
             gene_symbol != "---" &
             gene_symbol != "NA")
  probe_gene_map$gene_symbol <- sapply(
    probe_gene_map$gene_symbol,
    function(x) {
      first_entry <- trimws(strsplit(x, "///")[[1]][1])
      fields <- strsplit(first_entry, " // ")[[1]]
      if (length(fields) >= 2) trimws(fields[2]) else NA_character_
    }
  )

  valid_probes <- intersect(rownames(expr_data), probe_gene_map$probe_id)
  expr_filtered <- expr_data[valid_probes, , drop = FALSE]

  expr_df <- as.data.frame(expr_filtered)
  expr_df$probe_id <- rownames(expr_df)
  expr_df <- merge(expr_df, probe_gene_map, by = "probe_id", all.x = TRUE)

  sample_cols <- setdiff(colnames(expr_df), c("probe_id", "gene_symbol"))


  expr_by_gene <- expr_df %>%
    filter(!is.na(gene_symbol)) %>%
    dplyr::select(-probe_id) %>%
    mutate(row_var = apply(dplyr::select(., all_of(sample_cols)), 1, var, na.rm = TRUE)) %>%
    group_by(gene_symbol) %>%
    slice_max(order_by = row_var, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(-row_var)

  gene_symbols <- expr_by_gene$gene_symbol
  expr_matrix  <- as.matrix(expr_by_gene[, sample_cols])
  rownames(expr_matrix) <- gene_symbols

  cat("  探针数:", nrow(expr_data), " -> 基因数:", nrow(expr_matrix), "\n")
  cat("  重复基因处理策略: 保留方差最大的探针 (slice_max)\n")

  return(list(
    expr_matrix   = expr_matrix,
    probe_gene_map = probe_gene_map,
    gene_col_used  = gene_col
  ))
}


process_gpl96 <- function(expr_data, feature_data) {
  gene_vec <- sapply(feature_data[["Gene Symbol"]], function(x) {
    x <- as.character(x)
    if (is.na(x) || x == "" || x == "---") return(NA_character_)
    trimws(strsplit(x, "///")[[1]][1])
  })


  keep <- !is.na(gene_vec) & gene_vec != "" & gene_vec != "---"


  common_probes <- intersect(rownames(expr_data), rownames(feature_data)[keep])
  gene_labels   <- gene_vec[keep][match(common_probes, rownames(feature_data)[keep])]
  expr_sub      <- expr_data[common_probes, , drop = FALSE]


  expr_df <- as.data.frame(expr_sub)
  expr_df$gene_symbol <- gene_labels
  sample_cols <- colnames(expr_sub)

  expr_by_gene <- expr_df %>%
    filter(!is.na(gene_symbol)) %>%
    mutate(row_var = apply(dplyr::select(., all_of(sample_cols)), 1, var, na.rm = TRUE)) %>%
    group_by(gene_symbol) %>%
    slice_max(order_by = row_var, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(-row_var)

  gene_symbols <- expr_by_gene$gene_symbol
  expr_matrix  <- as.matrix(expr_by_gene[, sample_cols])
  rownames(expr_matrix) <- gene_symbols

  cat("  探针数:", nrow(expr_data), "-> 基因数:", nrow(expr_matrix), "\n")
  return(expr_matrix)
}




print_group_info <- function(pheno_data, dataset_name) {
  cat("\n--- ", dataset_name, " 临床信息详细 ---\n")
  cat("样本总数:", nrow(pheno_data), "\n")
  cat("所有metadata字段:\n")
  cat(paste("  -", colnames(pheno_data)), sep = "\n")


  group_cols <- c("characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2",
                  "characteristics_ch1.3", "characteristics_ch1.4",
                  "source_name_ch1", "title", "description",
                  grep("disease|condition|group|status|type|tissue|diagnosis|sample",
                       colnames(pheno_data), ignore.case = TRUE, value = TRUE))
  group_cols <- unique(group_cols[group_cols %in% colnames(pheno_data)])

  cat("\n可能包含分组信息的列:\n")
  for (col in group_cols) {
    cat("\n  列名: [", col, "]\n")
    tbl <- table(pheno_data[[col]], useNA = "ifany")
    if (length(tbl) <= 30) {
      for (i in seq_along(tbl)) {
        cat("    ", names(tbl)[i], ": ", tbl[i], " 个样本\n", sep = "")
      }
    } else {
      cat("    （唯一值过多:", length(tbl), "个，仅显示前10）\n")
      tbl_sorted <- sort(tbl, decreasing = TRUE)
      for (i in 1:min(10, length(tbl_sorted))) {
        cat("    ", names(tbl_sorted)[i], ": ", tbl_sorted[i], " 个样本\n", sep = "")
      }
    }
  }


  chars_cols <- grep("characteristics", colnames(pheno_data), ignore.case = TRUE, value = TRUE)
  if (length(chars_cols) > 0) {
    cat("\n\n--- characteristics 列完整汇总 ---\n")
    for (col in chars_cols) {
      cat("\n  [", col, "]:\n")
      tbl <- table(pheno_data[[col]], useNA = "ifany")
      for (i in seq_along(tbl)) {
        cat("    ", names(tbl)[i], ": ", tbl[i], "\n", sep = "")
      }
    }
  }
}






cat("########################################################################\n")
cat("#  1. 训练集：GSE57338 (GPL11532, 芯片数据)                            #\n")
cat("########################################################################\n\n")

dest_dir_1 <- file.path(ORIGINAL_DIR, "training_GSE57338")


cat(">>> Step 1: 下载 GSE57338 series matrix...\n")
gz_file_1 <- download_geo_matrix("GSE57338", dest_dir_1)


if (!is.null(gz_file_1) && file.exists(gz_file_1)) {
  cat("\n>>> Step 2: 加载数据...\n")
  GSE57338 <- getGEO(filename = gz_file_1, getGPL = FALSE)

  expr_53778 <- exprs(GSE57338)
  pheno_53778 <- pData(GSE57338)
  feature_53778 <- fData(GSE57338)

  cat(">>> 下载 GPL11532 注释文件...\n")
  gpl11532 <- getGEO("GPL11532", destdir = dest_dir_1)
  feature_53778 <- Table(gpl11532)
  rownames(feature_53778) <- feature_53778$ID

  cat("GPL11532 列名:\n")
  print(colnames(feature_53778))
  cat("gene_assignment前3行:\n")
  print(head(feature_53778$gene_assignment, 3))

  cat("平台:", annotation(GSE57338), "\n")
  cat("样本数:", ncol(expr_53778), "\n")
  cat("探针数:", nrow(expr_53778), "\n")


  print_group_info(pheno_53778, "GSE57338")


  cat("\n\n>>> Step 3: 探针ID转基因Symbol...\n")
  cat("  Feature Data 可用列:\n")
  cat("  ", paste(colnames(feature_53778), collapse = ", "), "\n")

  gene_result_53778 <- probe_to_gene(expr_53778, feature_53778)

  cat("feature_53778列名:\n")
  print(colnames(feature_53778))
  cat("\ngene_assignment前3行:\n")
  print(head(feature_53778$gene_assignment, 3))
  cat("\n行名前3个:\n")
  print(head(rownames(feature_53778), 3))
  cat("\n表达矩阵行名前3个:\n")
  print(head(rownames(expr_53778), 3))

  if (!is.null(gene_result_53778)) {
    expr_gene_53778 <- gene_result_53778$expr_matrix


    write.csv(expr_gene_53778, file = file.path(dest_dir_1, "GSE57338_gene_expression.csv"))
    write.csv(pheno_53778, file = file.path(dest_dir_1, "GSE57338_metadata.csv"))
    save(GSE57338, expr_gene_53778, pheno_53778, feature_53778,
         file = file.path(dest_dir_1, "GSE57338_processed.RData"))
    cat("  基因表达矩阵已保存\n")
  }


  write.csv(expr_53778, file = file.path(dest_dir_1, "GSE57338_probe_expression.csv"))

  cat("\n>>> GSE57338 处理完成!\n\n")
} else {
  cat("!!! GSE57338 文件下载失败，请手动下载后重新运行 !!!\n\n")
}
cat("feature列名:\n")
print(colnames(fData(GSE57338)))
cat("\ngene_assignment前3行:\n")
print(head(fData(GSE57338)$gene_assignment, 3))





cat("########################################################################\n")
cat("#  2. 验证集：GSE5406 (GPL96, 表达谱芯片)                              #\n")
cat("########################################################################\n\n")

dest_dir_2 <- file.path(ORIGINAL_DIR, "validation_GSE5406")


cat(">>> Step 1: 下载 GSE5406 series matrix...\n")
gz_file_2 <- download_geo_matrix("GSE5406", dest_dir_2)


if (!is.null(gz_file_2) && file.exists(gz_file_2)) {
  cat("\n>>> Step 2: 加载数据...\n")
  gse5406 <- getGEO(filename = gz_file_2, getGPL = FALSE)

  expr_5406  <- exprs(gse5406)
  pheno_5406 <- pData(gse5406)


  gse5406_family <- getGEO(
    filename = file.path(dest_dir_2, "GSE5406_family.soft.gz")
  )
  gpl_list     <- GPLList(gse5406_family)
  gpl96_obj    <- gpl_list[[1]]
  feature_5406 <- Table(gpl96_obj)
  rownames(feature_5406) <- feature_5406$ID

  cat("feature_5406列名:\n")
  print(colnames(feature_5406))
  cat("Gene Symbol前3行:\n")
  print(head(feature_5406[["Gene Symbol"]], 3))


  print_group_info(pheno_5406, "GSE5406")


  cat("\n\n>>> Step 3: 探针ID转基因Symbol...\n")
  cat("  Feature Data 可用列:\n")
  cat("  ", paste(colnames(feature_5406), collapse = ", "), "\n")

  expr_gene_5406 <- process_gpl96(expr_5406, feature_5406)
  cat("基因数:", nrow(expr_gene_5406), "\n")

  write.csv(expr_gene_5406,
            file = file.path(dest_dir_2, "GSE5406_gene_expression.csv"))
  write.csv(pheno_5406,
            file = file.path(dest_dir_2, "GSE5406_metadata.csv"))
  save(gse5406, expr_gene_5406, pheno_5406, feature_5406,
       file = file.path(dest_dir_2, "GSE5406_processed.RData"))
  cat("验证集基因表达矩阵已保存\n")

  cat("\n>>> GSE5406 处理完成!\n\n")
} else {
  cat("!!! GSE5406 文件下载失败，请手动下载后重新运行 !!!\n\n")
}






cat("########################################################################\n")
cat("#  3. 单细胞数据集：GSE183852 (GPL24676, scRNA-seq)                     #\n")
cat("########################################################################\n\n")

dest_dir_3 <- file.path(ORIGINAL_DIR, "singlecell_GSE183852")


cat(">>> Step 1: 下载 GSE183852 series matrix (获取样本metadata)...\n")
gz_file_3 <- download_geo_matrix("GSE183852", dest_dir_3)

if (!is.null(gz_file_3) && file.exists(gz_file_3)) {
  cat("\n>>> Step 2: 加载metadata...\n")
  gse183852 <- getGEO(filename = gz_file_3, getGPL = FALSE)

  pheno_183852 <- pData(gse183852)


  print_group_info(pheno_183852, "GSE183852")


  cat("\n\n--- GSE183852 每个样本的完整信息 ---\n")
  chars_cols <- grep("characteristics|source_name|title", colnames(pheno_183852),
                     ignore.case = TRUE, value = TRUE)
  if (length(chars_cols) > 0) {
    for (i in 1:nrow(pheno_183852)) {
      cat("\n  样本", i, ":", rownames(pheno_183852)[i], "\n")
      for (col in chars_cols) {
        cat("    ", col, ": ", as.character(pheno_183852[i, col]), "\n", sep = "")
      }
    }
  }

  write.csv(pheno_183852, file = file.path(dest_dir_3, "GSE183852_metadata.csv"))
  save(gse183852, pheno_183852, file = file.path(dest_dir_3, "GSE183852_metadata.RData"))
}


cat("\n\n>>> Step 3: 下载 GSE183852 supplementary 文件（单细胞原始数据）...\n")
cat("  注意：单细胞数据文件可能很大，请耐心等待...\n")
tryCatch({
  supp_files <- getGEOSuppFiles("GSE183852",
                                baseDir = dest_dir_3,
                                makeDirectory = FALSE)
  cat("\n  补充文件下载完成！文件列表：\n")
  for (f in rownames(supp_files)) {
    fsize <- file.info(f)$size
    cat("  -", basename(f), " (", round(fsize / 1024^2, 1), "MB)\n")
  }
}, error = function(e) {
  cat("  补充文件下载失败：", conditionMessage(e), "\n")
  cat("  请手动从 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183852 下载\n")
  cat("  supplementary files 并放置到:", dest_dir_3, "\n")
})


cat("\n  目录内所有文件：\n")
all_files <- list.files(dest_dir_3, recursive = TRUE, full.names = TRUE)
for (f in all_files) {
  fsize <- file.info(f)$size
  cat("  -", basename(f), " (", round(fsize / 1024^2, 1), "MB)\n")
}

cat("\n>>> GSE183852 处理完成!\n\n")






cat("########################################################################\n")
cat("#  4. 吞噬调节因子相关基因集 (MSigDB)                                  #\n")
cat("########################################################################\n\n")

dest_dir_4 <- file.path(ORIGINAL_DIR, "gene_sets")

cat(">>> 正在从 MSigDB 获取吞噬相关基因集...\n\n")


all_gene_sets <- msigdbr(species = "Homo sapiens")


target_sets <- c(
  "GOBP_OPSONIZATION",
  "GOBP_PHAGOCYTOSIS_RECOGNITION",
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",
  "REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS",
  "REACTOME_RESPONSE_OF_MTB_TO_PHAGOCYTOSIS",
  "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
  "WP_MICROGLIA_PATHOGEN_PHAGOCYTOSIS_PATHWAY"
)

all_phagocytosis_genes <- c()
gene_set_list <- list()

for (gs_name in target_sets) {
  gs_data <- all_gene_sets %>% filter(gs_name == !!gs_name)

  if (nrow(gs_data) > 0) {
    genes <- unique(gs_data$gene_symbol)
    gene_set_list[[gs_name]] <- genes
    all_phagocytosis_genes <- c(all_phagocytosis_genes, genes)

    cat("基因集:", gs_name, "\n")
    cat("  来源:", unique(gs_data$gs_cat), "/", unique(gs_data$gs_subcat), "\n")
    cat("  基因数量:", length(genes), "\n")
    cat("  基因列表:", paste(sort(genes), collapse = ", "), "\n\n")


    write.csv(data.frame(gene_symbol = sort(genes)),
              file = file.path(dest_dir_4, paste0(gs_name, ".csv")),
              row.names = FALSE)
  } else {
    cat("!!! 基因集未找到:", gs_name, "\n")
    cat("    尝试模糊搜索...\n")
    fuzzy <- all_gene_sets %>%
      filter(grepl(gsub("_", ".*", gs_name), gs_name, ignore.case = TRUE)) %>%
      pull(gs_name) %>% unique()
    if (length(fuzzy) > 0) {
      cat("    可能的匹配:", paste(head(fuzzy, 5), collapse = "\n                 "), "\n")
    }
    cat("\n")
  }
}


all_phagocytosis_genes_unique <- unique(all_phagocytosis_genes)

cat("\n========================================\n")
cat("吞噬因子基因集汇总\n")
cat("========================================\n\n")
cat("各基因集基因数量:\n")
for (gs_name in names(gene_set_list)) {
  cat("  ", gs_name, ": ", length(gene_set_list[[gs_name]]), " 个基因\n", sep = "")
}
cat("\n合并前总基因数（含重复）:", length(all_phagocytosis_genes), "\n")
cat("合并去重后总基因数:", length(all_phagocytosis_genes_unique), "\n")


write.csv(data.frame(gene_symbol = sort(all_phagocytosis_genes_unique)),
          file = file.path(dest_dir_4, "Phagocytosis_all_genes_unique.csv"),
          row.names = FALSE)
save(gene_set_list, all_phagocytosis_genes_unique,
     file = file.path(dest_dir_4, "phagocytosis_gene_sets.RData"))

cat("\n去重后完整基因列表:\n")
cat(paste(sort(all_phagocytosis_genes_unique), collapse = ", "), "\n")


cat("\n\n--- 基因集间重叠矩阵 ---\n")
set_names <- names(gene_set_list)
overlap_matrix <- matrix(0, nrow = length(set_names), ncol = length(set_names),
                         dimnames = list(set_names, set_names))
for (i in seq_along(set_names)) {
  for (j in seq_along(set_names)) {
    overlap_matrix[i, j] <- length(intersect(gene_set_list[[set_names[i]]],
                                             gene_set_list[[set_names[j]]]))
  }
}
print(overlap_matrix)

write.csv(overlap_matrix, file = file.path(dest_dir_4, "gene_set_overlap_matrix.csv"))






cat("\n\n")
cat("================================================================\n")
cat("                     数据下载完成 - 总结报告                     \n")
cat("================================================================\n\n")

cat("1. 训练集 GSE57338 (GPL11532):\n")
cat("   路径:", file.path(dest_dir_1), "\n")
if (exists("expr_53778")) {
  cat("   样本数:", ncol(expr_53778), "\n")
  cat("   探针数:", nrow(expr_53778), "\n")
  if (exists("expr_gene_53778")) {
    cat("   基因数:", nrow(expr_gene_53778), "\n")
  }
}
cat("   文件: GSE57338_gene_expression.csv, GSE57338_metadata.csv, GSE57338_processed.RData\n\n")

cat("2. 验证集 GSE5406 (GPL96):\n")
cat("   路径:", file.path(dest_dir_2), "\n")
if (exists("expr_5406")) {
  cat("   样本数:", ncol(expr_5406), "\n")
  cat("   探针数:", nrow(expr_5406), "\n")
  if (exists("expr_gene_5406")) {
    cat("   基因数:", nrow(expr_gene_5406), "\n")
  }
}
cat("   文件: GSE5406_gene_expression.csv, GSE5406_metadata.csv, GSE5406_processed.RData\n\n")

cat("3. 单细胞数据集 GSE183852 (GPL24676):\n")
cat("   路径:", file.path(dest_dir_3), "\n")
if (exists("pheno_183852")) {
  cat("   样本数:", nrow(pheno_183852), "\n")
}
cat("   文件: GSE183852_metadata.csv, GSE183852_metadata.RData, 及supplementary文件\n\n")

cat("4. 吞噬调节因子基因集:\n")
cat("   路径:", file.path(dest_dir_4), "\n")
cat("   基因集数量: ", length(gene_set_list), "\n")
cat("   合并去重后基因数:", length(all_phagocytosis_genes_unique), "\n")
cat("   文件: 各基因集CSV + Phagocytosis_all_genes_unique.csv + phagocytosis_gene_sets.RData\n\n")
