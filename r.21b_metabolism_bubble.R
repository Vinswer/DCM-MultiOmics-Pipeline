cat("[INFO] R:", R.version$version.string, "\n")


for (pkg in c("Seurat", "SeuratObject", "AUCell", "ggplot2", "dplyr", "tidyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[SKIP] 缺少必要包: %s，跳过\n", pkg))
    quit(status = 0)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(AUCell)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

WORK    <- ""
out_dir <- file.path(WORK, "15_metabolism")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


save_plot <- function(plot_obj, base_path, width, height, dpi = 200) {
  ggsave(paste0(base_path, ".png"),
         plot = plot_obj, width = width, height = height,
         dpi = dpi, bg = "white")
  cat("[INFO] 已保存", paste0(basename(base_path), ".png"), "\n")

  ggsave(paste0(base_path, ".pdf"),
         plot = plot_obj, width = width, height = height,
         bg = "white")
  cat("[INFO] 已保存", paste0(basename(base_path), ".pdf"), "\n")
}


black_text_theme <- theme(
  text             = element_text(color = "black"),
  axis.text        = element_text(color = "black"),
  axis.text.x      = element_text(color = "black", angle = 45, hjust = 1),
  axis.text.y      = element_text(color = "black"),
  axis.title       = element_text(color = "black"),
  plot.title       = element_text(color = "black", face = "bold"),
  plot.subtitle    = element_text(color = "black"),
  legend.text      = element_text(color = "black"),
  legend.title     = element_text(color = "black"),
  strip.text       = element_text(color = "black", face = "bold")
)


cat("[1/6] 读取 KEGG 代谢通路 GMT...\n")
gmt_path <- "/data/nas2/Software/miniconda3/envs/public_R/lib/R/library/scMetabolism/data/KEGG_metabolism_nc.gmt"

if (!file.exists(gmt_path)) {
  cat("[ERROR] GMT 文件不存在:", gmt_path, "\n")
  quit(status = 1)
}


parse_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)
  result <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) < 3) next
    pathway_name <- parts[1]
    genes <- parts[-(1:2)]
    result[[pathway_name]] <- genes
  }
  return(result)
}

gene_sets_list <- parse_gmt(gmt_path)
cat(sprintf("[INFO] 读取 %d 条代谢通路\n", length(gene_sets_list)))


cat("[2/6] 读取 seurat_annotated.rds...\n")
sc <- readRDS(file.path(WORK, "r.03_single_cell/seurat_annotated.rds"))
cat(sprintf("[INFO] Seurat 对象: %d 基因 × %d 细胞\n", nrow(sc), ncol(sc)))
cat("[INFO] 细胞类型:", paste(names(table(sc$cell_type)), collapse=", "), "\n")


cat("[3/6] 提取表达矩阵 (Seurat v5 兼容)...\n")
expr_mat <- tryCatch({
  SeuratObject::LayerData(sc, assay = "RNA", layer = "data")
}, error = function(e1) {
  tryCatch({
    sc@assays$RNA@data
  }, error = function(e2) {
    NULL
  })
})

if (is.null(expr_mat) || nrow(expr_mat) == 0) {
  cat("[ERROR] 无法提取表达矩阵\n")
  quit(status = 1)
}
cat(sprintf("[INFO] 表达矩阵: %d 基因 × %d 细胞\n", nrow(expr_mat), ncol(expr_mat)))


cat("[4/6] 运行 AUCell (buildRankings + calcAUC)...\n")


gene_sets_filtered <- lapply(gene_sets_list, function(gs) {
  intersect(gs, rownames(expr_mat))
})
gene_sets_filtered <- gene_sets_filtered[sapply(gene_sets_filtered, length) >= 5]
cat(sprintf("[INFO] 有效通路（≥5个基因）: %d 条\n", length(gene_sets_filtered)))

cat("[INFO] 构建基因排名矩阵（约需数分钟）...\n")
set.seed(42)
cell_rankings <- AUCell_buildRankings(
  expr_mat,
  nCores    = 4,
  plotStats = FALSE
)
cat("[INFO] 计算 AUC 分数...\n")
auc_obj    <- AUCell_calcAUC(gene_sets_filtered, cell_rankings, nCores = 4)
scores_mat <- t(getAUC(auc_obj))
cat(sprintf("[INFO] AUC 分数矩阵: %d 细胞 × %d 通路\n", nrow(scores_mat), ncol(scores_mat)))


cat("[5/6] 整合细胞类型信息并聚合...\n")
scores_df <- as.data.frame(scores_mat)
scores_df$cell_type <- sc$cell_type[rownames(scores_df)]


write.csv(scores_df, file.path(out_dir, "metabolism_scores_AUCell.csv"),
          row.names = TRUE)
cat("[INFO] 已保存 metabolism_scores_AUCell.csv\n")


agg_df <- scores_df %>%
  group_by(cell_type) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")


long_df <- agg_df %>%
  pivot_longer(-cell_type, names_to = "pathway", values_to = "mean_auc")

cat(sprintf("[INFO] 聚合后: %d 行（细胞类型 × 通路）\n", nrow(long_df)))


cat("[6/6] 绘制气泡图...\n")


pathway_var <- long_df %>%
  group_by(pathway) %>%
  summarise(var_auc = var(mean_auc, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(var_auc))

top_pathways <- head(pathway_var$pathway, 20)

plot_df <- long_df %>%
  filter(pathway %in% top_pathways) %>%
  mutate(
    pathway_short = ifelse(nchar(pathway) > 45,
                           paste0(substr(pathway, 1, 42), "..."),
                           pathway)
  ) %>%
  group_by(pathway) %>%
  mutate(
    pathway_mean = mean(mean_auc, na.rm = TRUE),
    pathway_sd   = sd(mean_auc, na.rm = TRUE),
    zscore       = ifelse(pathway_sd > 0, (mean_auc - pathway_mean) / pathway_sd, 0)
  ) %>%
  ungroup()

p <- ggplot(plot_df, aes(x = cell_type, y = pathway_short,
                         size = mean_auc, color = zscore)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(1, 8), name = "Mean AUC") +
  scale_color_gradient2(
    low      = "#3B4CC0",
    mid      = "#DDDDDD",
    high     = "#B40426",
    midpoint = 0,
    name     = "Z-score"
  ) +
  labs(
    title    = "Metabolic Activity Across Cell Types",
    subtitle = "Top 20 variable KEGG metabolic pathways (AUCell scores)",
    x        = "Cell Type",
    y        = "KEGG Metabolic Pathway"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9,  color = "black"),
    axis.text.y      = element_text(size = 8,  color = "black"),
    axis.title.x     = element_text(color = "black"),
    axis.title.y     = element_text(color = "black"),
    plot.title       = element_text(face = "bold", size = 13, color = "black"),
    plot.subtitle    = element_text(size = 9, color = "black"),
    legend.text      = element_text(color = "black"),
    legend.title     = element_text(color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    legend.position  = "right"
  )

save_plot(p,
          base_path = file.path(out_dir, "Fig28B_metabolism_bubble"),
          width = 14, height = 10)


if ("Condition" %in% colnames(sc@meta.data)) {
  cat("[INFO] 绘制 DCM vs Donor 分组图...\n")
  scores_df2 <- as.data.frame(scores_mat)
  scores_df2$cell_type <- sc$cell_type[rownames(scores_df2)]
  scores_df2$Condition <- sc$Condition[rownames(scores_df2)]

  agg_df2 <- scores_df2 %>%
    group_by(cell_type, Condition) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

  long_df2 <- agg_df2 %>%
    pivot_longer(-c(cell_type, Condition), names_to = "pathway", values_to = "mean_auc") %>%
    filter(pathway %in% top_pathways) %>%
    mutate(
      pathway_short = ifelse(nchar(pathway) > 45,
                             paste0(substr(pathway, 1, 42), "..."),
                             pathway)
    ) %>%
    group_by(pathway) %>%
    mutate(
      pathway_mean = mean(mean_auc, na.rm = TRUE),
      pathway_sd   = sd(mean_auc, na.rm = TRUE),
      zscore       = ifelse(pathway_sd > 0, (mean_auc - pathway_mean) / pathway_sd, 0)
    ) %>%
    ungroup()

  p2 <- ggplot(long_df2, aes(x = cell_type, y = pathway_short,
                             size = mean_auc, color = zscore)) +
    geom_point(alpha = 0.85) +
    facet_wrap(~Condition, ncol = 2) +
    scale_size_continuous(range = c(1, 7), name = "Mean AUC") +
    scale_color_gradient2(
      low = "#3B4CC0", mid = "#DDDDDD", high = "#B40426",
      midpoint = 0, name = "Z-score"
    ) +
    labs(
      title = "Metabolic Activity: DCM vs Donor",
      x = "Cell Type", y = "KEGG Metabolic Pathway"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8,  color = "black"),
      axis.text.y      = element_text(size = 7.5, color = "black"),
      axis.title.x     = element_text(color = "black"),
      axis.title.y     = element_text(color = "black"),
      plot.title       = element_text(face = "bold", size = 12, color = "black"),
      legend.text      = element_text(color = "black"),
      legend.title     = element_text(color = "black"),
      strip.text       = element_text(face = "bold", size = 11, color = "black"),
      strip.background = element_rect(fill = "#F0F4FF"),
      panel.grid.major = element_line(color = "grey90")
    )

  save_plot(p2,
            base_path = file.path(out_dir, "Fig28B_metabolism_DCMvsDonor"),
            width = 18, height = 10)
}

cat("\n[DONE] Figure 28B 代谢气泡图完成！\n")
