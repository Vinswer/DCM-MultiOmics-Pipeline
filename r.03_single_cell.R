rm(list = ls()); gc()
Sys.setenv(R_FUTURE_PLAN          = "sequential")
Sys.setenv(OMP_NUM_THREADS        = "1")
Sys.setenv(OPENBLAS_NUM_THREADS   = "1")
Sys.setenv(MKL_NUM_THREADS        = "1")
options(mc.cores                  = 1L)
options(Ncpus                     = 1L)
options(future.plan               = "sequential")
options(future.globals.maxSize    = 8000 * 1024^2)
options(future.globals.method.default = "ordered")
options(future.startup.script     = FALSE)
args <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4
cat(sprintf("[INFO] 使用核心数: %d\n", n_cores))




required_pkgs <- c("Seurat", "ggplot2", "dplyr", "patchwork", "tidyr",
                   "scales", "RColorBrewer", "ggrepel", "AUCell",
                   "GSEABase", "grDevices")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[SKIP] 缺少必需包: %s，退出。\n", pkg))
    quit(status = 0)
  }
}
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
  library(ggrepel)
  library(grDevices)
})


has_AUCell <- requireNamespace("AUCell", quietly = TRUE) &&
  requireNamespace("GSEABase", quietly = TRUE)
if (has_AUCell) {
  suppressPackageStartupMessages({
    library(AUCell)
    library(GSEABase)
  })
  cat("[INFO] AUCell 可用，将用于基因集富集分析。\n")
} else {
  cat("[WARN] AUCell/GSEABase 不可用，将用 AddModuleScore 替代。\n")
}

select  <- dplyr::select
filter  <- dplyr::filter
rename  <- dplyr::rename
mutate  <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
summarize <- dplyr::summarize
group_by  <- dplyr::group_by
ungroup   <- dplyr::ungroup
left_join <- dplyr::left_join
slice      <- dplyr::slice
slice_max  <- dplyr::slice_max
slice_min  <- dplyr::slice_min
pull       <- dplyr::pull
distinct   <- dplyr::distinct
n_distinct <- dplyr::n_distinct
cat("[INFO] dplyr 命名空间冲突已修复。\n")

cat(sprintf("[INFO] Seurat 版本: %s\n", packageVersion("Seurat")))


if (requireNamespace("future", quietly=TRUE)) {
  options(future.globals.method.default = "ordered")
  options(future.plan               = "sequential")
  options(future.startup.script     = FALSE)
  future::plan(future::sequential)
}
cat("[INFO] future 并行已禁用。\n")




SERVER_DIR  <- ""
INPUT_FILE  <- file.path(SERVER_DIR, "00_rawdata/singlecell_GSE183852/GSE183852_DCM_Cells.Robj.gz")
OUTPUT_DIR  <- file.path(SERVER_DIR, "03_single_cell")

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat(sprintf("[INFO] 已创建输出目录: %s\n", OUTPUT_DIR))
}




save_plot <- function(plot_obj, filename_base,
                      width = 10, height = 8, dpi = 150) {
  for (ext in c("png", "pdf")) {
    out_path <- file.path(OUTPUT_DIR, paste0(filename_base, ".", ext))
    if (ext == "png") {
      png(out_path, width = width, height = height,
          units = "in", res = dpi)
    } else {
      pdf(out_path, width = width, height = height)
    }
    print(plot_obj)
    dev.off()
    cat(sprintf("[INFO] 图形已保存: %s\n", out_path))
  }
}




cat("[INFO] 正在加载 HDCM Seurat 对象...\n")
if (!file.exists(INPUT_FILE)) {
  cat(sprintf("[ERROR] 输入文件不存在: %s\n", INPUT_FILE))
  quit(status = 1)
}

e <- new.env()
tmp_file <- tempfile(fileext = ".Robj")
system(paste("gunzip -c", shQuote(INPUT_FILE), ">", shQuote(tmp_file)))
load(tmp_file, envir = e)
file.remove(tmp_file)
sc <- e[["HDCM"]]
rm(e); gc()

cat(sprintf("[INFO] 对象加载完成: %d genes × %d cells\n",
            nrow(sc), ncol(sc)))


if (inherits(sc[["RNA"]], "Assay5")) {
  cat("[INFO] 检测到 Seurat v5 Assay5，执行 JoinLayers()...\n")
  sc[["RNA"]] <- JoinLayers(sc[["RNA"]])
}




n_genes_before <- nrow(sc)
n_cells_before <- ncol(sc)
cat(sprintf("[STAT] QC_before_genes: %d\n", n_genes_before))
cat(sprintf("[STAT] QC_before_cells: %d\n", n_cells_before))


if (!"percent.mito" %in% colnames(sc@meta.data)) {
  mito_genes <- grep("^MT-", rownames(sc), value = TRUE, ignore.case = TRUE)
  if (length(mito_genes) > 0) {
    sc[["percent.mito"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  } else {
    sc[["percent.mito"]] <- 0
  }
}




cat("[INFO] 绘制质控前 QC 小提琴图...\n")
qc_feats <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mito"),
                      colnames(sc@meta.data))

if (length(qc_feats) > 0) {
  qc_df_pre <- sc@meta.data[, qc_feats, drop = FALSE]
  qc_df_pre$cell    <- rownames(qc_df_pre)
  qc_df_pre$cluster <- "All cells"
  qc_long_pre <- tidyr::pivot_longer(qc_df_pre,
                                     cols      = all_of(qc_feats),
                                     names_to  = "feature",
                                     values_to = "value")
  p_qc_before <- ggplot(qc_long_pre,
                        aes(x = cluster, y = value, fill = feature)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.5) +
    facet_wrap(~ feature, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title   = element_text(size = 13, hjust = 0.5),
          legend.position = "none") +
    ggtitle(sprintf("QC Metrics Before Filtering  (cells=%d, genes=%d)",
                    n_cells_before, n_genes_before))

  save_plot(p_qc_before, "05_QC_violin_before",
            width = 4 * length(qc_feats), height = 5)
}










cat("[INFO] QC 标准(1): 过滤低覆盖基因（<3 个细胞表达）...\n")
counts_mat <- GetAssayData(sc, assay = "RNA", layer = "counts")
if (!inherits(counts_mat, "dgCMatrix") && !is.matrix(counts_mat)) {
  counts_mat <- tryCatch(
    GetAssayData(sc, assay = "RNA", slot = "counts"),
    error = function(e) NULL
  )
}
if (!is.null(counts_mat)) {
  genes_keep <- rowSums(counts_mat > 0) >= 3
  genes_keep_names <- rownames(counts_mat)[genes_keep]
  genes_keep_names <- intersect(genes_keep_names, rownames(sc))
  n_genes_removed  <- nrow(sc) - length(genes_keep_names)
  cat(sprintf("[INFO] 基因过滤：保留 %d / %d 基因（移除 %d 个低覆盖基因）\n",
              length(genes_keep_names), nrow(sc), n_genes_removed))
  sc <- sc[genes_keep_names, ]
} else {
  cat("[WARN] 无法获取 counts 矩阵，跳过基因层面过滤。\n")
}


cat("[INFO] QC 标准(2)(3)(4): 过滤低质量细胞...\n")
meta <- sc@meta.data


crit2 <- meta$nFeature_RNA >= 200 & meta$nFeature_RNA <= 4000

crit3 <- meta$percent.mito < 15

crit4 <- meta$nCount_RNA < 10000 & meta$nFeature_RNA > 200

cells_keep <- crit2 & crit3 & crit4
n_removed  <- sum(!cells_keep)

cat(sprintf("[INFO]   标准(2) 不通过（nFeature 超范围）: %d 细胞\n",
            sum(!crit2)))
cat(sprintf("[INFO]   标准(3) 不通过（mito ≥ 15%%）:      %d 细胞\n",
            sum(!crit3)))
cat(sprintf("[INFO]   标准(4) 不通过（count≥10000 或 feat≤200）: %d 细胞\n",
            sum(!crit4)))
cat(sprintf("[INFO] 共移除 %d 细胞，保留 %d 细胞。\n",
            n_removed, sum(cells_keep)))

sc <- sc[, cells_keep]

n_genes_after <- nrow(sc)
n_cells_after <- ncol(sc)
cat(sprintf("[STAT] QC_after_genes: %d\n", n_genes_after))
cat(sprintf("[STAT] QC_after_cells: %d\n", n_cells_after))
cat(sprintf("[STAT] QC_removed_cells: %d\n",
            n_cells_before - n_cells_after))




if (length(qc_feats) > 0) {
  qc_df_post <- sc@meta.data[, qc_feats, drop = FALSE]
  qc_df_post$cell    <- rownames(qc_df_post)
  qc_df_post$cluster <- "All cells"
  qc_long_post <- tidyr::pivot_longer(qc_df_post,
                                      cols      = all_of(qc_feats),
                                      names_to  = "feature",
                                      values_to = "value")
  p_qc_after <- ggplot(qc_long_post,
                       aes(x = cluster, y = value, fill = feature)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.5) +
    facet_wrap(~ feature, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title   = element_text(size = 13, hjust = 0.5),
          legend.position = "none") +
    ggtitle(sprintf("QC Metrics After Filtering  (cells=%d, genes=%d)",
                    n_cells_after, n_genes_after))

  save_plot(p_qc_after, "05_QC_violin_after",
            width = 4 * length(qc_feats), height = 5)


  qc_long_pre$phase  <- sprintf("Before (n=%d)", n_cells_before)
  qc_long_post$phase <- sprintf("After  (n=%d)", n_cells_after)
  qc_comb <- rbind(qc_long_pre[, c("phase","feature","value")],
                   qc_long_post[, c("phase","feature","value")])
  qc_comb$phase <- factor(qc_comb$phase)

  p_qc_comb <- ggplot(qc_comb,
                      aes(x = phase, y = value, fill = phase)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 0.1,
                 fill = "white", alpha = 0.5) +
    facet_wrap(~ feature, scales = "free_y") +
    scale_fill_manual(values = c("#E8836A", "#6AA8E8")) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          plot.title   = element_text(size = 13, hjust = 0.5),
          legend.position = "none") +
    ggtitle("QC Metrics: Before vs After Filtering")

  save_plot(p_qc_comb, "05_QC_violin_comparison",
            width = 4 * length(qc_feats), height = 5)
}






if (!"SCT" %in% names(sc@assays)) {
  cat("[INFO] 运行 SCTransform...\n")
  sc <- SCTransform(sc, vars.to.regress = "percent.mito",
                    verbose = FALSE)
}


if (!"pca" %in% names(sc@reductions)) {
  cat("[INFO] 运行 RunPCA...\n")
  sc <- RunPCA(sc, assay = "SCT", verbose = FALSE)
}




cat("[INFO] 绘制高变基因图...\n")


hvg_assay <- if ("SCT" %in% names(sc@assays)) "SCT" else "RNA"










































if (!exists("feat_info")) feat_info <- NULL
cat("[INFO] 绘制高变基因图...\n")

hvg_assay <- if ("SCT" %in% names(sc@assays)) "SCT" else "RNA"
feat_info  <- NULL


if (hvg_assay == "SCT") {
  feat_info <- tryCatch({
    sct_model <- sc@assays$SCT@SCTModel.list[[1]]
    gm        <- sct_model@feature.attributes
    data.frame(
      gene           = rownames(gm),
      avg_expression = gm$gmean,
      std_variance   = gm$residual_variance,
      is_variable    = rownames(gm) %in% VariableFeatures(sc, assay = "SCT"),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}


if (is.null(feat_info)) {
  tryCatch({
    sc_tmp    <- FindVariableFeatures(sc, assay = "RNA",
                                      selection.method = "vst",
                                      nfeatures = 3000, verbose = FALSE)
    hvg_df    <- HVFInfo(sc_tmp, assay = "RNA")
    feat_info <- data.frame(
      gene           = rownames(hvg_df),
      avg_expression = hvg_df$mean,
      std_variance   = hvg_df$variance.standardized,
      is_variable    = rownames(hvg_df) %in% VariableFeatures(sc_tmp, assay = "RNA"),
      stringsAsFactors = FALSE
    )
    rm(sc_tmp); gc()
  }, error = function(e) {
    cat(sprintf("[WARN] 无法提取高变基因信息: %s\n", conditionMessage(e)))
  })
}

if (!is.null(feat_info)) {
  feat_info <- feat_info[!is.na(feat_info$avg_expression) &
                           !is.na(feat_info$std_variance), ]

  n_var_count    <- sum(feat_info$is_variable)
  n_nonvar_count <- sum(!feat_info$is_variable)


  label_df <- feat_info %>%
    dplyr::arrange(desc(std_variance)) %>%
    dplyr::slice(1:10)

  p_hvg <- ggplot(feat_info,
                  aes(x     = avg_expression,
                      y     = std_variance,
                      color = is_variable)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "red"),
      labels = c(sprintf("Non-variable count: %d", n_nonvar_count),
                 sprintf("Variable count: %d",     n_var_count))
    ) +
    scale_x_log10(labels = scales::label_scientific()) +
    geom_text_repel(data = label_df,
                    aes(label = gene),
                    size = 3, color = "black",
                    max.overlaps = 20, seed = 42) +
    labs(x     = "Average Expression",
         y     = "Standardized Variance",
         color = NULL) +
    theme_classic() +
    theme(plot.title      = element_text(hjust = 0.5, size = 13),
          legend.position = c(0.85, 0.85),
          legend.text     = element_text(size = 9))

  save_plot(p_hvg, "05_HVG_plot", width = 10, height = 8)
}




cat("[INFO] 绘制 ElbowPlot（碎石图，展示前 20 个 PC）...\n")

pca_std         <- sc@reductions$pca@stdev
n_pcs_available <- length(pca_std)
N_SHOW          <- min(20L, n_pcs_available)


pct_var <- pca_std^2 / sum(pca_std^2) * 100
cum_var <- cumsum(pct_var)
delta   <- diff(pca_std[1:N_SHOW])

threshold     <- mean(abs(delta)) * 0.5
elbow_idx_raw <- which(abs(delta) < threshold)[1]
if (is.na(elbow_idx_raw) || elbow_idx_raw < 5)
  elbow_idx_raw <- min(15L, N_SHOW)
selected_dims <- min(elbow_idx_raw + 1L, N_SHOW)

cat(sprintf("[INFO] PCA 可用维度数: %d\n", n_pcs_available))
cat(sprintf("[INFO] 自动选择 PC 数量: %d (累计方差 %.1f%%)\n",
            selected_dims, cum_var[selected_dims]))
cat(sprintf("[STAT] Selected_PCs: %d\n", selected_dims))

elbow_df <- data.frame(PC  = seq_len(N_SHOW),
                       std = pca_std[seq_len(N_SHOW)])

p_elbow <- ggplot(elbow_df, aes(x = PC, y = std)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_vline(xintercept = selected_dims,
             linetype = "dashed", color = "red", linewidth = 0.9) +
  annotate("text",
           x     = selected_dims + 0.3,
           y     = max(elbow_df$std) * 0.92,
           label = sprintf("Selected: %d PCs\n(Cum. var. %.1f%%)",
                           selected_dims, cum_var[selected_dims]),
           hjust = 0, color = "red", size = 4) +
  scale_x_continuous(breaks = seq(1, N_SHOW, by = 2)) +
  labs(x     = "PC",
       y     = "Standard Deviation",
       title = sprintf("Elbow Plot — Top %d PCs", N_SHOW)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

save_plot(p_elbow, "05_ElbowPlot", width = 9, height = 6)


cat("[INFO] 绘制 PCA 散点图...\n")
tryCatch({
  n_pca_cols <- min(2, ncol(sc@reductions$pca@cell.embeddings))
  pca_emb <- as.data.frame(sc@reductions$pca@cell.embeddings[, 1:n_pca_cols])
  colnames(pca_emb) <- c("PC_1", "PC_2")

  color_by <- if ("orig.ident" %in% colnames(sc@meta.data)) "orig.ident" else "seurat_clusters"
  pca_emb$group <- sc@meta.data[[color_by]]

  n_grp   <- length(unique(pca_emb$group))
  grp_pal <- colorRampPalette(brewer.pal(min(n_grp, 12), "Paired"))(n_grp)

  p_pca <- ggplot(pca_emb, aes(x = PC_1, y = PC_2, color = group)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = grp_pal, name = color_by) +
    labs(title = "PCA Plot",
         x = "PC_1", y = "PC_2") +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0, size = 14),
      legend.title    = element_text(size = 10),
      legend.text     = element_text(size = 9),
      legend.key.size = unit(0.4, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

  save_plot(p_pca, "05_PCA_scatter", width = 10, height = 8)
}, error = function(e) {
  cat(sprintf("[WARN] PCA 散点图绘制失败: %s\n", conditionMessage(e)))
})




RESOLUTION  <- 0.5
cluster_col <- paste0("SCT_snn_res.", RESOLUTION)

if (!"umap" %in% names(sc@reductions)) {
  cat(sprintf("[INFO] 运行 FindNeighbors + FindClusters (res=%.1f) + RunUMAP...\n",
              RESOLUTION))
  sc <- FindNeighbors(sc, reduction = "pca",
                      dims = 1:selected_dims, verbose = FALSE)
  sc <- FindClusters(sc, resolution = RESOLUTION, verbose = FALSE)
  sc <- RunUMAP(sc, reduction = "pca",
                dims = 1:selected_dims,
                n.components = 2L,
                seed.use = 42, verbose = FALSE)
  cat("[INFO] UMAP 运行完成。\n")
} else {
  cat("[INFO] 已检测到 UMAP，跳过 RunUMAP。\n")
  if (!cluster_col %in% colnames(sc@meta.data)) {
    sc <- FindNeighbors(sc, reduction = "pca",
                        dims = 1:selected_dims, verbose = FALSE)
    sc <- FindClusters(sc, resolution = RESOLUTION, verbose = FALSE)
  }
}


if (!cluster_col %in% colnames(sc@meta.data)) {
  avail_cols <- grep("SCT_snn_res", colnames(sc@meta.data), value = TRUE)
  if (length(avail_cols) > 0) {
    cluster_col <- avail_cols[1]
    cat(sprintf("[INFO] 使用聚类列: %s\n", cluster_col))
  } else {
    cat("[ERROR] 找不到任何聚类结果列。\n")
    quit(status = 1)
  }
}

Idents(sc) <- cluster_col
sc$seurat_clusters <- sc@meta.data[[cluster_col]]

n_clusters    <- length(levels(Idents(sc)))
n_total_cells <- ncol(sc)

cat(sprintf("[STAT] SC_total_cells: %d\n",   n_total_cells))
cat(sprintf("[STAT] SC_cluster_count: %d\n", n_clusters))
cat(sprintf("[STAT] SC_resolution: %.1f\n",  RESOLUTION))




cat("[INFO] 绘制 UMAP 聚类图...\n")
p_umap_cluster <- DimPlot(sc, reduction = "umap",
                          label = TRUE, label.size = 4,
                          pt.size = 0.3, repel = TRUE) +
  ggtitle("UMAP — Clusters (SCT_snn_res.0.5)") +
  theme(plot.title   = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 10),
        legend.text  = element_text(size = 8))

save_plot(p_umap_cluster, "05_UMAP_clusters", width = 12, height = 9)




cat("[INFO] 开始细胞类型注释（参考 PMID:35959412）...\n")



cell_markers <- list(
  Cardiomyocytes = c("TNNT2","MYH7","MYL2","ACTC1","TNNI3","TTN","RYR2","PLN","MYBPC3","TNNC1"),
  Fibroblasts    = c("DCN","LUM","COL1A1","POSTN","TCF21","PDGFRA","MFAP4","CFD","NEGR1","C7"),
  Endothelial    = c("PECAM1","CDH5","VWF","CLDN5","FLT1","KDR","EMCN","ESAM","TIE1","ENG"),
  Pericytes_SMC  = c("ACTA2","MYH11","TAGLN","CALD1","CNN1","PDGFRB","RGS5","NOTCH3","MCAM","ABCC9"),
  Macrophages    = c("CD68","LYZ","CSF1R","C1QA","C1QB","C1QC","MRC1","CD163",
                     "S100A8","S100A9","VCAN","FCN1","CD14","VSIG4","MARCO","FCER1G","FCN3"),
  T_cells        = c("CD3D","CD3E","CD3G","CD4","CD8A","TRAC","TRBC1","IL7R","CCR7","TCF7"),
  NK_Bcells      = c("MS4A1","CD19","CD79A","GNLY","NKG7","KLRD1","GZMA","GZMB","NCAM1","FCGR3A"),
  Neural         = c("NRXN1","SNAP25","SYP","ENO2","GFAP","S100B","NCAM1","L1CAM","NEFM","NEFL")
)


cell_markers_filt <- lapply(cell_markers, function(g) {
  intersect(g, rownames(sc))
})

cell_markers_filt <- Filter(function(g) length(g) >= 2, cell_markers_filt)

cat("[INFO] 细胞类型对应有效基因数:\n")
for (ct in names(cell_markers_filt)) {
  cat(sprintf("  %-20s: %d genes\n", ct, length(cell_markers_filt[[ct]])))
}


cat("[INFO] 计算各细胞类型模块分数...\n")
score_cols <- c()
for (ct in names(cell_markers_filt)) {
  col_name <- paste0("score_", gsub(" |/", "_", ct))
  sc <- tryCatch(
    AddModuleScore(sc,
                   features = list(cell_markers_filt[[ct]]),
                   name     = col_name,
                   ctrl     = 50,
                   seed     = 42),
    error = function(e) {
      cat(sprintf("[WARN] AddModuleScore 出错 (%s): %s\n",
                  ct, conditionMessage(e)))
      sc
    }
  )

  actual_col <- paste0(col_name, "1")
  if (actual_col %in% colnames(sc@meta.data)) {
    score_cols <- c(score_cols, setNames(actual_col, ct))
  }
}


if (length(score_cols) > 0) {
  score_mat <- as.matrix(sc@meta.data[, score_cols, drop = FALSE])
  colnames(score_mat) <- names(score_cols)


  best_type      <- apply(score_mat, 1, which.max)
  best_type_name <- colnames(score_mat)[best_type]
  sc$cell_type   <- best_type_name


  cluster_annotation <- sc@meta.data %>%
    dplyr::group_by(seurat_clusters, cell_type) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(seurat_clusters, cell_type)

  cat("[INFO] Cluster → Cell Type 注释:\n")
  print(as.data.frame(cluster_annotation))


  anno_map <- setNames(as.character(cluster_annotation$cell_type),
                       as.character(cluster_annotation$seurat_clusters))
  cell_anno_vec <- anno_map[as.character(sc@meta.data$seurat_clusters)]
  names(cell_anno_vec) <- colnames(sc)
  sc <- AddMetaData(sc, metadata = cell_anno_vec, col.name = "annotated_cell_type")


  write.table(cluster_annotation,
              file.path(OUTPUT_DIR, "05_cluster_annotation.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("[INFO] 输出各细胞类型 Marker Gene 汇总...\n")

  marker_summary_lines <- c()

  for (ct in names(cell_markers_filt)) {

    ct_cells <- sum(sc$annotated_cell_type == ct, na.rm = TRUE)
    genes_str <- paste(cell_markers_filt[[ct]], collapse = "、")
    line <- sprintf("%s（%d个，%s）", ct, ct_cells, genes_str)
    marker_summary_lines <- c(marker_summary_lines, line)
    cat(sprintf("[MARKER] %s\n", line))
  }

  writeLines(marker_summary_lines,
             file.path(OUTPUT_DIR, "05_celltype_marker_summary.txt"))
  cat("[INFO] Marker Gene 汇总已保存: 05_celltype_marker_summary.txt\n")

  total_cells <- ncol(sc)

  celltype_summary <- sc@meta.data %>%
    dplyr::group_by(annotated_cell_type) %>%
    dplyr::summarise(cell_count = n(), .groups = "drop") %>%
    dplyr::arrange(desc(cell_count)) %>%
    dplyr::mutate(
      proportion = round(cell_count / total_cells * 100, 4),
      proportion_pct = paste0(round(cell_count / total_cells * 100, 2), "%")
    )

  write.csv(celltype_summary,
            file.path(OUTPUT_DIR, "05_celltype_summary.csv"),
            row.names = FALSE)
  cat("[INFO] 细胞类型汇总 CSV 已保存\n")
  cat("[INFO] 各细胞类型统计:\n")
  print(as.data.frame(celltype_summary))


  cat("[INFO] 绘制细胞类型注释 UMAP（图2类型）...\n")


  ct_types   <- unique(sc$annotated_cell_type)
  n_types    <- length(ct_types)
  ct_colors  <- setNames(
    colorRampPalette(brewer.pal(min(n_types, 12), "Set3"))(n_types),
    ct_types
  )

  Idents(sc) <- "annotated_cell_type"
  p_anno_umap <- DimPlot(sc, reduction = "umap",
                         group.by = "annotated_cell_type",
                         label = TRUE, label.size = 4,
                         pt.size = 0.3, repel = TRUE,
                         cols = ct_colors) +
    ggtitle("Cell Type Annotation") +
    theme(plot.title   = element_text(hjust = 0, size = 14),
          legend.title = element_text(size = 10),
          legend.text  = element_text(size = 9))

  save_plot(p_anno_umap, "05_UMAP_celltype_annotation",
            width = 13, height = 10)









  Idents(sc) <- cluster_col
  avg_expr <- AverageExpression(sc, assay = "SCT", layer = "data", verbose = FALSE)$SCT


  cluster_annotation <- lapply(colnames(avg_expr), function(cl) {
    cl_expr <- avg_expr[, cl, drop = TRUE]
    scores <- sapply(cell_markers_filt, function(genes) {
      g <- intersect(genes, names(cl_expr))
      if (length(g) == 0) return(0)
      mean(cl_expr[g])
    })
    data.frame(
      seurat_clusters = cl,
      cell_type = names(which.max(scores)),
      stringsAsFactors = FALSE
    )
  })
  cluster_annotation <- do.call(rbind, cluster_annotation)

} else {
  cat("[WARN] 未能计算任何细胞类型模块分数，跳过细胞注释图。\n")
  sc$annotated_cell_type <- "Unknown"
}










































































cat("[INFO] 绘制改进版细胞类型标志基因气泡图...\n")

DefaultAssay(sc) <- "SCT"
Idents(sc) <- "annotated_cell_type"




markers_anno_file <- file.path(OUTPUT_DIR, "05_markers_by_celltype.rds")

if (file.exists(markers_anno_file)) {
  cat("[INFO] 检测到已有 celltype markers 文件，直接加载。\n")
  markers_anno <- readRDS(markers_anno_file)
} else {
  cat("[INFO] 运行 FindAllMarkers（按细胞类型）...\n")
  markers_anno <- tryCatch({
    FindAllMarkers(
      sc,
      assay          = "SCT",
      only.pos       = TRUE,
      min.pct        = 0.25,
      logfc.threshold = 0.5,
      test.use       = "wilcox",
      verbose        = FALSE
    )
  }, error = function(e) {
    cat(sprintf("[WARN] FindAllMarkers 出错: %s\n", conditionMessage(e)))
    NULL
  })
  if (!is.null(markers_anno)) {
    saveRDS(markers_anno, markers_anno_file)
  }
}




if (!is.null(markers_anno) && nrow(markers_anno) > 0) {






  markers_filtered <- markers_anno %>%
    dplyr::filter(
      p_val_adj       < 0.05,
      avg_log2FC      > 0.5,
      pct.1           > 0.30,
      pct.2           < 0.50
    ) %>%
    dplyr::mutate(
      specificity = pct.1 / (pct.2 + 0.01)
    )


  top_markers_anno <- markers_filtered %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = specificity, n = 4, with_ties = FALSE) %>%
    dplyr::ungroup()

  cat("[INFO] 筛选后各细胞类型 Top Marker:\n")
  print(top_markers_anno %>%
          dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2, specificity))

  dot_genes_final <- unique(top_markers_anno$gene)
  dot_genes_final <- intersect(dot_genes_final, rownames(sc))

  cat(sprintf("[INFO] 最终用于气泡图的基因数: %d\n", length(dot_genes_final)))




  if (length(dot_genes_final) > 0) {
    tryCatch({

      expr_mat <- FetchData(sc, vars = dot_genes_final)
      expr_mat$cell_type <- sc$annotated_cell_type

      dot_long <- tidyr::pivot_longer(
        expr_mat,
        cols      = all_of(dot_genes_final),
        names_to  = "gene",
        values_to = "expr"
      )

      dot_sum <- dot_long %>%
        dplyr::group_by(cell_type, gene) %>%
        dplyr::summarise(
          avg_expr = mean(expm1(expr)),
          pct_expr = mean(expr > 0) * 100,
          .groups  = "drop"
        )




      gene_order <- top_markers_anno %>%
        dplyr::arrange(cluster) %>%
        dplyr::pull(gene) %>%
        unique()
      gene_order <- intersect(gene_order, dot_genes_final)

      ct_order <- names(sort(table(sc$annotated_cell_type), decreasing = TRUE))

      dot_sum$gene      <- factor(dot_sum$gene,      levels = gene_order)
      dot_sum$cell_type <- factor(dot_sum$cell_type, levels = rev(ct_order))






      gene_group_df <- top_markers_anno %>%
        dplyr::arrange(cluster) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(n_genes = n_distinct(gene), .groups = "drop") %>%
        dplyr::mutate(end_pos = cumsum(n_genes),
                      start_pos = end_pos - n_genes + 1)


      vline_pos <- gene_group_df$end_pos[-nrow(gene_group_df)] + 0.5


      q95 <- quantile(dot_sum$avg_expr, 0.95, na.rm = TRUE)
      dot_sum$avg_expr_capped <- pmin(dot_sum$avg_expr, q95)

      p_dot_improved <- ggplot(dot_sum,
                               aes(x     = gene,
                                   y     = cell_type,
                                   size  = pct_expr,
                                   color = avg_expr_capped)) +
        geom_point() +

        geom_vline(xintercept = vline_pos,
                   color      = "grey70",
                   linewidth  = 0.4,
                   linetype   = "dashed") +
        scale_size_continuous(
          range  = c(0.5, 8),
          breaks = c(0, 25, 50, 75, 100),
          name   = "% Expressed"
        ) +
        scale_color_gradientn(
          colors = c("lightgrey", "#FFD0B0", "#FF6B35", "#C0392B"),
          name   = "Avg\nExpression\n(capped 95%)"
        ) +

        annotate(
          "text",
          x     = (gene_group_df$start_pos + gene_group_df$end_pos) / 2,
          y     = length(unique(dot_sum$cell_type)) + 0.8,
          label = gene_group_df$cluster,
          size  = 2.8,
          angle = 30,
          hjust = 0.5,
          color = "grey30"
        ) +
        ggtitle("Cell Type Marker Genes — Improved DotPlot\n(Data-driven, specificity-filtered)") +
        theme_bw() +
        theme(
          axis.text.x     = element_text(size = 8, angle = 45,
                                         vjust = 1, hjust = 1),
          axis.text.y     = element_text(size = 10),
          axis.title      = element_blank(),
          plot.title      = element_text(hjust = 0.5, size = 12),
          panel.grid.major = element_line(color = "grey92"),
          legend.title    = element_text(size = 9),
          legend.text     = element_text(size = 8),
          plot.margin     = margin(t = 30, r = 10, b = 10, l = 10)
        )

      save_plot(p_dot_improved, "05_dotplot_celltype_improved",
                width  = max(14, length(dot_genes_final) * 0.55),
                height = 8)

      cat("[INFO] 改进版气泡图已保存。\n")

    }, error = function(e) {
      cat(sprintf("[WARN] 改进版气泡图绘制出错: %s\n", conditionMessage(e)))
    })
  }

} else {



  cat("[WARN] FindAllMarkers 未返回结果，使用筛选后的预设 marker。\n")

  representative_markers_v2 <- list(
    Cardiomyocytes = c("TNNT2", "MYH7",   "MYL2",  "MYBPC3"),
    Fibroblasts    = c("DCN",   "LUM",    "COL1A1","POSTN"),
    Endothelial    = c("PECAM1","CDH5",   "VWF",   "CLDN5"),
    Pericytes_SMC  = c("ACTA2", "MYH11",  "RGS5",  "PDGFRB"),
    Myeloid        = c("C1QA",  "C1QB",   "VSIG4", "MARCO"),
    Monocyte       = c("S100A8","VCAN",   "FCN1",  "CLEC12A"),
    T_cells        = c("CD3D",  "IL7R",   "TCF7",  "CCR7"),
    NK_Bcells      = c("GNLY",  "NKG7",   "CD79A", "MS4A1"),
    Neural         = c("NRXN1", "SNAP25", "SYP",   "NEFL")
  )


  dot_genes_fallback <- unique(unlist(representative_markers_v2))
  dot_genes_fallback <- intersect(dot_genes_fallback, rownames(sc))

  if (length(dot_genes_fallback) > 0) {
    p_dot_fallback <- DotPlot(
      sc,
      features    = dot_genes_fallback,
      group.by    = "annotated_cell_type",
      dot.scale   = 8,
      col.min     = 0
    ) +
      scale_color_gradientn(
        colors = c("lightgrey","#FFD0B0","#FF6B35","#C0392B")
      ) +
      RotatedAxis() +
      ggtitle("Cell Type Marker Genes — DotPlot (Fallback)") +
      theme(plot.title = element_text(hjust = 0.5, size = 12))

    save_plot(p_dot_fallback, "05_dotplot_celltype_fallback",
              width = max(14, length(dot_genes_fallback) * 0.5),
              height = 8)
  }
}


Idents(sc) <- cluster_col




cat("[INFO] 开始心力衰竭 vs 对照细胞丰度分析...\n")



group_candidates <- c("group", "Group", "condition", "Condition",
                      "disease", "Disease", "Status", "status",
                      "sample_type", "orig.ident", "SampleType")
group_col <- NULL
for (cand in group_candidates) {
  if (cand %in% colnames(sc@meta.data)) {
    vals <- unique(sc@meta.data[[cand]])

    if (any(grepl("HF|heart.fail|DCM|fail|disease|case",
                  vals, ignore.case = TRUE)) ||
        any(grepl("ctrl|Donor|normal|healthy|non",
                  vals, ignore.case = TRUE))) {
      group_col <- cand
      break
    }
  }
}

if (is.null(group_col)) {

  if ("orig.ident" %in% colnames(sc@meta.data)) {
    group_col <- "orig.ident"
    cat(sprintf("[INFO] 未找到标准分组列，使用: %s\n", group_col))
  } else {
    cat("[WARN] 找不到分组列，跳过心力衰竭 vs 对照分析。\n")
    group_col <- NULL
  }
}

if (!is.null(group_col)) {
  group_vals <- unique(sc@meta.data[[group_col]])
  cat(sprintf("[INFO] 分组列: %s，取值: %s\n",
              group_col, paste(group_vals, collapse = ", ")))




  hf_patterns  <- c("DCM","HF","heart.fail","fail","disease","case")
  ctrl_patterns <- c("Donor","donor","ctrl","Donor","normal","healthy","non")

  hf_levels   <- group_vals[grepl(paste(hf_patterns,  collapse="|"),
                                  group_vals, ignore.case=TRUE)]
  ctrl_levels <- group_vals[grepl(paste(ctrl_patterns, collapse="|"),
                                  group_vals, ignore.case=TRUE)]

  if (length(hf_levels) == 0 || length(ctrl_levels) == 0) {

    ident_counts <- table(sc@meta.data[[group_col]])
    med_n  <- median(ident_counts)
    hf_levels   <- names(ident_counts)[ident_counts >  med_n]
    ctrl_levels <- names(ident_counts)[ident_counts <= med_n]
    cat(sprintf("[INFO] 自动推断 HF 组: %s\n",
                paste(hf_levels, collapse=", ")))
    cat(sprintf("[INFO] 自动推断 Donor 组: %s\n",
                paste(ctrl_levels, collapse=", ")))
  }


  sc@meta.data$disease_group <- ifelse(
    sc@meta.data[[group_col]] %in% hf_levels,   "DCM",
    ifelse(sc@meta.data[[group_col]] %in% ctrl_levels, "Donor", NA)
  )


  cat(sprintf("[INFO] DCM cells: %d\n",
              sum(sc$disease_group == "DCM")))
  cat(sprintf("[INFO] Donor cells:      %d\n",
              sum(sc$disease_group == "Donor")))




  cat("[INFO] 绘制各细胞类型占比直方图（Figure21A）...\n")

  prop_df <- sc@meta.data %>%
    group_by(disease_group, annotated_cell_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(disease_group) %>%
    mutate(prop = n / sum(n) * 100) %>%
    ungroup()

  ct_pal <- colorRampPalette(
    brewer.pal(min(length(unique(prop_df$annotated_cell_type)), 12), "Set3")
  )(length(unique(prop_df$annotated_cell_type)))

  p_prop_bar <- ggplot(prop_df,
                       aes(x = disease_group, y = prop,
                           fill = annotated_cell_type)) +
    geom_bar(stat = "identity", position = "stack", color = "white",
             linewidth = 0.2) +
    scale_fill_manual(values = ct_pal, name = "Cell Type") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Group", y = "Proportion (%)",
         title = "Cell Type Proportions: DCM vs Donor") +
    theme_bw() +
    theme(plot.title   = element_text(hjust = 0.5, size = 13),
          axis.title   = element_text(size = 11),
          legend.text  = element_text(size = 9))

  save_plot(p_prop_bar, "05_Figure21A_celltype_proportion",
            width = 8, height = 7)




  cat("[INFO] 执行 Wilcoxon 秩和检验（细胞丰度差异）...\n")


  sample_col <- if ("sample" %in% colnames(sc@meta.data)) "sample"
  else if ("orig.ident" %in% colnames(sc@meta.data)) "orig.ident"
  else group_col

  cell_count_per_sample <- sc@meta.data %>%
    group_by(.data[[sample_col]], disease_group, annotated_cell_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[sample_col]], disease_group) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  cell_types_all <- unique(cell_count_per_sample$annotated_cell_type)

  wilcox_res <- lapply(cell_types_all, function(ct) {
    df_ct <- cell_count_per_sample %>%
      dplyr::filter(annotated_cell_type == ct)
    hf_vals   <- df_ct$prop[df_ct$disease_group == "DCM"]
    ctrl_vals <- df_ct$prop[df_ct$disease_group == "Donor"]
    if (length(hf_vals) < 2 || length(ctrl_vals) < 2)
      return(NULL)
    wt <- wilcox.test(hf_vals, ctrl_vals, exact = FALSE)
    data.frame(
      cell_type   = ct,
      p_value     = wt$p.value,
      mean_HF     = mean(hf_vals),
      mean_ctrl   = mean(ctrl_vals),
      log2FC_prop = log2((mean(hf_vals) + 1e-6) / (mean(ctrl_vals) + 1e-6)),
      stringsAsFactors = FALSE
    )
  })
  wilcox_df <- do.call(rbind, Filter(Negate(is.null), wilcox_res))
  wilcox_df$significant <- wilcox_df$p_value < 0.05

  write.table(wilcox_df,
              file.path(OUTPUT_DIR, "05_wilcoxon_celltype_abundance.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("[INFO] Wilcoxon 检验结果:\n")
  print(wilcox_df[order(wilcox_df$p_value), ])


  diff_cells <- wilcox_df$cell_type[wilcox_df$significant]
  cat(sprintf("[STAT] 差异显著细胞类型 (p<0.05): %s\n",
              if (length(diff_cells) > 0)
                paste(diff_cells, collapse = ", ")
              else "none"))


  p_volcano_abund <- ggplot(wilcox_df,
                            aes(x = log2FC_prop, y = -log10(p_value),
                                color = significant,
                                label = cell_type)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#E8534A"),
                       labels = c("ns", "p < 0.05")) +
    geom_text_repel(size = 3.5, max.overlaps = 20) +
    labs(x     = "log2 Fold Change (Proportion)",
         y     = "-log10(p-value)",
         title = "Cell Type Abundance: DCM vs Donor",
         color = "Significance") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 13))

  save_plot(p_volcano_abund, "05_Figure21A_abundance_volcano",
            width = 9, height = 7)




  cat("[INFO] 鉴定吞噬调节核心关键细胞类型...\n")

  phagocytic_genes <- c("FCN3", "FCER1G", "MAP2K1")
  phago_avail      <- intersect(phagocytic_genes, rownames(sc))

  if (length(phago_avail) < 1) {
    cat("[WARN] 吞噬调节基因均不在数据集中，跳过核心细胞鉴定。\n")
  } else {
    cat(sprintf("[INFO] 可用吞噬调节基因: %s\n",
                paste(phago_avail, collapse=", ")))


    cells_to_analyze <- if (length(diff_cells) > 0) diff_cells else cell_types_all

    sc_diff <- if (length(diff_cells) > 0) {
      subset(sc, annotated_cell_type %in% diff_cells)
    } else sc

    sc_diff <- AddModuleScore(sc_diff,
                              features = list(phago_avail),
                              name     = "phago_score",
                              ctrl     = 50, seed = 42)


    if (has_AUCell && length(phago_avail) >= 2) {
      cat("[INFO] 使用 AUCell 计算吞噬调节基因集富集分数...\n")
      if (inherits(sc_diff[["RNA"]], "Assay5")) {
        sc_diff[["RNA"]] <- JoinLayers(sc_diff[["RNA"]])
      }
      tryCatch({
        expr_raw <- GetAssayData(sc_diff, assay = "RNA", layer = "counts")
        if (inherits(expr_raw, "dgCMatrix") || is.matrix(expr_raw)) {
          gs   <- GeneSet(phago_avail, setName = "Phagocytic_regulation")
          gscl <- GeneSetCollection(list(gs))

          cells_rankings  <- AUCell_buildRankings(expr_raw,
                                                  plotStats = FALSE,
                                                  verbose   = FALSE)
          cells_AUC       <- AUCell_calcAUC(gscl, cells_rankings,
                                            aucMaxRank = ceiling(0.05 * nrow(expr_raw)),
                                            verbose    = FALSE)
          auc_vals        <- as.numeric(getAUC(cells_AUC)[1, ])
          sc_diff$phago_AUC <- auc_vals
          cat("[INFO] AUCell 富集分数计算完成。\n")
        }
      }, error = function(e) {
        cat(sprintf("[WARN] AUCell 计算出错: %s\n", conditionMessage(e)))
      })
    }


    score_use <- if ("phago_AUC" %in% colnames(sc_diff@meta.data))
      "phago_AUC"
    else "phago_score1"


    phago_score_by_ct <- sc_diff@meta.data %>%
      group_by(annotated_cell_type, disease_group) %>%
      summarise(mean_score = mean(.data[[score_use]], na.rm = TRUE),
                .groups = "drop")


    phago_wilcox <- lapply(cells_to_analyze, function(ct) {
      df_ct  <- sc_diff@meta.data[sc_diff$annotated_cell_type == ct, ]
      hf_s   <- df_ct[[score_use]][df_ct$disease_group == "DCM"]
      ctrl_s <- df_ct[[score_use]][df_ct$disease_group == "Donor"]
      if (length(hf_s) < 2 || length(ctrl_s) < 2) return(NULL)
      wt <- wilcox.test(hf_s, ctrl_s, exact = FALSE)
      data.frame(cell_type  = ct,
                 p_value    = wt$p.value,
                 mean_HF    = mean(hf_s,   na.rm = TRUE),
                 mean_ctrl  = mean(ctrl_s, na.rm = TRUE),
                 log2FC     = log2((mean(hf_s,   na.rm=TRUE) + 1e-8) /
                                     (mean(ctrl_s, na.rm=TRUE) + 1e-8)),
                 stringsAsFactors = FALSE)
    })
    phago_wilcox_df <- do.call(rbind, Filter(Negate(is.null), phago_wilcox))

    if (!is.null(phago_wilcox_df) && nrow(phago_wilcox_df) > 0) {
      phago_wilcox_df$sig <- phago_wilcox_df$p_value < 0.05
      phago_wilcox_df$n_phago_genes <- length(phago_avail)

      write.table(phago_wilcox_df,
                  file.path(OUTPUT_DIR, "05_phago_enrichment_wilcoxon.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE)


      sig_phago_cts <- phago_wilcox_df$cell_type[phago_wilcox_df$sig]



      hf_key_markers_2024 <- list(
        Monocyte          = c("S100A8","S100A9","LYZ","CD14","FCN1","VCAN"),
        Myeloid           = c("CD68","C1QA","C1QB","VSIG4","MARCO","MRC1"),
        Endothelial_cells = c("PECAM1","VWF","CDH5","FLT1","KDR","ENG"),
        Endothelial       = c("LYVE1","PROX1","FLT4","PDPN","CCL21"),
        Fibroblasts       = c("DCN","LUM","COL1A1","POSTN","THY1")
      )


      if (inherits(sc_diff[["RNA"]], "Assay5")) {
        sc_diff[["RNA"]] <- JoinLayers(sc_diff[["RNA"]])
      }
      expr_phago <- tryCatch(
        GetAssayData(sc_diff, assay="RNA", layer="data"),
        error=function(e) sc_diff@assays$RNA@data
      )
      avail_phago <- intersect(phago_avail, rownames(expr_phago))

      mean_by_ct <- sapply(unique(sc_diff$annotated_cell_type), function(ct) {
        cells_ct <- colnames(sc_diff)[sc_diff$annotated_cell_type == ct]
        if (length(cells_ct) < 5) return(0)
        avail_phago <- intersect(phago_avail, rownames(expr_phago))
        if (length(avail_phago) == 0) return(0)
        mean(rowMeans(as.matrix(expr_phago[avail_phago, cells_ct, drop=FALSE])))
      })
      names(mean_by_ct) <- unique(sc_diff$annotated_cell_type)

      cat("[INFO] 各细胞类型吞噬调节基因平均表达量:\n")
      print(round(sort(mean_by_ct, decreasing=TRUE), 4))

      key_cell_type <- names(which.max(mean_by_ct))
      cat(sprintf("\n[RESULT] ===== 核心关键细胞类型（基于实际表达量）: %s =====\n\n",
                  key_cell_type))
      cat(sprintf("[STAT] Key_cell_type: %s\n", key_cell_type))
    }




    cat("[INFO] 绘制生物标志物在差异细胞中的表达气泡图...\n")

    biomarkers <- c("FCN3", "MAP2K1", "FCER1G")
    bio_avail  <- intersect(biomarkers, rownames(sc_diff))

    if (length(bio_avail) > 0) {
      Idents(sc_diff) <- "annotated_cell_type"

      if (inherits(sc_diff[["RNA"]], "Assay5")) {
        sc_diff[["RNA"]] <- JoinLayers(sc_diff[["RNA"]])
      }

      bio_expr   <- FetchData(sc_diff, vars = bio_avail)
      bio_expr$cell_type    <- sc_diff$annotated_cell_type
      bio_expr$disease_group <- sc_diff$disease_group

      bio_long <- tidyr::pivot_longer(bio_expr,
                                      cols      = all_of(bio_avail),
                                      names_to  = "gene",
                                      values_to = "expr")

      bio_summary <- bio_long %>%
        group_by(cell_type, gene, disease_group) %>%
        summarise(avg_expr = mean(expm1(expr)),
                  pct_expr = mean(expr > 0) * 100,
                  .groups = "drop")

      bio_summary$cell_type <- factor(
        bio_summary$cell_type,
        levels = rev(sort(unique(bio_summary$cell_type)))
      )
      bio_summary$disease_group <- factor(bio_summary$disease_group,
                                          levels = c("Donor", "DCM"))

      p_bio_bubble <- ggplot(bio_summary,
                             aes(x = gene, y = cell_type,
                                 size = pct_expr, color = avg_expr)) +
        geom_point() +
        facet_grid(. ~ disease_group) +
        scale_size_continuous(range = c(0.5, 10),
                              breaks = c(20, 40, 60, 80),
                              name = "% Expressed") +
        scale_color_gradient(low = "white", high = "red",
                             name = "Avg Expression") +
        labs(x     = "Gene",
             y     = "Cell Type",
             title = "Biomarker Expression in Differential Cell Types") +
        theme_bw() +
        theme(axis.text.x  = element_text(size = 10, angle = 0),
              axis.text.y  = element_text(size = 9),
              axis.title   = element_text(size = 11),
              plot.title   = element_text(hjust = 0.5, size = 13),
              strip.text   = element_text(size = 11),
              legend.text  = element_text(size = 8))

      save_plot(p_bio_bubble, "05_Figure21B_biomarker_bubble",
                width = 10, height = 8)
    } else {
      cat("[WARN] 生物标志物基因均不在数据集中，跳过气泡图。\n")
    }




    cat("[INFO] 绘制 DCM vs Donor 差异细胞比较图...\n")

    ct_compare_df <- sc_diff@meta.data %>%
      group_by(disease_group, annotated_cell_type) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(disease_group) %>%
      mutate(prop = n / sum(n) * 100) %>%
      ungroup()


    if (exists("wilcox_df") && nrow(wilcox_df) > 0) {
      ct_compare_df <- ct_compare_df %>%
        dplyr::left_join(
          wilcox_df %>%
            dplyr::select(cell_type, significant) %>%
            dplyr::rename(annotated_cell_type = cell_type),
          by = "annotated_cell_type"
        )
      ct_compare_df$sig_label <- ifelse(
        !is.na(ct_compare_df$significant) & ct_compare_df$significant,
        "*", ""
      )
    } else {
      ct_compare_df$sig_label <- ""
    }

    p_compare <- ggplot(ct_compare_df,
                        aes(x = annotated_cell_type, y = prop,
                            fill = disease_group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8),
               width = 0.7) +
      scale_fill_manual(values = c("Donor" = "#6AA8E8",
                                   "DCM" = "#E8534A"),
                        name = "Group") +
      geom_text(aes(label = sig_label,
                    y     = prop + 0.5),
                position = position_dodge(width = 0.8),
                size = 5, color = "black", vjust = 0) +
      labs(x     = "Cell Type",
           y     = "Proportion (%)",
           title = "Cell Type Proportions: DCM vs Donor\n(* = p < 0.05, Wilcoxon)") +
      theme_bw() +
      theme(axis.text.x  = element_text(size = 9, angle = 45, hjust = 1),
            axis.title   = element_text(size = 11),
            plot.title   = element_text(hjust = 0.5, size = 12),
            legend.text  = element_text(size = 10))

    save_plot(p_compare, "05_Figure21_HF_vs_Donor_proportion",
              width = 12, height = 7)
  }
}




cat("[INFO] 绘制各 Cluster QC 小提琴图...\n")
qc_feats2 <- intersect(c("nFeature_RNA","nCount_RNA","percent.mito"),
                       colnames(sc@meta.data))
Idents(sc) <- cluster_col

if (length(qc_feats2) > 0) {
  qc_df_cl <- sc@meta.data[, c("seurat_clusters", qc_feats2), drop = FALSE]
  qc_df_cl$cell <- rownames(qc_df_cl)
  qc_long_cl <- tidyr::pivot_longer(qc_df_cl,
                                    cols      = all_of(qc_feats2),
                                    names_to  = "feature",
                                    values_to = "value")
  qc_long_cl$cluster <- factor(qc_long_cl$seurat_clusters)

  p_qc_cl <- ggplot(qc_long_cl, aes(x = cluster, y = value,
                                    fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_jitter(width = 0.15, size = 0.1, alpha = 0.3) +
    facet_wrap(~ feature, scales = "free_y",
               ncol = length(qc_feats2)) +
    theme_bw() +
    theme(axis.text.x  = element_text(size = 7, angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          plot.title   = element_text(size = 11, hjust = 0.5),
          legend.position = "none") +
    ggtitle("QC Metrics by Cluster")

  save_plot(p_qc_cl, "05_QC_violinplot",
            width  = 5 * length(qc_feats2),
            height = 6)
}




cat("[INFO] 绘制关键基因 FeaturePlot...\n")
key_genes       <- c("FCN3", "MAP2K1", "FCER1G")
key_genes_exist <- intersect(key_genes, rownames(sc))
missing_genes   <- setdiff(key_genes, rownames(sc))

if (length(missing_genes) > 0)
  cat(sprintf("[WARN] 以下基因不存在，跳过: %s\n",
              paste(missing_genes, collapse = ", ")))

if (length(key_genes_exist) > 0) {
  umap_emb <- sc@reductions[["umap"]]@cell.embeddings
  umap_1   <- umap_emb[, 1]
  umap_2   <- umap_emb[, 2]

  gene_expr <- tryCatch(
    FetchData(sc, vars = key_genes_exist),
    error = function(e) {
      mat <- tryCatch(
        as.matrix(GetAssayData(sc, assay="SCT", layer="data")[key_genes_exist,,drop=FALSE]),
        error = function(e2)
          GetAssayData(sc, assay="SCT", slot="data")[key_genes_exist,,drop=FALSE]
      )
      as.data.frame(t(as.matrix(mat)))
    }
  )

  n_cols <- min(3L, length(key_genes_exist))
  n_rows <- ceiling(length(key_genes_exist) / n_cols)

  plot_list <- lapply(key_genes_exist, function(gene) {
    expr_col  <- if (gene %in% colnames(gene_expr)) gene else tolower(gene)
    expr_vals <- if (expr_col %in% colnames(gene_expr))
      gene_expr[[expr_col]] else rep(0, length(umap_1))
    df_g <- data.frame(umap_1 = umap_1, umap_2 = umap_2,
                       expr   = as.numeric(expr_vals))
    df_g <- df_g[order(df_g$expr), ]
    ggplot(df_g, aes(x = umap_1, y = umap_2, color = expr)) +
      geom_point(size = 0.3, alpha = 0.8) +
      scale_color_gradient(low = "lightgrey", high = "red",
                           name = "expr") +
      ggtitle(gene) +
      theme_bw() +
      theme(plot.title  = element_text(size = 12, hjust = 0.5),
            legend.text = element_text(size = 7),
            axis.title  = element_text(size = 9),
            panel.grid  = element_blank())
  })

  p_feat <- patchwork::wrap_plots(plot_list, ncol = n_cols)
  save_plot(p_feat, "05_featureplot_key_genes",
            width  = 5 * n_cols,
            height = 4.5 * n_rows)
}




cat("[INFO] 计算各 Cluster Top Markers 并绘制 Dotplot...\n")
markers_file <- file.path(OUTPUT_DIR, "05_all_markers.rds")

if (file.exists(markers_file)) {
  cat("[INFO] 检测到已有 markers 文件，直接加载。\n")
  all_markers <- readRDS(markers_file)
} else {
  all_markers <- tryCatch({
    FindAllMarkers(sc, assay = "SCT", only.pos = TRUE,
                   min.pct = 0.25, logfc.threshold = 0.25,
                   verbose = FALSE)
  }, error = function(err) {
    cat(sprintf("[WARN] FindAllMarkers 出错: %s\n", conditionMessage(err)))
    NULL
  })
  if (!is.null(all_markers))
    saveRDS(all_markers, markers_file)
}

if (!is.null(all_markers) && nrow(all_markers) > 0) {
  top_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 3) %>%
    ungroup()

  dot_genes <- unique(top_markers$gene)
  dot_genes  <- intersect(dot_genes, rownames(sc))
  if (length(dot_genes) > 60) dot_genes <- dot_genes[1:60]

  if (length(dot_genes) > 0) {
    tryCatch({
      mat <- FetchData(sc, vars = dot_genes)
      mat$cluster <- as.character(Idents(sc))

      dot_long <- tidyr::pivot_longer(mat,
                                      cols      = all_of(dot_genes),
                                      names_to  = "gene",
                                      values_to = "expr")
      dot_summary <- dot_long %>%
        group_by(cluster, gene) %>%
        summarise(avg_expr = mean(expm1(expr)),
                  pct_expr = mean(expr > 0) * 100,
                  .groups  = "drop")

      dot_summary$gene    <- factor(dot_summary$gene, levels = dot_genes)
      dot_summary$cluster <- factor(
        dot_summary$cluster,
        levels = sort(unique(dot_summary$cluster))
      )

      p4 <- ggplot(dot_summary,
                   aes(x = gene, y = cluster,
                       size = pct_expr, color = avg_expr)) +
        geom_point() +
        scale_size_continuous(range = c(0.5, 6),
                              name = "% Expressed") +
        scale_color_gradient(low = "lightgrey", high = "blue",
                             name = "Avg Expr") +
        ggtitle("Top Cluster Markers — DotPlot") +
        theme_bw() +
        theme(axis.text.x  = element_text(size = 7, angle = 90,
                                          vjust = 0.5, hjust = 1),
              axis.text.y  = element_text(size = 8),
              axis.title   = element_blank(),
              plot.title   = element_text(hjust = 0.5, size = 13))

      save_plot(p4, "05_dotplot_clusters",
                width  = max(14, length(dot_genes) * 0.25),
                height = 8)
    }, error = function(err) {
      cat(sprintf("[WARN] DotPlot 绘制出错: %s\n", conditionMessage(err)))
    })
  }
}




cluster_stats <- table(Idents(sc))
cat("[INFO] 各 Cluster 细胞数:\n")
print(cluster_stats)
write.table(
  data.frame(Cluster   = names(cluster_stats),
             CellCount = as.integer(cluster_stats)),
  file.path(OUTPUT_DIR, "05_cluster_cellcounts.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)




rds_path <- file.path(OUTPUT_DIR, "seurat_clustered.rds")
cat(sprintf("[INFO] 保存 Seurat 对象至: %s\n", rds_path))
saveRDS(sc, rds_path)




cat("\n[INFO] ===== r.05_singlecell_v2.R 执行完成 =====\n")
cat(sprintf("[STAT] QC_before_cells:   %d\n", n_cells_before))
cat(sprintf("[STAT] QC_before_genes:   %d\n", n_genes_before))
cat(sprintf("[STAT] QC_after_cells:    %d\n", n_cells_after))
cat(sprintf("[STAT] QC_after_genes:    %d\n", n_genes_after))
cat(sprintf("[STAT] Selected_PCs:      %d\n", selected_dims))
cat(sprintf("[STAT] SC_total_cells:    %d\n", n_total_cells))
cat(sprintf("[STAT] SC_cluster_count:  %d\n", n_clusters))
cat(sprintf("[STAT] SC_resolution:     %.1f\n", RESOLUTION))
if (exists("key_cell_type"))
  cat(sprintf("[STAT] Key_cell_type:     %s\n", key_cell_type))
if (exists("diff_cells") && length(diff_cells) > 0)
  cat(sprintf("[STAT] Diff_cell_types:   %s\n",
              paste(diff_cells, collapse=", ")))
