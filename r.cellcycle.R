SERVER_DIR   <- ""
OUT_DIR      <- file.path(SERVER_DIR, "17_cellcycle_macrophage")
SC_PATH      <- file.path(SERVER_DIR, "03_single_cell/seurat_clustered.rds")


CELL_TYPE_COL <- "cell_type"
TARGET_CELL   <- "Macrophages"
CONDITION_COL <- "Condition"
DCM_LABEL     <- "DCM"
DONOR_LABEL   <- "Donor"
MARKER_GENES  <- c("FCN3", "FCER1G", "MAP2K1")


CONDITION_COLORS <- c("DCM" = "#E64B35", "Donor" = "#4DBBD5")

PHASE_COLORS     <- c("G1" = "#4DBBD5", "S" = "#F39B7F", "G2M" = "#00A087")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(paste0("[INFO] 输出目录: ", OUT_DIR, "\n"))




required_pkgs <- c("Seurat", "ggplot2", "dplyr", "tidyr", "ggpubr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("[SKIP] 缺少必要包: ", pkg, "，终止\n"))
    quit(status = 0)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
})
cat(paste0("[INFO] Seurat 版本: ", as.character(packageVersion("Seurat")), "\n"))




if (!file.exists(SC_PATH)) {
  cat("[SKIP] 单细胞 RDS 文件不存在:", SC_PATH, "\n")
  quit(status = 0)
}

cat("[INFO] 读取单细胞对象:", SC_PATH, "\n")
sc_full <- readRDS(SC_PATH)
cat(paste0("[INFO] 总细胞数: ", ncol(sc_full),
           " | 基因数: ", nrow(sc_full), "\n"))


for (col in c(CELL_TYPE_COL, CONDITION_COL)) {
  if (!col %in% colnames(sc_full@meta.data)) {
    cat(paste0("[ERROR] meta.data 中不存在列: ", col, "\n"))
    cat("[INFO] 当前可用列:", paste(colnames(sc_full@meta.data), collapse = ", "), "\n")
    quit(status = 1)
  }
}


cond_vals <- unique(sc_full@meta.data[[CONDITION_COL]])
cat("[INFO] Condition 取值:", paste(cond_vals, collapse = ", "), "\n")
donor_candidates <- setdiff(cond_vals, DCM_LABEL)
if (length(donor_candidates) == 1) {
  DONOR_LABEL <- donor_candidates
  cat(paste0("[INFO] 自动检测对照组标签: ", DONOR_LABEL, "\n"))
} else {
  cat(paste0("[INFO] 使用默认对照组标签: ", DONOR_LABEL, "\n"))
}
CONDITION_COLORS <- setNames(
  c("#E64B35", "#4DBBD5"),
  c(DCM_LABEL, DONOR_LABEL)
)




cat(paste0("[INFO] 提取 ", TARGET_CELL, " 亚群...\n"))


cell_mask <- grepl(TARGET_CELL, sc_full@meta.data[[CELL_TYPE_COL]],
                   ignore.case = TRUE)
if (sum(cell_mask) == 0) {
  cat(paste0("[ERROR] 未找到包含 '", TARGET_CELL, "' 的细胞\n"))
  cat("[INFO] 现有细胞类型:\n")
  print(table(sc_full@meta.data[[CELL_TYPE_COL]]))
  quit(status = 1)
}

sc <- sc_full[, cell_mask]
cat(paste0("[INFO] ", TARGET_CELL, " 细胞数: ", ncol(sc), "\n"))
cat("[INFO] 各 Condition 细胞数:\n")
print(table(sc@meta.data[[CONDITION_COL]]))




tryCatch({
  DefaultAssay(sc) <- "RNA"
}, error = function(e) {
  cat(paste0("[WARN] 设置 DefaultAssay 失败: ", e$message, "\n"))
})

tryCatch({
  sc <- JoinLayers(sc)
  cat("[INFO] JoinLayers 完成（Seurat v5）\n")
}, error = function(e) {
  cat("[INFO] JoinLayers 跳过（Seurat v4）\n")
})


has_data <- tryCatch({
  d <- GetAssayData(sc, layer = "data")
  nrow(d) > 0 && sum(d) > 0
}, error = function(e) FALSE)

if (!has_data) {
  cat("[INFO] 执行 NormalizeData...\n")
  tryCatch({
    sc <- NormalizeData(sc, verbose = FALSE)
  }, error = function(e) {
    cat(paste0("[WARN] NormalizeData 失败: ", e$message, "\n"))
  })
}




cat("[INFO] 获取细胞周期基因集...\n")
s_genes <- g2m_genes <- NULL

tryCatch({
  s_genes   <- cc.genes.updated.2019$s.genes
  g2m_genes <- cc.genes.updated.2019$g2m.genes
  cat(paste0("[INFO] 使用 cc.genes.updated.2019: S=",
             length(s_genes), " G2M=", length(g2m_genes), "\n"))
}, error = function(e) {
  tryCatch({
    s_genes   <<- cc.genes$s.genes
    g2m_genes <<- cc.genes$g2m.genes
    cat("[INFO] 使用 cc.genes\n")
  }, error = function(e2) {
    cat("[ERROR] 无法获取细胞周期基因集\n")
    quit(status = 1)
  })
})

s_genes_use   <- intersect(s_genes,   rownames(sc))
g2m_genes_use <- intersect(g2m_genes, rownames(sc))
cat(paste0("[INFO] 实际使用基因: S=", length(s_genes_use),
           " G2M=", length(g2m_genes_use), "\n"))


cc_success <- FALSE
for (assay_try in c("SCT", "RNA")) {
  if (!assay_try %in% names(sc@assays)) next
  tryCatch({
    DefaultAssay(sc) <- assay_try
    sc <- CellCycleScoring(sc,
                           s.features   = s_genes_use,
                           g2m.features = g2m_genes_use,
                           set.ident    = FALSE)
    cat(paste0("[INFO] CellCycleScoring 完成（assay: ", assay_try, "）\n"))
    cc_success <- TRUE
    break
  }, error = function(e) {
    cat(paste0("[WARN] CellCycleScoring(", assay_try, ") 失败: ", e$message, "\n"))
  })
}


if (!cc_success) {
  cat("[INFO] 使用手动均值评分...\n")
  data_mat <- tryCatch(
    GetAssayData(sc, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(sc, assay = "RNA", slot = "data")
  )
  calc_score <- function(genes, mat) {
    g_use <- intersect(genes, rownames(mat))
    if (length(g_use) == 0) return(rep(0, ncol(mat)))
    colMeans(as.matrix(mat[g_use, , drop = FALSE]), na.rm = TRUE)
  }
  sc@meta.data$S.Score   <- calc_score(s_genes_use, data_mat)
  sc@meta.data$G2M.Score <- calc_score(g2m_genes_use, data_mat)
  sc@meta.data$Phase <- with(sc@meta.data,
                             ifelse(S.Score > G2M.Score & S.Score > 0.01, "S",
                                    ifelse(G2M.Score > S.Score & G2M.Score > 0.01, "G2M", "G1")))
  cc_success <- TRUE
}

sc@meta.data$Phase <- factor(sc@meta.data$Phase, levels = c("G1", "S", "G2M"))

cat("[INFO] 各周期细胞数:\n")
print(table(sc@meta.data$Phase))
cat("[INFO] DCM vs Donor × Phase:\n")
print(table(sc@meta.data[[CONDITION_COL]], sc@meta.data$Phase))




meta_save <- sc@meta.data[, c("S.Score", "G2M.Score", "Phase",
                              CONDITION_COL, CELL_TYPE_COL)]
meta_save$cell_id <- rownames(meta_save)
write.csv(meta_save,
          file.path(OUT_DIR, "cellcycle_macrophage_scores.csv"),
          row.names = FALSE)
cat("[INFO] 评分 CSV 已保存\n")





markers_use <- intersect(MARKER_GENES, rownames(sc))
if (length(markers_use) == 0) {
  cat("[WARN] 标志物基因均不在对象中，图3/图4将跳过\n")
  cat("[INFO] 当前基因名示例（前20）:", paste(head(rownames(sc), 20), collapse = ", "), "\n")
} else {
  cat(paste0("[INFO] 找到标志物基因: ", paste(markers_use, collapse = ", "), "\n"))
  missing_markers <- setdiff(MARKER_GENES, rownames(sc))
  if (length(missing_markers) > 0)
    cat(paste0("[WARN] 以下基因未找到: ", paste(missing_markers, collapse = ", "), "\n"))
}

get_expr_df <- function(sc_obj, genes, condition_col, phase_col = "Phase") {
  if (length(genes) == 0) return(NULL)
  expr_mat <- tryCatch(
    GetAssayData(sc_obj, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(sc_obj, assay = "RNA", slot = "data")
  )
  expr_sub <- as.data.frame(t(as.matrix(expr_mat[genes, , drop = FALSE])))
  expr_sub$Phase     <- sc_obj@meta.data[[phase_col]]
  expr_sub$Condition <- sc_obj@meta.data[[condition_col]]
  expr_sub$cell_id   <- rownames(expr_sub)
  return(expr_sub)
}

expr_df <- get_expr_df(sc, markers_use, CONDITION_COL)




umap_png <- file.path(OUT_DIR, "cellcycle_umap_macrophage.png")

umap_df <- NULL
for (key in c("umap", "UMAP")) {
  if (key %in% names(sc@reductions)) {
    tryCatch({
      emb <- sc@reductions[[key]]@cell.embeddings
      umap_df <- as.data.frame(emb)
      colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
      break
    }, error = function(e) {})
  }
}


if (is.null(umap_df)) {
  cat("[INFO] Macrophages 子集无 UMAP，重新计算 PCA + UMAP...\n")
  tryCatch({
    DefaultAssay(sc) <- "RNA"
    sc <- FindVariableFeatures(sc, nfeatures = 2000, verbose = FALSE)
    sc <- ScaleData(sc, verbose = FALSE)
    sc <- RunPCA(sc, npcs = 20, verbose = FALSE)
    sc <- RunUMAP(sc, dims = 1:20, verbose = FALSE, seed.use = 42)
    emb <- sc@reductions[["umap"]]@cell.embeddings
    umap_df <- as.data.frame(emb)
    colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
    cat("[INFO] UMAP 计算完成\n")
  }, error = function(e) {
    cat(paste0("[WARN] UMAP 计算失败: ", e$message, "\n"))
  })
}

if (!is.null(umap_df)) {
  umap_df$Phase     <- sc@meta.data$Phase
  umap_df$Condition <- sc@meta.data[[CONDITION_COL]]
  umap_df$Phase     <- factor(umap_df$Phase, levels = c("G1", "S", "G2M"))

  set.seed(42)
  umap_df <- umap_df[sample(nrow(umap_df)), ]
  pt_size <- if (nrow(umap_df) > 5000) 0.3 else if (nrow(umap_df) > 1000) 0.5 else 1.0

  n_G1  <- sum(umap_df$Phase == "G1",  na.rm = TRUE)
  n_S   <- sum(umap_df$Phase == "S",   na.rm = TRUE)
  n_G2M <- sum(umap_df$Phase == "G2M", na.rm = TRUE)


  p1a <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Phase)) +
    geom_point(size = pt_size, alpha = 0.75) +
    scale_color_manual(values = PHASE_COLORS, name = "Phase") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = "Macrophages — Cell Cycle Phase",
         subtitle = paste0("G1: ", n_G1, "  S: ", n_S, "  G2M: ", n_G2M),
         x = "UMAP 1", y = "UMAP 2") +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          panel.grid    = element_blank(),
          legend.title  = element_text(face = "bold"))

  p1b <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Condition)) +
    geom_point(size = pt_size, alpha = 0.75) +
    scale_color_manual(values = CONDITION_COLORS, name = "Condition") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = "Macrophages — Condition",
         x = "UMAP 1", y = "UMAP 2") +
    theme_bw(base_size = 11) +
    theme(plot.title  = element_text(face = "bold", size = 12),
          panel.grid  = element_blank(),
          legend.title = element_text(face = "bold"))

  p_umap_combined <- ggarrange(p1a, p1b, ncol = 2, nrow = 1)

  png(umap_png, width = 2000, height = 900, res = 150)
  print(p_umap_combined)
  dev.off()
  cat(paste0("[INFO] 图1 UMAP PNG 已保存: ", umap_png, "\n"))

  umap_pdf <- sub("\\.png$", ".pdf", umap_png)
  cairo_pdf(umap_pdf, width = 2000/150, height = 900/150)
  print(p_umap_combined)
  dev.off()
  cat(paste0("[INFO] 图1 UMAP PDF 已保存: ", umap_pdf, "\n"))
} else {
  cat("[WARN] 无 UMAP 坐标，跳过图1\n")
}





prop_png <- file.path(OUT_DIR, "cellcycle_proportion_DCMvsDonor.png")

tryCatch({
  prop_df <- sc@meta.data %>%
    dplyr::select(Phase, Condition = !!sym(CONDITION_COL)) %>%
    dplyr::filter(Condition %in% c(DCM_LABEL, DONOR_LABEL)) %>%
    dplyr::mutate(Phase = factor(Phase, levels = c("G1", "S", "G2M")),
                  Condition = factor(Condition, levels = c(DCM_LABEL, DONOR_LABEL)))

  prop_summary <- prop_df %>%
    dplyr::group_by(Condition, Phase) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::group_by(Condition) %>%
    dplyr::mutate(
      total = sum(n),
      prop  = n / total * 100
    ) %>%
    dplyr::ungroup()


  fisher_results <- lapply(levels(prop_df$Phase), function(ph) {
    tbl <- table(prop_df$Condition, prop_df$Phase == ph)
    if (all(dim(tbl) == c(2, 2))) {
      p <- fisher.test(tbl)$p.value
    } else {
      p <- NA
    }
    data.frame(Phase = ph, p_value = p)
  })
  fisher_df <- do.call(rbind, fisher_results)
  fisher_df$label <- dplyr::case_when(
    fisher_df$p_value < 0.001 ~ "***",
    fisher_df$p_value < 0.01  ~ "**",
    fisher_df$p_value < 0.05  ~ "*",
    TRUE                       ~ "ns"
  )
  cat("[INFO] Fisher 检验结果:\n")
  print(fisher_df)

  p2 <- ggplot(prop_summary,
               aes(x = Condition, y = prop, fill = Phase)) +
    geom_bar(stat = "identity", position = "stack",
             width = 0.55, color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(prop >= 3,
                                 paste0(round(prop, 1), "%\n(n=", n, ")"),
                                 "")),
              position = position_stack(vjust = 0.5),
              size = 3.2, color = "white", fontface = "bold") +
    scale_fill_manual(values = PHASE_COLORS, name = "Phase") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 102)) +
    labs(title    = paste0("Macrophages Cell Cycle Distribution\n",
                           DCM_LABEL, " vs ", DONOR_LABEL),
         subtitle = "Percentage of cells in each phase (Fisher's exact test)",
         x = NULL, y = "Proportion (%)") +
    theme_bw(base_size = 12) +
    theme(plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
          legend.title  = element_text(face = "bold"),
          axis.text.x   = element_text(size = 12, face = "bold"),
          panel.grid    = element_blank())


  sig_annotations <- prop_summary %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(y_pos = 101, .groups = "drop")


  caption_text <- paste(
    sapply(1:nrow(fisher_df), function(i) {
      sprintf("Phase %s: %s (p=%.3f)",
              fisher_df$Phase[i], fisher_df$label[i],
              ifelse(is.na(fisher_df$p_value[i]), NA, fisher_df$p_value[i]))
    }),
    collapse = "\n"
  )

  p2 <- p2 + labs(caption = caption_text) +
    theme(plot.caption = element_text(size = 8, color = "grey40",
                                      hjust = 0, face = "italic"))

  png(prop_png, width = 900, height = 1000, res = 150)
  print(p2)
  dev.off()
  cat(paste0("[INFO] 图2 比例图 PNG 已保存: ", prop_png, "\n"))

  prop_pdf <- sub("\\.png$", ".pdf", prop_png)
  cairo_pdf(prop_pdf, width = 900/150, height = 1000/150)
  print(p2)
  dev.off()
  cat(paste0("[INFO] 图2 比例图 PDF 已保存: ", prop_pdf, "\n"))
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat(paste0("[WARN] 图2 绘制失败: ", e$message, "\n"))
})






violin_png <- file.path(OUT_DIR, "cellcycle_violin_markers.png")

if (!is.null(expr_df) && length(markers_use) > 0) {
  tryCatch({
    vln_long <- tidyr::pivot_longer(
      expr_df,
      cols      = all_of(markers_use),
      names_to  = "Gene",
      values_to = "Expression"
    )
    vln_long$Phase     <- factor(vln_long$Phase, levels = c("G1", "S", "G2M"))
    vln_long$Condition <- factor(vln_long$Condition,
                                 levels = c(DCM_LABEL, DONOR_LABEL))
    vln_long$Gene      <- factor(vln_long$Gene, levels = markers_use)

    p3 <- ggplot(vln_long,
                 aes(x = Phase, y = Expression, fill = Phase)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
      geom_boxplot(width = 0.12, outlier.size = 0.2,
                   fill = "white", alpha = 0.7) +
      facet_grid(Gene ~ Condition, scales = "free_y") +
      scale_fill_manual(values = PHASE_COLORS) +
      stat_compare_means(
        comparisons = list(c("G1", "S"), c("G1", "G2M"), c("S", "G2M")),
        method      = "wilcox.test",
        label       = "p.signif",
        size        = 3,
        tip.length  = 0.01
      ) +
      labs(
        title    = "Marker Gene Expression across Cell Cycle Phases",
        subtitle = paste0("Macrophages — ", DCM_LABEL, " vs ", DONOR_LABEL),
        x        = "Cell Cycle Phase",
        y        = "Normalized Expression"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        strip.text    = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "#F0F4FF"),
        legend.position  = "none",
        panel.grid.major.x = element_blank()
      )


    fig_h <- max(800, length(markers_use) * 320)
    png(violin_png, width = 1400, height = fig_h, res = 150)
    print(p3)
    dev.off()
    cat(paste0("[INFO] 图3 小提琴图 PNG 已保存: ", violin_png, "\n"))

    violin_pdf <- sub("\\.png$", ".pdf", violin_png)
    cairo_pdf(violin_pdf, width = 1400/150, height = fig_h/150)
    print(p3)
    dev.off()
    cat(paste0("[INFO] 图3 小提琴图 PDF 已保存: ", violin_pdf, "\n"))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(paste0("[WARN] 图3 绘制失败: ", e$message, "\n"))
  })
} else {
  cat("[WARN] 标志物基因不可用，跳过图3\n")
}







bubble_png <- file.path(OUT_DIR, "cellcycle_bubble_markers.png")

if (!is.null(expr_df) && length(markers_use) > 0) {
  tryCatch({
    bub_df <- tidyr::pivot_longer(
      expr_df,
      cols      = all_of(markers_use),
      names_to  = "Gene",
      values_to = "Expression"
    )
    bub_df$Phase     <- factor(bub_df$Phase, levels = c("G1", "S", "G2M"))
    bub_df$Condition <- factor(bub_df$Condition,
                               levels = c(DONOR_LABEL, DCM_LABEL))
    bub_df$Gene      <- factor(bub_df$Gene, levels = rev(markers_use))

    bub_summary <- bub_df %>%
      dplyr::group_by(Gene, Phase, Condition) %>%
      dplyr::summarise(
        mean_expr    = mean(Expression, na.rm = TRUE),
        pct_expr     = mean(Expression > 0, na.rm = TRUE) * 100,
        .groups      = "drop"
      )


    bub_summary <- bub_summary %>%
      dplyr::group_by(Gene) %>%
      dplyr::mutate(
        mean_zscore = scale(mean_expr)[, 1]
      ) %>%
      dplyr::ungroup()

    p4 <- ggplot(bub_summary,
                 aes(x = Phase, y = Gene,
                     size = pct_expr, color = mean_zscore)) +
      geom_point(alpha = 0.9) +
      facet_wrap(~ Condition, ncol = 2) +
      scale_size_continuous(
        range  = c(2, 12),
        name   = "% Expressed",
        breaks = c(10, 30, 50, 70, 90)
      ) +
      scale_color_gradient2(
        low      = "#3B4CC0",
        mid      = "#DDDDDD",
        high     = "#B40426",
        midpoint = 0,
        name     = "Mean Expr\n(Z-score)"
      ) +
      labs(
        title    = "Marker Gene Activity across Cell Cycle Phases",
        subtitle = paste0("Macrophages — ", DONOR_LABEL, " vs ", DCM_LABEL,
                          "\nBubble size = % expressed cells; Color = mean expression (z-score)"),
        x = "Cell Cycle Phase",
        y = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        strip.text       = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "#F0F4FF"),
        axis.text.y      = element_text(size = 11, face = "bold.italic"),
        axis.text.x      = element_text(size = 11),
        panel.grid.major = element_line(color = "grey90"),
        legend.position  = "right"
      )

    png(bubble_png, width = 1400, height = 700, res = 150)
    print(p4)
    dev.off()
    cat(paste0("[INFO] 图4 气泡图 PNG 已保存: ", bubble_png, "\n"))

    bubble_pdf <- sub("\\.png$", ".pdf", bubble_png)
    cairo_pdf(bubble_pdf, width = 1400/150, height = 700/150)
    print(p4)
    dev.off()
    cat(paste0("[INFO] 图4 气泡图 PDF 已保存: ", bubble_pdf, "\n"))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(paste0("[WARN] 图4 绘制失败: ", e$message, "\n"))
  })
} else {
  cat("[WARN] 标志物基因不可用，跳过图4\n")
}




cat("\n========================================\n")
cat("[STAT] Target_cell_type:", TARGET_CELL, "\n")
cat("[STAT] Total_Macrophages:", ncol(sc), "\n")
dcm_n   <- sum(sc@meta.data[[CONDITION_COL]] == DCM_LABEL,   na.rm = TRUE)
donor_n <- sum(sc@meta.data[[CONDITION_COL]] == DONOR_LABEL, na.rm = TRUE)
cat("[STAT] DCM_Macrophages:", dcm_n, "\n")
cat("[STAT] Donor_Macrophages:", donor_n, "\n")
phase_tbl <- table(sc@meta.data$Phase)
cat("[STAT] G1:", as.integer(phase_tbl["G1"]),
    "S:", as.integer(phase_tbl["S"]),
    "G2M:", as.integer(phase_tbl["G2M"]), "\n")
cat("[INFO] r.22_cellcycle_macrophage.R 运行完成\n")
cat("========================================\n")
