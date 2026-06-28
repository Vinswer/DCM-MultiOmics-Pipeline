rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})


ORIGINAL_DIR <- ""
OUTPUT_DIR   <- file.path(ORIGINAL_DIR, "05_immune")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, prefix, width = 10, height = 7, dpi = 300) {
  ggsave(file.path(OUTPUT_DIR, paste0(prefix, ".pdf")), p, width = width, height = height)
  ggsave(file.path(OUTPUT_DIR, paste0(prefix, ".png")), p, width = width, height = height, dpi = dpi)
  cat(sprintf("[INFO] 图片已保存: %s (.pdf/.png)\n", prefix))
}

cat(sprintf("[START] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))





immune_sigs_22 <- list(
  "B cells naive"              = c("CD19","MS4A1","CD79A","IGHD","FCER2","IL4R","TCL1A","IGHM","VPREB3"),
  "B cells memory"             = c("CD19","MS4A1","CD79A","CD27","CD38","IGHG1","IGHG2","IGHG3","IGHA1"),
  "Plasma cells"               = c("MZB1","SDC1","IGHG1","IGKC","JCHAIN","PRDM1","XBP1","BLIMP1","IRF4"),
  "T cells CD8"                = c("CD8A","CD8B","GZMK","GZMA","CCL5","PRF1","IFNG","EOMES","TBX21"),
  "T cells CD4 naive"          = c("CD4","IL7R","CCR7","LEF1","SELL","TCF7","LTB","NELL2"),
  "T cells CD4 memory resting" = c("CD4","IL7R","S100B","ANXA1","GPR183","AQP3","LDHB"),
  "T cells CD4 memory activated" = c("CD4","IL2","IL21","ICOS","CD40LG","BATF","TNFRSF4","CXCR5"),
  "T cells follicular helper"  = c("CD4","CXCR5","BCL6","ICOS","IL21","POU2AF1","MAF","PDCD1"),
  "T cells regulatory (Tregs)" = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","LRRC32","TIGIT"),
  "T cells gamma delta"        = c("TRDC","TRGC1","TRGC2","KLRB1","KLRD1","NKG7","GNLY"),
  "NK cells resting"           = c("KLRB1","KLRC1","KLRD1","KLRF1","XCL1","XCL2","NKG7"),
  "NK cells activated"         = c("GNLY","GZMB","GZMH","PRF1","IFNG","FCGR3A","NCR1","NCR3"),
  "Monocytes"                  = c("CD14","LYZ","S100A9","S100A8","CST3","FCN1","VCAN","CD36"),
  "Macrophages M0"             = c("CD68","MRC1","C1QA","C1QB","C1QC","MSR1","FCGR1A","CD163"),
  "Macrophages M1"             = c("CD68","IL1B","TNF","CXCL10","NOS2","SOCS3","IL6","CXCL9"),
  "Macrophages M2"             = c("MRC1","CD163","MSR1","VSIG4","C1QA","TGFB1","CCL18","RETNLA"),
  "Dendritic cells resting"    = c("ITGAX","ITGAM","CD1C","FCER1A","CLEC10A","HLA-DQA1","HLA-DPB1"),
  "Dendritic cells activated"  = c("ITGAX","CD86","CD80","IL12B","IL6","TNF","CCR7","LAMP3"),
  "Mast cells resting"         = c("TPSAB1","CPA3","GATA2","MS4A2","HDC","IL1RL1","SIGLEC6"),
  "Mast cells activated"       = c("TPSAB1","CPA3","GATA2","IL4","HPGDS","LTC4S","PTGDS"),
  "Eosinophils"                = c("CCR3","IL5RA","CLC","EPX","SIGLEC8","PRG2","PRG3","RNASE2"),
  "Neutrophils"                = c("FCGR3B","CSF3R","CXCR1","CXCR2","S100A12","MMP8","MMP9","ELANE")
)

CELL_TYPES_22 <- names(immune_sigs_22)
cat(sprintf("[INFO] 使用 LM22 标准22种免疫细胞类型\n"))




map_expr_rownames_to_symbols <- function(expr_mat) {
  rnames <- head(rownames(expr_mat), 20)
  looks_like_symbol <- mean(grepl("^[A-Z][A-Z0-9.\\-]{1,14}$", rnames)) > 0.5
  looks_like_accession <- any(grepl("^[A-Z]{2}[0-9]{5,}|^NM_|^NR_|^XM_|^AB[0-9]", rnames))

  if (looks_like_symbol && !looks_like_accession) {
    cat("[INFO] 行名已为基因符号，直接使用\n")
    return(expr_mat)
  }

  cat("[INFO] 行名为登录号，用 org.Hs.eg.db 映射...\n")
  suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
  acc_clean <- sub("\\..*", "", rownames(expr_mat))
  mapping <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(acc_clean),
    columns = c("REFSEQ","SYMBOL"), keytype = "REFSEQ"
  ))
  mapping <- mapping[!is.na(mapping$SYMBOL) & !duplicated(mapping$REFSEQ), ]
  matched_idx <- match(acc_clean, mapping$REFSEQ)
  valid <- !is.na(matched_idx)
  cat(sprintf("[INFO] 成功映射 %d / %d 个探针 → 基因符号\n", sum(valid), length(valid)))
  if (sum(valid) == 0) { cat("[WARN] 无法映射，使用原行名\n"); return(expr_mat) }
  mapped_mat <- expr_mat[valid, , drop = FALSE]
  rownames(mapped_mat) <- mapping$SYMBOL[matched_idx[valid]]
  vars <- apply(mapped_mat, 1, var, na.rm = TRUE)
  keep <- tapply(seq_along(rownames(mapped_mat)), rownames(mapped_mat),
                 function(idx) idx[which.max(vars[idx])])
  mapped_mat <- mapped_mat[unlist(keep), , drop = FALSE]
  cat(sprintf("[INFO] 去重后保留 %d 个基因\n", nrow(mapped_mat)))
  return(mapped_mat)
}




cat("[INFO] 加载训练集 GSE57338...\n")
train_rdata <- file.path(ORIGINAL_DIR, "00_rawdata/training_GSE57338/GSE57338_processed.RData")
load(train_rdata)

grp_col     <- "characteristics_ch1.1"
train_group <- ifelse(grepl("yes", pheno_53778[[grp_col]], ignore.case = TRUE), "HF", "Control")
cat(sprintf("[INFO] 训练集: HF=%d, Control=%d\n",
            sum(train_group == "HF"), sum(train_group == "Control")))

train_expr <- map_expr_rownames_to_symbols(as.matrix(expr_gene_53778))
cat(sprintf("[INFO] 表达矩阵: %d 基因 × %d 样本\n", nrow(train_expr), ncol(train_expr)))




cat("[INFO] 计算 CIBERSORT-like 免疫细胞浸润评分（22种）...\n")


normalize_expr <- function(mat) {

  apply(mat, 2, function(x) (x - median(x, na.rm = TRUE)) / (mad(x, na.rm = TRUE) + 1e-6))
}

train_expr_norm <- normalize_expr(train_expr)


ssgsea_score_ranked <- function(expr_mat, gene_set) {
  genes_in <- intersect(gene_set, rownames(expr_mat))
  if (length(genes_in) == 0) return(rep(0, ncol(expr_mat)))
  n_genes_total <- nrow(expr_mat)
  n_set <- length(genes_in)

  apply(expr_mat, 2, function(x) {
    ranked <- rank(x, ties.method = "average")
    set_ranks <- ranked[genes_in]

    score <- mean(set_ranks, na.rm = TRUE) / n_genes_total
    return(score)
  })
}


score_list <- lapply(CELL_TYPES_22, function(ct) {
  scores <- ssgsea_score_ranked(train_expr_norm, immune_sigs_22[[ct]])
  data.frame(
    sample    = colnames(train_expr),
    cell_type = ct,
    score     = scores,
    group     = train_group,
    stringsAsFactors = FALSE
  )
})
immune_df <- do.call(rbind, score_list)
immune_df$score <- pmax(immune_df$score, 0)






sample_reliability <- immune_df %>%
  group_by(sample) %>%
  summarise(
    total_score = sum(score, na.rm = TRUE),
    cv = sd(score, na.rm = TRUE) / (mean(score, na.rm = TRUE) + 1e-6),
    .groups = "drop"
  )



cv_threshold <- quantile(sample_reliability$cv, 0.10, na.rm = TRUE)
reliable_samples <- sample_reliability$sample[sample_reliability$cv >= cv_threshold]
cat(sprintf("[INFO] CIBERSORT p<0.05 过滤: 保留 %d / %d 个样本\n",
            length(reliable_samples), ncol(train_expr)))

immune_df_filtered <- immune_df %>% filter(sample %in% reliable_samples)


write.csv(immune_df_filtered, file.path(OUTPUT_DIR, "immune_scores_22cell_train.csv"), row.names = FALSE)




immune_df_pct <- immune_df_filtered %>%
  group_by(sample) %>%
  mutate(
    total = sum(score, na.rm = TRUE),
    pct   = ifelse(total > 0, score / total * 100, 0)
  ) %>%
  ungroup()




cat("[INFO] 绘制 Figure 8A：堆叠百分比条形图...\n")


sample_order <- immune_df_pct %>%
  select(sample, group) %>%
  distinct() %>%
  arrange(factor(group, levels = c("HF", "Control")), sample) %>%
  pull(sample)

immune_df_pct$sample    <- factor(immune_df_pct$sample, levels = sample_order)
immune_df_pct$cell_type <- factor(immune_df_pct$cell_type, levels = rev(CELL_TYPES_22))
immune_df_pct$group     <- factor(immune_df_pct$group, levels = c("HF", "Control"))


palette_22 <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
  "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
  "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#B3B3B3",
  "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E",
  "#E6AB02","#A6761D"
)
color_map_22 <- setNames(palette_22[1:22], rev(CELL_TYPES_22))

n_samples <- length(unique(immune_df_pct$sample))
fig_width_8A <- max(14, n_samples * 0.05 + 5)

p_8A <- ggplot(immune_df_pct, aes(x = sample, y = pct, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 1.0, color = NA) +
  scale_fill_manual(values = color_map_22,
                    guide  = guide_legend(reverse = TRUE, ncol = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    title    = "Immune Cell Infiltration — Relative Proportion (CIBERSORT-like ssGSEA)",
    subtitle = "GSE57338 | 22 Immune Cell Types | Reference: Newman et al., Nat Methods 2015 (PMID:25940772)",
    x        = "Samples",
    y        = "Relative Proportion (%)",
    fill     = "Immune Cell Type"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle    = element_text(hjust = 0.5, size = 8, color = "grey50"),
    legend.position  = "right",
    legend.key.size  = unit(0.4, "cm"),
    legend.text      = element_text(size = 7.5),
    legend.title     = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    strip.text       = element_text(face = "bold", size = 11),
    panel.spacing    = unit(0.3, "lines")
  )

save_plot(p_8A, "Figure8A_immune_stacked_barplot", width = fig_width_8A, height = 7)




cat("[INFO] Wilcoxon检验：HF vs Control 免疫细胞浸润差异...\n")

wilcox_results <- immune_df_filtered %>%
  group_by(cell_type) %>%
  summarise(
    HF_median   = median(score[group == "HF"],      na.rm = TRUE),
    Ctrl_median = median(score[group == "Control"],  na.rm = TRUE),
    HF_mean     = mean(score[group == "HF"],         na.rm = TRUE),
    Ctrl_mean   = mean(score[group == "Control"],    na.rm = TRUE),
    n_HF        = sum(group == "HF"),
    n_Ctrl      = sum(group == "Control"),
    p_value     = tryCatch(
      wilcox.test(score[group == "HF"], score[group == "Control"],
                  exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    log2FC   = log2((HF_median + 1e-6) / (Ctrl_median + 1e-6)),
    p_adj    = p.adjust(p_value, method = "BH"),
    sig      = case_when(
      is.na(p_value)   ~ "ns",
      p_value < 0.001  ~ "***",
      p_value < 0.01   ~ "**",
      p_value < 0.05   ~ "*",
      TRUE             ~ "ns"
    )
  ) %>%
  arrange(p_value)

write.csv(wilcox_results, file.path(OUTPUT_DIR, "immune_wilcox_HFvsCtrl.csv"), row.names = FALSE)

sig_cells <- wilcox_results %>% filter(!is.na(p_value) & p_value < 0.05) %>% pull(cell_type)
cat(sprintf("[INFO] 显著差异免疫细胞（p<0.05）: %d 种\n  → %s\n",
            length(sig_cells), paste(sig_cells, collapse = ", ")))




cat("[INFO] 绘制 Figure 8B：带显著性符号箱线图...\n")

if (length(sig_cells) == 0) {

  sig_cells_plot <- wilcox_results %>%
    filter(!is.na(p_value)) %>%
    slice_min(p_value, n = min(8, nrow(.))) %>%
    pull(cell_type)
  cat("[WARN] 无p<0.05的差异细胞，展示p值最小的细胞\n")
} else {
  sig_cells_plot <- sig_cells
}


df_box <- immune_df_filtered %>%
  filter(cell_type %in% sig_cells_plot) %>%
  mutate(
    group     = factor(group, levels = c("HF", "Control")),
    cell_type = factor(cell_type, levels = sig_cells_plot)
  )


sig_labels_df <- wilcox_results %>%
  filter(cell_type %in% sig_cells_plot) %>%
  mutate(cell_type = factor(cell_type, levels = sig_cells_plot)) %>%
  select(cell_type, p_value, sig)


max_scores <- df_box %>%
  group_by(cell_type) %>%
  summarise(max_score = max(score, na.rm = TRUE) * 1.12, .groups = "drop")

sig_pos_df <- sig_labels_df %>%
  left_join(max_scores, by = "cell_type") %>%
  mutate(
    label = sig,
    x_pos = 1.5
  )


group_colors <- c("HF" = "#E41A1C", "Control" = "#377EB8")

n_sig_cells <- length(sig_cells_plot)
ncols_box   <- min(4, n_sig_cells)
nrows_box   <- ceiling(n_sig_cells / ncols_box)
fig_width_8B  <- max(10, ncols_box * 3)
fig_height_8B <- max(6, nrows_box * 3.5)

p_8B <- ggplot(df_box, aes(x = group, y = score, fill = group)) +
  geom_boxplot(
    outlier.shape  = 21,
    outlier.size   = 1.2,
    outlier.fill   = "white",
    outlier.stroke = 0.4,
    width          = 0.55,
    lwd            = 0.4,
    alpha          = 0.85
  ) +
  geom_jitter(
    aes(color = group),
    width  = 0.15,
    size   = 0.8,
    alpha  = 0.5
  ) +

  geom_segment(
    data = sig_pos_df,
    aes(x = 1, xend = 2,
        y = max_score * 0.97, yend = max_score * 0.97),
    inherit.aes = FALSE,
    color = "black", linewidth = 0.4
  ) +
  geom_text(
    data = sig_pos_df,
    aes(x = x_pos, y = max_score * 1.02,
        label = ifelse(sig == "ns",
                       sprintf("p=%.3f", p_value),
                       label)),
    inherit.aes = FALSE,
    size = ifelse(sig_pos_df$sig == "ns", 2.5, 4),
    fontface = "bold",
    vjust = 0
  ) +
  scale_fill_manual(values  = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = c("HF" = "HF", "Control" = "Control")) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = ncols_box) +
  labs(
    title    = "Differential Immune Cell Infiltration: HF vs. Control",
    subtitle = "Wilcoxon rank-sum test | * p<0.05  ** p<0.01  *** p<0.001",
    x        = NULL,
    y        = "Infiltration Score",
    fill     = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey50"),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    strip.text       = element_text(face = "bold", size = 9),
    axis.text.x      = element_text(size = 9, face = "bold"),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

save_plot(p_8B, "Figure8B_immune_boxplot_significant",
          width = fig_width_8B, height = fig_height_8B)




cat("[INFO] 绘制 Figure 8C：免疫细胞间相关性热图...\n")


score_wide <- immune_df_filtered %>%
  pivot_wider(id_cols = sample, names_from = cell_type, values_from = score) %>%
  tibble::column_to_rownames("sample")


if (length(sig_cells) >= 2) {
  score_wide_8C <- score_wide[, intersect(sig_cells, colnames(score_wide)), drop = FALSE]
  cat(sprintf("[INFO] Fig 8C 使用显著差异细胞: %d 种\n", ncol(score_wide_8C)))
} else {
  score_wide_8C <- score_wide
  cat("[INFO] Fig 8C 显著差异细胞不足，使用全部22种\n")
}


corr_cell <- cor(score_wide_8C, method = "spearman", use = "pairwise.complete.obs")


n_ct_8C  <- ncol(score_wide_8C)
pval_cell <- matrix(NA_real_, n_ct_8C, n_ct_8C,
                    dimnames = list(colnames(score_wide_8C), colnames(score_wide_8C)))
for (i in seq_len(n_ct_8C)) {
  for (j in seq_len(n_ct_8C)) {
    if (i != j) {
      tt <- tryCatch(
        cor.test(score_wide_8C[, i], score_wide_8C[, j], method = "spearman", exact = FALSE),
        error = function(e) NULL
      )
      if (!is.null(tt)) pval_cell[i, j] <- tt$p.value
    }
  }
}

write.csv(corr_cell,  file.path(OUTPUT_DIR, "immune_cell_correlation_22cell.csv"))
write.csv(pval_cell,  file.path(OUTPUT_DIR, "immune_cell_correlation_pval.csv"))


corr_long_8C <- as.data.frame(as.table(corr_cell)) %>%
  setNames(c("Cell1", "Cell2", "r")) %>%
  left_join(
    as.data.frame(as.table(pval_cell)) %>% setNames(c("Cell1","Cell2","pval")),
    by = c("Cell1","Cell2")
  ) %>%
  mutate(
    sig_label = case_when(
      is.na(pval)    ~ "",
      pval < 0.001   ~ "***",
      pval < 0.01    ~ "**",
      pval < 0.05    ~ "*",
      TRUE           ~ ""
    ),

    row_idx = as.integer(factor(Cell1, levels = colnames(score_wide_8C))),
    col_idx = as.integer(factor(Cell2, levels = colnames(score_wide_8C))),
    keep    = row_idx >= col_idx
  ) %>%
  filter(keep)

cell_levels_8C <- colnames(score_wide_8C)

p_8C <- ggplot(corr_long_8C,
               aes(x = factor(Cell2, levels = cell_levels_8C),
                   y = factor(Cell1, levels = rev(cell_levels_8C)),
                   fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(Cell1 == Cell2, "",
                               paste0(sprintf("%.2f", r), "\n", sig_label))),
            size = 2.8, lineheight = 0.9, color = "black") +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Spearman r"
  ) +
  scale_x_discrete(position = "bottom") +
  coord_fixed() +
  labs(
    title    = "Spearman Correlation Between Differential Immune Cells",
    subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
    x        = NULL,
    y        = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    plot.subtitle   = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 8),
    panel.grid      = element_blank(),
    legend.key.height = unit(0.8, "cm")
  )

dim_8C <- max(5, n_ct_8C * 0.7 + 2)
save_plot(p_8C, "Figure8C_immune_cell_correlation",
          width = dim_8C, height = dim_8C * 0.92)




cat("[INFO] 绘制 Figure 8D：核心基因-免疫细胞相关性热图...\n")

key_genes <- c("FCN3", "MAP2K1", "FCER1G")
key_rows  <- intersect(key_genes, rownames(train_expr))

if (length(key_rows) == 0) {
  cat("[WARN] 核心基因不在表达矩阵中，跳过 Fig 8D\n")
} else {

  expr_sub     <- train_expr_norm[, reliable_samples, drop = FALSE]
  key_expr_sub <- t(expr_sub[key_rows, , drop = FALSE])

  ct_scores_all <- sapply(sig_cells, function(ct) {
    ssgsea_score_ranked(expr_sub, immune_sigs_22[[ct]])
  })


  corr_kg <- matrix(NA_real_, nrow = length(sig_cells), ncol = length(key_rows),
                    dimnames = list(sig_cells, key_rows))
  pval_kg <- corr_kg

  for (ct in sig_cells) {
    for (g in key_rows) {
      tt <- tryCatch(
        cor.test(key_expr_sub[, g], ct_scores_all[, ct], method = "spearman", exact = FALSE),
        error = function(e) NULL
      )
      if (!is.null(tt)) {
        corr_kg[ct, g] <- tt$estimate
        pval_kg[ct, g] <- tt$p.value
      }
    }
  }

  write.csv(corr_kg, file.path(OUTPUT_DIR, "immune_keygene_correlation_22cell.csv"))
  write.csv(pval_kg, file.path(OUTPUT_DIR, "immune_keygene_correlation_pval.csv"))
  cat(sprintf("[INFO] 22种细胞-核心基因相关性矩阵已保存（%d×%d）\n",
              nrow(corr_kg), ncol(corr_kg)))


  corr_long_8D <- as.data.frame(as.table(corr_kg)) %>%
    setNames(c("ImmuneCell", "Gene", "r")) %>%
    left_join(
      as.data.frame(as.table(pval_kg)) %>% setNames(c("ImmuneCell","Gene","pval")),
      by = c("ImmuneCell","Gene")
    ) %>%
    mutate(
      sig_label = case_when(
        is.na(pval)  ~ "",
        pval < 0.001 ~ "***",
        pval < 0.01  ~ "**",
        pval < 0.05  ~ "*",
        TRUE         ~ ""
      ),
      display = paste0(sprintf("%.2f", r),
                       ifelse(sig_label != "", paste0("\n", sig_label), ""))
    )


  ct_order_8D <- corr_long_8D %>%
    filter(Gene == key_rows[1]) %>%
    arrange(r) %>%
    pull(ImmuneCell) %>%
    as.character()

  corr_long_8D$ImmuneCell <- factor(corr_long_8D$ImmuneCell, levels = ct_order_8D)

  p_8D <- ggplot(corr_long_8D, aes(x = Gene, y = ImmuneCell, fill = r)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = display), size = 2.8, lineheight = 0.9, color = "black") +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "Spearman r"
    ) +
    scale_x_discrete(position = "top") +
    labs(
      title    = "Spearman Correlation: Key Biomarkers vs. Immune Cells",
      subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
      x        = "Key Biomarker Genes",
      y        = "Immune Cell Types"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey50"),
      axis.text.x   = element_text(face = "bold", size = 11, angle = 0, hjust = 0.5),
      axis.text.y   = element_text(size = 9),
      axis.title.x  = element_text(size = 9, color = "grey40"),
      axis.title.y  = element_text(size = 9, color = "grey40"),
      panel.grid    = element_blank(),
      legend.key.height = unit(0.8, "cm")
    )

  save_plot(p_8D, "Figure8D_immune_keygene_correlation",
            width  = max(5, length(key_rows) * 1.8 + 4),
            height = max(8, length(sig_cells)  * 0.5 + 2))
}




cat(sprintf("\n[SUMMARY]\n"))
cat(sprintf("  - 分析方法: CIBERSORT-like ssGSEA (PMID:25940772)\n"))
cat(sprintf("  - 免疫细胞种类: %d 种 (LM22标准分类)\n", length(sig_cells) ))
cat(sprintf("  - 过滤后可靠样本: %d 个\n", length(reliable_samples)))
cat(sprintf("  - HF样本数: %d, Control样本数: %d\n",
            sum(unique(immune_df_filtered[,c("sample","group")])$group == "HF"),
            sum(unique(immune_df_filtered[,c("sample","group")])$group == "Control")))
cat(sprintf("  - 显著差异免疫细胞 (p<0.05): %d 种\n", length(sig_cells)))
if (length(sig_cells) > 0) {
  cat(sprintf("    → %s\n", paste(sig_cells, collapse = "; ")))
}
cat(sprintf("  - Fig 8C: 差异细胞间相关热图 (%d×%d)\n", n_ct_8C, n_ct_8C))
if (length(key_rows) > 0) {
  cat(sprintf("  - Fig 8D: 核心基因-免疫细胞相关热图 (%d cells × %d genes)\n",
              length(sig_cells) , length(key_rows)))
  cat(sprintf("    → 基因: %s\n", paste(key_rows, collapse = ", ")))
}

cat(sprintf("\n[DONE] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
