rm(list = ls()); gc()

ORIGINAL_DIR <- ""
work_dir <- file.path(ORIGINAL_DIR, "00_rawdata")
if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

output_dir <- file.path(ORIGINAL_DIR, "01_deg")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir    <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(grid)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(stringr)

save_plot <- function(p, prefix, width = 7, height = 6, dpi = 300) {
  pdf_path <- file.path(fig_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(fig_dir, paste0(prefix, ".png"))
  ggsave(pdf_path, plot = p, width = width, height = height)
  ggsave(png_path, plot = p, width = width, height = height, dpi = dpi)
  cat("  已保存：", pdf_path, "\n")
  cat("  已保存：", png_path, "\n")
}

save_pheatmap <- function(ph, prefix, width = 10, height = 8, dpi = 300) {
  pdf_path <- file.path(fig_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(fig_dir, paste0(prefix, ".png"))
  pdf(pdf_path, width = width, height = height)
  grid.newpage(); grid.draw(ph$gtable)
  dev.off()
  png(png_path, width = width, height = height, units = "in", res = dpi)
  grid.newpage(); grid.draw(ph$gtable)
  dev.off()
  cat("  已保存：", pdf_path, "\n")
  cat("  已保存：", png_path, "\n")
}


cat(">>> 读取训练集 GSE57338 ...\n")
train_expr <- read.csv(
  file.path(work_dir, "training_GSE57338", "GSE57338_gene_expression.csv"),
  row.names = 1, check.names = FALSE
)
train_meta <- read.csv(
  file.path(work_dir, "training_GSE57338", "GSE57338_group.csv"),
  stringsAsFactors = FALSE
)
train_meta$group <- ifelse(train_meta$group == "yes", "HF",
                           ifelse(train_meta$group == "no", "Normal", NA))
rownames(train_meta) <- train_meta$sample_id

cat(">>> 读取吞噬调节因子基因集 ...\n")
phago_genes <- read.csv(
  file.path(work_dir, "gene_sets", "Phagocytosis_all_genes_unique.csv"),
  stringsAsFactors = FALSE
)$gene_symbol

train_group_col    <- "group"
train_hf_label     <- "HF"
train_normal_label <- "Normal"

align_data <- function(expr, meta) {
  common_samples <- intersect(colnames(expr), rownames(meta))
  expr <- expr[, common_samples, drop = FALSE]
  meta <- meta[common_samples, , drop = FALSE]
  list(expr = expr, meta = meta)
}


train <- align_data(train_expr, train_meta)
cat(sprintf("训练集：%d 基因 × %d 样本\n", nrow(train$expr), ncol(train$expr)))
cat(sprintf("  基因数：%d，示例行名：%s\n",
            nrow(train_expr),
            paste(head(rownames(train_expr), 5), collapse = ", ")))

run_limma <- function(expr, meta, grp_col, hf_lab, normal_lab) {
  group <- as.character(meta[[grp_col]])
  keep  <- group %in% c(hf_lab, normal_lab)
  cat(sprintf("  筛选后样本数：%d\n", sum(keep)))
  expr  <- expr[, keep, drop = FALSE]
  group <- group[keep]
  group_factor <- factor(group, levels = c(normal_lab, hf_lab))
  design <- model.matrix(~ group_factor)
  mat  <- data.matrix(expr)
  fit  <- lmFit(mat, design)
  fit  <- eBayes(fit)
  res  <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "P")
  res$gene <- rownames(res)
  res
}

deg_train <- run_limma(train$expr, train$meta, train_group_col, train_hf_label, train_normal_label)
write.csv(deg_train, file.path(output_dir, "DEG_GSE57338_full.csv"), row.names = FALSE)

filter_degs <- function(res, padj_cut = 0.05, lfc_cut = 0.5) {
  res %>%
    filter(adj.P.Val < padj_cut, abs(logFC) > lfc_cut) %>%
    mutate(direction = ifelse(logFC > 0, "Up", "Down"))
}

sig_train <- filter_degs(deg_train)
cat(sprintf("训练集显著 DEG：%d 个（Up: %d，Down: %d）\n",
            nrow(sig_train), sum(sig_train$direction == "Up"), sum(sig_train$direction == "Down")))

write.csv(sig_train, file.path(output_dir, "DEG_GSE57338_significant.csv"), row.names = FALSE)


plot_volcano <- function(res, title, prefix, top_n = 10,
                         padj_cut = 0.05, lfc_cut = 0.5) {
  res <- res %>%
    mutate(
      sig = case_when(
        adj.P.Val < padj_cut & logFC >  lfc_cut ~ "Up",
        adj.P.Val < padj_cut & logFC < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      ),
      neg_log10_padj = -log10(adj.P.Val)
    )

  top_up   <- res %>% filter(sig == "Up")   %>% arrange(desc(logFC)) %>% head(top_n)
  top_down <- res %>% filter(sig == "Down") %>% arrange(logFC)       %>% head(top_n)
  label_genes <- bind_rows(top_up, top_down)

  color_vals <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")

  p <- ggplot(res, aes(x = logFC, y = neg_log10_padj, color = sig)) +
    geom_point(size = 1.2, alpha = 0.7) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    geom_text_repel(
      data          = label_genes,
      aes(label     = gene),
      size          = 3,
      max.overlaps  = 30,
      segment.color = "grey50",
      show.legend   = FALSE
    ) +
    scale_color_manual(
      values = color_vals,
      labels = c("Up" = "Up-regulated", "Down" = "Down-regulated", "NS" = "Not significant")
    ) +
    labs(
      title = title,
      x     = expression(log[2]~Fold~Change),
      y     = expression(-log[10]~(adjusted~p~value)),
      color = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top"
    ) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("Up: ",   sum(res$sig == "Up"),
                            "\nDown: ", sum(res$sig == "Down")),
             size = 3.5, color = "grey30")

  save_plot(p, prefix, width = 7, height = 6)
  invisible(p)
}

cat(">>> 绘制火山图 ...\n")

plot_volcano(deg_train,
             title  = "Volcano Plot - Training Set (GSE57338)\nHF vs Normal",
             prefix = "01_volcano_GSE57338")

plot_heatmap <- function(expr, meta, deg_df, grp_col, hf_lab, normal_lab,
                         title, prefix,
                         show_top_n = 5) {
  group <- as.character(meta[[grp_col]])
  keep  <- group %in% c(hf_lab, normal_lab)
  expr  <- expr[, keep, drop = FALSE]
  group <- group[keep]


  genes       <- deg_df$gene
  genes_avail <- intersect(genes, rownames(expr))
  if (length(genes_avail) == 0) stop("热图基因集与表达矩阵无交集")

  mat   <- data.matrix(expr[genes_avail, ])
  mat_z <- t(scale(t(mat)))
  mat_z[is.nan(mat_z)] <- 0


  sample_order <- order(group)
  mat_z        <- mat_z[, sample_order]
  group_sorted <- group[sample_order]


  top_up   <- deg_df %>% filter(direction == "Up")   %>%
    arrange(desc(logFC)) %>% head(show_top_n) %>% pull(gene)
  top_down <- deg_df %>% filter(direction == "Down") %>%
    arrange(logFC)       %>% head(show_top_n) %>% pull(gene)
  label_genes <- intersect(c(top_up, top_down), genes_avail)


  row_labels        <- rep("", length(genes_avail))
  names(row_labels) <- genes_avail
  row_labels[label_genes] <- label_genes

  anno_col <- data.frame(Group = group_sorted, row.names = colnames(mat_z))
  anno_colors <- list(
    Group = setNames(c("#E41A1C", "#377EB8"), c(hf_lab, normal_lab))
  )

  n_genes      <- length(genes_avail)
  fontsize_row <- 7

  cellheight <- if (n_genes > 200) NA else if (n_genes > 100) 3 else 5

  ph <- pheatmap(
    mat_z,
    annotation_col    = anno_col,
    annotation_colors = anno_colors,
    cluster_cols      = FALSE,
    show_colnames     = FALSE,
    show_rownames     = TRUE,
    labels_row        = row_labels,
    fontsize_row      = fontsize_row,
    cellheight        = cellheight,
    color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    breaks            = seq(-2, 2, length.out = 101),
    main              = title,
    border_color      = NA,
    annotation_legend = TRUE,
    gaps_row          = NULL,
    silent            = TRUE
  )

  height <- max(8, min(n_genes * 0.06 + 4, 40))
  save_pheatmap(ph, prefix, width = 12, height = height)
  invisible(ph)
}
cat(">>> 绘制热图 ...\n")
plot_heatmap(train$expr, train$meta, sig_train,
             train_group_col, train_hf_label, train_normal_label,
             title  = "Significant DEGs Heatmap - Training Set (GSE57338)",
             prefix = "02_heatmap_GSE57338_DEGs")


cat(">>> 交集分析（训练集 DEG × 吞噬调节因子基因集）...\n")

candidate_genes <- intersect(sig_train$gene, phago_genes)
cat(sprintf("训练集显著 DEG 数目：%d\n",  length(sig_train$gene)))
cat(sprintf("吞噬基因集数目：%d\n",        length(phago_genes)))
cat(sprintf("交集候选基因数目：%d\n",      length(candidate_genes)))
cat("候选基因：", paste(candidate_genes, collapse = ", "), "\n")

write.csv(data.frame(gene = candidate_genes),
          file.path(output_dir, "Candidate_genes_DEG_x_Phagocytosis.csv"),
          row.names = FALSE)

n_deg        <- length(sig_train$gene)
n_phago      <- length(phago_genes)
n_intersect  <- length(candidate_genes)
n_only_deg   <- n_deg   - n_intersect
n_only_phago <- n_phago - n_intersect
total_all <- n_only_deg + n_only_phago + n_intersect

pct_only_deg   <- round(n_only_deg   / total_all * 100, 2)
pct_intersect  <- round(n_intersect  / total_all * 100, 2)
pct_only_phago <- round(n_only_phago / total_all * 100, 2)


cat(">>> 绘制 Venn 图 ...\n")
total_deg   <- n_deg
total_phago <- n_phago


make_circle <- function(cx, cy, r, n = 300) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

circle_left  <- make_circle(-0.55, 0, 1.0)
circle_right <- make_circle( 0.55, 0, 0.85)
circle_left$group  <- "left"
circle_right$group <- "right"

p_venn <- ggplot() +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = unit(c(20, 20, 20, 20), "pt")
  ) +
  geom_polygon(data = circle_left,  aes(x = x, y = y),
               fill = "#8B0000", alpha = 0.80, color = NA) +
  geom_polygon(data = circle_right, aes(x = x, y = y),
               fill = "#1B4F72", alpha = 0.80, color = NA) +
  annotate("text", x = -0.90, y = 0,
           label = sprintf("%d\n(%.2f%%)", n_only_deg, pct_only_deg),
           size = 8, fontface = "bold", color = "#FFAAAA", lineheight = 1.1) +
  annotate("text", x =  0.05, y = 0,
           label = sprintf("%d\n(%.2f%%)", n_intersect, pct_intersect),
           size = 8, fontface = "bold", color = "#FFE0E0", lineheight = 1.1) +
  annotate("text", x =  0.95, y = 0,
           label = sprintf("%d\n(%.2f%%)", n_only_phago, pct_only_phago),
           size = 8, fontface = "bold", color = "#AAD4E8", lineheight = 1.1) +
  annotate("text", x = -1.10, y = 1.25,
           label = "Training Set DEGs\n(HF vs Normal)",
           size = 4.5, color = "#E07070", fontface = "bold", hjust = 0.5,
           lineheight = 1.3) +
  annotate("text", x =  1.10, y = 1.25,
           label = "Phagocytosis\nRegulator Genes",
           size = 4.5, color = "#70B0D0", fontface = "bold", hjust = 0.5,
           lineheight = 1.3) +
  coord_fixed(xlim = c(-2.0, 2.0), ylim = c(-1.4, 1.4))

save_plot(p_venn, "03_Venn_DEG_Phagocytosis", width = 8, height = 7)


cat(">>> GO / KEGG 富集分析 ...\n")

entrez_ids <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = candidate_genes,
    columns = "ENTREZID",
    keytype = "SYMBOL"
  )
)
entrez_ids <- entrez_ids[!is.na(entrez_ids$ENTREZID), ]
entrez_vec <- unique(entrez_ids$ENTREZID)
cat(sprintf("候选基因转换为 Entrez ID：%d 个\n", length(entrez_vec)))

go_res <- enrichGO(
  gene          = entrez_vec,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

if (!is.null(go_res) && nrow(go_res@result) > 0) {
  write.csv(as.data.frame(go_res),
            file.path(output_dir, "GO_enrichment_candidate_genes.csv"),
            row.names = FALSE)

  go_df <- as.data.frame(go_res) %>%
    filter(p.adjust < 0.05) %>%
    group_by(ONTOLOGY) %>%
    slice_max(order_by = Count, n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      ONTOLOGY    = factor(ONTOLOGY, levels = c("BP", "CC", "MF")),
      Description = str_wrap(Description, width = 40)
    )

  ont_colors <- c("BP" = "#E74C3C", "CC" = "#2ECC71", "MF" = "#3498DB")
  n_terms    <- nrow(go_df)
  go_h <- max(12, n_terms * 0.5 + 5)
  go_w       <- 14

  p_go <- ggplot(go_df, aes(x = Count, y = Description, fill = ONTOLOGY)) +
    geom_bar(stat = "identity", width = 0.75) +
    scale_fill_manual(values = ont_colors, name = "Ontology") +
    facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = "GO Enrichment Analysis\n(Candidate Phagocytosis-related DEGs in HF)",
      x     = "Count",
      y     = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title       = element_text(size = 16, face = "bold"),
      axis.text.y      = element_text(size = 12),
      axis.text.x      = element_text(size = 12),
      strip.text       = element_text(size = 13, face = "bold"),
    )

  save_plot(p_go, "04_GO_barplot_candidate_genes", width = go_w, height = go_h)
}


kegg_res <- enrichKEGG(
  gene          = entrez_vec,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
  kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(kegg_res),
            file.path(output_dir, "KEGG_enrichment_candidate_genes.csv"),
            row.names = FALSE)

  kegg_df <- as.data.frame(kegg_res) %>%
    filter(p.adjust < 0.05) %>%
    slice_max(order_by = Count, n = 20, with_ties = FALSE) %>%
    mutate(
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      Description = factor(Description,
                           levels = Description[order(GeneRatio_num)])
    )

  p_kegg <- ggplot(kegg_df,
                   aes(x = GeneRatio_num, y = Description,
                       color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradientn(
      colors = c("#E41A1C", "#FF7F00", "#4DAF4A", "#377EB8"),
      name   = "q value"
    ) +
    scale_size_continuous(name = "Count", range = c(3, 10)) +
    labs(
      title = "KEGG Pathway Enrichment Analysis\n(Candidate Phagocytosis-related DEGs in HF)",
      x     = "GeneRatio",
      y     = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y     = element_text(size = 9),
      legend.position = "right"
    )

  n_kegg <- nrow(kegg_df)
  kegg_h <- max(5, min(n_kegg * 0.42 + 2, 18))
  save_plot(p_kegg, "05_KEGG_bubble_candidate_genes", width = 10, height = kegg_h)
} else {
  cat("  KEGG 富集无显著结果，跳过绘图。\n")
}

cat("\n════════════════════════════════════════\n")
cat("分析完成！输出目录：", output_dir, "\n")
cat("  DEG_GSE57338_full.csv\n")
cat("  DEG_GSE57338_significant.csv\n")
cat("  Candidate_genes_DEG_x_Phagocytosis.csv\n")
cat("  GO_enrichment_candidate_genes.csv\n")
cat("  KEGG_enrichment_candidate_genes.csv\n")
cat("  figures/01_volcano_GSE57338.pdf/.png\n")
cat("  figures/02_heatmap_GSE57338_DEGs.pdf/.png\n")
cat("  figures/03_Venn_DEG_Phagocytosis.pdf/.png\n")
cat("  figures/04_GO_barplot_candidate_genes.pdf/.png\n")
cat("  figures/05_KEGG_bubble_candidate_genes.pdf/.png\n")
cat("════════════════════════════════════════\n")
