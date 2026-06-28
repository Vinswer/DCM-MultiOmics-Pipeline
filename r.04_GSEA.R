rm(list = ls()); gc()

args    <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args) >= 1 && !is.na(as.integer(args[1]))) as.integer(args[1]) else 4L
cat(sprintf("[INFO] Using n_cores = %d\n", n_cores))

required_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "ggplot2",
                   "dplyr", "stringr", "enrichplot", "patchwork", "DOSE",
                   "hgu133plus2.db", "AnnotationDbi")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat("[SKIP] Required packages not available:", paste(missing_pkgs, collapse = ", "), "\n")
  quit(save = "no", status = 0, runLast = FALSE)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(enrichplot)
  library(patchwork)
  library(DOSE)
})




ORIGINAL_DIR <- ""
OUTPUT_DIR   <- file.path(ORIGINAL_DIR, "04_GSEA")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


TARGET_GENES <- c("FCN3", "FCER1G", "MAP2K1")


NES_CUTOFF  <- 1
FDR_CUTOFF  <- 0.25
PVAL_CUTOFF <- 0.05
TOP_N       <- 10

cat(sprintf("[INFO] 筛选阈值: NES > %.1f, FDR < %.2f, p < %.2f, TOP %d\n",
            NES_CUTOFF, FDR_CUTOFF, PVAL_CUTOFF, TOP_N))




save_plot <- function(p, prefix, width = 12, height = 8, dpi = 300) {
  ggsave(file.path(OUTPUT_DIR, paste0(prefix, ".pdf")),
         plot = p, width = width, height = height)
  ggsave(file.path(OUTPUT_DIR, paste0(prefix, ".png")),
         plot = p, width = width, height = height, dpi = dpi)
  cat(sprintf("[INFO] 图片已保存: %s\n", prefix))
}





expr_file <- file.path(ORIGINAL_DIR, "00_rawdata", "training_GSE57338", "GSE57338_gene_expression.csv")
cat(sprintf("[INFO] 读取表达矩阵: %s\n", expr_file))

if (!file.exists(expr_file)) {
  stop(sprintf("[ERROR] 表达矩阵文件不存在: %s\n请确认文件路径正确", expr_file))
}

expr_mat <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
cat(sprintf("[INFO] 表达矩阵维度: %d 基因 × %d 样本\n", nrow(expr_mat), ncol(expr_mat)))



























missing_targets <- TARGET_GENES[!TARGET_GENES %in% rownames(expr_mat)]
if (length(missing_targets) > 0) {
  stop(sprintf("[ERROR] 目标基因不在表达矩阵中: %s", paste(missing_targets, collapse = ", ")))
}
cat(sprintf("[INFO] 目标基因均存在: %s\n", paste(TARGET_GENES, collapse = ", ")))




run_gsea_corr <- function(target_gene, expr_mat,
                          nes_cutoff  = 1,
                          fdr_cutoff  = 0.25,
                          pval_cutoff = 0.05,
                          top_n       = 10) {

  cat(sprintf("\n══════════ [TARGET] %s ══════════\n", target_gene))




  target_expr <- expr_mat[target_gene, ]
  other_genes  <- rownames(expr_mat)[rownames(expr_mat) != target_gene]

  cat(sprintf("[INFO] 计算 %s 与 %d 个基因的Spearman相关性...\n",
              target_gene, length(other_genes)))

  cor_vals <- vapply(other_genes, function(g) {
    tryCatch(
      stats::cor(target_expr, expr_mat[g, ],
                 method = "spearman", use = "complete.obs"),
      error = function(e) NA_real_
    )
  }, numeric(1))

  cor_vals <- cor_vals[!is.na(cor_vals)]


  gene_list <- sort(cor_vals, decreasing = TRUE)
  cat(sprintf("[INFO] 排序基因数: %d（最大r=%.4f，最小r=%.4f）\n",
              length(gene_list), gene_list[1], tail(gene_list, 1)))


  cor_df <- data.frame(
    gene         = names(gene_list),
    spearman_r   = as.numeric(gene_list),
    target_gene  = target_gene,
    row.names    = NULL
  )
  write.csv(cor_df,
            file.path(OUTPUT_DIR, sprintf("correlation_spearman_%s.csv", target_gene)),
            row.names = FALSE)




  gene_df <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db, drop = TRUE)
  cat(sprintf("[INFO] ENTREZID转换成功: %d / %d\n", nrow(gene_df), length(gene_list)))

  gene_list_entrez <- gene_list[gene_df$SYMBOL]
  names(gene_list_entrez) <- gene_df$ENTREZID
  gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

  gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]

  results <- list()




  cat("[INFO] 运行 gseGO (BP)...\n")
  set.seed(42)
  gsea_go <- tryCatch(
    gseGO(geneList      = gene_list_entrez,
          OrgDb         = org.Hs.eg.db,
          ont           = "BP",
          minGSSize     = 15,
          maxGSSize     = 500,
          pvalueCutoff  = 1,
          pAdjustMethod = "BH",
          verbose       = FALSE,
          nPermSimple   = 1000),
    error = function(e) { cat("[WARN] gseGO失败:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    res_go <- as.data.frame(gsea_go) %>%
      mutate(Database = "GO_BP", target_gene = target_gene) %>%

      filter(abs(NES) > nes_cutoff,
             qvalue   < fdr_cutoff,
             pvalue   < pval_cutoff)
    cat(sprintf("[INFO] GO BP 筛选后: %d 条 (正向=%d, 负向=%d)\n",
                nrow(res_go), sum(res_go$NES > 0), sum(res_go$NES < 0)))
    if (nrow(res_go) > 0) {
      write.csv(res_go,
                file.path(OUTPUT_DIR, sprintf("GSEA_GO_%s.csv", target_gene)),
                row.names = FALSE)
      results[["GO_BP"]] <- list(obj = gsea_go, df = res_go)
    }
  } else {
    cat("[WARN] GO BP 无结果\n")
  }




  cat("[INFO] 运行 gseKEGG...\n")
  set.seed(42)
  gsea_kegg <- tryCatch(
    gseKEGG(geneList          = gene_list_entrez,
            organism          = "hsa",
            minGSSize         = 15,
            maxGSSize         = 500,
            pvalueCutoff      = 1,
            pAdjustMethod     = "BH",
            verbose           = FALSE,
            use_internal_data = TRUE,
            nPermSimple       = 1000),
    error = function(e) { cat("[WARN] gseKEGG失败:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
    res_kegg <- as.data.frame(gsea_kegg) %>%
      mutate(Database = "KEGG", target_gene = target_gene) %>%
      filter(abs(NES) > nes_cutoff,
             qvalue   < fdr_cutoff,
             pvalue   < pval_cutoff)
    cat(sprintf("[INFO] KEGG 筛选后: %d 条 (正向=%d, 负向=%d)\n",
                nrow(res_kegg), sum(res_kegg$NES > 0), sum(res_kegg$NES < 0)))
    if (nrow(res_kegg) > 0) {
      write.csv(res_kegg,
                file.path(OUTPUT_DIR, sprintf("GSEA_KEGG_%s.csv", target_gene)),
                row.names = FALSE)
      results[["KEGG"]] <- list(obj = gsea_kegg, df = res_kegg)
    }
  } else {
    cat("[WARN] KEGG 无显著结果\n")
  }





  if (length(results) > 0) {


    all_sig <- bind_rows(lapply(results, function(x) x$df)) %>%
      arrange(desc(abs(NES))) %>%
      slice_head(n = top_n)

    cat(sprintf("[INFO] TOP%d 通路（|NES|最大）:\n", top_n))
    print(all_sig[, c("Description", "Database", "NES", "pvalue", "qvalue")])



    for (db_name in names(results)) {
      res_df  <- results[[db_name]]$df
      gsea_obj <- results[[db_name]]$obj

      if (nrow(res_df) == 0) next


      top_ids <- res_df %>%
        arrange(desc(abs(NES))) %>%
        slice_head(n = min(top_n, nrow(res_df))) %>%
        pull(ID)


      n_paths  <- length(top_ids)
      path_colors <- if (n_paths <= 8) {
        RColorBrewer::brewer.pal(max(3, n_paths), "Dark2")[seq_len(n_paths)]
      } else {
        rainbow(n_paths)
      }


      if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
        path_colors <- rainbow(n_paths)
      }

      p_gsea <- tryCatch({
        enrichplot::gseaplot2(
          x          = gsea_obj,
          geneSetID  = top_ids,
          title      = sprintf("GSEA - %s\n(%s, TOP%d)", target_gene, db_name, n_paths),
          color      = path_colors,
          base_size  = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots   = 1:3,
          pvalue_table = FALSE
        )
      }, error = function(e) {
        cat(sprintf("[WARN] gseaplot2失败(%s): %s\n", db_name, conditionMessage(e)))
        NULL
      })

      if (!is.null(p_gsea)) {
        save_plot(p_gsea,
                  sprintf("GSEA_enrichplot_%s_%s", target_gene, db_name),
                  width = 10, height = 8)
      }
    }


    top10_bar <- all_sig %>%
      mutate(
        Label     = str_wrap(paste0("[", Database, "] ", Description), width = 48),
        Direction = ifelse(NES > 0,
                           sprintf("Pos (↑ corr. %s)", target_gene),
                           sprintf("Neg (↓ corr. %s)", target_gene))
      )

    p_bar <- ggplot(top10_bar,
                    aes(x = NES,
                        y = reorder(Label, NES),
                        fill = -log10(pvalue))) +
      geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey30", linetype = "dashed", linewidth = 0.5) +
      scale_fill_gradientn(
        colors = c("#4575B4", "#91BFDB", "#FEE090", "#FC8D59", "#D73027"),
        name   = "-log10(p-value)"
      ) +
      labs(
        title    = sprintf("GSEA TOP%d Pathways — %s", top_n, target_gene),
        subtitle = sprintf("Ranked by Spearman correlation | NES>%.1f, FDR<%.2f, p<%.2f",
                           nes_cutoff, fdr_cutoff, pval_cutoff),
        x = "Normalized Enrichment Score (NES)", y = NULL,
        caption = sprintf("Dataset: GSE57338 (training set) | Method: Spearman")
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.y      = element_text(size = 8),
        legend.position  = "right",
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        plot.caption     = element_text(size = 8, color = "grey50"),
        panel.grid.major.y = element_blank()
      )

    save_plot(p_bar,
              sprintf("GSEA_barplot_TOP10_%s", target_gene),
              width = 13, height = 7)


    p_dot <- ggplot(all_sig,
                    aes(x      = NES,
                        y      = reorder(str_wrap(Description, 40), NES),
                        size   = setSize,
                        color  = qvalue)) +
      geom_point(alpha = 0.85) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_gradientn(
        colors = c("#D73027", "#FC8D59", "#FEE090", "#91BFDB", "#4575B4"),
        name   = "FDR (q-value)",
        limits = c(0, fdr_cutoff)
      ) +
      scale_size_continuous(name = "Gene Set Size", range = c(3, 10)) +
      labs(
        title   = sprintf("%s — GSEA Dot Plot (TOP%d)", target_gene, top_n),
        x       = "NES", y = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.y     = element_text(size = 8),
        legend.position = "right",
        plot.title      = element_text(face = "bold")
      )

    save_plot(p_dot,
              sprintf("GSEA_dotplot_TOP10_%s", target_gene),
              width = 13, height = 7)


    write.csv(all_sig,
              file.path(OUTPUT_DIR, sprintf("GSEA_TOP10_%s.csv", target_gene)),
              row.names = FALSE)
  }

  return(results)
}




cat(sprintf("[START] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
all_results <- list()
for (tg in TARGET_GENES) {
  all_results[[tg]] <- run_gsea_corr(
    target_gene = tg,
    expr_mat    = expr_mat,
    nes_cutoff  = NES_CUTOFF,
    fdr_cutoff  = FDR_CUTOFF,
    pval_cutoff = PVAL_CUTOFF,
    top_n       = TOP_N
  )
}





cat("\n[INFO] 绘制整合GSEA图谱...\n")

all_top10_list <- lapply(TARGET_GENES, function(tg) {
  all_sig <- bind_rows(lapply(all_results[[tg]], function(x) x$df))
  if (nrow(all_sig) == 0) return(NULL)
  all_sig %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = TOP_N) %>%
    mutate(target_gene = tg)
})
all_top10_list <- Filter(Negate(is.null), all_top10_list)

if (length(all_top10_list) > 0) {
  combined_top <- bind_rows(all_top10_list)


  combined_top <- combined_top %>%
    mutate(
      Short_desc = str_wrap(Description, width = 40),
      Direction  = ifelse(NES > 0, "Positive", "Negative"),
      Gene_label = factor(target_gene, levels = TARGET_GENES)
    )

  p_combined <- ggplot(combined_top,
                       aes(x     = Gene_label,
                           y     = reorder(Short_desc, NES),
                           size  = abs(NES),
                           color = qvalue,
                           shape = Direction)) +
    geom_point(alpha = 0.85) +
    scale_color_gradientn(
      colors = c("#D73027", "#FC8D59", "#FEE090", "#91BFDB", "#4575B4"),
      name   = "FDR",
      limits = c(0, FDR_CUTOFF)
    ) +
    scale_size_continuous(name = "|NES|", range = c(3, 10)) +
    scale_shape_manual(values = c("Positive" = 16, "Negative" = 17),
                       name   = "Direction") +
    labs(
      title    = "Integrated GSEA Landscape: HF Phagocytosis-Immunity & Myocardial Remodeling",
      subtitle = sprintf("Core biomarkers: %s | GSE57338 (training set) | Spearman correlation-based ranking",
                         paste(TARGET_GENES, collapse = ", ")),
      x = "Core Biomarker", y = "Pathway",
      caption = sprintf("NES>%.1f, FDR<%.2f, p<%.2f | GO-BP + KEGG",
                        NES_CUTOFF, FDR_CUTOFF, PVAL_CUTOFF)
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y      = element_text(size = 7.5),
      axis.text.x      = element_text(face = "bold", size = 11),
      legend.position  = "right",
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 8.5, color = "grey40"),
      plot.caption     = element_text(size = 8, color = "grey50"),
      panel.grid.major.x = element_blank()
    )

  save_plot(p_combined, "GSEA_integrated_landscape",
            width = 15, height = max(10, nrow(combined_top) * 0.35 + 3))


  p_facet <- ggplot(combined_top,
                    aes(x = NES,
                        y = reorder(Short_desc, NES),
                        fill = -log10(pvalue))) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_gradientn(
      colors = c("#4575B4", "#91BFDB", "#FEE090", "#FC8D59", "#D73027"),
      name   = "-log10(p)"
    ) +
    facet_wrap(~ target_gene, scales = "free_y", ncol = 3) +
    labs(
      title    = "GSEA TOP10 Pathways per Core Biomarker",
      subtitle = "HF Phagocytosis-Immunity & Myocardial Remodeling | GSE57338",
      x = "NES", y = NULL,
      caption = sprintf("NES>%.1f, FDR<%.2f, p<%.2f | Spearman ranking | GO-BP + KEGG",
                        NES_CUTOFF, FDR_CUTOFF, PVAL_CUTOFF)
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.text       = element_text(face = "bold", size = 11),
      axis.text.y      = element_text(size = 7.5),
      legend.position  = "bottom",
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(color = "grey40", size = 9)
    )

  save_plot(p_facet, "GSEA_facet_TOP10_all_genes",
            width = 18, height = 10)


  write.csv(combined_top,
            file.path(OUTPUT_DIR, "GSEA_integrated_all_top10.csv"),
            row.names = FALSE)
  cat(sprintf("[INFO] 整合结果已保存（共 %d 条通路）\n", nrow(combined_top)))
}




cat("\n══════════ [SUMMARY] ══════════\n")
for (tg in TARGET_GENES) {
  res_list <- all_results[[tg]]
  n_go_pos   <- if (!is.null(res_list[["GO_BP"]])) sum(res_list[["GO_BP"]]$df$NES > 0) else 0L
  n_go_neg   <- if (!is.null(res_list[["GO_BP"]])) sum(res_list[["GO_BP"]]$df$NES < 0) else 0L
  n_kegg_pos <- if (!is.null(res_list[["KEGG"]])) sum(res_list[["KEGG"]]$df$NES > 0) else 0L
  n_kegg_neg <- if (!is.null(res_list[["KEGG"]])) sum(res_list[["KEGG"]]$df$NES < 0) else 0L
  cat(sprintf("[STAT] %s | GO_BP pos=%d neg=%d | KEGG pos=%d neg=%d | 合计=%d\n",
              tg, n_go_pos, n_go_neg, n_kegg_pos, n_kegg_neg,
              n_go_pos + n_go_neg + n_kegg_pos + n_kegg_neg))
}

cat(sprintf("[DONE] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
