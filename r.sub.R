args    <- commandArgs(trailingOnly = TRUE)
n_cores <- ifelse(length(args) >= 1, as.integer(args[1]), 4)
cat("[INFO] 使用核心数:", n_cores, "\n")

SERVER_DIR   <- ""
KEY_GENES    <- c("FCN3", "FCER1G", "MAP2K1")
KEY_CELL     <- "Macrophages"

cat("[INFO] 关键基因:", paste(KEY_GENES, collapse = ", "), "\n")
cat("[INFO] 关键细胞:", KEY_CELL, "\n")

calc_plot_height <- function(n_items,
                             base_height   = 5,
                             per_item_inch = 0.35,
                             min_height    = 6,
                             max_height    = 24) {
  h <- base_height + n_items * per_item_inch
  max(min_height, min(max_height, h))
}

wrap_pathway_name <- function(x, width = 55) {
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3), "..."), x)
}

required_pkgs_ucell <- c("UCell", "Seurat", "ggplot2", "dplyr", "tidyr")
skip_ucell <- FALSE
for (pkg in required_pkgs_ucell) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("[SKIP] 缺少必要包:", pkg, "，跳过 UCell 分析\n")
    skip_ucell <- TRUE; break
  }
}

if (!skip_ucell) {
  suppressPackageStartupMessages({
    library(UCell); library(Seurat)
    library(ggplot2); library(dplyr); library(tidyr)
  })
  cat("[INFO] Seurat:", as.character(packageVersion("Seurat")),
      "| UCell:", as.character(packageVersion("UCell")), "\n")

  input_rds_ucell   <- file.path(SERVER_DIR, "03_single_cell/seurat_annotated.rds")
  key_genes_file    <- file.path(SERVER_DIR, "02_ml/key_genes.txt")
  output_dir_ucell  <- file.path(SERVER_DIR, "16_ssgsea")
  dir.create(output_dir_ucell, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(key_genes_file)) {
    key_genes_ucell <- readLines(key_genes_file)
    key_genes_ucell <- key_genes_ucell[nchar(trimws(key_genes_ucell)) > 0]
  } else {
    key_genes_ucell <- KEY_GENES
    cat("[WARN] key_genes.txt 不存在，使用全局 KEY_GENES\n")
  }
  cat("[INFO] key_genes:", paste(key_genes_ucell, collapse = ", "), "\n")

  if (!file.exists(input_rds_ucell)) {
    cat("[SKIP] 输入文件不存在:", input_rds_ucell, "\n")
  } else {
    cat("[INFO] 读取:", input_rds_ucell, "\n")
    sc_u <- readRDS(input_rds_ucell)
    cat("[INFO] 细胞数:", ncol(sc_u), "| 基因数:", nrow(sc_u), "\n")

    tryCatch({
      if (inherits(sc_u[["RNA"]], "Assay5"))
        sc_u[["RNA"]] <- JoinLayers(sc_u[["RNA"]])
    }, error = function(e)
      cat("[WARN] JoinLayers 失败:", conditionMessage(e), "\n"))

    meta_u <- sc_u@meta.data
    group_col_u <- if ("celltype_annotation" %in% colnames(meta_u)) "celltype_annotation" else
      if ("cell_type"           %in% colnames(meta_u)) "cell_type" else "seurat_clusters"
    cat("[INFO] 分组列:", group_col_u, "\n")

    key_genes_exist_u <- intersect(key_genes_ucell, rownames(sc_u))
    missing_u <- setdiff(key_genes_ucell, rownames(sc_u))
    if (length(missing_u) > 0)
      cat("[WARN] 基因不在数据集中，已忽略:", paste(missing_u, collapse = ", "), "\n")

    if (length(key_genes_exist_u) == 0) {
      cat("[SKIP] key_genes 均不存在于数据集，跳过 UCell\n")
    } else {
      counts_mat_u <- tryCatch(
        GetAssayData(sc_u, assay = "RNA", layer = "counts"),
        error = function(e) tryCatch(
          GetAssayData(sc_u, assay = "RNA", slot  = "counts"),
          error = function(e2) { cat("[SKIP] 无法获取 counts\n"); NULL }
        )
      )

      if (!is.null(counts_mat_u)) {
        gene_sets_u <- list(HF_biomarkers = key_genes_exist_u)

        bpparam_u <- tryCatch({
          if (requireNamespace("BiocParallel", quietly = TRUE)) {
            if (n_cores > 1) BiocParallel::MulticoreParam(n_cores)
            else             BiocParallel::SerialParam()
          } else NULL
        }, error = function(e) NULL)

        sc_u <- tryCatch({
          if (!is.null(bpparam_u))
            AddModuleScore_UCell(sc_u, features = gene_sets_u, assay = "RNA", BPPARAM = bpparam_u)
          else
            AddModuleScore_UCell(sc_u, features = gene_sets_u, assay = "RNA")
        }, error = function(e) {
          cat("[WARN] AddModuleScore_UCell 失败:", conditionMessage(e), "\n")
          tryCatch({
            ucell_s <- UCell::ScoreSignatures_UCell(
              matrix = counts_mat_u, features = gene_sets_u,
              maxRank = 1500, ncores = 1)
            sc_u@meta.data$HF_biomarkers_UCell <-
              ucell_s[rownames(sc_u@meta.data), "HF_biomarkers_UCell"]
            sc_u
          }, error = function(e2) {
            cat("[SKIP] UCell 全部失败\n"); sc_u
          })
        })

        ucell_col_u <- "HF_biomarkers_UCell"
        if (ucell_col_u %in% colnames(sc_u@meta.data)) {
          ucell_vec_u <- sc_u@meta.data[[ucell_col_u]]
          ucell_mean_u  <- mean(ucell_vec_u, na.rm = TRUE)
          ucell_high_u  <- sum(ucell_vec_u > 0.1, na.rm = TRUE)
          cat(sprintf("[STAT] UCell_mean_score: %.3f\n", ucell_mean_u))
          cat(sprintf("[STAT] UCell_high_score_cells: %d\n", ucell_high_u))

          score_df_u <- data.frame(
            cell      = rownames(sc_u@meta.data),
            cell_type = sc_u@meta.data[[group_col_u]],
            UCell     = ucell_vec_u, stringsAsFactors = FALSE)

          celltype_scores_u <- score_df_u %>%
            group_by(cell_type) %>%
            summarise(mean_UCell = mean(UCell, na.rm = TRUE),
                      median_UCell = median(UCell, na.rm = TRUE),
                      n_cells = n(), n_high = sum(UCell > 0.1, na.rm = TRUE),
                      .groups = "drop") %>%
            arrange(desc(mean_UCell))

          write.csv(celltype_scores_u,
                    file.path(output_dir_ucell, "ssgsea_celltype_scores.csv"),
                    row.names = FALSE)
          cat("[INFO] ssgsea_celltype_scores.csv 已保存\n")

          if (!"umap" %in% names(sc_u@reductions)) {
            tryCatch({
              if (!"pca" %in% names(sc_u@reductions)) {
                sc_u <- NormalizeData(sc_u, verbose = FALSE)
                sc_u <- FindVariableFeatures(sc_u, nfeatures = 2000, verbose = FALSE)
                sc_u <- ScaleData(sc_u, verbose = FALSE)
                sc_u <- RunPCA(sc_u, npcs = 30, verbose = FALSE)
              }
              sc_u <- RunUMAP(sc_u, reduction = "pca", dims = 1:20,
                              n.components = 2L, seed.use = 42, verbose = FALSE)
            }, error = function(e) cat("[WARN] RunUMAP 失败\n"))
          }

          tryCatch({
            if (!"umap" %in% names(sc_u@reductions)) stop("无 UMAP")
            umap_emb <- as.data.frame(sc_u@reductions[["umap"]]@cell.embeddings)
            colnames(umap_emb) <- c("UMAP_1", "UMAP_2")
            umap_emb[[ucell_col_u]] <- sc_u@meta.data[rownames(umap_emb), ucell_col_u]
            umap_emb <- umap_emb[order(umap_emb[[ucell_col_u]]), ]

            p_umap_u <- ggplot(umap_emb, aes(UMAP_1, UMAP_2, color = .data[[ucell_col_u]])) +
              geom_point(size = 0.4, alpha = 0.85) +
              scale_color_gradientn(
                colours = c("lightgrey", "#FEE08B", "#F46D43", "#A50026"),
                name = "UCell\nScore") +
              ggtitle(paste0("UCell Score — HF Biomarkers (",
                             paste(key_genes_exist_u, collapse = " / "), ")")) +
              theme_bw() +
              theme(plot.title = element_text(hjust = 0.5, size = 12),
                    panel.grid = element_blank())
            ggsave(file.path(output_dir_ucell, "ssgsea_umap.png"),
                   p_umap_u, width = 9, height = 7, dpi = 150)
            cat("[INFO] ssgsea_umap.png 已保存\n")
          }, error = function(e) cat("[WARN] ssgsea_umap.png 失败\n"))

          tryCatch({
            ct_order_u <- celltype_scores_u %>% arrange(desc(mean_UCell)) %>% pull(cell_type)
            score_df_u$cell_type <- factor(score_df_u$cell_type, levels = ct_order_u)
            n_ct_u <- length(ct_order_u)
            ct_colors_u <- scales::hue_pal()(n_ct_u)
            names(ct_colors_u) <- ct_order_u

            p_vln_u <- ggplot(score_df_u, aes(cell_type, UCell, fill = cell_type)) +
              geom_violin(scale = "width", trim = TRUE, linewidth = 0.3) +
              geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.2,
                           outlier.alpha = 0.3, linewidth = 0.3) +
              geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", linewidth = 0.5) +
              scale_fill_manual(values = ct_colors_u) +
              ggtitle("UCell Score Distribution by Cell Type") +
              ylab("UCell Score") +
              annotate("text", x = n_ct_u * 0.85, y = 0.105,
                       label = "threshold = 0.1", color = "red", size = 3) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                    axis.title.x = element_blank(),
                    plot.title = element_text(hjust = 0.5, size = 12),
                    legend.position = "none")
            vln_w <- max(10, n_ct_u * 0.9)
            ggsave(file.path(output_dir_ucell, "ssgsea_violin.png"),
                   p_vln_u, width = vln_w, height = 7, dpi = 150)
            cat("[INFO] ssgsea_violin.png 已保存\n")
          }, error = function(e) cat("[WARN] ssgsea_violin.png 失败\n"))
        }
      }
    }
  }
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
})

HAS_REACTOME_23b <- requireNamespace("ReactomePA", quietly = TRUE)
cat("[INFO] ReactomePA 可用:", HAS_REACTOME_23b, "\n")

out_dir_23b <- file.path(SERVER_DIR, "r.03_single_cell/monocyte_subtype")
markers_file_23b <- file.path(out_dir_23b, "subtype_markers.csv")

if (!file.exists(markers_file_23b)) {
  cat("[SKIP] subtype_markers.csv 不存在，跳过 Figure 24A\n")
} else {
  markers_23b <- read.csv(markers_file_23b, stringsAsFactors = FALSE)
  cat(sprintf("[INFO] markers: %d 行\n", nrow(markers_23b)))

  sig_markers_23b <- markers_23b %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
    arrange(cluster, p_val_adj)
  clusters_23b <- sort(unique(sig_markers_23b$cluster))
  cat(sprintf("[INFO] 亚型数量: %d\n", length(clusters_23b)))

  cat("[INFO] SYMBOL → ENTREZID...\n")
  gene_list_23b <- lapply(clusters_23b, function(cl) {
    genes <- unique(sig_markers_23b %>% filter(cluster == cl) %>% pull(gene))
    genes <- genes[nchar(genes) > 0]
    if (length(genes) < 5) { cat(sprintf("[WARN] Cluster %s 基因 < 5，跳过\n", cl)); return(NULL) }
    mapped <- suppressMessages(bitr(genes, fromType = "SYMBOL",
                                    toType = "ENTREZID", OrgDb = org.Hs.eg.db))
    if (nrow(mapped) == 0) { cat(sprintf("[WARN] Cluster %s 映射为空\n", cl)); return(NULL) }
    mapped$ENTREZID
  })

  names(gene_list_23b) <- paste0("Subtype_", clusters_23b)
  gene_list_23b <- Filter(Negate(is.null), gene_list_23b)


  MAC_SUBTYPE_MARKERS <- list(

    "FOLR2+ Res-Mac"     = c("FOLR2", "LYVE1", "MRC1", "CD163", "STAB1"),
    "C1Q+ Res-Mac"       = c("C1QA", "C1QB", "C1QC", "RNASE1", "APOE"),
    "SPP1+ Mac"          = c("SPP1", "TREM2", "GPNMB", "FABP5", "LPL"),
    "LYVE1+ Mac"         = c("LYVE1", "FOLR2", "SELENOP", "STAB1"),

    "IL1B+ Inflam-Mac"   = c("IL1B", "IL6", "CCL3", "CCL4", "CXCL8", "TNF"),
    "ISG15+ IFN-Mac"     = c("ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "IFI44L"),
    "S100A8+ Mono-Mac"   = c("S100A8", "S100A9", "CD14", "FCN1", "VCAN"),

    "FCN3+ Mac"          = c("FCN3", "FCER1G", "MS4A7", "CX3CR1"),

    "MKI67+ Prolif-Mac"  = c("MKI67", "TOP2A", "PCNA", "STMN1"),

    "CD68+ Mac"          = c("CD68", "CD14", "CSF1R", "ITGAM", "MRC1")
  )


  annotate_subtype <- function(cluster_id, markers_df, marker_dict,
                               cell_type_suffix = "Macrophages") {

    cluster_genes <- markers_df %>%
      filter(cluster == cluster_id, p_val_adj < 0.05, avg_log2FC > 0.25) %>%
      arrange(p_val_adj) %>%
      pull(gene) %>%
      unique()

    if (length(cluster_genes) == 0) return(paste0("Unknown ", cell_type_suffix))

    hit_counts <- sapply(marker_dict, function(markers) {
      sum(markers %in% cluster_genes)
    })

    best_idx <- which.max(hit_counts)
    if (hit_counts[best_idx] == 0) return(paste0("Unknown ", cell_type_suffix))

    subtype_name <- names(marker_dict)[best_idx]
    top_marker   <- marker_dict[[best_idx]][marker_dict[[best_idx]] %in% cluster_genes][1]
    cat(sprintf("  Cluster %s → %s (top hit: %s, matched: %d genes)\n",
                cluster_id, subtype_name, top_marker, hit_counts[best_idx]))
    return(subtype_name)
  }


  cat("[INFO] 对巨噬细胞亚型进行 marker 注释...\n")
  subtype_annotation <- setNames(
    sapply(clusters_23b, function(cl)
      annotate_subtype(cl, sig_markers_23b, MAC_SUBTYPE_MARKERS, "Macrophages")),
    paste0("Subtype_", clusters_23b)
  )

  cat("[INFO] 注释结果：\n")
  for (nm in names(subtype_annotation))
    cat(sprintf("  %s → %s\n", nm, subtype_annotation[[nm]]))


  subtype_annotation <- setNames(
    sapply(clusters_23b, function(cl)
      annotate_subtype(cl, sig_markers_23b, MAC_SUBTYPE_MARKERS, "Macrophages")),
    paste0("Subtype_", clusters_23b)
  )


  seen <- list()
  for (i in seq_along(subtype_annotation)) {
    nm <- subtype_annotation[i]
    if (!is.null(seen[[nm]])) {
      seen[[nm]] <- seen[[nm]] + 1L

      first_pos <- which(subtype_annotation == nm)[1]
      if (!grepl("-\\d+$", subtype_annotation[first_pos]))
        subtype_annotation[first_pos] <- paste0(subtype_annotation[first_pos], "-1")
      subtype_annotation[i] <- paste0(nm, "-", seen[[nm]])
    } else {
      seen[[nm]] <- 1L
    }
  }


  cat("[INFO] 注释结果（去重后）：\n")
  for (nm in names(subtype_annotation))
    cat(sprintf("  %s → %s\n", nm, subtype_annotation[[nm]]))


  annotation_df <- data.frame(
    cluster_id  = names(subtype_annotation),
    annotation  = unname(subtype_annotation),
    stringsAsFactors = FALSE
  )
  write.csv(annotation_df,
            file.path(out_dir_23b, "subtype_annotation_map.csv"),
            row.names = FALSE)
  cat("[INFO] subtype_annotation_map.csv 已保存\n")

  names(gene_list_23b) <- subtype_annotation[names(gene_list_23b)]

  if (length(gene_list_23b) > 0) {

    cat("[INFO] compareCluster GO BP...\n")
    ck_go_23b <- tryCatch(
      compareCluster(gene_list_23b, fun = "enrichGO",
                     OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05,
                     qvalueCutoff  = 0.2,  readable = TRUE),
      error = function(e) { cat("[WARN] GO 失败:", conditionMessage(e), "\n"); NULL }
    )

    if (!is.null(ck_go_23b) && nrow(ck_go_23b) > 0) {

      go_df <- as.data.frame(ck_go_23b)
      go_df$Description <- wrap_pathway_name(go_df$Description, width = 55)
      ck_go_23b@compareClusterResult$Description <- go_df$Description

      n_go_terms <- length(unique(go_df$Description))
      go_h <- calc_plot_height(n_go_terms, base_height = 5,
                               per_item_inch = 0.32, min_height = 7, max_height = 28)

      p_go_23b <- dotplot(ck_go_23b, showCategory = 5,
                          title = "Fibroblast Subtype — GO Biological Process Enrichment",
                          font.size = 9) +
        theme(
          plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text.x     = element_text(angle = 45, hjust = 1),
          axis.text.y     = element_text(size = 8, lineheight = 0.85),
          plot.margin = ggplot2::margin(t = 10, r = 20, b = 10, l = 160)
        )
      ggsave(file.path(out_dir_23b, "subtype_GO_enrichment.png"),
             p_go_23b, width = 12, height = go_h, dpi = 300, bg = "white")
      ggsave(file.path(out_dir_23b, "subtype_GO_enrichment.pdf"),
             p_go_23b, width = 12, height = go_h)
      write.csv(as.data.frame(ck_go_23b),
                file.path(out_dir_23b, "subtype_GO_enrichment.csv"), row.names = FALSE)
      cat(sprintf("[OUTPUT] %s\n", file.path(out_dir_23b, "subtype_GO_enrichment.png")))
      cat(sprintf("[OUTPUT] %s\n", file.path(out_dir_23b, "subtype_GO_enrichment.pdf")))
      cat(sprintf("[OUTPUT] %s\n", file.path(out_dir_23b, "subtype_GO_enrichment.csv")))
      cat(sprintf("[INFO] GO 图高: %.1f in\n", go_h))
    } else {
      cat("[WARN] GO 富集无显著结果\n")
    }


    cat("[INFO] compareCluster KEGG...\n")
    ck_kegg_23b <- tryCatch(
      compareCluster(gene_list_23b, fun = "enrichKEGG",
                     organism = "hsa", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                     use_internal_data = TRUE),
      error = function(e) { cat("[WARN] KEGG 失败:", conditionMessage(e), "\n"); NULL }
    )

    if (!is.null(ck_kegg_23b) && nrow(ck_kegg_23b) > 0) {

      ck_kegg_23b_r <- tryCatch(
        setReadable(ck_kegg_23b, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
        error = function(e) ck_kegg_23b
      )

      kegg_df <- as.data.frame(ck_kegg_23b_r)

      kegg_df$Description <- ifelse(
        grepl("^hsa:", kegg_df$Description),
        kegg_df$ID,
        kegg_df$Description
      )
      kegg_df$Description <- wrap_pathway_name(kegg_df$Description, width = 55)
      ck_kegg_23b_r@compareClusterResult$Description <- kegg_df$Description

      n_kegg_terms <- length(unique(kegg_df$Description))
      kegg_h <- calc_plot_height(n_kegg_terms, base_height = 5,
                                 per_item_inch = 0.32, min_height = 7, max_height = 28)

      p_kegg_23b <- dotplot(ck_kegg_23b_r, showCategory = 5,
                            title = "Fibroblast Subtype — KEGG Pathway Enrichment",
                            font.size = 9) +
        theme(
          plot.title  = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8, lineheight = 0.85),
          plot.margin = ggplot2::margin(t = 10, r = 20, b = 10, l = 160)
        )
      ggsave(file.path(out_dir_23b, "subtype_KEGG_enrichment.png"),
             p_kegg_23b, width = 12, height = kegg_h, dpi = 300, bg = "white")
      ggsave(file.path(out_dir_23b, "subtype_KEGG_enrichment.pdf"),
             p_kegg_23b, width = 12, height = kegg_h)
      write.csv(as.data.frame(ck_kegg_23b),
                file.path(out_dir_23b, "subtype_KEGG_enrichment.csv"), row.names = FALSE)
      cat(sprintf("[OUTPUT] %s\n", file.path(out_dir_23b, "subtype_KEGG_enrichment.png")))
      cat(sprintf("[OUTPUT] %s\n", file.path(out_dir_23b, "subtype_KEGG_enrichment.pdf")))
      cat(sprintf("[OUTPUT] %s\n", file.path(out_dir_23b, "subtype_KEGG_enrichment.csv")))
      cat(sprintf("[INFO] KEGG 图高: %.1f in\n", kegg_h))
    } else {
      cat("[WARN] KEGG 富集无显著结果\n")
    }


    if (HAS_REACTOME_23b) {
      cat("[INFO] compareCluster Reactome...\n")
      suppressPackageStartupMessages(library(ReactomePA))
      ck_reactome_23b <- tryCatch(
        compareCluster(gene_list_23b, fun = "enrichPathway",
                       organism = "human", pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2, readable = TRUE),
        error = function(e) { cat("[WARN] Reactome 失败:", conditionMessage(e), "\n"); NULL }
      )

      if (!is.null(ck_reactome_23b) && nrow(ck_reactome_23b) > 0) {

        rct_df <- as.data.frame(ck_reactome_23b)
        rct_df$Description <- wrap_pathway_name(rct_df$Description, width = 55)
        ck_reactome_23b@compareClusterResult$Description <- rct_df$Description

        n_rct_terms <- length(unique(rct_df$Description))
        rct_h <- calc_plot_height(n_rct_terms, base_height = 5,
                                  per_item_inch = 0.32, min_height = 7, max_height = 28)

        p_reactome_23b <- dotplot(ck_reactome_23b, showCategory = 5,
                                  title = "Fibroblast Subtype — Reactome Pathway Enrichment",
                                  font.size = 9) +
          theme(
            plot.title  = element_text(hjust = 0.5, face = "bold", size = 11),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8, lineheight = 0.85),
            plot.margin = ggplot2::margin(t = 10, r = 20, b = 10, l = 160)
          )
        ggsave(file.path(out_dir_23b, "subtype_Reactome_enrichment.png"),
               p_reactome_23b, width = 12, height = rct_h, dpi = 300, bg = "white")
        ggsave(file.path(out_dir_23b, "subtype_Reactome_enrichment.pdf"),
               p_reactome_23b, width = 12, height = rct_h)
        write.csv(as.data.frame(ck_reactome_23b),
                  file.path(out_dir_23b, "subtype_Reactome_enrichment.csv"), row.names = FALSE)
        cat(sprintf("[OUTPUT] %s\n", file.path(out_dir_23b, "subtype_Reactome_enrichment.png")))
        cat(sprintf("[OUTPUT] %s\n", file.path(out_dir_23b, "subtype_Reactome_enrichment.pdf")))
        cat(sprintf("[OUTPUT] %s\n", file.path(out_dir_23b, "subtype_Reactome_enrichment.csv")))
        cat(sprintf("[INFO] Reactome 图高: %.1f in\n", rct_h))
      } else {
        cat("[WARN] Reactome 富集无显著结果\n")
      }
    } else {
      cat("[INFO] ReactomePA 未安装，跳过 Reactome\n")
    }
  }
}

cat("[DONE] PART 2 完成\n")

required_pkgs_mono <- c("monocle", "Seurat", "ggplot2", "BiocGenerics",
                        "VGAM", "igraph", "DDRTree", "irlba")
skip_mono <- FALSE
for (pkg in required_pkgs_mono) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("[SKIP] 缺少包:", pkg, "\n"); skip_mono <- TRUE; break
  }
}

if (!skip_mono) {
  suppressPackageStartupMessages({
    library(monocle); library(Seurat); library(ggplot2)
  })

  input_rds_primary_mono  <- file.path(SERVER_DIR, "r.03_single_cell/seurat_keycells.rds")
  input_rds_fallback_mono <- file.path(SERVER_DIR, "r.03_single_cell/seurat_annotated.rds")
  output_dir_mono         <- file.path(SERVER_DIR, "14_trajectory")
  dir.create(output_dir_mono, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(input_rds_primary_mono)) {
    cat("[INFO] 读取:", input_rds_primary_mono, "\n")
    sc_m <- readRDS(input_rds_primary_mono)
  } else if (file.exists(input_rds_fallback_mono)) {
    cat("[WARN] 使用 fallback:", input_rds_fallback_mono, "\n")
    sc_m <- readRDS(input_rds_fallback_mono)
  } else {
    cat("[SKIP] 无输入文件，跳过 Monocle2\n")
    sc_m <- NULL
  }

  if (!is.null(sc_m)) {
    cat("[INFO] 原始细胞数:", ncol(sc_m), "\n")

    group_col_m <- if ("celltype_annotation" %in% colnames(sc_m@meta.data)) "celltype_annotation" else
      if ("cell_type"           %in% colnames(sc_m@meta.data)) "cell_type" else NULL

    if (!is.null(group_col_m)) {
      mono_cells <- colnames(sc_m)[sc_m@meta.data[[group_col_m]] %in% KEY_CELL]
      if (length(mono_cells) >= 50) {
        cat(sprintf("[INFO] 筛选 %s 细胞: %d 个\n", KEY_CELL, length(mono_cells)))
        sc_m <- sc_m[, mono_cells]
      } else {
        cat(sprintf("[WARN] %s 细胞数 < 50（%d），使用全部细胞\n",
                    KEY_CELL, length(mono_cells)))
      }
    }

    set.seed(42)
    if (ncol(sc_m) > 20000) {
      sc_m <- sc_m[, sample(colnames(sc_m), 10000)]
      cat("[INFO] 抽样后:", ncol(sc_m), "\n")
    }

    group_col_m2 <- if (!is.null(group_col_m)) group_col_m else "seurat_clusters"

    expr_mat_m <- tryCatch(
      as.matrix(GetAssayData(sc_m, assay = "RNA", layer = "counts")),
      error = function(e) tryCatch(
        as.matrix(GetAssayData(sc_m, assay = "RNA", slot  = "counts")),
        error = function(e2) { cat("[SKIP] 无法获取 counts\n"); NULL }
      )
    )

    if (!is.null(expr_mat_m)) {
      pd_m  <- new("AnnotatedDataFrame", data = sc_m@meta.data)
      fd_m  <- new("AnnotatedDataFrame",
                   data = data.frame(gene_short_name = rownames(expr_mat_m),
                                     row.names = rownames(expr_mat_m)))
      cds_m <- tryCatch(
        newCellDataSet(expr_mat_m, phenoData = pd_m, featureData = fd_m,
                       expressionFamily = negbinomial.size()),
        error = function(e) { cat("[SKIP] CellDataSet 失败\n"); NULL }
      )

      if (!is.null(cds_m)) {
        cat("[INFO] estimateSizeFactors...\n")
        cds_m <- estimateSizeFactors(cds_m)
        cds_m <- tryCatch(estimateDispersions(cds_m),
                          error = function(e) { cat("[WARN] estimateDispersions 跳过\n"); cds_m })

        valid_key_m <- KEY_GENES[KEY_GENES %in% rownames(cds_m)]
        hvg_m <- tryCatch(VariableFeatures(sc_m), error = function(e) character(0))
        hvg_m <- head(hvg_m[hvg_m %in% rownames(cds_m)], 500)
        ordering_m <- unique(c(valid_key_m, hvg_m))
        ordering_m <- head(ordering_m, 500)
        cat(sprintf("[INFO] ordering genes: %d\n", length(ordering_m)))

        if (length(ordering_m) > 0) {
          cds_m <- setOrderingFilter(cds_m, ordering_m)
          cds_m <- tryCatch(
            reduceDimension(cds_m, max_components = 2, method = "DDRTree"),
            error = function(e) { cat("[SKIP] reduceDimension 失败\n"); NULL }
          )

          if (!is.null(cds_m)) {
            cds_m <- tryCatch(orderCells(cds_m),
                              error = function(e) { cat("[SKIP] orderCells 失败\n"); NULL })

            if (!is.null(cds_m)) {
              cat("[INFO] 拟时序完成\n")
              pseudo_r <- range(pData(cds_m)$Pseudotime, na.rm = TRUE)
              cat(sprintf("[STAT] Trajectory_cell_count: %d\n", ncol(cds_m)))
              cat(sprintf("[STAT] Pseudotime_range: %.2f-%.2f\n", pseudo_r[1], pseudo_r[2]))

              write.csv(
                data.frame(cell_id = rownames(pData(cds_m)),
                           Pseudotime = pData(cds_m)$Pseudotime,
                           State = pData(cds_m)$State),
                file.path(output_dir_mono, "monocle2_pseudotime.csv"), row.names = FALSE
              )

              tryCatch({
                p_traj <- if (group_col_m2 %in% colnames(pData(cds_m)))
                  plot_cell_trajectory(cds_m, color_by = group_col_m2) +
                  ggtitle(paste0("Monocle2 Trajectory — ", KEY_CELL))
                else
                  plot_cell_trajectory(cds_m, color_by = "State") + ggtitle("State Trajectory")
                ggsave(file.path(output_dir_mono, "trajectory_umap.png"),
                       p_traj, width = 10, height = 8, dpi = 150)
                cat("[INFO] trajectory_umap.png 已保存\n")
              }, error = function(e) cat("[WARN] trajectory_umap.png 失败\n"))

              tryCatch({
                p_pt <- plot_cell_trajectory(cds_m, color_by = "Pseudotime") +
                  ggtitle("Pseudotime Trajectory") +
                  scale_color_viridis_c(option = "plasma")
                ggsave(file.path(output_dir_mono, "pseudotime_umap.png"),
                       p_pt, width = 10, height = 8, dpi = 150)
                cat("[INFO] pseudotime_umap.png 已保存\n")
              }, error = function(e) cat("[WARN] pseudotime_umap.png 失败\n"))

              tryCatch({
                genes_plot_m <- valid_key_m[valid_key_m %in% rownames(cds_m)]
                if (length(genes_plot_m) > 0) {


                  HARDCODE_ANNO_M <- c(
                    "0" = "HP+ Mac",
                    "1" = "S100A8+ Mono-Mac",
                    "2" = "VSIG4+ Res-Mac",
                    "3" = "MMP19+ Inflam-Mac",
                    "4" = "VMO1+ Mac",
                    "5" = "IFIT+ IFN-Mac",
                    "6" = "TSPAN18+ Mac",
                    "7" = "PTH1R+ Mac",
                    "8" = "MYOC+ Mac"
                  )


                  cluster_col <- NULL
                  for (cand in c("seurat_clusters", "RNA_snn_res.0.8", "cluster")) {
                    if (cand %in% colnames(pData(cds_m))) { cluster_col <- cand; break }
                  }

                  if (!is.null(cluster_col)) {
                    cl_vec <- as.character(pData(cds_m)[[cluster_col]])
                    anno_vec <- ifelse(
                      cl_vec %in% names(HARDCODE_ANNO_M),
                      HARDCODE_ANNO_M[cl_vec],
                      paste0("Mac_", cl_vec)
                    )
                    pData(cds_m)$mac_subtype <- anno_vec
                    color_col <- "mac_subtype"
                    cat("[INFO] 拟时序图使用巨噬细胞亚型注释分色\n")
                  } else {

                    color_col <- "State"
                    cat("[WARN] 找不到簇列，使用 State 分色\n")
                  }

                  p_gn <- plot_genes_in_pseudotime(
                    cds_m[genes_plot_m, ],
                    color_by = color_col,
                    ncol     = min(3, length(genes_plot_m))) +
                    ggtitle(paste0("Key Genes along Pseudotime — ", KEY_CELL, " Subtypes"))

                  ggsave(file.path(output_dir_mono, "key_genes_trajectory.png"), p_gn,
                         width  = 4 * min(3, length(genes_plot_m)),
                         height = 4 * ceiling(length(genes_plot_m) / 3),
                         dpi    = 150)
                  cat("[INFO] key_genes_trajectory.png 已保存\n")
                }
              }, error = function(e) cat("[WARN] key_genes_trajectory.png 失败:", conditionMessage(e), "\n"))

              saveRDS(cds_m, file.path(output_dir_mono, "monocle2_cds.rds"))
              cat("[INFO] monocle2_cds.rds 已保存\n")
            }
          }
        }
      }
    }
  }
}


traj_dir_fix <- file.path(SERVER_DIR, "14_trajectory")
cds_file_fix <- file.path(traj_dir_fix, "monocle2_cds.rds")

if (!file.exists(cds_file_fix)) {
  cat("[SKIP] monocle2_cds.rds 不存在\n")
} else {
  suppressPackageStartupMessages({ library(ggplot2); library(monocle) })
  cat("[INFO] 读取 CDS...\n")
  cds_fix <- readRDS(cds_file_fix)

  reducedDim_fix <- t(reducedDimS(cds_fix))
  pdata_fix <- pData(cds_fix)

  plot_df_fix <- data.frame(
    Component_1 = reducedDim_fix[, 1],
    Component_2 = reducedDim_fix[, 2],
    Pseudotime  = pdata_fix$Pseudotime,
    State       = as.factor(pdata_fix$State),
    CellType    = if ("cell_type"       %in% colnames(pdata_fix)) pdata_fix$cell_type else
      if ("celltype_use"    %in% colnames(pdata_fix)) pdata_fix$celltype_use else
        rep("Unknown", nrow(pdata_fix))
  )

  tree_coords_fix <- tryCatch({
    mst_fix <- minSpanningTree(cds_fix)
    el_fix  <- igraph::get.edgelist(mst_fix)
    dp_fix  <- t(reducedDimK(cds_fix))
    data.frame(x = dp_fix[el_fix[,1], 1], y = dp_fix[el_fix[,1], 2],
               xend = dp_fix[el_fix[,2], 1], yend = dp_fix[el_fix[,2], 2])
  }, error = function(e) { cat("[WARN] 树坐标提取失败\n"); NULL })

  p1_fix <- ggplot(plot_df_fix, aes(Component_1, Component_2)) +
    geom_point(aes(color = State), size = 0.3, alpha = 0.5) +
    scale_color_brewer(palette = "Set1") +
    labs(title = paste0("Monocle2 Trajectory — ", KEY_CELL, " (State)"),
         x = "Component 1", y = "Component 2", color = "State") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  if (!is.null(tree_coords_fix))
    p1_fix <- p1_fix + geom_segment(data = tree_coords_fix,
                                    aes(x, y, xend = xend, yend = yend),
                                    color = "black", linewidth = 1, inherit.aes = FALSE)
  ggsave(file.path(traj_dir_fix, "trajectory_umap.png"), p1_fix,
         width = 8, height = 6, dpi = 150)
  cat("[INFO] trajectory_umap.png (fix) 已保存\n")

  p2_fix <- ggplot(plot_df_fix, aes(Component_1, Component_2)) +
    geom_point(aes(color = Pseudotime), size = 0.3, alpha = 0.5) +
    scale_color_viridis_c(option = "plasma") +
    labs(title = paste0("Monocle2 Pseudotime — ", KEY_CELL),
         x = "Component 1", y = "Component 2", color = "Pseudotime") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  if (!is.null(tree_coords_fix))
    p2_fix <- p2_fix + geom_segment(data = tree_coords_fix,
                                    aes(x, y, xend = xend, yend = yend),
                                    color = "white", linewidth = 1, inherit.aes = FALSE)
  ggsave(file.path(traj_dir_fix, "pseudotime_umap.png"), p2_fix,
         width = 8, height = 6, dpi = 150)
  cat("[INFO] pseudotime_umap.png (fix) 已保存\n")
}


suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

OUT_DIR_RCT <- file.path(SERVER_DIR, "16_reactome")
dir.create(OUT_DIR_RCT, showWarnings = FALSE, recursive = TRUE)

deg_path_rct <- file.path(SERVER_DIR, "r.01_deg/DEG_GSE57338_full.csv")
if (!file.exists(deg_path_rct)) {
  cat("[SKIP] DEG 文件不存在\n")
} else {
  cat("[INFO] 读取 DEG:", deg_path_rct, "\n")
  deg_rct <- read.csv(deg_path_rct, stringsAsFactors = FALSE)
  cat("[INFO] DEG 行数:", nrow(deg_rct), "\n")

  gene_col_rct <- NULL
  for (cn in c("gene","Gene","GENE","gene_name","symbol","Symbol","hgnc_symbol"))
    if (cn %in% colnames(deg_rct)) { gene_col_rct <- cn; break }
  if (is.null(gene_col_rct)) gene_col_rct <- colnames(deg_rct)[1]

  fc_col_rct <- NULL
  for (cn in c("logFC","log2FoldChange","log2FC","LogFC","avg_log2FC","log2_fold_change"))
    if (cn %in% colnames(deg_rct)) { fc_col_rct <- cn; break }
  if (is.null(fc_col_rct)) fc_col_rct <- colnames(deg_rct)[2]

  pval_col_rct <- NULL
  for (cn in c("adj.P.Val","padj","FDR","p_val_adj","P.Value","pvalue","p.value"))
    if (cn %in% colnames(deg_rct)) { pval_col_rct <- cn; break }

  cat(sprintf("[INFO] 基因列: %s | FC列: %s | pval列: %s\n",
              gene_col_rct, fc_col_rct,
              ifelse(is.null(pval_col_rct), "未找到", pval_col_rct)))

  sc_path_rct <- NULL
  for (p in c(file.path(SERVER_DIR, "r.03_single_cell/seurat_keycells.rds"),
              file.path(SERVER_DIR, "r.03_single_cell/seurat_annotated.rds"))) {
    if (file.exists(p)) { sc_path_rct <- p; break }
  }

  sc_obj_rct <- NULL
  if (!is.null(sc_path_rct) && requireNamespace("Seurat", quietly = TRUE)) {
    suppressPackageStartupMessages(library(Seurat))
    tryCatch({
      sc_obj_rct <- readRDS(sc_path_rct)
      cat("[INFO] 单细胞对象:", ncol(sc_obj_rct), "cells\n")
    }, error = function(e) cat("[WARN] 单细胞读取失败\n"))
  }

  deg_clean_rct <- deg_rct[!is.na(deg_rct[[gene_col_rct]]) &
                             deg_rct[[gene_col_rct]] != "" &
                             !is.na(deg_rct[[fc_col_rct]]), ]
  deg_clean_rct <- deg_clean_rct[!duplicated(deg_clean_rct[[gene_col_rct]]), ]
  gene_list_rct <- sort(setNames(as.numeric(deg_clean_rct[[fc_col_rct]]),
                                 deg_clean_rct[[gene_col_rct]]), decreasing = TRUE)

  sig_genes_rct <- if (!is.null(pval_col_rct))
    deg_clean_rct[[gene_col_rct]][!is.na(deg_clean_rct[[pval_col_rct]]) &
                                    deg_clean_rct[[pval_col_rct]] < 0.05]
  else
    deg_clean_rct[[gene_col_rct]][abs(deg_clean_rct[[fc_col_rct]]) > 1]
  cat("[INFO] 显著 DEG:", length(sig_genes_rct), "\n")

  gene_entrez_rct <- NULL; sig_entrez_rct <- NULL
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
    tryCatch({
      emap <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = names(gene_list_rct),
                                    column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
      valid_e <- !is.na(emap)
      gene_entrez_rct <- gene_list_rct[valid_e]
      names(gene_entrez_rct) <- emap[valid_e]
      gene_entrez_rct <- sort(gene_entrez_rct, decreasing = TRUE)
      cat("[INFO] Entrez IDs:", length(gene_entrez_rct), "\n")

      if (length(sig_genes_rct) > 0) {
        semap <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = sig_genes_rct,
                                       column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
        sig_entrez_rct <- semap[!is.na(semap)]
      }
    }, error = function(e) cat("[WARN] Entrez 转换失败\n"))
  }


  reactomegsa_ok_rct <- FALSE
  reactomegsa_res_rct <- NULL
  gsea_res_rct <- NULL
  ora_res_rct  <- NULL
  method_rct   <- "none"
  result_df_rct <- NULL
  save_plot_both <- function(plot_obj, basename, width, height, dpi = 150) {
    png_path <- file.path(OUT_DIR_RCT, paste0(basename, ".png"))
    pdf_path <- file.path(OUT_DIR_RCT, paste0(basename, ".pdf"))
    ggsave(png_path, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
    ggsave(pdf_path, plot_obj, width = width, height = height)
    cat(sprintf("[OUTPUT] %s\n", png_path))
    cat(sprintf("[OUTPUT] %s\n", pdf_path))
  }


  if (requireNamespace("ReactomeGSA", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ReactomeGSA))
    tryCatch({
      bulk_expr_rct <- matrix(deg_clean_rct[[fc_col_rct]], nrow = nrow(deg_clean_rct),
                              ncol = 1, dimnames = list(deg_clean_rct[[gene_col_rct]], "logFC"))
      req_rct <- ReactomeAnalysisRequest(method = "ssGSEA")
      req_rct <- add_dataset(req_rct, expression_values = bulk_expr_rct,
                             name = "Bulk_RNA", type = "rnaseq_counts",
                             comparison_factor = "logFC")
      if (!is.null(sc_obj_rct)) {
        tryCatch({
          DefaultAssay(sc_obj_rct) <- "RNA"
          sc_cnt_rct <- GetAssayData(sc_obj_rct, layer = "counts")
          if (ncol(sc_cnt_rct) > 500)
            sc_cnt_rct <- sc_cnt_rct[, sample(ncol(sc_cnt_rct), 500)]
          req_rct <- add_dataset(req_rct, expression_values = as.matrix(sc_cnt_rct),
                                 name = "scRNA", type = "rnaseq_counts")
        }, error = function(e) cat("[WARN] scRNA 添加失败:", conditionMessage(e), "\n"))
      }
      cat("[INFO] 提交 ReactomeGSA 请求（需外网）...\n")
      reactomegsa_res_rct <- performAnalysis(req_rct)
      reactomegsa_ok_rct  <- TRUE
      method_rct <- "ReactomeGSA"
      cat("[INFO] ReactomeGSA 请求成功\n")
    }, error = function(e) {
      cat("[WARN] ReactomeGSA 失败，完整报错:\n")
      cat("  ", conditionMessage(e), "\n")
      cat("[INFO] 自动切换至 ReactomePA（离线模式）\n")
    })
  } else {
    cat("[INFO] ReactomeGSA 未安装，直接使用 ReactomePA\n")
  }


  if (!reactomegsa_ok_rct && !is.null(gene_entrez_rct) &&
      length(gene_entrez_rct) >= 10 &&
      requireNamespace("ReactomePA", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ReactomePA))
    tryCatch({
      cat("[INFO] 运行 ReactomePA::gsePathway...\n")
      gsea_res_rct <- gsePathway(gene_entrez_rct, organism = "human",
                                 pvalueCutoff = 0.2, pAdjustMethod = "BH",
                                 verbose = FALSE, seed = 42)
      cat("[INFO] GSEA 通路数:", nrow(gsea_res_rct@result), "\n")
      method_rct <- "ReactomePA_GSEA"
    }, error = function(e) cat("[WARN] GSEA 失败:", conditionMessage(e), "\n"))
  }

  if (!reactomegsa_ok_rct && is.null(gsea_res_rct) &&
      !is.null(sig_entrez_rct) && length(sig_entrez_rct) >= 5 &&
      requireNamespace("ReactomePA", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ReactomePA))
    tryCatch({
      cat("[INFO] 运行 ReactomePA::enrichPathway (ORA)...\n")
      ora_res_rct <- enrichPathway(
        gene          = sig_entrez_rct,
        universe      = if (!is.null(gene_entrez_rct)) names(gene_entrez_rct) else NULL,
        organism      = "human", pvalueCutoff  = 0.05,
        pAdjustMethod = "BH",   qvalueCutoff  = 0.2, readable = TRUE)
      cat("[INFO] ORA 通路数:", nrow(ora_res_rct@result), "\n")
      method_rct <- "ReactomePA_ORA"
    }, error = function(e) cat("[WARN] ORA 失败:", conditionMessage(e), "\n"))
  }

  pathway_count_rct <- 0
  if (reactomegsa_ok_rct && !is.null(reactomegsa_res_rct)) {
    tryCatch({
      result_df_rct <- result_summary(reactomegsa_res_rct)
      pathway_count_rct <- nrow(result_df_rct)
    }, error = function(e) cat("[WARN] ReactomeGSA result_summary 失败:", conditionMessage(e), "\n"))
  } else if (!is.null(gsea_res_rct)) {
    result_df_rct <- as.data.frame(gsea_res_rct)
    pathway_count_rct <- nrow(result_df_rct)
  } else if (!is.null(ora_res_rct)) {
    result_df_rct <- as.data.frame(ora_res_rct)
    pathway_count_rct <- nrow(result_df_rct)
  }

  csv_path_rct <- file.path(OUT_DIR_RCT, "reactome_results.csv")
  if (!is.null(result_df_rct) && nrow(result_df_rct) > 0) {
    write.csv(result_df_rct, csv_path_rct, row.names = FALSE)
    cat(sprintf("[OUTPUT] %s\n", csv_path_rct))
  } else {
    write.csv(data.frame(pathway = character(), pvalue = numeric(), padj = numeric()),
              csv_path_rct, row.names = FALSE)
    cat(sprintf("[OUTPUT] %s (空结果占位)\n", csv_path_rct))
  }


  if (reactomegsa_ok_rct && !is.null(reactomegsa_res_rct)) {
    tryCatch({
      png_p <- file.path(OUT_DIR_RCT, "reactome_heatmap.png")
      pdf_p <- file.path(OUT_DIR_RCT, "reactome_heatmap.pdf")
      png(png_p, width = 1400, height = 1000, res = 150)
      plot_overview(reactomegsa_res_rct)
      dev.off()
      pdf(pdf_p, width = 1400 / 150, height = 1000 / 150)
      plot_overview(reactomegsa_res_rct)
      dev.off()
      cat(sprintf("[OUTPUT] %s\n", png_p))
      cat(sprintf("[OUTPUT] %s\n", pdf_p))
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat("[WARN] ReactomeGSA 热图失败:", conditionMessage(e), "\n")
    })
  }


  if (!is.null(gsea_res_rct)) {
    tryCatch({
      plot_df_g <- as.data.frame(gsea_res_rct)
      plot_df_g <- head(plot_df_g[order(plot_df_g$p.adjust), ], 20)
      if (nrow(plot_df_g) > 0) {
        plot_df_g$Description <- wrap_pathway_name(plot_df_g$Description, 60)
        plot_df_g$Description <- factor(plot_df_g$Description,
                                        levels = rev(plot_df_g$Description))
        nes_col_g <- if ("NES" %in% colnames(plot_df_g)) "NES" else colnames(plot_df_g)[3]
        rct_h_g   <- calc_plot_height(nrow(plot_df_g), per_item_inch = 0.38,
                                      min_height = 7, max_height = 24)
        p_g <- ggplot(plot_df_g,
                      aes(x = .data[[nes_col_g]], y = Description,
                          color = p.adjust, size = setSize)) +
          geom_point() +
          scale_color_gradient(low = "red", high = "blue", name = "Adj. P-value") +
          scale_size_continuous(name = "Gene Set Size") +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
          labs(title    = paste0("Reactome Pathway GSEA (", method_rct, ")"),
               subtitle = paste0("Top ", nrow(plot_df_g), " pathways | pAdj < 0.2"),
               x = "Normalized Enrichment Score (NES)", y = NULL) +
          theme_bw(base_size = 11) +
          theme(plot.title    = element_text(face = "bold", size = 13),
                plot.subtitle = element_text(size = 9, color = "grey40"),
                axis.text.y   = element_text(size = 8, lineheight = 0.85),
                plot.margin   = ggplot2::margin(t = 10, r = 20, b = 10, l = 160))
        save_plot_both(p_g, "reactome_GSEA_dotplot", width = 13, height = rct_h_g)
      }
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat("[WARN] GSEA 绘图失败:", conditionMessage(e), "\n")
    })
  }


  if (!is.null(ora_res_rct)) {
    tryCatch({
      plot_df_o <- as.data.frame(ora_res_rct)
      plot_df_o <- head(plot_df_o[order(plot_df_o$p.adjust), ], 20)
      if (nrow(plot_df_o) > 0) {
        plot_df_o$Description    <- wrap_pathway_name(plot_df_o$Description, 60)
        plot_df_o$Description    <- factor(plot_df_o$Description,
                                           levels = rev(plot_df_o$Description))
        plot_df_o$neg_log10_padj <- -log10(plot_df_o$p.adjust + 1e-300)
        gr_parts_o               <- strsplit(as.character(plot_df_o$GeneRatio), "/")
        plot_df_o$GeneRatioNum   <- sapply(gr_parts_o, function(x)
          if (length(x) == 2) as.numeric(x[1]) / as.numeric(x[2]) else NA)
        rct_h_o <- calc_plot_height(nrow(plot_df_o), per_item_inch = 0.38,
                                    min_height = 7, max_height = 24)
        p_o <- ggplot(plot_df_o, aes(GeneRatioNum, Description,
                                     color = neg_log10_padj, size = Count)) +
          geom_point() +
          scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
          scale_size_continuous(name = "Gene Count") +
          labs(title    = paste0("Reactome Pathway ORA (", method_rct, ")"),
               subtitle = paste0("Top ", nrow(plot_df_o), " pathways | pAdj < 0.05"),
               x = "Gene Ratio", y = NULL) +
          theme_bw(base_size = 11) +
          theme(plot.title    = element_text(face = "bold", size = 13),
                plot.subtitle = element_text(size = 9, color = "grey40"),
                axis.text.y   = element_text(size = 8, lineheight = 0.85),
                plot.margin   = ggplot2::margin(t = 10, r = 20, b = 10, l = 160))
        save_plot_both(p_o, "reactome_ORA_dotplot", width = 13, height = rct_h_o)
      }
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat("[WARN] ORA 绘图失败:", conditionMessage(e), "\n")
    })
  }


  if (is.null(gsea_res_rct) && is.null(ora_res_rct) && !reactomegsa_ok_rct) {
    tryCatch({
      p_empty_rct <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, size = 6, color = "grey50",
                 label = paste0("Reactome 通路分析未获得显著结果\n方法: ", method_rct)) +
        xlim(0, 1) + ylim(0, 1) +
        labs(title = "Reactome Pathway Analysis") + theme_void() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      save_plot_both(p_empty_rct, "reactome_heatmap", width = 10, height = 6)
    }, error = function(e) { if (dev.cur() > 1) dev.off() })
  }

  cat("[STAT] Reactome_pathway_count:", pathway_count_rct, "\n")
  cat("[INFO] 分析方法:", method_rct, "\n")
}
