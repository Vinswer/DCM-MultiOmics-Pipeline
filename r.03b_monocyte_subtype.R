rm(list = ls()); gc()

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

wrap_pathway_name <- function(x, width = 45) {

  sapply(x, function(nm) {
    if (nchar(nm) <= width) return(nm)
    wrapped <- strwrap(nm, width = width)
    paste(wrapped, collapse = "\n")
  }, USE.NAMES = FALSE)
}

resolve_kegg_names <- function(df) {
  need_fix <- grepl("^hsa[0-9]{5}$|^hsa:[0-9]+$", df$Description)
  if (!any(need_fix)) return(df)

  fixed_desc     <- df$Description
  kegg_ids_raw   <- df$ID[need_fix]
  kegg_ids_clean <- gsub("^hsa:", "", kegg_ids_raw)
  kegg_ids_clean <- ifelse(grepl("^hsa", kegg_ids_clean),
                           kegg_ids_clean, paste0("hsa", kegg_ids_clean))

  resolved <- rep(NA_character_, length(kegg_ids_clean))
  if (requireNamespace("KEGGREST", quietly = TRUE)) {
    tryCatch({
      batch_size <- 10
      for (i in seq(1, length(kegg_ids_clean), by = batch_size)) {
        idx   <- i:min(i + batch_size - 1, length(kegg_ids_clean))
        info  <- KEGGREST::keggGet(kegg_ids_clean[idx])
        for (j in seq_along(info)) {
          nm <- info[[j]]$NAME
          if (!is.null(nm) && length(nm) > 0) {
            nm_clean <- gsub(" - Homo sapiens.*$", "", nm[1])
            nm_clean <- gsub(" - .*$", "", nm_clean)
            resolved[idx[j]] <- nm_clean
          }
        }
        Sys.sleep(0.3)
      }
    }, error = function(e) cat("[WARN] KEGGREST 查询失败:", conditionMessage(e), "\n"))
  }

  for (i in seq_along(kegg_ids_raw)) {
    orig_row <- which(need_fix)[i]
    fixed_desc[orig_row] <- if (!is.na(resolved[i]) && nchar(resolved[i]) > 0)
      resolved[i] else kegg_ids_clean[i]
  }
  df$Description <- fixed_desc
  df
}

plot_top10_global <- function(result_df,
                              db_name,
                              total_count,
                              out_path_png,
                              out_path_pdf = NULL,
                              top_n        = 10,
                              label_width  = 40) {

  library(ggplot2); library(dplyr)

  if (!"Cluster" %in% colnames(result_df)) {
    cat(sprintf("[WARN] %s top%d: 缺少 Cluster 列，跳过\n", db_name, top_n))
    return(invisible(NULL))
  }

  top_pathways <- result_df %>%
    group_by(Description) %>%
    summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(min_padj) %>%
    slice_head(n = top_n) %>%
    pull(Description)

  if (length(top_pathways) == 0) {
    cat(sprintf("[WARN] %s top%d: 无数据\n", db_name, top_n)); return(invisible(NULL))
  }


  plot_df <- result_df %>% filter(Description %in% top_pathways)

  if ("GeneRatio" %in% colnames(plot_df) && is.character(plot_df$GeneRatio)) {
    parts <- strsplit(as.character(plot_df$GeneRatio), "/")
    plot_df$GeneRatioNum <- sapply(parts, function(x)
      if (length(x) == 2) as.numeric(x[1]) / as.numeric(x[2]) else NA_real_)
    x_var <- "GeneRatioNum"; x_lab <- "Gene Ratio"
  } else if ("NES" %in% colnames(plot_df)) {
    x_var <- "NES"; x_lab <- "NES"
  } else {
    x_var <- "GeneRatioNum"
    plot_df$GeneRatioNum <- NA_real_; x_lab <- "Gene Ratio"
  }

  size_var <- if ("Count"   %in% colnames(plot_df)) "Count" else
    if ("setSize" %in% colnames(plot_df)) "setSize" else NULL


  pathway_order <- result_df %>%
    filter(Description %in% top_pathways) %>%
    group_by(Description) %>%
    summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(min_padj) %>%
    pull(Description)

  plot_df$Description <- factor(plot_df$Description, levels = rev(pathway_order))


  levels(plot_df$Description) <- wrap_pathway_name(levels(plot_df$Description), label_width)


  n_rows  <- length(unique(plot_df$Description))
  plot_h  <- max(6, n_rows * 0.60 + 3)

  n_clust <- length(unique(plot_df$Cluster))
  plot_w  <- max(10, n_clust * 1.8 + 5)


  main_title <- sprintf("Subtype %s Enrichment \u2014 Top %d Pathways (Total: %d)",
                        db_name, top_n, total_count)

  p <- ggplot(plot_df,
              aes_string(x = "Cluster", y = "Description",
                         color = "p.adjust",
                         size  = if (!is.null(size_var)) size_var else "p.adjust")) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue",
                         name = "Adj. P-value") +
    { if (!is.null(size_var) && size_var == "Count")
      scale_size_continuous(name = "Gene Count",   range = c(2, 12))
      else if (!is.null(size_var) && size_var == "setSize")
        scale_size_continuous(name = "Set Size",     range = c(2, 12))
      else
        scale_size_continuous(name = "Adj. P-value", range = c(2, 8))
    } +
    labs(title = main_title, x = "Cluster", y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y     = element_text(size = 10, lineheight = 0.9),
      axis.title.x    = element_text(size = 11),
      panel.grid.major = element_line(color = "grey90"),
      legend.position = "right",
      plot.margin     = ggplot2::margin(t = 15, r = 20, b = 15, l = 10)
    )

  ggsave(out_path_png, p,
         width = plot_w, height = plot_h, dpi = 300, bg = "white")
  cat(sprintf("[INFO] %s top%d 图已保存: %s（%.0fw × %.1fh in）\n",
              db_name, top_n, basename(out_path_png), plot_w, plot_h))

  if (!is.null(out_path_pdf))
    ggsave(out_path_pdf, p, width = plot_w, height = plot_h)

  invisible(p)
}

wrap_short <- function(x, width = 55) {
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3), "..."), x)
}


calc_full_height <- function(n_items, min_height = 8, max_height = 48) {
  max(min_height, min(max_height, 4 + n_items * 0.28))
}


required_pkgs_ucell <- c("UCell", "Seurat", "ggplot2", "dplyr", "tidyr")
skip_ucell <- FALSE
for (pkg in required_pkgs_ucell) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("[SKIP] 缺少必要包:", pkg, "\n"); skip_ucell <- TRUE; break
  }
}

if (!skip_ucell) {
  suppressPackageStartupMessages({
    library(UCell); library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  })
  cat("[INFO] Seurat:", as.character(packageVersion("Seurat")),
      "| UCell:", as.character(packageVersion("UCell")), "\n")

  input_rds_ucell  <- file.path(SERVER_DIR, "r.03_single_cell/seurat_annotated.rds")
  key_genes_file   <- file.path(SERVER_DIR, "r.02_ml/key_genes.txt")
  output_dir_ucell <- file.path(SERVER_DIR, "20_ssgsea")
  dir.create(output_dir_ucell, recursive = TRUE, showWarnings = FALSE)

  key_genes_ucell <- if (file.exists(key_genes_file)) {
    kg <- readLines(key_genes_file); kg[nchar(trimws(kg)) > 0]
  } else { cat("[WARN] key_genes.txt 不存在\n"); KEY_GENES }

  if (!file.exists(input_rds_ucell)) {
    cat("[SKIP] 输入文件不存在\n")
  } else {
    sc_u <- readRDS(input_rds_ucell)
    tryCatch({ if (inherits(sc_u[["RNA"]], "Assay5")) sc_u[["RNA"]] <- JoinLayers(sc_u[["RNA"]]) },
             error = function(e) NULL)

    meta_u      <- sc_u@meta.data
    group_col_u <- if ("celltype_annotation" %in% colnames(meta_u)) "celltype_annotation" else
      if ("cell_type"           %in% colnames(meta_u)) "cell_type" else "seurat_clusters"

    key_genes_exist_u <- intersect(key_genes_ucell, rownames(sc_u))
    missing_u         <- setdiff(key_genes_ucell, rownames(sc_u))
    if (length(missing_u) > 0)
      cat("[WARN] 缺失基因:", paste(missing_u, collapse = ", "), "\n")

    if (length(key_genes_exist_u) > 0) {
      counts_mat_u <- tryCatch(GetAssayData(sc_u, assay="RNA", layer="counts"),
                               error = function(e) tryCatch(GetAssayData(sc_u, assay="RNA", slot="counts"),
                                                            error = function(e2) NULL))

      if (!is.null(counts_mat_u)) {
        gene_sets_u <- list(HF_biomarkers = key_genes_exist_u)
        bpparam_u <- tryCatch(
          if (requireNamespace("BiocParallel", quietly=TRUE))
            if (n_cores > 1) BiocParallel::MulticoreParam(n_cores)
          else BiocParallel::SerialParam()
          else NULL, error = function(e) NULL)

        sc_u <- tryCatch(
          if (!is.null(bpparam_u))
            AddModuleScore_UCell(sc_u, features=gene_sets_u, assay="RNA", BPPARAM=bpparam_u)
          else AddModuleScore_UCell(sc_u, features=gene_sets_u, assay="RNA"),
          error = function(e) sc_u)

        ucell_col_u <- "HF_biomarkers_UCell"
        if (ucell_col_u %in% colnames(sc_u@meta.data)) {
          ucell_vec_u <- sc_u@meta.data[[ucell_col_u]]
          cat(sprintf("[STAT] UCell_mean_score: %.3f\n", mean(ucell_vec_u, na.rm=TRUE)))

          score_df_u <- data.frame(cell=rownames(sc_u@meta.data),
                                   cell_type=sc_u@meta.data[[group_col_u]],
                                   UCell=ucell_vec_u, stringsAsFactors=FALSE)
          celltype_scores_u <- score_df_u %>%
            group_by(cell_type) %>%
            summarise(mean_UCell=mean(UCell,na.rm=TRUE),
                      median_UCell=median(UCell,na.rm=TRUE),
                      n_cells=n(), n_high=sum(UCell>0.1,na.rm=TRUE), .groups="drop") %>%
            arrange(desc(mean_UCell))
          write.csv(celltype_scores_u,
                    file.path(output_dir_ucell,"ssgsea_celltype_scores.csv"), row.names=FALSE)


          if (!"umap" %in% names(sc_u@reductions)) {
            tryCatch({
              if (!"pca" %in% names(sc_u@reductions)) {
                sc_u <- NormalizeData(sc_u, verbose=FALSE)
                sc_u <- FindVariableFeatures(sc_u, nfeatures=2000, verbose=FALSE)
                sc_u <- ScaleData(sc_u, verbose=FALSE)
                sc_u <- RunPCA(sc_u, npcs=30, verbose=FALSE)
              }
              sc_u <- RunUMAP(sc_u, reduction="pca", dims=1:20,
                              n.components=2L, seed.use=42, verbose=FALSE)
            }, error=function(e) NULL)
          }

          tryCatch({
            umap_emb <- as.data.frame(sc_u@reductions[["umap"]]@cell.embeddings)
            colnames(umap_emb) <- c("UMAP_1","UMAP_2")
            umap_emb[[ucell_col_u]] <- sc_u@meta.data[rownames(umap_emb), ucell_col_u]
            umap_emb <- umap_emb[order(umap_emb[[ucell_col_u]]),]
            p_umap_u <- ggplot(umap_emb, aes(UMAP_1, UMAP_2, color=.data[[ucell_col_u]])) +
              geom_point(size=0.4, alpha=0.85) +
              scale_color_gradientn(colours=c("lightgrey","#FEE08B","#F46D43","#A50026"),
                                    name="UCell\nScore") +
              ggtitle(paste0("UCell Score — HF Biomarkers (",
                             paste(key_genes_exist_u, collapse=" / "),")")) +
              theme_bw() + theme(plot.title=element_text(hjust=0.5,size=12),
                                 panel.grid=element_blank())
            ggsave(file.path(output_dir_ucell,"ssgsea_umap.png"), p_umap_u, width=9, height=7, dpi=150)
            ggsave(file.path(output_dir_ucell,"ssgsea_umap.pdf"), p_umap_u, width=9, height=7)
          }, error=function(e) NULL)

          tryCatch({
            ct_order_u <- celltype_scores_u %>% arrange(desc(mean_UCell)) %>% pull(cell_type)
            score_df_u$cell_type <- factor(score_df_u$cell_type, levels=ct_order_u)
            n_ct_u <- length(ct_order_u)
            ct_colors_u <- scales::hue_pal()(n_ct_u); names(ct_colors_u) <- ct_order_u
            p_vln_u <- ggplot(score_df_u, aes(cell_type, UCell, fill=cell_type)) +
              geom_violin(scale="width", trim=TRUE, linewidth=0.3) +
              geom_boxplot(width=0.1, fill="white", outlier.size=0.2,
                           outlier.alpha=0.3, linewidth=0.3) +
              geom_hline(yintercept=0.1, linetype="dashed", color="red", linewidth=0.5) +
              scale_fill_manual(values=ct_colors_u) +
              ggtitle("UCell Score Distribution by Cell Type") + ylab("UCell Score") +
              theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
                                 axis.title.x=element_blank(),
                                 plot.title=element_text(hjust=0.5,size=12),
                                 legend.position="none")
            vln_w <- max(10, n_ct_u * 0.9)
            ggsave(file.path(output_dir_ucell,"ssgsea_violin.png"), p_vln_u, width=vln_w, height=7, dpi=150)
            ggsave(file.path(output_dir_ucell,"ssgsea_violin.pdf"), p_vln_u, width=vln_w, height=7)
          }, error=function(e) NULL)
        }
      }
    }
  }
}
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
  library(AnnotationDbi);   library(ggplot2); library(dplyr)
})
HAS_REACTOME <- requireNamespace("ReactomePA", quietly=TRUE)
cat("[INFO] ReactomePA 可用:", HAS_REACTOME, "\n")

out_dir      <- file.path(SERVER_DIR, "r.03_single_cell/monocyte_subtype")
markers_file <- file.path(out_dir, "subtype_markers.csv")

if (!file.exists(markers_file)) {
  cat("[SKIP] subtype_markers.csv 不存在\n")
} else {
  markers   <- read.csv(markers_file, stringsAsFactors=FALSE)
  sig_mk    <- markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
    arrange(cluster, p_val_adj)
  clusters  <- sort(unique(sig_mk$cluster))

  gene_list <- lapply(clusters, function(cl) {
    genes <- unique(sig_mk %>% filter(cluster==cl) %>% pull(gene))
    genes <- genes[nchar(genes) > 0]
    if (length(genes) < 5) return(NULL)
    mapped <- suppressMessages(bitr(genes, fromType="SYMBOL",
                                    toType="ENTREZID", OrgDb=org.Hs.eg.db))
    if (nrow(mapped)==0) return(NULL)
    mapped$ENTREZID
  })

  HARDCODE_ANNO <- c(
    "0" = "S100A8+ Mono-Mac",
    "1" = "S100A8+ Mono-Mac-2",
    "2" = "C1Q+ Res-Mac",
    "3" = "IL1B+ Inflam-Mac",
    "4" = "ISG15+ IFN-Mac",
    "5" = "IL1B+ Inflam-Mac-2",
    "6" = "C1Q+ Res-Mac-2",
    "7" = "SPP1+ Mac",
    "8" = "C1Q+ Res-Mac-3"
  )

  clusters_chr   <- as.character(clusters)
  cluster_labels <- sapply(clusters_chr, function(cl) {
    anno <- HARDCODE_ANNO[cl]
    if (!is.na(anno)) anno else paste0("Subtype_", cl)
  })
  names(cluster_labels) <- clusters_chr

  cat("[INFO] 最终亚型标签：\n")
  for (i in seq_along(cluster_labels))
    cat(sprintf("  %s → %s\n", names(cluster_labels)[i], cluster_labels[i]))
  names(gene_list) <- cluster_labels[clusters_chr]
  gene_list <- Filter(Negate(is.null), gene_list)

  if (length(gene_list) > 0) {

    ck_go <- tryCatch(
      compareCluster(gene_list, fun="enrichGO", OrgDb=org.Hs.eg.db, ont="BP",
                     pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
                     readable=TRUE),
      error=function(e) { cat("[WARN] GO:", conditionMessage(e),"\n"); NULL })

    if (!is.null(ck_go) && nrow(ck_go) > 0) {
      go_df     <- as.data.frame(ck_go)
      total_go  <- length(unique(go_df$Description))
      cat(sprintf("[INFO] GO 总通路: %d\n", total_go))


      go_df$Description <- wrap_short(go_df$Description, 55)
      ck_go@compareClusterResult$Description <- go_df$Description
      n_go  <- length(unique(go_df$Description))
      go_h  <- calc_full_height(n_go, min_height=8, max_height=48)
      p_go  <- dotplot(ck_go, showCategory=5,
                       title=sprintf("Subtype GO BP Enrichment (Total: %d pathways)", total_go),
                       font.size=9) +
        theme(plot.title=element_text(hjust=0.5,face="bold",size=11),
              axis.text.x=element_text(angle=45,hjust=1),
              axis.text.y=element_text(size=8,lineheight=0.85),
              plot.margin=ggplot2::margin(10,20,10,160))
      ggsave(file.path(out_dir,"subtype_GO_enrichment.png"), p_go,
             width=12, height=go_h, dpi=300, bg="white", limitsize=FALSE)
      ggsave(file.path(out_dir,"subtype_GO_enrichment.pdf"), p_go,
             width=12, height=go_h, limitsize=FALSE)
      write.csv(go_df, file.path(out_dir,"subtype_GO_enrichment.csv"), row.names=FALSE)
      cat(sprintf("[INFO] subtype_GO_enrichment.png 已保存 (%.1f in)\n", go_h))


      plot_top10_global(
        result_df    = go_df,
        db_name      = "GO BP",
        total_count  = total_go,
        out_path_png = file.path(out_dir, "subtype_GO_top10.png"),
        out_path_pdf = file.path(out_dir, "subtype_GO_top10.pdf"),
        top_n        = 10,
        label_width  = 40
      )
    } else cat("[WARN] GO 无显著结果\n")


    ck_kegg <- tryCatch(
      compareCluster(gene_list, fun="enrichKEGG", organism="hsa",
                     pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
                     use_internal_data=FALSE),
      error=function(e) {
        cat("[WARN] KEGG online 失败，用内置数据\n")
        tryCatch(compareCluster(gene_list, fun="enrichKEGG", organism="hsa",
                                pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
                                use_internal_data=TRUE),
                 error=function(e2) NULL)
      })

    if (!is.null(ck_kegg) && nrow(ck_kegg) > 0) {

      ck_kegg_r <- tryCatch(setReadable(ck_kegg, OrgDb=org.Hs.eg.db, keyType="ENTREZID"),
                            error=function(e) ck_kegg)
      kegg_df <- as.data.frame(ck_kegg_r)

      is_id_like <- grepl("^hsa[0-9]{5}$|^hsa:[0-9]+$", kegg_df$Description)
      if (any(is_id_like)) {
        cat(sprintf("[INFO] 修复 %d 条 KEGG ID 描述...\n", sum(is_id_like)))
        kegg_df <- resolve_kegg_names(kegg_df)
      }
      total_kegg <- length(unique(kegg_df$Description))
      cat(sprintf("[INFO] KEGG 总通路: %d\n", total_kegg))

      kegg_df$Description <- wrap_short(kegg_df$Description, 55)
      ck_kegg_r@compareClusterResult$Description <- kegg_df$Description
      n_kegg  <- length(unique(kegg_df$Description))
      kegg_h  <- calc_full_height(n_kegg, min_height=7, max_height=48)
      p_kegg  <- dotplot(ck_kegg_r, showCategory=5,
                         title=sprintf("Subtype KEGG Enrichment (Total: %d pathways)", total_kegg),
                         font.size=9) +
        theme(plot.title=element_text(hjust=0.5,face="bold",size=11),
              axis.text.x=element_text(angle=45,hjust=1),
              axis.text.y=element_text(size=8,lineheight=0.85),
              plot.margin=ggplot2::margin(10,20,10,160))
      ggsave(file.path(out_dir,"subtype_KEGG_enrichment.png"), p_kegg,
             width=12, height=kegg_h, dpi=300, bg="white", limitsize=FALSE)
      ggsave(file.path(out_dir,"subtype_KEGG_enrichment.pdf"), p_kegg,
             width=12, height=kegg_h, limitsize=FALSE)
      write.csv(kegg_df, file.path(out_dir,"subtype_KEGG_enrichment.csv"), row.names=FALSE)
      cat(sprintf("[INFO] subtype_KEGG_enrichment.png 已保存 (%.1f in)\n", kegg_h))


      plot_top10_global(
        result_df    = kegg_df,
        db_name      = "KEGG",
        total_count  = total_kegg,
        out_path_png = file.path(out_dir, "subtype_KEGG_top10.png"),
        out_path_pdf = file.path(out_dir, "subtype_KEGG_top10.pdf"),
        top_n        = 10,
        label_width  = 40
      )
    } else cat("[WARN] KEGG 无显著结果\n")

    if (HAS_REACTOME) {
      cat("[INFO] Reactome enrichment...\n")
      suppressPackageStartupMessages(library(ReactomePA))
      ck_rct <- tryCatch(
        compareCluster(gene_list, fun="enrichPathway", organism="human",
                       pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE),
        error=function(e) { cat("[WARN] Reactome:", conditionMessage(e),"\n"); NULL })

      if (!is.null(ck_rct) && nrow(ck_rct) > 0) {
        rct_df     <- as.data.frame(ck_rct)
        total_rct  <- length(unique(rct_df$Description))
        cat(sprintf("[INFO] Reactome 总通路: %d\n", total_rct))


        rct_df$Description <- wrap_short(rct_df$Description, 55)
        ck_rct@compareClusterResult$Description <- rct_df$Description
        n_rct  <- length(unique(rct_df$Description))
        rct_h  <- calc_full_height(n_rct, min_height=7, max_height=48)
        p_rct  <- dotplot(ck_rct, showCategory=5,
                          title=sprintf("Subtype Reactome Enrichment (Total: %d pathways)", total_rct),
                          font.size=9) +
          theme(plot.title=element_text(hjust=0.5,face="bold",size=11),
                axis.text.x=element_text(angle=45,hjust=1),
                axis.text.y=element_text(size=8,lineheight=0.85),
                plot.margin=ggplot2::margin(10,20,10,160))
        ggsave(file.path(out_dir,"subtype_Reactome_enrichment.png"), p_rct,
               width=12, height=rct_h, dpi=300, bg="white", limitsize=FALSE)
        ggsave(file.path(out_dir,"subtype_Reactome_enrichment.pdf"), p_rct,
               width=12, height=rct_h, limitsize=FALSE)
        write.csv(rct_df, file.path(out_dir,"subtype_Reactome_enrichment.csv"), row.names=FALSE)
        cat(sprintf("[INFO] subtype_Reactome_enrichment.png 已保存 (%.1f in)\n", rct_h))


        plot_top10_global(
          result_df    = rct_df,
          db_name      = "Reactome",
          total_count  = total_rct,
          out_path_png = file.path(out_dir, "subtype_Reactome_top10.png"),
          out_path_pdf = file.path(out_dir, "subtype_Reactome_top10.pdf"),
          top_n        = 10,
          label_width  = 40
        )
      } else cat("[WARN] Reactome 无显著结果\n")
    } else {
      cat("[INFO] ReactomePA 未安装，跳过\n")
    }
  }
}


required_pkgs_mono <- c("monocle","Seurat","ggplot2","BiocGenerics",
                        "VGAM","igraph","DDRTree","irlba")
skip_mono <- FALSE
for (pkg in required_pkgs_mono) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat("[SKIP] 缺少包:", pkg, "\n"); skip_mono <- TRUE; break
  }
}

if (!skip_mono) {
  suppressPackageStartupMessages({ library(monocle); library(Seurat); library(ggplot2) })

  input_primary  <- file.path(SERVER_DIR, "r.03_single_cell/seurat_keycells.rds")
  input_fallback <- file.path(SERVER_DIR, "r.03_single_cell/seurat_annotated.rds")
  output_dir_mono <- file.path(SERVER_DIR, "14_trajectory")
  dir.create(output_dir_mono, recursive=TRUE, showWarnings=FALSE)

  sc_m <- if (file.exists(input_primary)) {
    cat("[INFO] 读取:", input_primary, "\n"); readRDS(input_primary)
  } else if (file.exists(input_fallback)) {
    cat("[WARN] 使用 fallback\n"); readRDS(input_fallback)
  } else { cat("[SKIP] 无输入文件\n"); NULL }

  if (!is.null(sc_m)) {
    group_col_m <- if ("celltype_annotation" %in% colnames(sc_m@meta.data)) "celltype_annotation" else
      if ("cell_type"           %in% colnames(sc_m@meta.data)) "cell_type" else NULL

    if (!is.null(group_col_m)) {
      mono_cells <- colnames(sc_m)[sc_m@meta.data[[group_col_m]] %in% KEY_CELL]
      if (length(mono_cells) >= 50) sc_m <- sc_m[, mono_cells]
    }
    set.seed(42)
    if (ncol(sc_m) > 20000) sc_m <- sc_m[, sample(colnames(sc_m), 10000)]
    group_col_m2 <- if (!is.null(group_col_m)) group_col_m else "seurat_clusters"

    expr_mat_m <- tryCatch(as.matrix(GetAssayData(sc_m, assay="RNA", layer="counts")),
                           error=function(e) tryCatch(as.matrix(GetAssayData(sc_m, assay="RNA", slot="counts")),
                                                      error=function(e2) NULL))

    if (!is.null(expr_mat_m)) {
      pd_m  <- new("AnnotatedDataFrame", data=sc_m@meta.data)
      fd_m  <- new("AnnotatedDataFrame",
                   data=data.frame(gene_short_name=rownames(expr_mat_m),
                                   row.names=rownames(expr_mat_m)))
      cds_m <- tryCatch(
        newCellDataSet(expr_mat_m, phenoData=pd_m, featureData=fd_m,
                       expressionFamily=negbinomial.size()),
        error=function(e) NULL)

      if (!is.null(cds_m)) {
        cds_m <- estimateSizeFactors(cds_m)
        cds_m <- tryCatch(estimateDispersions(cds_m), error=function(e) cds_m)

        valid_key_m <- KEY_GENES[KEY_GENES %in% rownames(cds_m)]
        hvg_m       <- head(tryCatch(VariableFeatures(sc_m), error=function(e) character(0)), 500)
        hvg_m       <- hvg_m[hvg_m %in% rownames(cds_m)]
        ordering_m  <- head(unique(c(valid_key_m, hvg_m)), 500)

        if (length(ordering_m) > 0) {
          cds_m <- setOrderingFilter(cds_m, ordering_m)
          cds_m <- tryCatch(reduceDimension(cds_m, max_components=2, method="DDRTree"),
                            error=function(e) NULL)
          if (!is.null(cds_m)) {
            cds_m <- tryCatch(orderCells(cds_m), error=function(e) NULL)
            if (!is.null(cds_m)) {
              pseudo_r <- range(pData(cds_m)$Pseudotime, na.rm=TRUE)
              cat(sprintf("[STAT] Trajectory_cells: %d | Pseudotime: %.2f-%.2f\n",
                          ncol(cds_m), pseudo_r[1], pseudo_r[2]))
              write.csv(data.frame(cell_id=rownames(pData(cds_m)),
                                   Pseudotime=pData(cds_m)$Pseudotime,
                                   State=pData(cds_m)$State),
                        file.path(output_dir_mono,"monocle2_pseudotime.csv"), row.names=FALSE)
              tryCatch({
                p_t <- if (group_col_m2 %in% colnames(pData(cds_m)))
                  plot_cell_trajectory(cds_m, color_by=group_col_m2) +
                  ggtitle(paste0("Monocle2 — ",KEY_CELL))
                else plot_cell_trajectory(cds_m, color_by="State")
                ggsave(file.path(output_dir_mono,"trajectory_umap.png"), p_t, width=10, height=8, dpi=150)
              }, error=function(e) NULL)
              tryCatch({
                p_pt <- plot_cell_trajectory(cds_m, color_by="Pseudotime") +
                  scale_color_viridis_c(option="plasma")
                ggsave(file.path(output_dir_mono,"pseudotime_umap.png"), p_pt, width=10, height=8, dpi=150)
              }, error=function(e) NULL)
              saveRDS(cds_m, file.path(output_dir_mono,"monocle2_cds.rds"))
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
  cds_fix <- readRDS(cds_file_fix)
  rd_fix  <- t(reducedDimS(cds_fix))
  pd_fix  <- pData(cds_fix)
  plot_df_fix <- data.frame(
    Component_1=rd_fix[,1], Component_2=rd_fix[,2],
    Pseudotime=pd_fix$Pseudotime, State=as.factor(pd_fix$State),
    CellType=if ("cell_type" %in% colnames(pd_fix)) pd_fix$cell_type else
      rep("Unknown", nrow(pd_fix)))

  tree_coords <- tryCatch({
    mst <- minSpanningTree(cds_fix); el <- igraph::get.edgelist(mst)
    dp  <- t(reducedDimK(cds_fix))
    data.frame(x=dp[el[,1],1], y=dp[el[,1],2], xend=dp[el[,2],1], yend=dp[el[,2],2])
  }, error=function(e) NULL)

  p1 <- ggplot(plot_df_fix, aes(Component_1, Component_2)) +
    geom_point(aes(color=State), size=0.3, alpha=0.5) +
    scale_color_brewer(palette="Set1") +
    labs(title=paste0("Monocle2 — ",KEY_CELL," (State)"),
         x="Component 1", y="Component 2") +
    theme_bw(12) + theme(plot.title=element_text(hjust=0.5,face="bold"))
  if (!is.null(tree_coords))
    p1 <- p1 + geom_segment(data=tree_coords, aes(x,y,xend=xend,yend=yend),
                            color="black", linewidth=1, inherit.aes=FALSE)
  ggsave(file.path(traj_dir_fix,"trajectory_umap.png"), p1, width=8, height=6, dpi=150)

  p2 <- ggplot(plot_df_fix, aes(Component_1, Component_2)) +
    geom_point(aes(color=Pseudotime), size=0.3, alpha=0.5) +
    scale_color_viridis_c(option="plasma") +
    labs(title=paste0("Monocle2 Pseudotime — ",KEY_CELL),
         x="Component 1", y="Component 2") +
    theme_bw(12) + theme(plot.title=element_text(hjust=0.5,face="bold"))
  if (!is.null(tree_coords))
    p2 <- p2 + geom_segment(data=tree_coords, aes(x,y,xend=xend,yend=yend),
                            color="white", linewidth=1, inherit.aes=FALSE)
  ggsave(file.path(traj_dir_fix,"pseudotime_umap.png"), p2, width=8, height=6, dpi=150)
  cat("[INFO] 轨迹图 (fix) 已保存\n")
}
cat("[DONE] PART 4 完成\n")

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
OUT_DIR_RCT <- file.path(SERVER_DIR, "16_reactome")
dir.create(OUT_DIR_RCT, showWarnings=FALSE, recursive=TRUE)

deg_path_rct <- file.path(SERVER_DIR, "r.01_deg/DEG_GSE57338_full.csv")
if (!file.exists(deg_path_rct)) {
  cat("[SKIP] DEG 文件不存在\n")
} else {
  deg_rct <- read.csv(deg_path_rct, stringsAsFactors=FALSE)

  find_col <- function(df, candidates) {
    for (cn in candidates) if (cn %in% colnames(df)) return(cn)
    return(colnames(df)[1])
  }
  gene_col_rct <- find_col(deg_rct, c("gene","Gene","GENE","gene_name","symbol","Symbol","hgnc_symbol"))
  fc_col_rct   <- find_col(deg_rct, c("logFC","log2FoldChange","log2FC","LogFC","avg_log2FC"))
  pval_col_rct <- find_col(deg_rct, c("adj.P.Val","padj","FDR","p_val_adj","P.Value","pvalue"))
  cat(sprintf("[INFO] 基因列:%s FC列:%s pval列:%s\n", gene_col_rct, fc_col_rct, pval_col_rct))

  sc_obj_rct <- NULL
  for (p in c(file.path(SERVER_DIR,"r.03_single_cell/seurat_keycells.rds"),
              file.path(SERVER_DIR,"r.03_single_cell/seurat_annotated.rds"))) {
    if (file.exists(p) && requireNamespace("Seurat",quietly=TRUE)) {
      suppressPackageStartupMessages(library(Seurat))
      tryCatch({ sc_obj_rct <- readRDS(p); break }, error=function(e) NULL)
    }
  }

  deg_clean <- deg_rct[!is.na(deg_rct[[gene_col_rct]]) & deg_rct[[gene_col_rct]]!="" &
                         !is.na(deg_rct[[fc_col_rct]]),]
  deg_clean <- deg_clean[!duplicated(deg_clean[[gene_col_rct]]),]
  gene_list_rct <- sort(setNames(as.numeric(deg_clean[[fc_col_rct]]),
                                 deg_clean[[gene_col_rct]]), decreasing=TRUE)
  sig_genes_rct <- deg_clean[[gene_col_rct]][!is.na(deg_clean[[pval_col_rct]]) &
                                               deg_clean[[pval_col_rct]] < 0.05]

  gene_entrez_rct <- NULL; sig_entrez_rct <- NULL
  if (requireNamespace("org.Hs.eg.db",quietly=TRUE)) {
    suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
    tryCatch({
      emap <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=names(gene_list_rct),
                                    column="ENTREZID", keytype="SYMBOL", multiVals="first")
      ve <- !is.na(emap)
      gene_entrez_rct <- sort(gene_list_rct[ve], decreasing=TRUE)
      names(gene_entrez_rct) <- emap[ve]
      if (length(sig_genes_rct)>0) {
        semap <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=sig_genes_rct,
                                       column="ENTREZID", keytype="SYMBOL", multiVals="first")
        sig_entrez_rct <- semap[!is.na(semap)]
      }
    }, error=function(e) NULL)
  }

  reactomegsa_ok <- FALSE; reactomegsa_res <- NULL
  gsea_res <- NULL; ora_res <- NULL; method_rct <- "none"; result_df_rct <- NULL

  if (requireNamespace("ReactomeGSA",quietly=TRUE)) {
    suppressPackageStartupMessages(library(ReactomeGSA))
    tryCatch({
      bulk_expr <- matrix(deg_clean[[fc_col_rct]], nrow=nrow(deg_clean), ncol=1,
                          dimnames=list(deg_clean[[gene_col_rct]],"logFC"))
      req <- ReactomeAnalysisRequest(method="ssGSEA")
      req <- add_dataset(req, expression_values=bulk_expr, name="Bulk_RNA",
                         type="rnaseq_counts", comparison_factor="logFC")
      reactomegsa_res <- performAnalysis(req)
      reactomegsa_ok  <- TRUE; method_rct <- "ReactomeGSA"
    }, error=function(e) NULL)
  }

  if (!reactomegsa_ok && !is.null(gene_entrez_rct) && length(gene_entrez_rct)>=10 &&
      requireNamespace("ReactomePA",quietly=TRUE)) {
    suppressPackageStartupMessages(library(ReactomePA))
    tryCatch({
      gsea_res   <- gsePathway(gene_entrez_rct, organism="human", pvalueCutoff=0.2,
                               pAdjustMethod="BH", verbose=FALSE, seed=42)
      method_rct <- "ReactomePA_GSEA"
    }, error=function(e) NULL)
  }

  if (!reactomegsa_ok && is.null(gsea_res) && !is.null(sig_entrez_rct) &&
      length(sig_entrez_rct)>=5 && requireNamespace("ReactomePA",quietly=TRUE)) {
    suppressPackageStartupMessages(library(ReactomePA))
    tryCatch({
      ora_res    <- enrichPathway(gene=sig_entrez_rct,
                                  universe=if(!is.null(gene_entrez_rct)) names(gene_entrez_rct) else NULL,
                                  organism="human", pvalueCutoff=0.05, pAdjustMethod="BH",
                                  qvalueCutoff=0.2, readable=TRUE)
      method_rct <- "ReactomePA_ORA"
    }, error=function(e) NULL)
  }

  pathway_count_rct <- 0
  if (reactomegsa_ok && !is.null(reactomegsa_res)) {
    tryCatch({ result_df_rct <- result_summary(reactomegsa_res)
    pathway_count_rct <- nrow(result_df_rct) }, error=function(e) NULL)
  } else if (!is.null(gsea_res)) {
    result_df_rct <- as.data.frame(gsea_res); pathway_count_rct <- nrow(result_df_rct)
  } else if (!is.null(ora_res)) {
    result_df_rct <- as.data.frame(ora_res); pathway_count_rct <- nrow(result_df_rct)
  }
  write.csv(if (!is.null(result_df_rct) && nrow(result_df_rct)>0) result_df_rct else
    data.frame(pathway=character(),pvalue=numeric(),padj=numeric()),
    file.path(OUT_DIR_RCT,"reactome_results.csv"), row.names=FALSE)

  heatmap_path <- file.path(OUT_DIR_RCT, "reactome_heatmap.png")
  top10_path   <- file.path(OUT_DIR_RCT, "reactome_top10.png")

  if (reactomegsa_ok && !is.null(reactomegsa_res)) {
    tryCatch({ png(heatmap_path, width=1400, height=1000, res=150)
      plot_overview(reactomegsa_res); dev.off() },
      error=function(e) { if(dev.cur()>1) dev.off() })
  }


  if (!file.exists(heatmap_path) && !is.null(gsea_res)) {
    tryCatch({
      df_g <- as.data.frame(gsea_res); total_g <- nrow(df_g)
      df_g <- head(df_g[order(df_g$p.adjust),], 20)
      df_g$Description <- wrap_short(df_g$Description, 60)
      df_g$Description <- factor(df_g$Description, levels=rev(df_g$Description))
      nes_col <- if ("NES" %in% colnames(df_g)) "NES" else colnames(df_g)[3]
      h_g <- calc_plot_height(nrow(df_g), per_item_inch=0.38, min_height=7, max_height=24)
      p_g <- ggplot(df_g, aes_string(nes_col,"Description",color="p.adjust",size="setSize")) +
        geom_point() +
        scale_color_gradient(low="red", high="blue", name="Adj. P") +
        geom_vline(xintercept=0, linetype="dashed", color="grey50") +
        labs(title=sprintf("Reactome GSEA — Top 20 (Total: %d) [%s]", total_g, method_rct),
             x="NES", y=NULL) +
        theme_bw(11) + theme(plot.title=element_text(face="bold",hjust=0.5),
                             axis.text.y=element_text(size=8,lineheight=0.85),
                             plot.margin=ggplot2::margin(10,20,10,160))
      ggsave(heatmap_path, p_g, width=13, height=h_g, dpi=150, bg="white")

      df_top <- head(as.data.frame(gsea_res)[order(as.data.frame(gsea_res)$p.adjust),], 10)
      df_top$Description <- wrap_short(df_top$Description, 60)
      df_top$Description <- factor(df_top$Description, levels=rev(df_top$Description))
      h_t <- calc_plot_height(10, per_item_inch=0.45, min_height=6, max_height=12)
      p_t <- ggplot(df_top, aes_string(nes_col,"Description",color="p.adjust",size="setSize")) +
        geom_point() +
        scale_color_gradient(low="red", high="blue", name="Adj. P") +
        geom_vline(xintercept=0, linetype="dashed", color="grey50") +
        labs(title=sprintf("Reactome GSEA — Top 10 Pathways (Total: %d)", total_g),
             x="NES", y=NULL, size="Set Size") +
        theme_bw(12) + theme(plot.title=element_text(face="bold",hjust=0.5),
                             axis.text.y=element_text(size=10,lineheight=0.9),
                             plot.margin=ggplot2::margin(15,20,15,10))
      ggsave(top10_path, p_t, width=11, height=h_t, dpi=300, bg="white")
      cat("[INFO] Reactome GSEA top10 已保存\n")
    }, error=function(e) { if(dev.cur()>1) dev.off() })
  }


  if (!file.exists(heatmap_path) && !is.null(ora_res)) {
    tryCatch({
      df_o <- as.data.frame(ora_res); total_o <- nrow(df_o)
      df_o <- head(df_o[order(df_o$p.adjust),], 20)
      df_o$Description <- wrap_short(df_o$Description, 60)
      df_o$Description <- factor(df_o$Description, levels=rev(df_o$Description))
      parts <- strsplit(as.character(df_o$GeneRatio),"/")
      df_o$GeneRatioNum <- sapply(parts, function(x)
        if(length(x)==2) as.numeric(x[1])/as.numeric(x[2]) else NA)
      h_o <- calc_plot_height(nrow(df_o), per_item_inch=0.38, min_height=7, max_height=24)
      p_o <- ggplot(df_o, aes(GeneRatioNum,Description,color=p.adjust,size=Count)) +
        geom_point() +
        scale_color_gradient(low="red", high="blue", name="Adj. P") +
        labs(title=sprintf("Reactome ORA — Top 20 (Total: %d) [%s]", total_o, method_rct),
             x="Gene Ratio", y=NULL) +
        theme_bw(11) + theme(plot.title=element_text(face="bold",hjust=0.5),
                             axis.text.y=element_text(size=8,lineheight=0.85),
                             plot.margin=ggplot2::margin(10,20,10,160))
      ggsave(heatmap_path, p_o, width=13, height=h_o, dpi=150, bg="white")

      df_top <- head(as.data.frame(ora_res)[order(as.data.frame(ora_res)$p.adjust),], 10)
      df_top$Description <- wrap_short(df_top$Description, 60)
      df_top$Description <- factor(df_top$Description, levels=rev(df_top$Description))
      parts2 <- strsplit(as.character(df_top$GeneRatio),"/")
      df_top$GeneRatioNum <- sapply(parts2, function(x)
        if(length(x)==2) as.numeric(x[1])/as.numeric(x[2]) else NA)
      h_t <- calc_plot_height(10, per_item_inch=0.45, min_height=6, max_height=12)
      p_t <- ggplot(df_top, aes(GeneRatioNum,Description,color=p.adjust,size=Count)) +
        geom_point() +
        scale_color_gradient(low="red", high="blue", name="Adj. P") +
        labs(title=sprintf("Reactome ORA — Top 10 Pathways (Total: %d)", total_o),
             x="Gene Ratio", y=NULL, size="Gene Count") +
        theme_bw(12) + theme(plot.title=element_text(face="bold",hjust=0.5),
                             axis.text.y=element_text(size=10,lineheight=0.9),
                             plot.margin=ggplot2::margin(15,20,15,10))
      ggsave(top10_path, p_t, width=11, height=h_t, dpi=300, bg="white")
      cat("[INFO] Reactome ORA top10 已保存\n")
    }, error=function(e) { if(dev.cur()>1) dev.off() })
  }


  if (!file.exists(heatmap_path)) {
    p_e <- ggplot() +
      annotate("text", x=0.5, y=0.5, size=6, color="grey50",
               label=paste0("No significant Reactome pathways\nMethod: ", method_rct)) +
      xlim(0,1) + ylim(0,1) + theme_void() +
      labs(title="Reactome Pathway Analysis") +
      theme(plot.title=element_text(face="bold",hjust=0.5))
    ggsave(heatmap_path, p_e, width=10, height=6, dpi=150)
  }

  cat("[STAT] Reactome_pathway_count:", pathway_count_rct, "\n")
  cat("[INFO] 方法:", method_rct, "\n")
}
cat("[DONE] PART 5 完成\n")

cat("\n========== 全部分析完成 ==========\n")
cat(sprintf("[DONE] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
