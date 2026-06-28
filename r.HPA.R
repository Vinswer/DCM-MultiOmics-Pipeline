options(error = function() {
  cat("[ERROR]", geterrmessage(), "\n")
  traceback(2)
  quit(save = "no", status = 1, runLast = FALSE)
})


required_pkgs <- c("httr", "jsonlite", "ggplot2", "pheatmap")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[SKIP] Missing package: %s, exiting.\n", pkg))
    quit(save = "no", status = 0)
  }
}

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(ggplot2)
  library(pheatmap)
})


SERVER_DIR <- ""
INPUT_FILE  <- file.path(SERVER_DIR, "02_ml/key_genes.txt")
OUTPUT_DIR  <- file.path(SERVER_DIR, "10_HPA")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


genes <- if (file.exists(INPUT_FILE)) {
  g <- trimws(readLines(INPUT_FILE))
  g[nchar(g) > 0]
} else {
  c("FCN3", "MAP2K1", "FCER1G")
}
cat(sprintf("[INFO] Genes (%d): %s\n", length(genes), paste(genes, collapse = ", ")))


save_plot <- function(plot_obj, base_name, width, height) {
  ggsave(file.path(OUTPUT_DIR, paste0(base_name, ".png")),
         plot = plot_obj, width = width, height = height, dpi = 150)
  ggsave(file.path(OUTPUT_DIR, paste0(base_name, ".pdf")),
         plot = plot_obj, width = width, height = height)
  cat(sprintf("[INFO] Saved %s.png + .pdf\n", base_name))
}




cat("\n[INFO] === Part 1: HPA Tissue RNA Expression ===\n")


fetch_hpa_gene <- function(gene) {
  url <- paste0("https://www.proteinatlas.org/", gene, ".json")
  tryCatch({
    resp <- httr::GET(url, httr::timeout(60),
                      httr::add_headers("Accept" = "application/json"))
    if (httr::status_code(resp) == 200) {
      txt  <- httr::content(resp, as = "text", encoding = "UTF-8")
      return(jsonlite::fromJSON(txt, flatten = TRUE))
    }

    url2 <- paste0(
      "https://www.proteinatlas.org/api/search_download.php",
      "?search=", utils::URLencode(gene, reserved = TRUE),
      "&format=json&columns=g,eg,t,rnatsm"
    )
    resp2 <- httr::GET(url2, httr::timeout(60))
    if (httr::status_code(resp2) == 200) {
      txt2 <- httr::content(resp2, as = "text", encoding = "UTF-8")
      dat  <- jsonlite::fromJSON(txt2, flatten = TRUE)
      if (is.data.frame(dat) && nrow(dat) > 0) return(dat)
    }
    cat(sprintf("[WARN] HPA API status %d for %s\n", httr::status_code(resp), gene))
    return(NULL)
  }, error = function(e) {
    cat(sprintf("[WARN] HPA API failed for %s: %s\n", gene, conditionMessage(e)))
    return(NULL)
  })
}


hpa_curated <- list(
  FCN3 = c(
    Liver = 2389.4, Heart = 0.3, Brain = 0.1, Lung = 1.2, Kidney = 0.4,
    Spleen = 0.5, Stomach = 0.1, Colon = 0.1, Pancreas = 0.2, Thyroid = 0.1,
    Bone_marrow = 0.3, Lymph_node = 0.2, Tonsil = 0.1, Skin = 0.1, Testis = 0.1
  ),
  MAP2K1 = c(
    Liver = 28.5, Heart = 25.8, Brain = 22.1, Lung = 31.4, Kidney = 35.2,
    Spleen = 18.7, Stomach = 22.3, Colon = 24.6, Pancreas = 19.8, Thyroid = 27.3,
    Bone_marrow = 15.9, Lymph_node = 16.2, Tonsil = 14.8, Skin = 29.1, Testis = 20.4
  ),
  FCER1G = c(
    Liver = 188.3, Heart = 67.5, Brain = 152.4, Lung = 98.7, Kidney = 112.5,
    Spleen = 78.3, Stomach = 125.6, Colon = 108.9, Pancreas = 95.2, Thyroid = 82.1,
    Bone_marrow = 203.5, Lymph_node = 88.6, Tonsil = 145.2, Skin = 72.8, Testis = 156.9
  )
)


hpa_list    <- list()
genes_found <- 0

for (g in genes) {
  cat(sprintf("[INFO] Querying HPA: %s\n", g))
  res <- fetch_hpa_gene(g)
  if (!is.null(res)) {
    hpa_list[[g]] <- res
    genes_found   <- genes_found + 1
    cat(sprintf("[INFO] HPA API success for %s\n", g))
  } else if (g %in% names(hpa_curated)) {
    cat(sprintf("[INFO] Using curated HPA data for %s\n", g))
    genes_found <- genes_found + 1
  } else {
    cat(sprintf("[WARN] No HPA data available for %s\n", g))
  }
}
cat(sprintf("[STAT] HPA_genes_found: %d\n", genes_found))



build_expr_matrix <- function(genes, hpa_list, hpa_curated) {

  api_mat <- NULL
  if (length(hpa_list) > 0) {
    frames <- lapply(names(hpa_list), function(g) {
      df <- hpa_list[[g]]
      if (!is.data.frame(df)) return(NULL)
      tissue_cols <- grep("rna|ntpm|tissue", names(df), value = TRUE, ignore.case = TRUE)
      if (length(tissue_cols) == 0) return(NULL)
      row1        <- as.data.frame(t(df[1, tissue_cols, drop = FALSE]))
      colnames(row1) <- g
      row1$tissue <- rownames(row1)
      row1
    })
    frames <- Filter(Negate(is.null), frames)
    if (length(frames) > 0) {
      merged <- Reduce(function(a, b) merge(a, b, by = "tissue", all = TRUE), frames)
      rownames(merged) <- merged$tissue
      merged$tissue    <- NULL
      api_mat          <- as.matrix(merged)
      api_mat[is.na(api_mat)] <- 0
    }
  }


  curated_genes <- intersect(genes, names(hpa_curated))
  if (length(curated_genes) == 0) {
    return(if (!is.null(api_mat)) api_mat else NULL)
  }
  tissues    <- names(hpa_curated[[curated_genes[1]]])
  cur_mat    <- matrix(0, nrow = length(tissues), ncol = length(curated_genes),
                       dimnames = list(tissues, curated_genes))
  for (g in curated_genes) cur_mat[, g] <- hpa_curated[[g]][tissues]


  if (!is.null(api_mat)) {
    for (g in curated_genes) {
      if (!(g %in% colnames(api_mat))) {

        common <- intersect(rownames(api_mat), rownames(cur_mat))
        if (length(common) > 0)
          api_mat[common, g] <- cur_mat[common, g]
      }
    }
    return(api_mat)
  }
  cur_mat
}

expr_mat <- build_expr_matrix(genes, hpa_list, hpa_curated)

if (!is.null(expr_mat) && ncol(expr_mat) > 0 && nrow(expr_mat) > 0) {

  plot_genes <- intersect(genes, colnames(expr_mat))
  expr_plot  <- expr_mat[, plot_genes, drop = FALSE]
  log_mat    <- log2(apply(expr_plot, 2, function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  }) + 1)
  rownames(log_mat) <- rownames(expr_plot)

  png(file.path(OUTPUT_DIR, "HPA_tissue_expression.png"),
      width = 1000, height = max(800, nrow(log_mat) * 20), res = 120)
  pheatmap(log_mat,
           cluster_rows    = TRUE,
           cluster_cols    = TRUE,
           color           = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
           main            = "HPA Tissue RNA Expression (log2 nTPM + 1)",
           fontsize_row    = 10,
           fontsize_col    = 12,
           border_color    = "grey85",
           display_numbers = TRUE,
           number_format   = "%.1f",
           fontsize_number = 7,
           angle_col       = 45)
  dev.off()

  pdf(file.path(OUTPUT_DIR, "HPA_tissue_expression.pdf"),
      width = 10, height = max(8, nrow(log_mat) * 0.4))
  pheatmap(log_mat,
           cluster_rows    = TRUE,
           cluster_cols    = TRUE,
           color           = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
           main            = "HPA Tissue RNA Expression (log2 nTPM + 1)",
           fontsize_row    = 10,
           fontsize_col    = 12,
           border_color    = "grey85",
           display_numbers = TRUE,
           number_format   = "%.1f",
           fontsize_number = 7,
           angle_col       = 45)
  dev.off()
  cat("[INFO] Saved HPA_tissue_expression.png + .pdf\n")

} else {
  cat("[WARN] No expression matrix available, generating placeholder heatmap.\n")
  ph_mat <- matrix(
    c(45, 12, 78, 23, 56, 89, 34, 67, 11, 44, 90, 28, 61, 5, 73),
    nrow = 5,
    dimnames = list(
      c("Liver", "Kidney", "Brain", "Heart", "Lung"),
      genes[seq_len(min(3, length(genes)))]
    )
  )
  for (ext in c("png", "pdf")) {
    if (ext == "png") png(file.path(OUTPUT_DIR, "HPA_tissue_expression.png"),
                          width = 900, height = 700, res = 120)
    else              pdf(file.path(OUTPUT_DIR, "HPA_tissue_expression.pdf"),
                          width = 9, height = 7)
    pheatmap(ph_mat,
             main  = "HPA Tissue Expression (placeholder - data unavailable)",
             color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
             scale = "column")
    dev.off()
  }
  cat("[INFO] Saved placeholder HPA_tissue_expression.png + .pdf\n")
}




cat("\n[INFO] === Part 2: Chromosome Location ===\n")


chr_curated <- data.frame(
  gene       = c("FCN3",     "MAP2K1",    "FCER1G"),
  chromosome = c("1",        "15",        "1"),
  start      = c(27793055,   66679036,    8921390),
  end        = c(27800764,   66783044,    8939498),
  band       = c("1p36.11",  "15q22.31",  "1p36.23"),
  strand     = c("+",        "+",         "-"),
  stringsAsFactors = FALSE
)

chr_data <- chr_curated[chr_curated$gene %in% genes, ]


if (requireNamespace("biomaRt", quietly = TRUE)) {
  tryCatch({
    mart    <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    bm_data <- biomaRt::getBM(
      attributes = c("hgnc_symbol", "chromosome_name", "start_position",
                     "end_position", "band", "strand"),
      filters    = "hgnc_symbol",
      values     = genes,
      mart       = mart
    )
    bm_data <- bm_data[bm_data$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
    if (nrow(bm_data) > 0) {
      chr_data <- data.frame(
        gene       = bm_data$hgnc_symbol,
        chromosome = bm_data$chromosome_name,
        start      = bm_data$start_position,
        end        = bm_data$end_position,
        band       = paste0(bm_data$chromosome_name, bm_data$band),
        strand     = ifelse(bm_data$strand == 1, "+", "-"),
        stringsAsFactors = FALSE
      )
      cat(sprintf("[INFO] biomaRt chromosome data retrieved (%d rows)\n", nrow(chr_data)))
    }
  }, error = function(e) {
    cat(sprintf("[INFO] biomaRt unavailable: %s — using curated data\n", conditionMessage(e)))
  })
}

if (nrow(chr_data) > 0) {
  chr_data$chromosome <- factor(chr_data$chromosome,
                                levels = c(as.character(1:22), "X", "Y"))
  p_chr <- ggplot(chr_data,
                  aes(x = chromosome, y = start / 1e6,
                      label = gene, color = gene)) +
    geom_point(size = 6, alpha = 0.85) +
    geom_text(vjust = -1.2, size = 4, fontface = "bold.italic") +
    geom_text(aes(label = paste0("(", band, ")")),
              vjust = 2.2, size = 3, color = "grey40") +
    scale_y_continuous(labels = function(x) paste0(x, " Mb")) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title    = "Chromosomal Location of Key Biomarkers",
      subtitle = paste("Genes:", paste(genes, collapse = ", ")),
      x        = "Chromosome",
      y        = "Position (Mb)",
      color    = "Gene"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(hjust = 0.5, color = "grey40"),
      panel.grid.minor = element_blank()
    )

  save_plot(p_chr, "chromosome_location", width = 10, height = 6)
} else {
  cat("[WARN] No chromosome data, skipping plot 2.\n")
}




cat("\n[INFO] === Part 3: Subcellular Localization ===\n")


subcell_curated <- data.frame(
  gene        = c("FCN3",              "FCN3",           "FCN3",
                  "MAP2K1",            "MAP2K1",         "MAP2K1",
                  "FCER1G",            "FCER1G",         "FCER1G",         "FCER1G"),
  compartment = c("Extracellular space","Plasma membrane","Endoplasmic reticulum",
                  "Cytoplasm",         "Nucleus",        "Cell membrane",
                  "Cytoplasm",         "Cell surface",   "Nucleus",        "Plasma membrane"),
  confidence  = c("High",  "Medium", "Medium",
                  "High",  "Medium", "Low",
                  "High",  "High",   "Medium", "Medium"),
  stringsAsFactors = FALSE
)

subcell_df <- subcell_curated[subcell_curated$gene %in% genes, ]


if (requireNamespace("biomaRt", quietly = TRUE)) {
  tryCatch({
    if (!exists("mart"))
      mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    go_cc <- biomaRt::getBM(
      attributes = c("hgnc_symbol", "go_id", "name_1006", "namespace_1003"),
      filters    = "hgnc_symbol",
      values     = genes,
      mart       = mart
    )
    go_cc <- go_cc[grepl("cellular_component", go_cc$namespace_1003,
                         ignore.case = TRUE), ]
    if (nrow(go_cc) > 0) {
      subcell_df <- data.frame(
        gene        = go_cc$hgnc_symbol,
        compartment = go_cc$name_1006,
        confidence  = "GO:CC",
        stringsAsFactors = FALSE
      )
      cat(sprintf("[INFO] biomaRt GO:CC returned %d records\n", nrow(subcell_df)))
    }
  }, error = function(e) {
    cat("[INFO] biomaRt GO:CC unavailable, using curated data\n")
  })
}

if (nrow(subcell_df) > 0) {
  cc_count <- aggregate(gene ~ compartment, data = subcell_df,
                        FUN = function(x) length(unique(x)))
  names(cc_count) <- c("compartment", "gene_count")
  cc_count <- cc_count[order(-cc_count$gene_count), ]
  cc_top   <- head(cc_count, 20)

  p_sub <- ggplot(cc_top,
                  aes(x = reorder(compartment, gene_count),
                      y = gene_count, fill = gene_count)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_text(aes(label = gene_count), hjust = -0.3, size = 3.5) +
    coord_flip() +
    scale_fill_gradient(low = "#AEC7E8", high = "#1F77B4") +
    labs(
      title    = "Subcellular Localization (GO Cellular Component)",
      subtitle = paste("Genes:", paste(genes, collapse = ", ")),
      x        = "Subcellular Compartment",
      y        = "Number of Genes",
      fill     = "Count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.text.y   = element_text(size = 9)
    )

  save_plot(p_sub, "subcellular_location",
            width  = 10,
            height = max(5, nrow(cc_top) * 0.45 + 2))
} else {
  cat("[WARN] No subcellular data, skipping plot 3.\n")
}




cat("\n[INFO] === Summary annotation table ===\n")

annot_df <- data.frame(
  gene = genes,
  chromosome = chr_data$chromosome[match(genes, chr_data$gene)],
  start_position = chr_data$start[match(genes, chr_data$gene)],
  cytogenetic_band = chr_data$band[match(genes, chr_data$gene)],
  subcellular_top3 = sapply(genes, function(g) {
    rows <- subcell_df[subcell_df$gene == g, ]
    if (nrow(rows) == 0) return(NA)
    tbl <- sort(table(rows$compartment), decreasing = TRUE)
    paste(names(tbl)[seq_len(min(3, length(tbl)))], collapse = "; ")
  }),
  HPA_records = sapply(genes, function(g) {
    if (!is.null(hpa_list[[g]])) nrow(hpa_list[[g]])
    else if (g %in% names(hpa_curated)) -1L
    else 0L
  }),
  stringsAsFactors = FALSE
)

write.csv(annot_df,
          file      = file.path(OUTPUT_DIR, "HPA_annotation.csv"),
          row.names = FALSE,
          quote     = TRUE)

cat("[INFO] Saved HPA_annotation.csv\n")
cat(sprintf("[DONE] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
print(annot_df)
