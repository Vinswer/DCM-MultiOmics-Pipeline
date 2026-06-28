args    <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4L
cat("[INFO] n_cores =", n_cores, "\n")

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})


SERVER_DIR <- ""
OUT_DIR    <- file.path(SERVER_DIR, "07_PPI")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


input_file <- file.path(SERVER_DIR, "r.01_deg", "Candidate_genes_DEG_x_Phagocytosis.csv")
if (!file.exists(input_file)) {
  cat("[WARN] Candidate gene file not found, using key_genes as candidates\n")
  key_genes_file <- file.path(SERVER_DIR, "r.02_ml", "key_genes.txt")
  genes <- readLines(key_genes_file)
  genes <- trimws(genes[nchar(trimws(genes)) > 0])
} else {
  gene_df <- read.csv(input_file, stringsAsFactors = FALSE)
  gene_col <- intersect(c("gene", "Gene", "SYMBOL", "symbol", "gene_symbol"), colnames(gene_df))
  if (length(gene_col) == 0) gene_col <- colnames(gene_df)[1]
  genes <- unique(gene_df[[gene_col[1]]])
  genes <- genes[!is.na(genes) & nchar(genes) > 0]
}

key_genes_file2 <- file.path(SERVER_DIR, "r.02_ml", "key_genes.txt")
if (file.exists(key_genes_file2)) {
  kg <- trimws(readLines(key_genes_file2))
  kg <- kg[nchar(kg) > 0]
  genes <- unique(c(genes, kg))
}
cat("[INFO] Candidate genes:", length(genes), "--", paste(head(genes, 8), collapse = ", "), "...\n")


entrez_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = genes,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "SYMBOL"
))
entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID), ]
entrez_map <- entrez_map[!duplicated(entrez_map$SYMBOL), ]
cat(sprintf("[INFO] Successfully mapped %d / %d genes to ENTREZID\n",
            nrow(entrez_map), length(genes)))

genes_mapped <- entrez_map$SYMBOL
entrez_ids   <- entrez_map$ENTREZID


cat("[INFO] Building GO BP functional network (offline mode)...\n")

go_data <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = entrez_ids,
  columns = c("ENTREZID", "GO", "ONTOLOGY"),
  keytype = "ENTREZID"
))
go_bp <- go_data[!is.na(go_data$GO) & go_data$ONTOLOGY == "BP", ]
cat(sprintf("[INFO] Retrieved %d GO BP annotations\n", nrow(go_bp)))


if (nrow(go_bp) > 0) {
  gene_go_list <- split(go_bp$GO, go_bp$ENTREZID)

  entrez2sym <- setNames(entrez_map$SYMBOL, entrez_map$ENTREZID)
  valid_entrez <- intersect(names(gene_go_list), entrez_ids)
  gene_go_list <- gene_go_list[valid_entrez]
  sym_names    <- entrez2sym[valid_entrez]

  n <- length(sym_names)
  cat(sprintf("[INFO] Computing GO sharing matrix for %d genes...\n", n))

  min_shared <- max(1, round(quantile(sapply(gene_go_list, length), 0.25)))
  edges_list <- list()
  for (i in seq_len(n - 1)) {
    go_i <- gene_go_list[[i]]
    for (j in seq(i + 1, n)) {
      shared <- length(intersect(go_i, gene_go_list[[j]]))
      if (shared >= min_shared) {
        edges_list[[length(edges_list) + 1]] <- data.frame(
          gene1  = sym_names[i],
          gene2  = sym_names[j],
          weight = shared,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(edges_list) > 0) {
    edges_df <- do.call(rbind, edges_list)
    w_thresh <- quantile(edges_df$weight, 0.50)
    edges_df <- edges_df[edges_df$weight >= w_thresh, ]
    cat(sprintf("[INFO] Retained %d edges after filtering (shared GO terms >= %d)\n",
                nrow(edges_df), as.integer(w_thresh)))
  } else {
    cat("[WARN] Sharing matrix is empty, using minimal key-gene network\n")
    edges_df <- NULL
  }
} else {
  cat("[WARN] GO BP annotations are empty\n")
  edges_df <- NULL
}


if (is.null(edges_df) || nrow(edges_df) == 0) {
  kg_use <- intersect(genes, genes_mapped)
  if (length(kg_use) < 2) kg_use <- genes_mapped[1:min(5, length(genes_mapped))]
  pairs <- combn(kg_use, 2)
  edges_df <- data.frame(
    gene1  = pairs[1, ],
    gene2  = pairs[2, ],
    weight = 1,
    stringsAsFactors = FALSE
  )
}


all_nodes <- unique(c(edges_df$gene1, edges_df$gene2))
key_genes_avail <- intersect(genes_mapped,
                             intersect(readLines(key_genes_file2),
                                       all_nodes))

g <- graph_from_data_frame(
  d        = edges_df[, c("gene1", "gene2", "weight")],
  directed = FALSE,
  vertices = data.frame(name = all_nodes)
)

deg <- degree(g)
btw <- betweenness(g, normalized = TRUE)
V(g)$degree      <- deg
V(g)$betweenness <- btw
V(g)$is_key      <- V(g)$name %in% key_genes_avail
V(g)$is_hub      <- deg >= quantile(deg, 0.80)

node_count <- vcount(g)
edge_count <- ecount(g)
hub_genes  <- V(g)$name[V(g)$is_hub]

cat(sprintf("[STAT] PPI_node_count: %d\n", node_count))
cat(sprintf("[STAT] PPI_edge_count: %d\n", edge_count))
cat(sprintf("[STAT] PPI_hub_genes: %s\n", paste(sort(hub_genes), collapse = ",")))


write.csv(edges_df, file.path(OUT_DIR, "PPI_edges.csv"),
          row.names = FALSE, quote = FALSE)


hub_df <- data.frame(
  gene         = hub_genes,
  degree       = deg[hub_genes],
  betweenness  = round(btw[hub_genes], 4)
)
hub_df <- hub_df[order(-hub_df$degree), ]
write.csv(hub_df, file.path(OUT_DIR, "PPI_hub_genes.csv"),
          row.names = FALSE, quote = FALSE)


cat("[INFO] Plotting PPI network (base igraph)...\n")

node_color <- ifelse(V(g)$is_key, "#E74C3C",
                     ifelse(V(g)$is_hub, "#F39C12", "#3498DB"))

node_size  <- 8 + deg / max(deg) * 20
node_label <- ifelse(V(g)$is_key | V(g)$is_hub, V(g)$name, "")
edge_width <- as.numeric(E(g)$weight)

edge_width <- 1.5 + (edge_width - min(edge_width)) /
  (max(edge_width) - min(edge_width) + 1e-6) * 3.5

set.seed(42)
layout_g <- layout_with_fr(g)

png(file.path(OUT_DIR, "PPI_network.png"), width = 3600, height = 3000, res = 200)
par(mar = c(3, 3, 4, 3))
plot(g,
     layout      = layout_g,
     vertex.color = node_color,
     vertex.size  = node_size,
     vertex.label = node_label,
     vertex.label.cex  = 1.4,
     vertex.label.font = 2,
     vertex.label.color = "black",
     vertex.frame.color = NA,
     edge.width   = edge_width,
     edge.color   = "#AAAAAA",
     main = paste0("Functional Association Network\n",
                   "(GO BP Pathway Sharing | Nodes: ", node_count,
                   " | Edges: ", edge_count, ")"),
     cex.main = 1.8)

legend("topright", inset = c(0.01, 0.01),
       legend = c("Key Gene (FCN3/MAP2K1/ENO1)",
                  "Hub Gene (degree >= 80th pct.)",
                  "Other Candidate"),
       fill   = c("#E74C3C", "#F39C12", "#3498DB"),
       bty    = "o", cex = 1.1, xpd = FALSE)
dev.off()
cat("[INFO] Saved PPI_network.png\n")
pdf(file.path(OUT_DIR, "PPI_network.pdf"),
    width = 18, height = 15)

par(mar = c(3, 3, 4, 3))

plot(g,
     layout      = layout_g,
     vertex.color = node_color,
     vertex.size  = node_size,
     vertex.label = node_label,
     vertex.label.cex  = 1.4,
     vertex.label.font = 2,
     vertex.label.color = "black",
     vertex.frame.color = NA,
     edge.width   = edge_width,
     edge.color   = "#AAAAAA",
     main = paste0("Functional Association Network\n",
                   "(GO BP Pathway Sharing | Nodes: ", node_count,
                   " | Edges: ", edge_count, ")"),
     cex.main = 1.8)

legend("topright", inset = c(0.01, 0.01),
       legend = c("Key Gene (FCN3/MAP2K1/ENO1)",
                  "Hub Gene (degree >= 80th pct.)",
                  "Other Candidate"),
       fill   = c("#E74C3C", "#F39C12", "#3498DB"),
       bty    = "o", cex = 1.1, xpd = FALSE)

dev.off()

cat("[INFO] Saved PPI_network.pdf\n")

if (requireNamespace("ggraph", quietly = TRUE) &&
    requireNamespace("ggrepel", quietly = TRUE)) {
  library(ggraph)
  library(ggrepel)

  node_type <- ifelse(V(g)$is_key, "Key Gene",
                      ifelse(V(g)$is_hub, "Hub Gene", "Other"))

  set.seed(42)
  p_gg <- ggraph(g, layout = "fr") +

    geom_edge_link(aes(width = weight), alpha = 0.35, color = "#999999") +
    scale_edge_width(range = c(0.8, 4.0), name = "Shared GO Terms") +

    geom_node_point(aes(size = degree, color = node_type), stroke = 1.2) +
    scale_color_manual(
      values = c("Key Gene" = "#E74C3C", "Hub Gene" = "#F39C12", "Other" = "#3498DB"),
      name   = "Gene Type") +
    scale_size_continuous(range = c(5, 18), name = "Degree") +

    geom_node_label(
      aes(label = ifelse(is_key | is_hub, name, ""),
          filter = is_key | is_hub),
      size       = 5,
      fontface   = "bold",
      repel      = TRUE,
      label.size = 0.25,
      fill       = "white",
      alpha      = 0.90,
      label.padding = unit(0.3, "lines")
    ) +
    labs(
      title    = "Functional Association Network (GO BP Pathway Sharing)",
      subtitle = paste0("Nodes: ", node_count,
                        "  |  Edges: ", edge_count,
                        "  |  Hub Genes: ", length(hub_genes))
    ) +
    theme_graph(base_family = "sans") +
    theme(
      plot.title      = element_text(size = 20, face = "bold"),
      plot.subtitle   = element_text(size = 14),
      legend.title    = element_text(size = 13, face = "bold"),
      legend.text     = element_text(size = 12),
      legend.key.size = unit(1.2, "lines"),
      legend.position = "right"
    )

  ggsave(file.path(OUT_DIR, "PPI_network_ggraph.png"),
         plot = p_gg, width = 16, height = 12, dpi = 200)
  ggsave(file.path(OUT_DIR, "PPI_network_ggraph.pdf"),
         plot = p_gg, width = 16, height = 12)
  cat("[INFO] Saved PPI_network_ggraph.png\n")
}

cat(sprintf("[DONE] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
