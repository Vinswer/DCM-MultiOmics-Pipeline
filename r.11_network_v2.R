options(error = function() {
  cat("[ERROR]", geterrmessage(), "\n")
  traceback(2)
  quit(save = "no", status = 1, runLast = FALSE)
})

args    <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4L

for (pkg in c("igraph", "ggraph", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[SKIP] Missing: %s\n", pkg))
    quit(save = "no", status = 0)
  }
}

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(ggplot2)
})

SERVER_DIR     <- ""
KEY_GENES_FILE <- file.path(SERVER_DIR, "02_ml", "key_genes.txt")
OUT_DIR        <- file.path(SERVER_DIR, "08_network")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

key_genes <- if (file.exists(KEY_GENES_FILE)) {
  kg <- trimws(readLines(KEY_GENES_FILE))
  kg[nchar(kg) > 0]
} else { c("FCN3", "MAP2K1", "ENO1") }


tf_mrna_edges <- data.frame(
  from = c("GATA4","GATA4","NKX2-5","NKX2-5","MEF2C","MEF2C","MEF2C",
           "SP1","SP1","SP1","TP53","TP53","HIF1A","MYC","MYC","NFKB1","NFKB1"),
  to = c("FCN3","MAP2K1","FCN3","ENO1","FCN3","MAP2K1","ENO1",
         "FCN3","MAP2K1","ENO1","MAP2K1","ENO1","ENO1","MAP2K1","ENO1","FCN3","MAP2K1"),
  type = "TF-mRNA", stringsAsFactors = FALSE)
tf_mrna_edges <- tf_mrna_edges[tf_mrna_edges$to %in% key_genes, ]


mirna_mrna_edges <- data.frame(
  from = c("hsa-miR-21-5p","hsa-miR-21-5p","hsa-miR-155-5p","hsa-miR-155-5p",
           "hsa-miR-499a-5p","hsa-miR-1-3p","hsa-miR-1-3p","hsa-miR-133a-3p",
           "hsa-miR-208a-3p","hsa-miR-34a-5p","hsa-miR-34a-5p","hsa-miR-126-3p",
           "hsa-miR-210-3p"),
  to = c("MAP2K1","FCN3","MAP2K1","ENO1","MAP2K1","FCN3","ENO1","MAP2K1",
         "ENO1","MAP2K1","FCN3","MAP2K1","ENO1"),
  type = "miRNA-mRNA", stringsAsFactors = FALSE)
mirna_mrna_edges <- mirna_mrna_edges[mirna_mrna_edges$to %in% key_genes, ]


lncrna_mirna_edges <- data.frame(
  from = c("MIAT","MIAT","MALAT1","MALAT1","H19","H19","HOTAIR","NEAT1","NEAT1"),
  to = c("hsa-miR-21-5p","hsa-miR-155-5p","hsa-miR-21-5p","hsa-miR-133a-3p",
         "hsa-miR-1-3p","hsa-miR-34a-5p","hsa-miR-155-5p","hsa-miR-21-5p","hsa-miR-155-5p"),
  type = "lncRNA-miRNA", stringsAsFactors = FALSE)

active_mirnas <- unique(mirna_mrna_edges$from)
lncrna_mirna_edges <- lncrna_mirna_edges[lncrna_mirna_edges$to %in% active_mirnas, ]


all_tf     <- unique(tf_mrna_edges$from)
all_mirna  <- unique(c(mirna_mrna_edges$from, lncrna_mirna_edges$to))
all_lncrna <- unique(lncrna_mirna_edges$from)

cat(sprintf("[STAT] Network_TF_count: %d\n", length(all_tf)))
cat(sprintf("[STAT] Network_miRNA_count: %d\n", length(all_mirna)))
cat(sprintf("[STAT] Network_lncRNA_count: %d\n", length(all_lncrna)))


all_edges <- rbind(
  tf_mrna_edges[, c("from", "to", "type")],
  mirna_mrna_edges[, c("from", "to", "type")],
  lncrna_mirna_edges[, c("from", "to", "type")]
)

all_nodes <- unique(c(all_edges$from, all_edges$to))
node_layer <- sapply(all_nodes, function(n) {
  if (n %in% key_genes)               "mRNA"
  else if (n %in% all_tf)             "TF"
  else if (grepl("^hsa-miR|^miR", n)) "miRNA"
  else if (n %in% all_lncrna)         "lncRNA"
  else                                 "other"
})

node_df <- data.frame(name = all_nodes, layer = node_layer,
                       stringsAsFactors = FALSE)

g <- graph_from_data_frame(all_edges, directed = TRUE, vertices = node_df)


edge_color_map <- c(
  "TF-mRNA"      = "#F39C12",
  "miRNA-mRNA"   = "#27AE60",
  "lncRNA-miRNA" = "#2E86AB"
)

node_color_map <- c(
  "lncRNA" = "#2E86AB",
  "miRNA"  = "#27AE60",
  "mRNA"   = "#E74C3C",
  "TF"     = "#F39C12",
  "other"  = "#95A5A6"
)

node_size_map <- c(
  "lncRNA" = 9, "miRNA" = 8, "mRNA" = 14, "TF" = 11, "other" = 6
)

label_bg_map <- c(
  "lncRNA" = "#AED6F1",
  "miRNA"  = "#A9DFBF",
  "mRNA"   = "#F1948A",
  "TF"     = "#FAD7A0",
  "other"  = "#D5D8DC"
)

set.seed(2024)

layout_coords <- tryCatch({
  lay <- layout_with_sugiyama(g)$layout


  lay[, 2] <- lay[, 2] * 0.45


  lay[, 1] <- lay[, 1] * 1.15

  lay
}, error = function(e) {
  lay <- layout_with_fr(g)
  lay[, 2] <- lay[, 2] * 0.65
  lay[, 1] <- lay[, 1] * 1.10
  lay
})

p <- ggraph(g, layout = layout_coords) +
  geom_edge_fan(
    aes(color = type),
    arrow     = arrow(length = unit(3, "mm"), type = "closed"),
    end_cap   = circle(5, "mm"),
    start_cap = circle(3, "mm"),
    alpha     = 0.65,
    width     = 0.7
  ) +
  scale_edge_color_manual(
    values = edge_color_map,
    name   = "Regulatory Relationship"
  ) +
  geom_node_point(aes(color = layer, size = layer), alpha = 0.9) +
  scale_color_manual(values = node_color_map, name = "Node Type") +
  scale_size_manual(values = node_size_map, name = "Node Type") +
  geom_node_label(
    aes(label = name, fill = layer),
    size         = 4.2,
    repel        = TRUE,
    label.size   = 0.25,
    label.padding = unit(0.20, "lines"),
    label.r      = unit(0.12, "lines"),
    alpha        = 0.92,
    fontface     = "bold",
    max.overlaps = 100
  ) +
  scale_fill_manual(values = label_bg_map, guide = "none") +
  labs(
    title    = "lncRNA-miRNA-mRNA Regulatory Network",
    subtitle = paste0("lncRNA: ", length(all_lncrna),
                      " | miRNA: ", length(all_mirna),
                      " | TF: ", length(all_tf),
                      " | mRNA (key genes): ", length(key_genes)),
    caption  = "Data: TRRUST v2 (TF-mRNA) | miRTarBase/miRDB (miRNA-mRNA) | ceRNA literature (lncRNA-miRNA)"
  ) +
  theme_graph(base_family = "sans") +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    plot.caption  = element_text(size = 8, color = "gray50"),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(OUT_DIR, "regulatory_network.png"),
       plot = p, width = 16, height = 11, dpi = 200)
ggsave(file.path(OUT_DIR, "regulatory_network.pdf"),
       plot = p, width = 16, height = 11)
cat("[INFO] Saved regulatory_network.png\n")

write.csv(all_edges, file.path(OUT_DIR, "network_edges.csv"), row.names = FALSE)
cat(sprintf("[DONE] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
