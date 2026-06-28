

args <- commandArgs(trailingOnly = TRUE)
n_cores <- ifelse(length(args) >= 1, as.integer(args[1]), 4)
cat("[INFO] 使用核心数:", n_cores, "\n")


required_pkgs <- c("scMetabolism", "AUCell", "Seurat", "ggplot2", "dplyr", "tidyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("[SKIP] 缺少必要包:", pkg, "，跳过代谢分析\n")
    quit(status = 0)
  }
}


suppressPackageStartupMessages({
  library(scMetabolism)
  library(AUCell)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

cat("[INFO] Seurat 版本:", as.character(packageVersion("Seurat")), "\n")


SERVER_DIR         <- ""
ORIGINAL_DIR       <- SERVER_DIR
input_rds_primary  <- file.path(SERVER_DIR, "03_single_cell/seurat_keycells.rds")
input_rds_fallback <- file.path(SERVER_DIR, "03_single_cell/seurat_annotated.rds")
output_dir         <- file.path(SERVER_DIR, "15_metabolism")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


if (file.exists(input_rds_primary)) {
  cat("[INFO] 读取文件:", input_rds_primary, "\n")
  sc <- readRDS(input_rds_primary)
} else if (file.exists(input_rds_fallback)) {
  cat("[WARN] 主文件不存在，使用 fallback:", input_rds_fallback, "\n")
  sc <- readRDS(input_rds_fallback)
} else {
  cat("[SKIP] 输入文件均不存在，跳过代谢分析\n")
  quit(status = 0)
}

cat("[INFO] 细胞数:", ncol(sc), "| 基因数:", nrow(sc), "\n")


tryCatch({
  if (inherits(sc[["RNA"]], "Assay5")) {
    cat("[INFO] 检测到 Seurat v5 Assay5，执行 JoinLayers()...\n")
    sc[["RNA"]] <- JoinLayers(sc[["RNA"]])
  }
}, error = function(e) {
  cat("[WARN] JoinLayers 失败（可跳过）:", conditionMessage(e), "\n")
})


meta <- sc@meta.data
if ("celltype_annotation" %in% colnames(meta)) {
  group_col <- "celltype_annotation"
  cat("[INFO] 使用分组列: celltype_annotation\n")
} else if ("cell_type" %in% colnames(meta)) {
  group_col <- "cell_type"
  cat("[INFO] 使用分组列: cell_type\n")
} else {
  group_col <- "seurat_clusters"
  cat("[WARN] 未找到细胞类型列，使用 seurat_clusters\n")
}
cat("[INFO] 细胞类型分布:\n")
print(table(meta[[group_col]]))


cat("[INFO] 运行 scMetabolism（AUCell + KEGG）...\n")
sc <- tryCatch({
  sc.metabolism.Seurat(
    obj              = sc,
    method           = "AUCell",
    metabolism.type  = "KEGG",
    ncores           = n_cores
  )
}, error = function(e) {
  cat("[WARN] sc.metabolism.Seurat 报错:", conditionMessage(e), "\n")
  cat("[WARN] 尝试单核心运行...\n")
  tryCatch({
    sc.metabolism.Seurat(
      obj             = sc,
      method          = "AUCell",
      metabolism.type = "KEGG",
      ncores          = 1
    )
  }, error = function(e2) {
    cat("[SKIP] scMetabolism 分析失败:", conditionMessage(e2), "\n")
    quit(status = 0)
  })
})


if (!"METABOLISM" %in% names(sc@assays)) {
  cat("[SKIP] METABOLISM assay 未生成，跳过后续分析\n")
  quit(status = 0)
}



get_metab_matrix <- function(sc) {
  is_valid_mat <- function(m) {
    !is.null(m) && (is.matrix(m) || inherits(m, "sparseMatrix") ||
                    inherits(m, "dgeMatrix") || is.data.frame(m)) &&
      prod(dim(m)) > 0
  }


  for (layer_name in c("counts", "data", "scale.data")) {
    for (getter_arg in c("layer", "slot")) {
      result <- tryCatch({
        args <- list(object = sc, assay = "METABOLISM")
        args[[getter_arg]] <- layer_name
        mat <- do.call(GetAssayData, args)
        if (is_valid_mat(mat)) {
          cat("[INFO] 成功：GetAssayData(", getter_arg, "='", layer_name, "') 维度:",
              dim(mat)[1], "x", dim(mat)[2], "\n")
          return(mat)
        }
        NULL
      }, error = function(e) NULL)
      if (!is.null(result)) return(result)
    }
  }


  for (slot_name in c("counts", "data", "scale.data")) {
    result <- tryCatch({
      mat <- slot(sc@assays[["METABOLISM"]], slot_name)
      if (is_valid_mat(mat)) {
        cat("[INFO] 成功：@assays[['METABOLISM']]@", slot_name, " 维度:",
            dim(mat)[1], "x", dim(mat)[2], "\n")
        return(mat)
      }
      NULL
    }, error = function(e) NULL)
    if (!is.null(result)) return(result)
  }


  result <- tryCatch({
    layers <- sc@assays[["METABOLISM"]]@layers
    for (lname in names(layers)) {
      mat <- layers[[lname]]
      if (is_valid_mat(mat)) {
        cat("[INFO] 成功：@layers[['", lname, "']] 维度:",
            dim(mat)[1], "x", dim(mat)[2], "\n")
        return(mat)
      }
    }
    NULL
  }, error = function(e) NULL)
  if (!is.null(result)) return(result)


  metab_obj <- sc@assays[["METABOLISM"]]
  if (is.list(metab_obj)) {
    cat("[INFO] METABOLISM assay 是 list，names:", paste(names(metab_obj), collapse = ", "), "\n")
    for (key in c("data", "counts", "scale.data", names(metab_obj))) {
      mat <- metab_obj[[key]]
      if (is_valid_mat(mat)) {
        cat("[INFO] 成功：list[['", key, "']] 维度:", dim(mat)[1], "x", dim(mat)[2], "\n")
        return(mat)
      }
    }
  }


  for (misc_key in c("METABOLISM", "metabolism", "AUCell")) {
    result <- tryCatch({
      mat <- sc@misc[[misc_key]]
      if (is_valid_mat(mat)) {
        cat("[INFO] 成功：@misc[['", misc_key, "']] 维度:", dim(mat)[1], "x", dim(mat)[2], "\n")
        return(mat)
      }
      NULL
    }, error = function(e) NULL)
    if (!is.null(result)) return(result)
  }


  cat("[DEBUG] sc@assays names:", paste(names(sc@assays), collapse = ", "), "\n")
  if ("METABOLISM" %in% names(sc@assays)) {
    obj <- sc@assays[["METABOLISM"]]
    cat("[DEBUG] METABOLISM class:", paste(class(obj), collapse = "/"), "\n")
    cat("[DEBUG] typeof:", typeof(obj), "\n")
    if (is.list(obj)) cat("[DEBUG] list names:", paste(names(obj), collapse = ", "), "\n")
    tryCatch(cat("[DEBUG] slotNames:", paste(slotNames(obj), collapse = ", "), "\n"),
             error = function(e) NULL)
  }
  cat("[DEBUG] sc@misc names:", paste(names(sc@misc), collapse = ", "), "\n")
  return(NULL)
}

metab_data_raw <- get_metab_matrix(sc)
if (is.null(metab_data_raw)) {
  cat("[SKIP] 无法获取 METABOLISM 数据矩阵，跳过\n")
  quit(status = 0)
}


cat("[INFO] 转换代谢矩阵为稠密格式...\n")
metab_dense  <- as.matrix(metab_data_raw)
metab_mat    <- as.data.frame(t(metab_dense))
n_pathways   <- ncol(metab_mat)
cat("[INFO] 代谢通路数:", n_pathways, "\n")
cat("[INFO] 代谢矩阵维度:", nrow(metab_mat), "cells ×", n_pathways, "pathways\n")


metab_out <- metab_mat
metab_out$cell        <- rownames(metab_out)
metab_out$cell_type   <- meta[rownames(metab_out), group_col]
write.csv(metab_out,
          file.path(output_dir, "metabolism_scores.csv"),
          row.names = FALSE)
cat("[INFO] metabolism_scores.csv 已保存\n")


metab_mat$cell_type <- meta[rownames(metab_mat), group_col]
metab_summary <- metab_mat %>%
  group_by(cell_type) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
mat_mean <- as.data.frame(metab_summary[, -1])
rownames(mat_mean) <- metab_summary$cell_type


pathway_var <- apply(mat_mean, 2, var, na.rm = TRUE)
top20_paths <- names(sort(pathway_var, decreasing = TRUE))[1:min(20, n_pathways)]
cat("[INFO] Top20 代谢通路（按方差）:\n")
print(top20_paths)

top_pathway <- top20_paths[1]
cat(sprintf("[STAT] Metabolism_pathway_count: %d\n", n_pathways))
cat(sprintf("[STAT] Metabolism_top_pathway: %s\n", top_pathway))


cat("[INFO] 绘制 metabolism_heatmap.png ...\n")
tryCatch({
  hm_data <- mat_mean[, top20_paths, drop = FALSE]

  hm_scaled <- as.data.frame(scale(hm_data))
  hm_scaled$cell_type <- rownames(hm_scaled)

  hm_long <- tidyr::pivot_longer(
    hm_scaled,
    cols      = -cell_type,
    names_to  = "pathway",
    values_to = "zscore"
  )
  hm_long$cell_type <- factor(hm_long$cell_type, levels = rownames(hm_scaled))
  hm_long$pathway   <- factor(hm_long$pathway,   levels = rev(top20_paths))

  p_hm <- ggplot(hm_long, aes(x = cell_type, y = pathway, fill = zscore)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#D6604D",
      midpoint = 0,
      name     = "Z-score"
    ) +
    ggtitle("Metabolic Pathway Activity by Cell Type (Top 20)") +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y  = element_text(size = 8),
      axis.title   = element_blank(),
      plot.title   = element_text(hjust = 0.5, size = 12),
      legend.title = element_text(size = 9),
      panel.grid   = element_blank()
    )

  ggsave(file.path(output_dir, "metabolism_heatmap.png"),
         plot = p_hm, width = max(10, nrow(mat_mean) * 1.2), height = 10, dpi = 150)
  ggsave(file.path(output_dir, "metabolism_heatmap.pdf"),
         plot = p_hm, width = max(10, nrow(mat_mean) * 1.2), height = 10)
  cat("[INFO] metabolism_heatmap.png 已保存\n")
}, error = function(e) {
  cat("[WARN] metabolism_heatmap.png 绘制失败:", conditionMessage(e), "\n")
  try(dev.off(), silent = TRUE)
})


cat("[INFO] 绘制 metabolism_umap.png ...\n")
tryCatch({

  if (!"umap" %in% names(sc@reductions)) {
    cat("[INFO] 无 UMAP，先运行 RunPCA + RunUMAP...\n")
    if (!"pca" %in% names(sc@reductions)) {
      sc <- RunPCA(sc, assay = "RNA", verbose = FALSE)
    }
    sc <- RunUMAP(sc, reduction = "pca", dims = 1:20,
                  n.components = 2L, seed.use = 42, verbose = FALSE)
  }
  umap_embed <- as.data.frame(sc@reductions[["umap"]]@cell.embeddings)
  colnames(umap_embed) <- c("UMAP_1", "UMAP_2")


  top1_score <- metab_dense[top_pathway, ]
  umap_embed[[top_pathway]] <- top1_score[rownames(umap_embed)]
  umap_embed$cell_type <- meta[rownames(umap_embed), group_col]


  umap_embed <- umap_embed[order(umap_embed[[top_pathway]]), ]

  p_umap <- ggplot(umap_embed,
                   aes(x = UMAP_1, y = UMAP_2,
                       color = .data[[top_pathway]])) +
    geom_point(size = 0.4, alpha = 0.8) +
    scale_color_gradient(low = "lightgrey", high = "#D6604D",
                         name = "AUCell\nScore") +
    ggtitle(paste0("UMAP — ", top_pathway)) +
    theme_bw() +
    theme(
      plot.title  = element_text(hjust = 0.5, size = 12),
      panel.grid  = element_blank(),
      axis.title  = element_text(size = 10),
      legend.text = element_text(size = 8)
    )

  ggsave(file.path(output_dir, "metabolism_umap.png"),
         plot = p_umap, width = 9, height = 7, dpi = 150)
  ggsave(file.path(output_dir, "metabolism_umap.pdf"),
         plot = p_umap, width = 9, height = 7)
  cat("[INFO] metabolism_umap.png 已保存\n")
}, error = function(e) {
  cat("[WARN] metabolism_umap.png 绘制失败:", conditionMessage(e), "\n")
  try(dev.off(), silent = TRUE)
})

cat("[INFO] 代谢分析完成！输出目录:", output_dir, "\n")
cat(sprintf("[STAT] Metabolism_pathway_count: %d\n", n_pathways))
cat(sprintf("[STAT] Metabolism_top_pathway: %s\n", top_pathway))
