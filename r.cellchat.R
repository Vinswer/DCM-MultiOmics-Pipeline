rm(list = ls()); gc()
cat("[INFO] r.16c_cellchat_cardio.R 开始\n")
cat(sprintf("[INFO] R: %s\n", R.version$version.string))

required <- c("CellChat", "Seurat", "ggplot2", "patchwork",
              "igraph", "dplyr", "tidyr", "RColorBrewer",
              "ggrepel", "scales")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[SKIP] 缺少必要包: %s\n", pkg)); quit(status = 0)
  }
}
suppressPackageStartupMessages({
  library(CellChat); library(Seurat)
  library(ggplot2);  library(patchwork); library(igraph)
  library(dplyr);    library(tidyr);    library(RColorBrewer)
  library(ggrepel);  library(scales)
})

select    <- dplyr::select;   filter    <- dplyr::filter
rename    <- dplyr::rename;   mutate    <- dplyr::mutate
arrange   <- dplyr::arrange;  summarise <- dplyr::summarise
group_by  <- dplyr::group_by; ungroup   <- dplyr::ungroup
left_join <- dplyr::left_join; slice    <- dplyr::slice
slice_max <- dplyr::slice_max; pull     <- dplyr::pull

WORK    <- ""
SC_RDS  <- file.path(WORK, "03_single_cell/seurat_clustered.rds")
OUT_CC  <- file.path(WORK, "11_cellchat");  dir.create(OUT_CC,  recursive=TRUE, showWarnings=FALSE)
OUT_SUB <- file.path(WORK, "12_cardio_sub"); dir.create(OUT_SUB, recursive=TRUE, showWarnings=FALSE)

KEY_CELL <- "Macrophages"
KEY_GENES <- c("FCN3","FCER1G","MAP2K1")

save_plot <- function(p, base, width=10, height=8, dpi=150, dir=OUT_CC) {
  for (ext in c("png","pdf")) {
    fp <- file.path(dir, paste0(base, ".", ext))
    if (ext=="png") png(fp, width=width, height=height, units="in", res=dpi)
    else            pdf(fp, width=width, height=height)
    print(p); dev.off()
    cat(sprintf("[INFO] 保存: %s\n", fp))
  }
}

sc <- readRDS(SC_RDS)
cat(sprintf("[INFO] 加载完成: %d genes × %d cells\n", nrow(sc), ncol(sc)))

group_col <- NULL
for (cand in c("Condition","condition","disease_group","Group","group","orig.ident")) {
  if (cand %in% colnames(sc@meta.data)) { group_col <- cand; break }
}
if (is.null(group_col)) {
  cat("[ERROR] 找不到分组列，退出。\n"); quit(status=1)
}
cat(sprintf("[INFO] 分组列: %s，取值: %s\n", group_col,
            paste(unique(sc@meta.data[[group_col]]), collapse=", ")))

group_vals <- unique(sc@meta.data[[group_col]])
hf_pat   <- "DCM|HF|heart.fail|fail|disease|case"
ctrl_pat <- "Donor|ctrl|control|normal|healthy|non"
hf_vals  <- group_vals[grepl(hf_pat,   group_vals, ignore.case=TRUE)]
ct_vals  <- group_vals[grepl(ctrl_pat, group_vals, ignore.case=TRUE)]

sc@meta.data$cc_group <- ifelse(
  sc@meta.data[[group_col]] %in% hf_vals,  "DCM",
  ifelse(sc@meta.data[[group_col]] %in% ct_vals, "Donor", NA)
)
sc <- subset(sc, !is.na(cc_group))
cat(sprintf("[INFO] DCM: %d cells | Donor: %d cells\n",
            sum(sc$cc_group=="DCM"), sum(sc$cc_group=="Donor")))

ct_col <- NULL
for (cand in c("annotated_cell_type","cell_type","celltype","Celltype")) {
  if (cand %in% colnames(sc@meta.data)) { ct_col <- cand; break }
}
if (is.null(ct_col)) { cat("[ERROR] 找不到细胞类型列\n"); quit(status=1) }
sc@meta.data$cc_celltype <- as.character(sc@meta.data[[ct_col]])

db <- CellChatDB.human
db_genes <- unique(unlist(strsplit(unique(c(
  unlist(db$interaction[,"ligand"]),
  unlist(db$interaction[,"receptor"]))), "_")))
db_genes <- db_genes[nchar(db_genes)>0 & !is.na(db_genes)]

full_mat <- tryCatch(
  SeuratObject::LayerData(sc, assay="RNA", layer="data"),
  error=function(e) tryCatch(sc@assays$RNA@data, error=function(e2) NULL)
)
if (is.null(full_mat)) { cat("[ERROR] 无法提取表达矩阵\n"); quit(status=1) }

keep_genes <- intersect(rownames(full_mat), db_genes)
sub_mat    <- full_mat[keep_genes, ]
cat(sprintf("[INFO] CellChat 子矩阵: %d 基因 × %d 细胞\n",
            nrow(sub_mat), ncol(sub_mat)))

meta_cc <- data.frame(
  cc_group    = sc$cc_group,
  cc_celltype = sc$cc_celltype,
  row.names   = colnames(sc)
)
rm(full_mat, sc); gc()

MAX_PER_TYPE <- 150
make_subset_cc <- function(grp, min_cells=5) {
  cells <- rownames(meta_cc)[meta_cc$cc_group == grp]
  ct    <- meta_cc[cells, "cc_celltype"]
  names(ct) <- cells
  ct_tab     <- table(ct)
  valid_cts  <- names(ct_tab[ct_tab >= min_cells])
  if (length(valid_cts) < 2) {
    cat(sprintf("[SKIP] %s 有效细胞类型不足\n", grp)); return(NULL)
  }
  keep <- cells[ct %in% valid_cts]
  ct2  <- ct[keep]
  set.seed(42)
  sampled <- unlist(lapply(valid_cts, function(c_) {
    idx <- names(ct2)[ct2==c_]
    if (length(idx)<=MAX_PER_TYPE) idx else sample(idx, MAX_PER_TYPE)
  }))
  mat_out  <- sub_mat[, sampled, drop=FALSE]
  meta_out <- data.frame(cell_type=ct2[sampled], row.names=sampled)
  cat(sprintf("[INFO] %s: %d 细胞，%d 类型\n",
              grp, ncol(mat_out), length(valid_cts)))
  list(mat=mat_out, meta=meta_out, valid_cts=valid_cts)
}

dcm_dat   <- make_subset_cc("DCM")
donor_dat <- make_subset_cc("Donor")

run_cellchat <- function(dat, label) {
  if (is.null(dat)) return(NULL)
  cat(sprintf("\n[INFO] === CellChat: %s ===\n", label))

  cc <- tryCatch(
    createCellChat(object=dat$mat, meta=dat$meta, group.by="cell_type"),
    error=function(e) { cat("[WARN]", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(cc)) return(NULL)
  cc@DB <- db

  cc <- tryCatch(subsetData(cc),                     error=function(e) cc)
  cc <- tryCatch(identifyOverExpressedGenes(cc),      error=function(e) cc)
  cc <- tryCatch(identifyOverExpressedInteractions(cc), error=function(e) cc)

  cc <- tryCatch(
    computeCommunProb(cc, type="triMean", population.size=FALSE),
    error=function(e) { cat("[WARN] computeCommunProb:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(cc)) return(NULL)

  cc <- tryCatch(filterCommunication(cc, min.cells=5), error=function(e) cc)
  cc <- tryCatch(computeCommunProbPathway(cc),          error=function(e) cc)
  cc <- tryCatch(aggregateNet(cc),                      error=function(e) cc)

  saveRDS(cc, file.path(OUT_CC, sprintf("cellchat_%s.rds", label)))
  cat(sprintf("[INFO] 已保存 cellchat_%s.rds\n", label))
  cc
}

cc_dcm   <- run_cellchat(dcm_dat,   "DCM")
rm(dcm_dat); gc()
cc_donor <- run_cellchat(donor_dat, "Donor")
rm(donor_dat); gc()

all_cts <- unique(c(
  if (!is.null(cc_dcm))   levels(cc_dcm@idents)   else NULL,
  if (!is.null(cc_donor)) levels(cc_donor@idents) else NULL
))
n_ct <- length(all_cts)
pal  <- setNames(
  colorRampPalette(brewer.pal(min(n_ct,12),"Set3"))(n_ct),
  all_cts
)

draw_circle_net <- function(cc, label, metric=c("count","weight"), color_use=NULL) {
  metric <- match.arg(metric)
  if (is.null(cc)) { cat(sprintf("[SKIP] %s %s 圈图\n", label, metric)); return(invisible(NULL)) }
  mat <- if (metric=="count") cc@net$count else cc@net$weight
  gS  <- as.numeric(table(cc@idents))
  clr <- if (!is.null(color_use)) color_use[match(rownames(mat), names(color_use))] else NULL
  title_str <- sprintf("%s\n(%s)",
                       if (metric=="count") "Number of Interactions" else "Interaction Strength",
                       label)
  par(mar=c(1,1,3,1))
  netVisual_circle(mat,
                   vertex.weight = gS,
                   weight.scale  = TRUE,
                   label.edge    = FALSE,
                   color.use     = clr,
                   title.name    = title_str)
}


if (!is.null(cc_dcm)) {
  for (ext in c("png","pdf")) {
    fp <- file.path(OUT_CC, paste0("Fig22A_circle_count_DCM.", ext))
    if (ext=="png") png(fp, width=9, height=9, units="in", res=150)
    else pdf(fp, width=9, height=9)
    draw_circle_net(cc_dcm, "DCM (Heart Failure)", "count", pal)
    dev.off(); cat(sprintf("[INFO] 保存: %s\n", fp))
  }
}

if (!is.null(cc_donor)) {
  for (ext in c("png","pdf")) {
    fp <- file.path(OUT_CC, paste0("Fig22A_circle_count_Donor.", ext))
    if (ext=="png") png(fp, width=9, height=9, units="in", res=150)
    else pdf(fp, width=9, height=9)
    draw_circle_net(cc_donor, "Donor (Normal Control)", "count", pal)
    dev.off(); cat(sprintf("[INFO] 保存: %s\n", fp))
  }
}


if (!is.null(cc_dcm)) {
  for (ext in c("png","pdf")) {
    fp <- file.path(OUT_CC, paste0("Fig22A_circle_weight_DCM.", ext))
    if (ext=="png") png(fp, width=9, height=9, units="in", res=150)
    else pdf(fp, width=9, height=9)
    draw_circle_net(cc_dcm, "DCM (Heart Failure)", "weight", pal)
    dev.off(); cat(sprintf("[INFO] 保存: %s\n", fp))
  }
}

if (!is.null(cc_donor)) {
  for (ext in c("png","pdf")) {
    fp <- file.path(OUT_CC, paste0("Fig22A_circle_weight_Donor.", ext))
    if (ext=="png") png(fp, width=9, height=9, units="in", res=150)
    else pdf(fp, width=9, height=9)
    draw_circle_net(cc_donor, "Donor (Normal Control)", "weight", pal)
    dev.off(); cat(sprintf("[INFO] 保存: %s\n", fp))
  }
}

for (ext in c("png","pdf")) {
  fp <- file.path(OUT_CC, paste0("Fig22A_circle_combined_2x2.", ext))
  if (ext=="png") png(fp, width=16, height=16, units="in", res=150)
  else pdf(fp, width=16, height=16)
  par(mfrow=c(2,2), mar=c(1,1,3,1))
  draw_circle_net(cc_dcm,   "DCM (Heart Failure)",   "count",  pal)
  draw_circle_net(cc_donor, "Donor (Normal Control)", "count",  pal)
  draw_circle_net(cc_dcm,   "DCM (Heart Failure)",   "weight", pal)
  draw_circle_net(cc_donor, "Donor (Normal Control)", "weight", pal)
  dev.off(); cat(sprintf("[INFO] 保存: %s\n", fp))
}
cat("[INFO] 圈图说明: 上行=通讯次数(Number of Interactions), 下行=通讯强度(Interaction Strength)\n")
cat("[INFO]           左列=DCM(心力衰竭), 右列=Donor(正常对照)\n")

for (obj in list(list(cc=cc_dcm, lbl="DCM"), list(cc=cc_donor, lbl="Donor"))) {
  if (!is.null(obj$cc)) {
    tryCatch({
      p_h <- netVisual_heatmap(obj$cc,
                               color.heatmap = "Reds",
                               title.name    = sprintf("Interaction Strength (%s)", obj$lbl))
      save_plot(p_h, sprintf("Fig22A_heatmap_%s", obj$lbl),
                width=9, height=7, dir=OUT_CC)
    }, error=function(e) cat(sprintf("[WARN] 热图 %s: %s\n", obj$lbl, conditionMessage(e))))
  }
}


get_lr_df <- function(cc, group_label) {
  if (is.null(cc)) return(NULL)
  tryCatch({
    df <- subsetCommunication(cc)
    if (nrow(df) == 0) return(NULL)

    mol_cols <- intersect(c("ligand","receptor"), colnames(df))
    if ("ligand" %in% colnames(df) && "receptor" %in% colnames(df)) {

      df$interaction_name_2 <- paste(df$ligand, df$receptor, sep=" - ")
    }

    pval_col <- intersect(c("pval","p.value","pvalue","prob"), colnames(df))[1]
    prob_col <- intersect(c("prob","probability","Probability"), colnames(df))[1]

    if (!is.na(pval_col))  df$pval_use <- df[[pval_col]]  else df$pval_use <- 1
    if (!is.na(prob_col))  df$prob_use <- df[[prob_col]]  else df$prob_use <- 0

    df$log2mean <- log2(df$prob_use * 1e4 + 1)
    df_filt <- df[df$pval_use <= 0.05, ]
    if (nrow(df_filt) == 0) {
      cat(sprintf("[WARN] %s: 筛选后无互作行，降低阈值重试\n", group_label))
      df_filt <- df[df$pval_use <= 0.05, ]
    }
    df_filt$group <- group_label

    cat(sprintf("[DEBUG] subsetCommunication 列名: %s\n",
                paste(colnames(df_filt), collapse=", ")))

    char_cols <- names(df_filt)[sapply(df_filt, is.character) | sapply(df_filt, is.factor)]
    df_filt$involves_key <- apply(
      df_filt[, char_cols, drop=FALSE], 1,
      function(row) any(grepl(KEY_CELL, row, ignore.case=TRUE))
    )

    cat(sprintf("[INFO] %s: 筛选后 %d 条互作（Cardiomyocytes 相关: %d 条）\n",
                group_label, nrow(df_filt), sum(df_filt$involves_key)))
    df_filt
  }, error=function(e) {
    cat(sprintf("[WARN] get_lr_df %s: %s\n", group_label, conditionMessage(e)))
    NULL
  })
}

lr_dcm   <- get_lr_df(cc_dcm,   "DCM")
lr_donor <- get_lr_df(cc_donor, "Donor")


plot_lr_bubble <- function(df, title_str, top_n=30, highlight_key=TRUE) {
  if (is.null(df) || nrow(df)==0) return(NULL)


  df <- df[order(-df$prob_use), ]
  if (nrow(df) > top_n) df <- df[1:top_n, ]


  df$interaction_label <- paste0(df$source, " → ", df$target,
                                 "\n(", df$ligand, " - ", df$receptor, ")")


  df$interaction_label <- ifelse(nchar(df$interaction_label)>60,
                                 paste0(substr(df$interaction_label,1,58),"…"),
                                 df$interaction_label)


  df$neg_log_p <- -log10(df$pval_use + 1e-300)

  color_vals <- if (highlight_key) {
    c("FALSE"="grey70", "TRUE"="#E8534A")
  } else {
    c("FALSE"="steelblue", "TRUE"="steelblue")
  }

  ggplot(df, aes(x=log2mean, y=reorder(interaction_label, log2mean),
                 size=neg_log_p, color=as.character(involves_key))) +
    geom_point(alpha=0.85) +
    scale_size_continuous(range=c(2,8), name="-log10(p-value)") +
    scale_color_manual(values=color_vals,
                       labels=c("Other interactions",
                                sprintf("%s-related", KEY_CELL)),
                       name="Interaction type") +
    geom_vline(xintercept=0.1, linetype="dashed", color="grey50") +
    labs(x="log2 mean (communication probability)",
         y="Ligand - Receptor Interaction",
         title=title_str) +
    theme_bw(base_size=10) +
    theme(plot.title  = element_text(hjust=0.5, size=12),
          axis.text.y = element_text(size=7),
          legend.position = "right")
}


if (!is.null(lr_dcm)) {
  p_bubble_dcm <- plot_lr_bubble(lr_dcm,
                                 sprintf("Ligand-Receptor Interactions (DCM)\n[%s highlighted]", KEY_CELL))
  if (!is.null(p_bubble_dcm))
    save_plot(p_bubble_dcm, "Fig22B_LR_bubble_DCM",
              width=14, height=10, dir=OUT_CC)
}

if (!is.null(lr_donor)) {
  p_bubble_donor <- plot_lr_bubble(lr_donor,
                                   sprintf("Ligand-Receptor Interactions (Donor)\n[%s highlighted]", KEY_CELL))
  if (!is.null(p_bubble_donor))
    save_plot(p_bubble_donor, "Fig22B_LR_bubble_Donor",
              width=14, height=10, dir=OUT_CC)
}


if (!is.null(lr_dcm) && !is.null(lr_donor)) {
  tryCatch({

    combined_lr <- rbind(
      head(lr_dcm[order(-lr_dcm$prob_use),], 20),
      head(lr_donor[order(-lr_donor$prob_use),], 20)
    )
    combined_lr$pathway_lr <- paste0(combined_lr$ligand," - ",combined_lr$receptor)
    combined_lr$cell_pair  <- paste0(combined_lr$source," → ",combined_lr$target)
    combined_lr$neg_log_p  <- -log10(combined_lr$pval_use + 1e-300)
    combined_lr$group      <- factor(combined_lr$group, levels=c("Donor","DCM"))
    combined_lr$involves_key <- as.character(combined_lr$involves_key)

    p_compare <- ggplot(combined_lr,
                        aes(x=group, y=reorder(pathway_lr, prob_use),
                            size=prob_use, color=involves_key)) +
      geom_point(alpha=0.85) +
      scale_size_continuous(range=c(2,9), name="Communication\nProbability") +
      scale_color_manual(values=c("FALSE"="grey70","TRUE"="#E8534A"),
                         labels=c("Other",sprintf("%s-related",KEY_CELL)),
                         name="Interaction") +
      facet_wrap(~cell_pair, scales="free_y", ncol=4) +
      labs(x=NULL, y="Ligand - Receptor",
           title=sprintf("LR Interactions: DCM vs Donor\n(%s highlighted)", KEY_CELL)) +
      theme_bw(base_size=9) +
      theme(plot.title  = element_text(hjust=0.5, size=12),
            axis.text.x = element_text(size=9),
            strip.text  = element_text(size=7))

    save_plot(p_compare, "Fig22B_LR_bubble_comparison",
              width=16, height=10, dir=OUT_CC)
  }, error=function(e)
    cat(sprintf("[WARN] 合并气泡图: %s\n", conditionMessage(e))))
}

for (obj in list(list(cc=cc_dcm, lbl="DCM"), list(cc=cc_donor, lbl="Donor"))) {
  if (!is.null(obj$cc)) {
    tryCatch({
      df <- subsetCommunication(obj$cc)
      write.csv(df, file.path(OUT_CC, sprintf("cellchat_LR_%s.csv", obj$lbl)),
                row.names=FALSE)
      cat(sprintf("[INFO] 保存 cellchat_LR_%s.csv (%d 行)\n", obj$lbl, nrow(df)))
    }, error=function(e) cat(sprintf("[WARN] CSV %s: %s\n", obj$lbl, conditionMessage(e))))
  }
}


sc_full <- readRDS(SC_RDS)

ct_col2 <- NULL
for (cand in c("annotated_cell_type","cell_type","celltype","Celltype")) {
  if (cand %in% colnames(sc_full@meta.data)) { ct_col2 <- cand; break }
}
if (is.null(ct_col2)) { cat("[ERROR] 找不到细胞类型列\n"); quit(status=1) }

cat(sprintf("[INFO] 可用细胞类型: %s\n",
            paste(sort(unique(sc_full@meta.data[[ct_col2]])), collapse=", ")))

cardio_cells <- colnames(sc_full)[sc_full@meta.data[[ct_col2]] == KEY_CELL]
cat(sprintf("[INFO] %s 细胞数: %d\n", KEY_CELL, length(cardio_cells)))

if (length(cardio_cells) < 10) {
  cat(sprintf("[WARN] 细胞数过少（%d），跳过 PART B。\n", length(cardio_cells)))
} else {

  sc_cardio <- subset(sc_full, cells = cardio_cells)

  group_col_b <- NULL
  for (cand in c("cc_group","Condition","condition","disease_group","Group","group")) {
    if (cand %in% colnames(sc_full@meta.data)) { group_col_b <- cand; break }
  }
  rm(sc_full); gc()
  cat(sprintf("[INFO] 子集: %d genes × %d cells\n", nrow(sc_cardio), ncol(sc_cardio)))

  if (!"percent.mito" %in% colnames(sc_cardio@meta.data))
    sc_cardio[["percent.mito"]] <- PercentageFeatureSet(sc_cardio, pattern="^MT-")

  n_before_qc2 <- ncol(sc_cardio)
  sc_cardio_filt <- tryCatch(
    subset(sc_cardio,
           subset = nFeature_RNA >= 200 &
             nFeature_RNA <= 5000 &
             percent.mito < 20  &
             nCount_RNA   < 15000),
    error = function(e) sc_cardio
  )
  if (ncol(sc_cardio_filt) >= 10) {
    sc_cardio <- sc_cardio_filt
  }
  n_after_qc2 <- ncol(sc_cardio)
  cat(sprintf("[INFO] 二次质控: %d → %d 细胞（移除 %d）\n",
              n_before_qc2, n_after_qc2, n_before_qc2 - n_after_qc2))

  sc_cardio <- NormalizeData(sc_cardio,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000,
                             verbose = FALSE)
  sc_cardio <- FindVariableFeatures(sc_cardio,
                                    selection.method = "vst",
                                    nfeatures = 2000,
                                    verbose = FALSE)
  sc_cardio <- ScaleData(sc_cardio, verbose = FALSE)


  cat("[B-4] PCA + ElbowPlot...\n")
  n_cells_c <- ncol(sc_cardio)
  max_npcs  <- max(2L, min(n_cells_c - 1L, 50L))
  sc_cardio <- RunPCA(sc_cardio, npcs = max_npcs, verbose = FALSE)

  pca_std_c <- sc_cardio@reductions$pca@stdev
  N_SHOW_C  <- min(20L, length(pca_std_c))
  pct_var_c <- pca_std_c^2 / sum(pca_std_c^2) * 100
  cum_var_c <- cumsum(pct_var_c)
  delta_c   <- diff(pca_std_c[seq_len(N_SHOW_C)])
  thr_c     <- mean(abs(delta_c)) * 0.5
  elb_c     <- which(abs(delta_c) < thr_c)[1]
  if (is.na(elb_c) || elb_c < 2) elb_c <- min(5L, N_SHOW_C)
  n_dims_c  <- min(elb_c + 1L, N_SHOW_C, max_npcs)
  cat(sprintf("[INFO] 选择 PC 数: %d (累计方差 %.1f%%)\n",
              n_dims_c, cum_var_c[n_dims_c]))

  elbow_df_c <- data.frame(PC  = seq_len(N_SHOW_C),
                           std = pca_std_c[seq_len(N_SHOW_C)])
  p_elbow_c <- ggplot(elbow_df_c, aes(x=PC, y=std)) +
    geom_point(size=2.5, color="steelblue") +
    geom_line(color="steelblue", linewidth=0.8) +
    geom_vline(xintercept=n_dims_c, linetype="dashed",
               color="red", linewidth=0.9) +
    annotate("text", x=n_dims_c+0.3, y=max(elbow_df_c$std)*0.92,
             label=sprintf("Selected: %d PCs\n(Cum. var. %.1f%%)",
                           n_dims_c, cum_var_c[n_dims_c]),
             hjust=0, color="red", size=3.5) +
    scale_x_continuous(breaks=seq(1, N_SHOW_C, by=2)) +
    labs(x="PC", y="Standard Deviation",
         title=sprintf("Elbow Plot — %s Subclustering", KEY_CELL)) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=12))
  save_plot(p_elbow_c, "Fig23_cardio_ElbowPlot",
            width=8, height=5, dir=OUT_SUB)


  n_neighbors_c <- max(2L, min(10L, n_cells_c - 1L))
  sc_cardio <- FindNeighbors(sc_cardio, dims=1:n_dims_c,
                             k.param=n_neighbors_c, verbose=FALSE)
  sc_cardio <- FindClusters(sc_cardio, resolution=0.8, verbose=FALSE)
  sc_cardio <- RunUMAP(sc_cardio, dims=1:n_dims_c,
                       n.neighbors=n_neighbors_c,
                       seed.use=42, verbose=FALSE)
  n_sub_clusters <- length(levels(sc_cardio$seurat_clusters))
  cat(sprintf("[INFO] 亚型聚类数（res=0.8）: %d\n", n_sub_clusters))


  Idents(sc_cardio) <- "seurat_clusters"
  markers_file_c <- file.path(OUT_SUB, "cardio_subtype_markers.rds")
  if (file.exists(markers_file_c)) {
    markers_c <- readRDS(markers_file_c)
    cat("[INFO] 加载已有 marker 文件\n")
  } else {
    markers_c <- tryCatch(
      FindAllMarkers(sc_cardio,
                     only.pos        = TRUE,
                     min.pct         = 0.25,
                     logfc.threshold = 0.25,
                     verbose         = FALSE),
      error=function(e) {
        cat(sprintf("[WARN] FindAllMarkers: %s\n", conditionMessage(e))); NULL
      }
    )
    if (!is.null(markers_c) && nrow(markers_c)>0)
      saveRDS(markers_c, markers_file_c)
  }

  has_singler <- requireNamespace("SingleR",    quietly=TRUE) &&
    requireNamespace("celldex",    quietly=TRUE) &&
    requireNamespace("SummarizedExperiment", quietly=TRUE)
  singler_labels <- NULL
  if (has_singler) {
    suppressPackageStartupMessages({
      library(SingleR); library(celldex); library(SummarizedExperiment)
    })
    tryCatch({
      ref_data <- celldex::HumanPrimaryCellAtlasData()
      expr_singler <- tryCatch(
        GetAssayData(sc_cardio, assay="RNA", layer="data"),
        error=function(e) sc_cardio@assays$RNA@data
      )
      sr_res <- SingleR(test      = expr_singler,
                        ref       = ref_data,
                        labels    = ref_data$label.main,
                        de.method = "wilcox")
      singler_labels <- sr_res$pruned.labels
      singler_labels[is.na(singler_labels)] <- "Unknown"
      names(singler_labels) <- rownames(sr_res)
      sc_cardio <- AddMetaData(sc_cardio,
                               metadata = singler_labels,
                               col.name = "singler_label")
      cat("[INFO] SingleR 注释完成\n")
      print(table(singler_labels))
    }, error=function(e)
      cat(sprintf("[WARN] SingleR 失败: %s\n", conditionMessage(e))))
  } else {
    cat("[WARN] SingleR/celldex 未安装，跳过\n")
  }

  cardio_subtype_markers <- list(
    CM_stressed      = c("HSPA1A","HSPA1B","HSPB1","HSP90AA1","DDIT3","ATF3","HSPA6"),
    CM_hypertrophic  = c("MYH7","NPPA","NPPB","ACTA1","MYL4","ANKRD1","XIRP2"),
    CM_ischemic      = c("LDHA","PKM","ENO1","HK2","VEGFA","SLC2A1","BNIP3"),
    CM_mature        = c("TNNT2","MYL2","MYBPC3","ACTC1","TNNC1","TTN","MYL3"),
    CM_proliferating = c("MKI67","TOP2A","PCNA","CDK1","CCNB1","CENPF","STMN1"),
    CM_apoptotic     = c("CASP3","BAX","BBC3","FAS","CYCS","APAF1","PMAIP1"),
    CM_fibrotic      = c("COL1A1","COL3A1","FN1","POSTN","ACTA2","TGM2","CTGF")
  )
  for (sub_ct in names(cardio_subtype_markers)) {
    genes_use <- intersect(cardio_subtype_markers[[sub_ct]], rownames(sc_cardio))
    if (length(genes_use) < 2) next
    sc_cardio <- tryCatch(
      AddModuleScore(sc_cardio,
                     features = list(genes_use),
                     name     = paste0("score_", sub_ct),
                     ctrl     = min(50L, nrow(sc_cardio)-1L),
                     seed     = 42),
      error=function(e) sc_cardio
    )
  }
  score_cols_c <- paste0("score_", names(cardio_subtype_markers), "1")
  score_cols_c <- score_cols_c[score_cols_c %in% colnames(sc_cardio@meta.data)]

  if (length(score_cols_c) > 0) {
    score_mat_c <- as.matrix(sc_cardio@meta.data[, score_cols_c, drop=FALSE])
    ct_names_c  <- gsub("score_(.+)1$","\\1", score_cols_c)
    colnames(score_mat_c) <- ct_names_c
    best_c <- apply(score_mat_c, 1, which.max)
    sc_cardio$subtype_module <- ct_names_c[best_c]

    meta_b <- sc_cardio@meta.data
    meta_b$seurat_clusters <- as.character(meta_b$seurat_clusters)
    sub_anno <- meta_b %>%
      dplyr::group_by(seurat_clusters, subtype_module) %>%
      dplyr::summarise(n=n(), .groups="drop") %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::slice_max(order_by=n, n=1, with_ties=FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(seurat_clusters, subtype_module)

    sub_map_b <- setNames(as.character(sub_anno$subtype_module),
                          as.character(sub_anno$seurat_clusters))

    HARDCODE_ANNO <- list(
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
    HARDCODE_ANNO <- unlist(HARDCODE_ANNO)

    clusters_chr <- as.character(sc_cardio$seurat_clusters)
    hardcode_vec <- ifelse(
      clusters_chr %in% names(HARDCODE_ANNO),
      HARDCODE_ANNO[clusters_chr],
      paste0("Mac_", clusters_chr)
    )
    names(hardcode_vec) <- colnames(sc_cardio)
    sc_cardio <- AddMetaData(sc_cardio, metadata=hardcode_vec,
                             col.name="annotated_subtype")
    cat("[INFO] HARDCODE_ANNO 已覆盖 annotated_subtype\n")
    cat("[INFO] 亚型注释结果:\n"); print(as.data.frame(sub_anno))
    write.table(
      sc_cardio@meta.data[, intersect(
        c("seurat_clusters","annotated_subtype","subtype_module","singler_label"),
        colnames(sc_cardio@meta.data)), drop=FALSE],
      file.path(OUT_SUB, "cardio_subtype_annotation.txt"),
      sep="\t", quote=FALSE, row.names=TRUE)
  } else {
    sc_cardio$annotated_subtype <- paste0("Sub_", sc_cardio$seurat_clusters)
    sub_map_b <- setNames(paste0("Sub_", levels(sc_cardio$seurat_clusters)),
                          levels(sc_cardio$seurat_clusters))
  }


  Idents(sc_cardio) <- "annotated_subtype"
  sub_cts_all <- sort(unique(sc_cardio$annotated_subtype))
  n_sub       <- length(sub_cts_all)
  sub_pal     <- setNames(
    colorRampPalette(brewer.pal(min(n_sub,9),"Set1"))(n_sub),
    sub_cts_all)

  p23a_sub <- DimPlot(sc_cardio, reduction="umap",
                      group.by="annotated_subtype",
                      label=TRUE, label.size=4,
                      pt.size=0.8, repel=TRUE, cols=sub_pal) +
    ggtitle(sprintf("%s Subtypes (res=0.8)", KEY_CELL)) +
    theme(plot.title=element_text(hjust=0.5, size=13))
  save_plot(p23a_sub, "Fig23A_UMAP_subtype", width=11, height=8, dir=OUT_SUB)

  if (!is.null(group_col_b) && group_col_b %in% colnames(sc_cardio@meta.data)) {
    grp_pal <- c("DCM"="#E8534A","Donor"="#6AA8E8",
                 "HeartFailure"="#E8534A","Control"="#6AA8E8")
    grp_use <- unique(sc_cardio@meta.data[[group_col_b]])
    grp_col <- grp_pal[grp_use]; grp_col[is.na(grp_col)] <- "grey70"

    p23a_grp <- DimPlot(sc_cardio, reduction="umap",
                        group.by=group_col_b, pt.size=0.8,
                        cols=grp_col) +
      ggtitle(sprintf("%s by Group", KEY_CELL)) +
      theme(plot.title=element_text(hjust=0.5, size=13))
    save_plot(p23a_grp, "Fig23A_UMAP_group", width=10, height=7, dir=OUT_SUB)

    p23a_split <- DimPlot(sc_cardio, reduction="umap",
                          group.by="annotated_subtype",
                          split.by=group_col_b,
                          pt.size=0.8, cols=sub_pal) +
      ggtitle(sprintf(" %s Split by Group", KEY_CELL)) +
      theme(plot.title=element_text(hjust=0.5, size=13))
    save_plot(p23a_split, "Fig23A_UMAP_split", width=16, height=7, dir=OUT_SUB)
  }

  if (!is.null(markers_c) && nrow(markers_c) > 0) {
    top3_vln <- markers_c %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(order_by=avg_log2FC, n=3) %>%
      dplyr::ungroup()
    vln_genes <- unique(top3_vln$gene)
    vln_genes  <- intersect(vln_genes, rownames(sc_cardio))
    if (length(vln_genes) > 24) vln_genes <- vln_genes[1:24]

    if (length(vln_genes) > 0) {
      Idents(sc_cardio) <- "annotated_subtype"
      tryCatch({
        p23b <- VlnPlot(sc_cardio,
                        features  = vln_genes,
                        group.by  = "annotated_subtype",
                        pt.size   = 0,
                        ncol      = min(4L, length(vln_genes)),
                        cols      = sub_pal) &
          theme(axis.text.x  = element_text(size=7, angle=45, hjust=1),
                axis.title.x = element_blank(),
                plot.title   = element_text(size=9))
        save_plot(p23b, "Fig23B_VlnPlot_subtype_markers",
                  width  = min(4L, length(vln_genes)) * 4,
                  height = ceiling(length(vln_genes)/4) * 3,
                  dir    = OUT_SUB)
      }, error=function(e)
        cat(sprintf("[WARN] VlnPlot: %s\n", conditionMessage(e))))
    }
  } else {
    cat("[WARN] 无 marker，跳过 VlnPlot\n")
  }

  if (!is.null(markers_c) && nrow(markers_c) > 0) {
    top5_dot <- markers_c %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(order_by=avg_log2FC, n=5) %>%
      dplyr::ungroup()
    dot_genes_c <- unique(top5_dot$gene)
    dot_genes_c <- intersect(dot_genes_c, rownames(sc_cardio))
    if (length(dot_genes_c) > 60) dot_genes_c <- dot_genes_c[1:60]

    if (length(dot_genes_c) > 0) {
      Idents(sc_cardio) <- "annotated_subtype"
      tryCatch({
        expr_dot <- FetchData(sc_cardio, vars=dot_genes_c)
        expr_dot$subtype <- sc_cardio$annotated_subtype
        dot_long_c <- tidyr::pivot_longer(expr_dot,
                                          cols=all_of(dot_genes_c),
                                          names_to="gene", values_to="expr")
        dot_sum_c <- dot_long_c %>%
          dplyr::group_by(subtype, gene) %>%
          dplyr::summarise(avg_expr = mean(expm1(expr)),
                           pct_expr = mean(expr>0)*100,
                           .groups  = "drop")

        sub_order <- unique(sub_map_b[as.character(
          sort(as.integer(names(sub_map_b))))])
        sub_order <- sub_order[!is.na(sub_order)]
        dot_sum_c$gene    <- factor(dot_sum_c$gene, levels=dot_genes_c)
        dot_sum_c$subtype <- factor(dot_sum_c$subtype, levels=rev(unique(sub_order)))

        p23c <- ggplot(dot_sum_c,
                       aes(x=gene, y=subtype,
                           size=pct_expr, color=avg_expr)) +
          geom_point() +
          scale_size_continuous(range=c(0.3,7), breaks=c(0,25,50,75,100),
                                name="% Expressed") +
          scale_color_gradient2(low="lightgrey", mid="#FFCCBB", high="red",
                                midpoint=median(dot_sum_c$avg_expr, na.rm=TRUE),
                                name="Avg Expr") +
          ggtitle("Top 5 Marker Genes per Subtype") +
          theme_bw(base_size=10) +
          theme(axis.text.x = element_text(size=7, angle=90, vjust=0.5, hjust=1),
                axis.text.y = element_text(size=9),
                axis.title  = element_blank(),
                plot.title  = element_text(hjust=0.5, size=13))
        save_plot(p23c, "Fig23C_DotPlot_subtype_markers",
                  width  = max(14, length(dot_genes_c)*0.25),
                  height = max(6,  n_sub*0.5+2),
                  dir    = OUT_SUB)
      }, error=function(e)
        cat(sprintf("[WARN] DotPlot: %s\n", conditionMessage(e))))
    }
  } else {
    cat("[WARN] 无 marker，跳过 DotPlot\n")
  }


  HARDCODE_ANNO <- list(
    "0"  = "HP+ Mac",
    "1"  = "S100A8+ Mono-Mac",
    "2"  = "VSIG4+ Res-Mac",
    "3"  = "MMP19+ Inflam-Mac",
    "4"  = "VMO1+ Mac",
    "5"  = "IFIT+ IFN-Mac",
    "6"  = "TSPAN18+ Mac",
    "7"  = "PTH1R+ Mac",
    "8"  = "MYOC+ Mac"
  )
  HARDCODE_ANNO <- unlist(HARDCODE_ANNO)

  key_avail <- intersect(KEY_GENES, rownames(sc_cardio))
  if (length(key_avail) > 0) {
    Idents(sc_cardio) <- "seurat_clusters"
    expr_key <- FetchData(sc_cardio, vars=key_avail)

    cells_keep <- clusters_num %in% names(HARDCODE_ANNO)
    cat(sprintf("[INFO] 过滤无注释簇后保留细胞: %d / %d\n",
                sum(cells_keep), length(cells_keep)))


    sc_sub      <- sc_cardio[, cells_keep]
    expr_key    <- FetchData(sc_sub, vars = key_avail)

    annotated_labels <- unname(HARDCODE_ANNO[
      as.character(sc_sub$seurat_clusters)])

    cluster_order_lab <- unname(HARDCODE_ANNO[
      as.character(sort(as.integer(names(HARDCODE_ANNO))))])

    expr_key$Cluster <- factor(annotated_labels, levels = cluster_order_lab)
    key_long <- tidyr::pivot_longer(expr_key,
                                    cols=all_of(key_avail),
                                    names_to="Gene", values_to="Expression")
    key_long$Gene <- factor(key_long$Gene, levels=KEY_GENES)
    n_cl   <- length(unique(key_long$Cluster))
    cl_pal <- setNames(
      colorRampPalette(brewer.pal(min(n_cl,10),"Paired"))(n_cl),
      levels(key_long$Cluster))

    p_box <- ggplot(key_long, aes(x=Cluster, y=Expression, fill=Cluster)) +
      geom_boxplot(outlier.size=0.5, outlier.alpha=0.4,
                   linewidth=0.4, width=0.7) +
      facet_grid(Gene~., scales="free_y") +
      scale_fill_manual(values=cl_pal, guide="none") +
      labs(x=NULL, y="Normalized Expression",
           title="Key Gene Expression by Subtype") +
      theme_bw(base_size=11) +
      theme(axis.text.x  = element_text(size=8, angle=45, hjust=1, vjust=1),
            strip.text   = element_text(size=11, face="bold"),
            plot.title   = element_text(hjust=0.5, size=13),
            plot.margin  = ggplot2::margin(t=10, r=20, b=40, l=10),
            panel.grid.major.x = element_blank())

    save_plot(p_box, "Fig_cardio_keygene_boxplot",
              width  = max(12, n_cl * 1.4 + 2),
              height = 4 * length(key_avail) + 1.5,
              dir    = OUT_SUB)

    p_vln2 <- ggplot(key_long, aes(x=Cluster, y=Expression, fill=Cluster)) +
      geom_violin(scale="width", trim=TRUE, alpha=0.75) +
      geom_boxplot(width=0.12, fill="white", outlier.size=0.3,
                   linewidth=0.4, alpha=0.8) +
      facet_grid(Gene~., scales="free_y") +
      scale_fill_manual(values=cl_pal, guide="none") +
      labs(x=NULL, y="Normalized Expression",
           title="Key Gene Expression by Subtype") +
      theme_bw(base_size=11) +
      theme(axis.text.x  = element_text(size=8, angle=45, hjust=1, vjust=1),
            strip.text   = element_text(size=11, face="bold"),
            plot.title   = element_text(hjust=0.5, size=13),
            plot.margin  = ggplot2::margin(t=10, r=20, b=40, l=10),
            panel.grid.major.x = element_blank())

    save_plot(p_vln2, "Fig_cardio_keygene_violin",
              width  = max(12, n_cl * 1.4 + 2),
              height = 4 * length(key_avail) + 1.5,
              dir    = OUT_SUB)
  } else {
    cat(sprintf("[WARN] 关键基因 %s 均不在数据集中\n",
                paste(KEY_GENES, collapse=",")))
  }


  saveRDS(sc_cardio, file.path(OUT_SUB, "seurat_cardiomyocytes.rds"))
  cat(sprintf("[INFO] 保存: %s\n", file.path(OUT_SUB,"seurat_cardiomyocytes.rds")))

}
