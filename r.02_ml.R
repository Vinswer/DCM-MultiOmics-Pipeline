rm(list = ls()); gc()

library(GEOquery)
library(limma)
library(e1071)
library(caret)
library(glmnet)
library(randomForest)
library(VennDiagram)
library(pROC)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(grid)
library(RColorBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)

ORIGINAL_DIR <- "/data/nas1/chengyuzhen_OD/project/02_project_1624"
OUTPUT_DIR <- file.path(ORIGINAL_DIR, "02_ml")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

candidate_file <- file.path(ORIGINAL_DIR, "project_1624_results_20260417/01_deg/Candidate_genes_DEG_x_Phagocytosis.csv")
train_geo_file <- file.path(ORIGINAL_DIR, "project_1624_results_20260417/00_rawdata/training_GSE57338/GSE57338_series_matrix.txt.gz")
valid_geo_file <- file.path(ORIGINAL_DIR, "project_1624_results_20260417/00_rawdata/validation_GSE5406/GSE5406_series_matrix.txt.gz")

candidate_df <- read.csv(candidate_file, stringsAsFactors = FALSE)
candidate_genes <- unique(trimws(candidate_df[, 1]))
cat("Candidate genes:", paste(candidate_genes, collapse = ", "), "\n")

train_gse <- getGEO(filename = train_geo_file, GSEMatrix = TRUE, getGPL = FALSE)
valid_gse  <- getGEO(filename = valid_geo_file, GSEMatrix = TRUE, getGPL = FALSE)

train_pheno <- pData(train_gse)
valid_pheno  <- pData(valid_gse)


get_labels_train <- function(pheno) {
  hf_col <- grep("heart.failure", colnames(pheno), ignore.case = TRUE, value = TRUE)
  if (length(hf_col) > 0) {
    vals <- tolower(trimws(as.character(pheno[[hf_col[1]]])))
    vals <- gsub(".*:\\s*", "", vals)
    cat("Train HF col values (sample):", head(vals, 5), "\n")
    if (any(vals == "yes") && any(vals == "no")) {
      return(ifelse(vals == "yes", "HF", "Control"))
    }
  }
  char_cols <- grep("characteristics", colnames(pheno), ignore.case = TRUE, value = TRUE)
  for (col in char_cols) {
    vals <- tolower(trimws(as.character(pheno[[col]])))
    vals_clean <- gsub(".*:\\s*", "", vals)
    if (any(vals_clean == "yes") && any(vals_clean == "no")) {
      return(ifelse(vals_clean == "yes", "HF", "Control"))
    }
  }
  all_text <- apply(pheno, 1, function(x) paste(tolower(x), collapse = " "))
  return(ifelse(grepl("non.failing|normal heart|healthy|donor", all_text), "Control", "HF"))
}

get_labels_valid <- function(pheno) {
  char_cols <- grep("characteristics", colnames(pheno), ignore.case = TRUE, value = TRUE)
  for (col in char_cols) {
    vals <- tolower(trimws(as.character(pheno[[col]])))
    if (any(grepl("heart failure|cardiomyopathy|failing", vals)) &&
        any(grepl("non.fail|normal|healthy|donor|unused", vals))) {
      return(ifelse(grepl("non.fail|normal|healthy|donor|unused", vals), "Control", "HF"))
    }
  }
  all_text <- apply(pheno, 1, function(x) paste(tolower(x), collapse = " "))
  cat("Valid label text sample:\n"); print(head(all_text, 3))
  return(ifelse(grepl("non.fail|normal|healthy|donor|unused", all_text), "Control", "HF"))
}

train_labels <- get_labels_train(train_pheno)
valid_labels  <- get_labels_valid(valid_pheno)
cat("Train label table:\n"); print(table(train_labels))
cat("Valid label table:\n");  print(table(valid_labels))

train_expr_file <- file.path(ORIGINAL_DIR,
                             "project_1624_results_20260417/00_rawdata/training_GSE57338/GSE57338_gene_expression.csv")
train_expr_full <- read.csv(train_expr_file, row.names = 1, check.names = FALSE)

train_expr_sym <- train_expr_full
train_key_genes <- intersect(candidate_genes, rownames(train_expr_full))

train_mat <- t(train_expr_full[train_key_genes, , drop = FALSE])


sample_names <- colnames(train_expr_full)
names(train_labels) <- colnames(train_gse)


common_samples <- intersect(sample_names, names(train_labels))
cat("共同样本数:", length(common_samples), "\n")

train_mat    <- train_mat[common_samples, , drop = FALSE]
train_labels <- train_labels[common_samples]

train_genes <- colnames(train_mat)
cat("Train genes for ML:", train_genes, "\n")
stopifnot(length(train_genes) > 0)

train_df <- data.frame(train_mat, label = factor(train_labels), check.names = FALSE)


save_plot <- function(filename_base, plot_fn, width = 6, height = 5) {
  pdf(paste0(filename_base, ".pdf"), width = width, height = height)
  plot_fn()
  dev.off()
  png(paste0(filename_base, ".png"), width = width * 100, height = height * 100, res = 100)
  plot_fn()
  dev.off()
}






library(xgboost)

set.seed(42)
x_train_xgb <- as.matrix(train_df[, train_genes])
y_train_xgb <- ifelse(train_df$label == "HF", 1, 0)

xgb_model <- xgboost(
  data        = x_train_xgb,
  label       = y_train_xgb,
  max_depth   = 2,
  eta         = 0.1,
  nrounds     = 25,
  objective   = "binary:logistic",
  eval_metric = "error",
  verbose     = 0
)

xgb_imp <- xgb.importance(feature_names = train_genes, model = xgb_model)
cat("XGBoost feature importance:\n"); print(xgb_imp)



xgb_genes_final <- as.character(xgb_imp$Feature[xgb_imp$Gain > 0])
cat("XGBoost selected genes (rel_gain > 0.05):", xgb_genes_final, "\n")



xgb_plot_fn <- function() {
  imp_sorted <- xgb_imp[order(xgb_imp$Gain, decreasing = FALSE), ]
  imp_sorted$Feature <- factor(imp_sorted$Feature, levels = imp_sorted$Feature)

  p <- ggplot(imp_sorted, aes(x = Feature, y = Gain)) +
    geom_bar(stat = "identity", fill = "#E41A1C", width = 0.6) +
    coord_flip() +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text        = element_text(size = 10)
    ) +
    labs(title = "XGBoost Feature Importance",
         x     = "",
         y     = "Gain")
  print(p)
}
save_plot(file.path(OUTPUT_DIR, "01_XGBoost"), xgb_plot_fn, width = 7, height = 5)


x_train <- as.matrix(train_df[, train_genes])
y_train <- ifelse(train_df$label == "HF", 1, 0)

set.seed(42)
lasso_cv <- cv.glmnet(
  x       = x_train,
  y       = y_train,
  family  = "binomial",
  alpha   = 1,
  nfolds  = 10
)



















lambda_strict <- lasso_cv$lambda.1se * 1.5
lasso_plot_fn <- function() {
  par(mfrow = c(1, 2), oma = c(0, 0, 0, 0))


  plot(lasso_cv, main = "")
  mtext("LASSO 10-fold CV", side = 3, line = 3, cex = 1.1, font = 2)
  abline(v = log(lasso_cv$lambda.1se), lty = 2, col = "grey50")
  abline(v = log(lambda_strict),       lty = 2, col = "red")
  text(x = log(lambda_strict), y = par("usr")[4] * 0.95,
       labels = expression(paste("1.5x ", lambda["1SE"])),
       col = "red", adj = c(1.1, 0), cex = 0.9)


  plot(lasso_cv$glmnet.fit, xvar = "lambda", label = FALSE, main = "")
  mtext("Coefficient Path", side = 3, line = 2.5, cex = 1.1, font = 2)

  abline(v = -log(lambda_strict), lty = 2, col = "red")
  text(x   = -log(lambda_strict),
       y   = par("usr")[3] + diff(par("usr")[3:4]) * 0.15,
       labels = expression(paste("1.5x ", lambda["1SE"])),
       col = "red", adj = c(1.1, 0), cex = 0.9)



  all_coef   <- as.matrix(coef(lasso_cv$glmnet.fit))
  gene_rows  <- rownames(all_coef)[rownames(all_coef) != "(Intercept)"]
  ever_nonzero <- gene_rows[apply(all_coef[gene_rows, , drop = FALSE], 1,
                                  function(r) any(r != 0))]


  glmnet_colors <- c("#BF0000", "#0000BF", "#00BF00", "#BFBF00",
                     "#BF00BF", "#00BFBF", "gray40", "orange",
                     "brown",   "pink")
  n_nonzero <- length(ever_nonzero)
  leg_cols  <- rep_len(glmnet_colors, n_nonzero)

  legend("topright",
         legend  = ever_nonzero,
         col     = leg_cols,
         lwd     = 1.8,
         cex     = 0.72,
         bty     = "n",
         title   = "Gene",
         title.col = "black",
         inset   = c(0.01, 0.01))


  par(mfrow = c(1, 1))
}

save_plot(file.path(OUTPUT_DIR, "02_LASSO"), lasso_plot_fn, width = 14, height = 5)

lambda_strict <- lasso_cv$lambda.1se
lasso_coef    <- coef(lasso_cv, s = lambda_strict)
lasso_genes   <- rownames(lasso_coef)[which(abs(as.numeric(lasso_coef)) > 0)]
lasso_genes   <- lasso_genes[lasso_genes != "(Intercept)"]
lasso_genes   <- lasso_genes[lasso_genes %in% train_genes]

cat("LASSO selected genes (1.5x lambda.1se):", lasso_genes, "\n")
cat("Lambda used:", lambda_strict, "| lambda.1se was:", lasso_cv$lambda.1se, "\n")


set.seed(42)
rf_model <- randomForest(
  x         = train_df[, train_genes],
  y         = train_df$label,
  ntree     = 1000,
  mtry      = max(1, floor(sqrt(length(train_genes)))),
  importance = TRUE
)

cat("Random Forest OOB error rate:", rf_model$err.rate[1000, "OOB"], "\n")

rf_plot_fn <- function() {
  err_data  <- rf_model$err.rate
  trees     <- seq_len(nrow(err_data))
  col_names <- colnames(err_data)
  hf_col    <- col_names[grep("HF|heart|fail",           col_names, ignore.case = TRUE)[1]]
  ctrl_col  <- col_names[grep("Control|normal|healthy",  col_names, ignore.case = TRUE)[1]]
  if (is.na(hf_col))   hf_col   <- col_names[2]
  if (is.na(ctrl_col)) ctrl_col <- col_names[3]

  plot(trees, err_data[, "OOB"], type = "l", col = "black", lwd = 2,
       xlab = "Number of Trees", ylab = "OOB Error Rate",
       ylim = c(0, min(1, max(err_data, na.rm = TRUE) + 0.05)),
       main = "Random Forest OOB Error")
  lines(trees, err_data[, hf_col],   col = "darkred",   lwd = 1.5, lty = 2)
  lines(trees, err_data[, ctrl_col], col = "darkgreen", lwd = 1.5, lty = 3)
  legend("topright",
         legend = c("OOB", "HF", "Control"),
         col    = c("black", "darkred", "darkgreen"),
         lty    = c(1, 2, 3), lwd = 1.5, bty = "n", cex = 0.9)
}
save_plot(file.path(OUTPUT_DIR, "03_RandomForest"), rf_plot_fn, width = 5, height = 5)


rf_imp     <- importance(rf_model)
rf_imp_df  <- data.frame(
  gene             = rownames(rf_imp),
  MeanDecreaseGini = rf_imp[, "MeanDecreaseGini"]
)
rf_imp_df$rel_imp <- rf_imp_df$MeanDecreaseGini / max(rf_imp_df$MeanDecreaseGini)
rf_genes          <- as.character(rf_imp_df$gene[rf_imp_df$rel_imp > 0.25])

rf_imp_plot_fn <- function() {
  imp_sorted <- rf_imp_df[order(rf_imp_df$MeanDecreaseGini, decreasing = FALSE), ]
  imp_sorted$gene    <- factor(imp_sorted$gene, levels = imp_sorted$gene)
  imp_sorted$selected <- ifelse(imp_sorted$rel_imp > 0.25, "Selected", "Not Selected")

  p <- ggplot(imp_sorted, aes(x = gene, y = MeanDecreaseGini, fill = selected)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = 0.25 * max(rf_imp_df$MeanDecreaseGini),
               linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text",
             x    = 1,
             y    = 0.25 * max(rf_imp_df$MeanDecreaseGini),
             label = "threshold = 0.25",
             color = "red", vjust = -0.5, hjust = 0, size = 3.5) +
    coord_flip() +
    scale_fill_manual(values = c("Selected" = "#E41A1C", "Not Selected" = "grey70")) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text        = element_text(size = 10),
      legend.position  = "top",
      legend.title     = element_blank()
    ) +
    labs(title = "Random Forest Feature Importance",
         x     = "",
         y     = "Mean Decrease Gini")
  print(p)
}
save_plot(file.path(OUTPUT_DIR, "03b_RF_Importance"), rf_imp_plot_fn, width = 7, height = 5)
cat("RF selected genes (relative importance > 0.25):", rf_genes, "\n")
cat("RF all importance:\n"); print(rf_imp_df[order(-rf_imp_df$rel_imp), ])




venn_list <- list(
  XGBoost      = xgb_genes_final,
  LASSO        = lasso_genes,
  RandomForest = rf_genes
)

draw_venn <- function() {
  grid.newpage()
  vp <- venn.diagram(
    x              = venn_list,
    filename       = NULL,
    category.names = c("XGBoost", "LASSO", "Random Forest"),
    col            = c("#E41A1C", "#377EB8", "#4DAF4A"),
    fill           = c(alpha("#E41A1C", 0.3), alpha("#377EB8", 0.3), alpha("#4DAF4A", 0.3)),
    cex            = 1.2, cat.cex = 1.0,
    cat.col        = c("#E41A1C", "#377EB8", "#4DAF4A"),
    margin         = 0.1
  )
  grid.draw(vp)
}
pdf(file.path(OUTPUT_DIR, "04_Venn_diagram.pdf"), width = 6, height = 6)
draw_venn()
dev.off()
png(file.path(OUTPUT_DIR, "04_Venn_diagram.png"), width = 600, height = 600, res = 100)
draw_venn()
dev.off()




key_genes <- Reduce(intersect, venn_list)
cat("Key genes (3-way intersection):", key_genes, "\n")

if (length(key_genes) == 0) {
  cat("No 3-way intersection. Using pairwise intersections.\n")
  key_genes <- unique(c(
    intersect(xgb_genes_final, lasso_genes),
    intersect(xgb_genes_final, rf_genes),
    intersect(lasso_genes,     rf_genes)
  ))
  cat("Pairwise intersection genes:", key_genes, "\n")
}
if (length(key_genes) == 0) {
  cat("No pairwise intersection. Using LASSO genes as fallback.\n")
  key_genes <- if (length(lasso_genes) > 0) lasso_genes else train_genes
}

key_genes <- key_genes[key_genes %in% train_genes]
writeLines(key_genes, file.path(OUTPUT_DIR, "key_genes.txt"))
cat("Final key genes:", paste(key_genes, collapse = ", "), "\n")


valid_expr_full <- read.csv(
  file.path(ORIGINAL_DIR, "project_1624_results_20260417/00_rawdata/validation_GSE5406/GSE5406_gene_expression.csv"),
  row.names = 1, check.names = FALSE
)


valid_key_genes <- intersect(key_genes, rownames(valid_expr_full))
cat("验证集中找到 key genes:", valid_key_genes, "\n")

valid_key <- t(valid_expr_full[valid_key_genes, , drop = FALSE])






train_key <- train_mat[, key_genes, drop = FALSE]


if (length(valid_key_genes) < length(key_genes)) {
  cat("Warning: the following key genes are absent in validation set and will be skipped for validation:\n")
  cat(setdiff(key_genes, valid_key_genes), "\n")
}
stopifnot(length(valid_key_genes) > 0)


plot_roc_multi <- function(expr_mat, labels, title_str, filename_base) {
  label_bin <- ifelse(labels == "HF", 1, 0)
  n_genes   <- ncol(expr_mat)
  colors    <- if (n_genes <= 8) {
    brewer.pal(max(3, n_genes), "Set1")[1:n_genes]
  } else {
    rainbow(n_genes)
  }

  roc_list <- lapply(seq_len(n_genes), function(i) {
    tryCatch(roc(label_bin, expr_mat[, i], quiet = TRUE), error = function(e) NULL)
  })
  names(roc_list) <- colnames(expr_mat)
  roc_list        <- Filter(Negate(is.null), roc_list)

  auc_vals    <- sapply(roc_list, function(r) round(as.numeric(auc(r)), 3))
  legend_text <- paste0(names(roc_list), " (AUC=", auc_vals, ")")

  do_plot <- function() {
    plot(roc_list[[1]], col = colors[1], lwd = 2, main = title_str,
         xlab = "1 - Specificity", ylab = "Sensitivity", legacy.axes = TRUE)
    if (length(roc_list) > 1) {
      for (i in 2:length(roc_list)) {
        lines(roc_list[[i]], col = colors[i], lwd = 2)
      }
    }
    abline(a = 0, b = 1, lty = 2, col = "gray60")
    legend("bottomright", legend = legend_text,
           col = colors[seq_along(roc_list)],
           lwd = 2, cex = 0.8, bty = "n")
  }

  pdf(paste0(filename_base, ".pdf"), width = 7, height = 6)
  do_plot()
  dev.off()
  png(paste0(filename_base, ".png"), width = 700, height = 600, res = 100)
  do_plot()
  dev.off()

  return(auc_vals)
}

auc_train <- plot_roc_multi(train_key, train_labels,
                            "Training Set (GSE57338) ROC",
                            file.path(OUTPUT_DIR, "05_ROC_training"))
auc_valid  <- plot_roc_multi(valid_key,  valid_labels,
                             "Validation Set (GSE5406) ROC",
                             file.path(OUTPUT_DIR, "06_ROC_validation"))

cat("Training AUC:\n");   print(auc_train)
cat("Validation AUC:\n"); print(auc_valid)




plot_boxplots <- function(expr_mat, labels, dataset_name, filename_base) {
  df_long       <- data.frame(expr_mat, check.names = FALSE)
  df_long$label <- factor(labels, levels = c("Control", "HF"))
  df_melt       <- reshape2::melt(df_long, id.vars = "label",
                                  variable.name = "Gene", value.name = "Expression")

  n_genes      <- ncol(expr_mat)
  ncols_facet  <- min(4, n_genes)
  nrows_facet  <- ceiling(n_genes / ncols_facet)

  p <- ggplot(df_melt, aes(x = label, y = Expression, fill = label)) +
    geom_boxplot(outlier.size = 0.8, width = 0.6, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
    facet_wrap(~ Gene, scales = "free_y", ncol = ncols_facet) +
    scale_fill_manual(values = c("Control" = "#4DAF4A", "HF" = "#E41A1C")) +
    stat_compare_means(method = "wilcox.test",
                       comparisons = list(c("Control", "HF")),
                       label = "p.signif", tip.length = 0.01) +
    theme_bw(base_size = 11) +
    theme(
      strip.text       = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "grey90"),
      axis.text.x      = element_text(size = 9),
      legend.position  = "none",
      panel.grid.minor = element_blank()
    ) +
    labs(title = paste0(dataset_name, " - Key Gene Expression"),
         x = "", y = "Expression Level")

  plot_w <- ncols_facet * 2.8
  plot_h <- nrows_facet * 3.0 + 0.5

  ggsave(paste0(filename_base, ".pdf"), plot = p, width = plot_w, height = plot_h)
  ggsave(paste0(filename_base, ".png"), plot = p, width = plot_w, height = plot_h, dpi = 100)
}

plot_boxplots(train_key, train_labels, "Training Set (GSE57338)",
              file.path(OUTPUT_DIR, "07_Boxplot_training"))
plot_boxplots(valid_key, valid_labels,  "Validation Set (GSE5406)",
              file.path(OUTPUT_DIR, "08_Boxplot_validation"))

cat("\n=== Analysis Complete ===\n")
cat("Key genes:", paste(key_genes, collapse = ", "), "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
