rm(list = ls()); gc()


suppressPackageStartupMessages({
  library(rms)
  library(ggplot2)
  library(dplyr)
  library(pROC)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ResourceSelection)
})


SERVER_DIR <- ""
setwd(SERVER_DIR)

out_dir <- "06_nomogram"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


key_genes <- readLines("02_ml/key_genes.txt")
key_genes <- trimws(key_genes[nchar(trimws(key_genes)) > 0])
cat("[INFO] key_genes:", paste(key_genes, collapse = ", "), "\n")


map_expr_to_symbol_df <- function(expr_mat) {
  rnames <- head(rownames(expr_mat), 10)

  looks_like_symbol <- all(grepl("^[A-Z][A-Z0-9.-]{1,10}$", rnames))
  looks_like_accession <- any(grepl("^[A-Z]{2}[0-9]{5,}|^NM_|^NR_|^XM_|^XR_|^AB[0-9]", rnames))

  if (looks_like_symbol && !looks_like_accession) {
    cat("[INFO] 行名已为基因符号，直接转置使用\n")
    expr_df <- as.data.frame(t(expr_mat))

    sym_dup <- colnames(expr_df)
    if (any(duplicated(sym_dup))) {
      var_v <- apply(expr_mat, 1, var, na.rm = TRUE)
      keep  <- tapply(seq_along(sym_dup), sym_dup, function(i) i[which.max(var_v[i])])
      expr_df <- expr_df[, unlist(keep), drop = FALSE]
      cat(sprintf("[INFO] 去重后 %d 个基因\n", ncol(expr_df)))
    }
    return(expr_df)
  }

  cat("[INFO] 行名为登录号，用 org.Hs.eg.db 映射...\n")
  acc_clean <- sub("\\..*", "", rownames(expr_mat))
  mapping <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = unique(acc_clean),
    columns = c("REFSEQ", "SYMBOL"),
    keytype = "REFSEQ"
  ))
  mapping <- mapping[!is.na(mapping$SYMBOL), ]
  mapping <- mapping[!duplicated(mapping$REFSEQ), ]

  matched_idx <- match(acc_clean, mapping$REFSEQ)
  valid       <- !is.na(matched_idx)

  cat(sprintf("[INFO] 成功映射 %d / %d 个登录号\n", sum(valid), length(valid)))
  if (sum(valid) == 0) stop("org.Hs.eg.db 未能映射任何行名，请检查数据格式。")

  expr_mapped        <- as.data.frame(t(expr_mat[valid, , drop = FALSE]))
  colnames(expr_mapped) <- mapping$SYMBOL[matched_idx[valid]]


  sym_dup <- colnames(expr_mapped)
  if (any(duplicated(sym_dup))) {
    var_v    <- apply(expr_mat[valid, ], 1, var, na.rm = TRUE)
    keep_idx <- tapply(seq_along(sym_dup), sym_dup, function(i) i[which.max(var_v[i])])
    expr_mapped <- expr_mapped[, unlist(keep_idx), drop = FALSE]
    cat(sprintf("[INFO] 去重后基因数: %d\n", ncol(expr_mapped)))
  }
  return(expr_mapped)
}


build_model_df <- function(expr_mat, pheno_df,
                           genes, group_col, hf_pattern) {

  expr_sym <- map_expr_to_symbol_df(expr_mat)


  avail <- intersect(genes, colnames(expr_sym))
  missing_g <- setdiff(genes, colnames(expr_sym))
  if (length(missing_g) > 0)
    cat("[WARN] 以下基因未映射到，跳过:", paste(missing_g, collapse = ", "), "\n")
  if (length(avail) == 0) stop("所有目标基因均未映射，无法建模。")
  cat("[INFO] 可用基因:", paste(avail, collapse = ", "), "\n")

  df <- expr_sym[, avail, drop = FALSE]


  sample_names <- rownames(df)
  pheno_match  <- pheno_df[match(sample_names, rownames(pheno_df)), , drop = FALSE]
  group_raw    <- pheno_match[[group_col]]
  df$outcome   <- as.integer(grepl(hf_pattern, group_raw, ignore.case = TRUE))

  cat(sprintf("[INFO] HF=%d, Normal=%d, 样本总数=%d\n",
              sum(df$outcome), sum(df$outcome == 0), nrow(df)))
  return(df)
}




cat("[INFO] 加载训练集 GSE57338...\n")
load("00_rawdata/training_GSE57338/GSE57338_processed.RData")



df_train <- build_model_df(
  expr_mat  = expr_gene_53778,
  pheno_df  = pheno_53778,
  genes     = key_genes,
  group_col = "characteristics_ch1.1",
  hf_pattern = "yes"
)

avail_genes <- setdiff(colnames(df_train), "outcome")
cat("[INFO] 训练集可用基因:", paste(avail_genes, collapse = ", "), "\n")
cat("[INFO] 训练集:", nrow(df_train), "样本 | HF:", sum(df_train$outcome),
    "| Normal:", sum(df_train$outcome == 0), "\n")




cat("[INFO] 加载验证集 GSE5406...\n")
load("00_rawdata/validation_GSE5406/GSE5406_processed.RData")


df_valid <- build_model_df(
  expr_mat  = expr_gene_5406,
  pheno_df  = pheno_5406,
  genes     = avail_genes,
  group_col = "characteristics_ch1",
  hf_pattern = "heart failure"
)

avail_genes_valid <- setdiff(colnames(df_valid), "outcome")
common_genes <- intersect(avail_genes, avail_genes_valid)
if (length(common_genes) == 0) stop("训练集与验证集无共同可用基因，终止。")
if (length(common_genes) < length(avail_genes)) {
  cat("[WARN] 验证集缺失基因，移除:",
      paste(setdiff(avail_genes, common_genes), collapse = ", "), "\n")
}

avail_genes <- common_genes
df_train <- df_train[, c(avail_genes, "outcome")]
df_valid  <- df_valid[,  c(avail_genes, "outcome")]




dd <- datadist(df_train)
options(datadist = "dd")

formula_str <- paste("outcome ~", paste(avail_genes, collapse = " + "))
lrm_fit <- lrm(as.formula(formula_str), data = df_train, x = TRUE, y = TRUE)

cat("[INFO] LRM 模型摘要（训练集）:\n")
print(lrm_fit)

auc_train <- as.numeric(lrm_fit$stats["C"])
cat(sprintf("[STAT] Nomogram_AUC_train: %.3f\n", auc_train))

pred_valid <- predict(lrm_fit, newdata = df_valid, type = "fitted")
roc_valid  <- roc(df_valid$outcome, pred_valid, quiet = TRUE)
auc_valid  <- as.numeric(auc(roc_valid))
cat(sprintf("[STAT] Nomogram_AUC_valid: %.3f\n", auc_valid))
cat(sprintf("[STAT] Nomogram_key_genes: %s\n", paste(avail_genes, collapse = ",")))














prob_ticks <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95)
nom <- nomogram(lrm_fit,
                fun         = plogis,
                fun.at      = prob_ticks,
                funlabel    = "Heart Failure Probability",
                lp          = FALSE,
                maxscale    = 100)


png(file.path(out_dir, "nomogram.png"), width = 3200, height = 2800, res = 300)
pdf(file.path(out_dir, "nomogram.pdf"), width = 3200, height = 2800)


par(mar = c(10, 4, 4, 2) + 0.1)

plot(nom,
     xfrac    = 0.35,
     cex.axis = 0.6,
     cex.var  = 0.9,
     tcl      = -0.2,
     main     = "Nomogram for Heart Failure Prediction")


coef_all <- coef(lrm_fit)
coef_all <- coef_all[names(coef_all) != "Intercept"]
or_all   <- exp(coef_all)
coef_text <- paste(sprintf("%s: beta=%.3f, OR=%.3f", names(coef_all), coef_all, or_all), collapse = ";  ")


mtext(coef_text, side = 1, line = 7, adj = 0.5, cex = 0.8)

dev.off()
cat("[INFO] 已保存nomogram.png\n")


mar_settings <- c(6, 6, 4, 3) + 0.1


pred_prob_train <- predict(lrm_fit, type = "fitted")
pred_prob_valid <- predict(lrm_fit, newdata = df_valid, type = "fitted")

hl_train <- ResourceSelection::hoslem.test(df_train$outcome, pred_prob_train, g = 10)
hl_valid <- ResourceSelection::hoslem.test(df_valid$outcome, pred_prob_valid, g = 10)

cat(sprintf("[STAT] HL_train: X2 = %.3f, df = %d, P = %.4g\n",
            as.numeric(hl_train$statistic),
            as.integer(hl_train$parameter),
            hl_train$p.value))
cat(sprintf("[STAT] HL_valid: X2 = %.3f, df = %d, P = %.4g\n",
            as.numeric(hl_valid$statistic),
            as.integer(hl_valid$parameter),
            hl_valid$p.value))


set.seed(42)
cal <- calibrate(lrm_fit, method = "boot", B = 500)
















hl_text <- sprintf("H-L test (train) P = %.3f", hl_train$p.value)

png(file.path(out_dir, "calibration.png"), width = 2800, height = 2400, res = 300)
par(mar = mar_settings)

plot(cal,
     xlab      = "Predicted Probability",
     ylab      = "Observed Probability",
     main      = "Bootstrap Calibration Curve (B=500)",
     subtitles = FALSE,
     lwd       = 2,
     xlim      = c(0, 1),
     ylim      = c(0, 1),
     legend    = FALSE)

abline(0, 1, lty = 2, col = "gray40", lwd = 1.5)






text(0.05, 0.95, labels = hl_text, adj = c(0, 1), cex = 0.9)

legend("bottomright",
       legend = c("Apparent", "Bias-corrected (Train)", "Ideal"),
       lty    = c(2, 1, 2),
       lwd    = c(1.5, 2, 1.5),
       col    = c("black", "black", "gray40"),
       bty    = "n", cex = 0.9)
dev.off()


pdf(file.path(out_dir, "calibration.pdf"), width = 7, height = 6)
par(mar = mar_settings)

plot(cal,
     xlab      = "Predicted Probability",
     ylab      = "Observed Probability",
     main      = "Bootstrap Calibration Curve (B=500)",
     subtitles = FALSE,
     lwd       = 2,
     xlim      = c(0, 1),
     ylim      = c(0, 1),
     legend    = FALSE)

abline(0, 1, lty = 2, col = "gray40", lwd = 1.5)






text(0.05, 0.95, labels = hl_text, adj = c(0, 1), cex = 0.9)

legend("bottomright",
       legend = c("Apparent", "Bias-corrected (Train)", "Ideal"),
       lty    = c(2, 1, 2),
       lwd    = c(1.5, 2, 1.5),
       col    = c("black", "black", "gray40"),
       bty    = "n", cex = 0.9)
dev.off()

cat("[INFO] 已保存校准曲线 PNG 和 PDF\n")


pred_prob_train <- predict(lrm_fit, type = "fitted")
thresholds <- seq(0, 0.99, by = 0.01)
n    <- nrow(df_train)
y    <- df_train$outcome
prev <- mean(y)

calc_nb <- function(thresh, probs, labels) {
  tp <- sum(probs >= thresh & labels == 1)
  fp <- sum(probs >= thresh & labels == 0)
  tp / n - fp / n * (thresh / (1 - thresh + 1e-9))
}

nb_model <- sapply(thresholds, calc_nb, probs = pred_prob_train, labels = y)
nb_all   <- prev - (1 - prev) * thresholds / (1 - thresholds + 1e-9)
nb_none  <- rep(0, length(thresholds))

dca_df <- data.frame(
  threshold   = rep(thresholds, 3),
  net_benefit = c(nb_model, nb_all, nb_none),
  Strategy    = rep(c("Nomogram", "Treat All", "Treat None"),
                    each = length(thresholds))
)
dca_df$net_benefit <- pmax(dca_df$net_benefit, -0.05)
single_gene_nb_list <- lapply(avail_genes, function(g) {
  probs_g <- predict(
    lrm(as.formula(paste("outcome ~", g)), data = df_train),
    type = "fitted"
  )
  nb_g <- sapply(thresholds, calc_nb, probs = probs_g, labels = y)
  data.frame(
    threshold   = thresholds,
    net_benefit = pmax(nb_g, -0.05),
    Strategy    = g
  )
})
single_gene_nb_df <- do.call(rbind, single_gene_nb_list)
dca_df <- rbind(dca_df, single_gene_nb_df)

p_dca <- ggplot(dca_df, aes(x = threshold, y = net_benefit,
                            color = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c("Nomogram"   = "blue",
                                "Treat All"  = "red",
                                "Treat None" = "black",
                                setNames(
                                  c("darkorange", "purple", "darkgreen"),
                                  avail_genes
                                ))) +
  scale_linetype_manual(values = c("Nomogram"   = "solid",
                                   "Treat All"  = "dashed",
                                   "Treat None" = "dotted",
                                   setNames(
                                     rep("solid", length(avail_genes)),
                                     avail_genes
                                   ))) +
  coord_cartesian(xlim = c(0, 1),
                  ylim = c(-0.05, max(nb_model, na.rm = TRUE) + 0.05)) +
  labs(title = "Decision Curve Analysis",
       x     = "Threshold Probability",
       y     = "Net Benefit") +
  theme_classic(base_size = 13) +
  theme(legend.position  = "right",
        plot.title       = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(out_dir, "DCA_curve.png"), plot = p_dca,
       width = 7, height = 6, dpi = 300)
ggsave(file.path(out_dir, "DCA_curve.pdf"), plot = p_dca,
       width = 7, height = 6)
cat("[INFO] 已保存 DCA_curve.png\n")






































































































library(pROC)

pred_train_roc <- predict(lrm_fit, type = "fitted")
roc_train <- roc(df_train$outcome, pred_train_roc, quiet = TRUE)
ci_train  <- ci.auc(roc_train, conf.level = 0.95)

n_train <- nrow(df_train)
label_train <- sprintf("Training (GSE57338, n=%d)\n(AUC=%.3f, 95%%CI: %.3f-%.3f)",
                       n_train,
                       as.numeric(auc(roc_train)),
                       ci_train[1], ci_train[3])

make_roc_df <- function(roc_obj, label) {
  data.frame(
    FPR   = 1 - roc_obj$specificities,
    TPR   = roc_obj$sensitivities,
    group = label
  )
}

roc_df <- make_roc_df(roc_train, label_train)
roc_df$group <- factor(roc_df$group, levels = label_train)


gene_subtitle <- paste("Predictors:", paste(avail_genes, collapse = " + "),
                       "| Training: GSE57338")

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR,
                            color = group, linetype = group)) +
  geom_line(linewidth = 1.0) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50", linewidth = 0.6) +
  scale_color_manual(values = c("red")) +
  scale_linetype_manual(values = c("solid")) +
  scale_x_continuous(breaks = seq(0, 1, 0.25),
                     labels = c("0.00","0.25","0.50","0.75","1.00")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title    = "ROC Curve for Nomogram Prediction Model",
       subtitle = gene_subtitle,
       x        = "1 - Specificity  (False Positive Rate)",
       y        = "Sensitivity (True Positive Rate)",
       color    = NULL,
       linetype = NULL) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
    legend.position      = c(0.72, 0.25),
    legend.background    = element_rect(fill = "white", color = "grey70"),
    legend.key.width     = unit(1.5, "cm"),
    legend.text          = element_text(size = 9),
    axis.title           = element_text(size = 12)
  )

ggsave(file.path(out_dir, "nomogram_ROC.png"), plot = p_roc,
       width = 7, height = 7, dpi = 300)
ggsave(file.path(out_dir, "nomogram_ROC.pdf"), plot = p_roc,
       width = 7, height = 7)
cat("[INFO] 已保存 nomogram_ROC.png 和 nomogram_ROC.pdf\n")

coef_df <- data.frame(
  Gene        = names(coef(lrm_fit)),
  Coefficient = as.numeric(coef(lrm_fit)),
  OR          = exp(as.numeric(coef(lrm_fit)))
)
write.csv(coef_df, file.path(out_dir, "nomogram_coef.csv"),
          row.names = FALSE, quote = FALSE)
cat("[INFO] 已保存 nomogram_coef.csv\n")

cat("[DONE] 06_nomogram 全部完成。\n")
