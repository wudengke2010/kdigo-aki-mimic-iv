###################################################
# AKI Paradox — R重绘全部统计图表 (v4)
# UNIFIED with Python statistics — all numbers verified
# 输出: submission_figures_R/ 目录下 PDF + PNG (300 DPI)
###################################################

suppressWarnings(suppressMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(survival)
  library(survminer)
  library(rms)
  library(forestplot)
  library(cmprsk)
  library(pROC)
}))

# ---- Global theme ----
theme_pub <- theme_bw(base_size = 12, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = "grey50"),
    strip.text = element_text(face = "bold", size = 11)
  )
theme_set(theme_pub)

# ---- Colors ----
C_PROT  <- "#2166AC"   # blue = protective
C_HARM  <- "#B2182B"   # red  = harmful
C_NEUT  <- "#4D4D4D"   # grey
C_AKI   <- "#E08214"   # orange
C_CKD   <- "#8073AC"   # purple
C_STAGE1 <- "#66C2A5"
C_STAGE2 <- "#FC8D62"
C_STAGE3 <- "#8DA0CB"

# ---- Load data ----
dat <- fread("aki_paradox_cohort_with_sofa.csv")
cat("Total rows:", nrow(dat), "\n")

# AKI cohort only
aki <- dat[max_aki_stage > 0]
aki[, `:=`(
  aki_stage = factor(max_aki_stage, levels = 1:3, labels = c("Stage 1","Stage 2","Stage 3")),
  aki_stage_num = max_aki_stage,
  ckd_f = factor(ckd, levels = c(0,1), labels = c("Pure AKI","AKI+CKD")),
  male = as.integer(gender == "M"),
  died = hospital_expire_flag,
  stage2 = as.integer(max_aki_stage == 2),
  stage3 = as.integer(max_aki_stage == 3)
)]
aki[, ckd_s2 := ckd * stage2]
aki[, ckd_s3 := ckd * stage3]
cat("AKI cohort:", nrow(aki), "\n")

# Survival data
sv <- fread("aki_survival_data.csv")
sv[, `:=`(
  aki_stage = factor(max_aki_stage, levels = 1:3, labels = c("Stage 1","Stage 2","Stage 3")),
  group = factor(ifelse(ckd == 1, "AKI+CKD", "Pure AKI"), levels = c("Pure AKI","AKI+CKD")),
  died = hospital_expire_flag
)]

# Output dir
outdir <- "submission_figures_R"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# VERIFIED STATISTICS (from Python compute_all_stats_v4.py)
# These are the SINGLE SOURCE OF TRUTH for all annotations
# ============================================================
# Unadjusted ORs (by stage, CKD vs Pure)
V_OR <- list(
  s1_unadj = list(or=1.323, lo=1.183, hi=1.480, p=0.0000),
  s2_unadj = list(or=1.231, lo=1.060, hi=1.428, p=0.0063),
  s3_unadj = list(or=0.755, lo=0.646, hi=0.883, p=0.0004),
  s1_sofa  = list(or=1.208, lo=1.064, hi=1.372, p=0.0035),
  s2_sofa  = list(or=1.107, lo=0.941, hi=1.303, p=0.2211),
  s3_sofa  = list(or=0.643, lo=0.537, hi=0.769, p=0.0000),
  # PSM without SOFA (full cohort, caliper=0.2, 5042 pairs)
  s1_psm   = list(or=1.232, lo=1.063, hi=1.429, p=0.0056),
  s2_psm   = list(or=1.241, lo=1.029, hi=1.495, p=0.0224),
  s3_psm   = list(or=0.623, lo=0.499, hi=0.779, p=0.0000),
  # PSM with SOFA (full cohort, caliper=0.2, 5072 pairs)
  s1_psm_s = list(or=1.084, lo=0.937, hi=1.254, p=0.2770),
  s2_psm_s = list(or=1.126, lo=0.936, hi=1.355, p=0.2064),
  s3_psm_s = list(or=0.587, lo=0.473, hi=0.728, p=0.0000)
)

# AUC values (verified from Python)
V_AUC <- c(A=0.717, B=0.717, C=0.718, D=0.747)
V_AUC_LO <- c(A=0.709, B=0.709, C=0.710, D=0.738)
V_AUC_HI <- c(A=0.725, B=0.726, C=0.727, D=0.755)

# SOFA stratified (cutoff <=5, verified)
V_SOFA <- list(
  low  = list(n=483, or=0.740, lo=0.480, hi=1.140, p=0.172, pure_mort=29.3, ckd_mort=23.5),
  med  = list(n=869, or=0.957, lo=0.713, hi=1.284, p=0.770, pure_mort=35.1, ckd_mort=34.1),
  high = list(n=1666, or=0.679, lo=0.552, hi=0.836, p=0.0003, pure_mort=49.2, ckd_mort=39.7)
)

# ============================================================
# FIGURE 2: Forest Plot — 3 stages x 4 methods
# Uses VERIFIED ORs from Python (hardcoded for consistency)
# ============================================================
cat("\n=== Figure 2: Forest Plot ===\n")

fd <- data.table(
  stage = rep(c("Stage 1","Stage 2","Stage 3"), each=4),
  method_short = rep(c("Unadjusted","SOFA-Adjusted","PSM","PSM+SOFA"), 3),
  method = rep(c("Unadjusted","SOFA-\nAdjusted","PSM","PSM+\nSOFA"), 3),
  or = c(1.323, 1.208, 1.232, 1.084,
         1.231, 1.107, 1.241, 1.126,
         0.755, 0.643, 0.623, 0.587),
  lo = c(1.183, 1.064, 1.063, 0.937,
         1.060, 0.941, 1.029, 0.936,
         0.646, 0.537, 0.499, 0.473),
  hi = c(1.480, 1.372, 1.429, 1.254,
         1.428, 1.303, 1.495, 1.355,
         0.883, 0.769, 0.779, 0.728),
  p = c(0.0000, 0.0035, 0.0056, 0.2770,
        0.0063, 0.2211, 0.0224, 0.2064,
        0.0004, 0.0000, 0.0000, 0.0000)
)

fd[, stage := factor(stage, levels = c("Stage 1","Stage 2","Stage 3"))]
fd[, method_f := factor(method_short, levels = c("Unadjusted","SOFA-Adjusted","PSM","PSM+SOFA"))]
fd[, label := sprintf("%.2f (%.2f-%.2f)", or, lo, hi)]
fd[, sig := ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))]

# Shapes per method
method_shapes <- c("Unadjusted"=1, "SOFA-Adjusted"=15, "PSM"=17, "PSM+SOFA"=18)
method_colors <- c("Unadjusted"="#999999", "SOFA-Adjusted"="#2166AC", "PSM"="#E08214", "PSM+SOFA"="#B2182B")

p2 <- ggplot(fd, aes(x = or, y = method_f, color = method_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.6) +
  geom_point(aes(shape = method_f), size = 3.5) +
  geom_text(aes(x = max(hi) + 0.15, label = label), hjust = 0, size = 3.2, color = "black", family = "sans") +
  facet_wrap(~ stage, ncol = 1, scales = "free_y") +
  scale_x_log10(limits = c(0.3, 2.8), breaks = c(0.5, 1, 1.5, 2)) +
  scale_shape_manual(values = method_shapes) +
  scale_color_manual(values = method_colors) +
  labs(
    title = "Figure 2. Forest Plot: AKI+CKD vs Pure AKI - OR for Hospital Mortality",
    x = "Odds Ratio (95% CI), log scale",
    y = "",
    caption = "KEY FINDING: Stage 3 AKI paradox STRENGTHENS after SOFA adjustment (OR 0.76 -> 0.64) and PSM+SOFA (OR 0.59)"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    plot.caption = element_text(face = "bold", color = "#B2182B", size = 10, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(15, 120, 15, 10)
  )

ggsave(file.path(outdir, "Figure_2_Forest_Plot.pdf"), p2, width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(outdir, "Figure_2_Forest_Plot.png"), p2, width = 10, height = 8, dpi = 300)
cat("  Saved Figure 2\n")

# ============================================================
# FIGURE 3: SOFA Distribution (A: boxplot, B: stratified forest)
# ============================================================
cat("\n=== Figure 3: SOFA Distribution ===\n")

# Panel A: SOFA boxplot by AKI stage x CKD
sofa_plot <- aki[, .(sofa_total, aki_stage, ckd_f)]
p3a <- ggplot(sofa_plot, aes(x = aki_stage, y = sofa_total, fill = ckd_f)) +
  geom_violin(outlier.size = 0.3, outlier.alpha = 0.2, width = 0.6, linewidth = 0.4, alpha = 0.7) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, width = 0.3, linewidth = 0.4, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Pure AKI" = C_AKI, "AKI+CKD" = C_CKD)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white", color = "black",
               position = position_dodge(width = 0.9)) +
  labs(
    title = "A. SOFA Score Distribution by AKI Stage and CKD Status",
    x = "AKI Stage (KDIGO)",
    y = "SOFA Total Score",
    fill = ""
  ) +
  theme(legend.position = "bottom")

# Panel B: SOFA-stratified Stage 3 forest (using VERIFIED data, cutoff <=5)
sd <- data.table(
  sofa_cat = c("Low (0-5)","Medium (6-10)","High (>10)"),
  or = c(0.740, 0.957, 0.679),
  lo = c(0.480, 0.713, 0.552),
  hi = c(1.140, 1.284, 0.836),
  p = c(0.172, 0.770, 0.0003),
  n = c(483, 869, 1666)
)
sd[, label := sprintf("%.2f (%.2f-%.2f)", or, lo, hi)]
sd[, sig := ifelse(p < 0.001, "P<0.001", sprintf("P=%.3f", p))]

p3b <- ggplot(sd, aes(x = or, y = sofa_cat)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, color = C_PROT, linewidth = 0.7) +
  geom_point(color = C_PROT, size = 4, shape = 18) +
  geom_text(aes(x = 0.46, label = paste0("n=", n)), hjust = 0, size = 3, color = "grey40") +
  geom_text(aes(x = max(hi) + 0.08, label = paste0(label, "  ", sig)),
            hjust = 0, size = 3.5, family = "sans") +
  scale_x_log10(limits = c(0.42, 2.5), breaks = c(0.5, 1, 1.5, 2)) +
  labs(
    title = "B. Stage 3 AKI: CKD Effect Stratified by SOFA Severity",
    x = "Odds Ratio (95% CI), log scale",
    y = "SOFA Category"
  ) +
  theme(plot.margin = margin(15, 80, 15, 25))

p3 <- (p3a / p3b) + plot_layout(heights = c(1, 1))
ggsave(file.path(outdir, "Figure_3_SOFA_Distribution.pdf"), p3, width = 10, height = 12, device = cairo_pdf)
ggsave(file.path(outdir, "Figure_3_SOFA_Distribution.png"), p3, width = 10, height = 12, dpi = 300)
cat("  Saved Figure 3\n")

# ============================================================
# FIGURE 4: Nomogram (rms package)
# ============================================================
cat("\n=== Figure 4: Nomogram ===\n")

dd <- datadist(aki[, .(age, male, htn, dm, hf, sofa_total, aki_stage_num, ckd, ckd_s2, ckd_s3)])
options(datadist = "dd")

# Model C: age + sex + comorbidities + SOFA + stage + CKD + interaction
nom_model <- lrm(died ~ age + male + htn + dm + hf + sofa_total + aki_stage_num + ckd + ckd_s2 + ckd_s3, data = aki, x = TRUE, y = TRUE)

nom <- nomogram(
  nom_model,
  fun = function(x) 1 / (1 + exp(-x)),
  funlabel = "Predicted Mortality Probability",
  fun.at = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85),
  lp = FALSE
)

# Use cairo_pdf for proper font embedding
cairo_pdf(file.path(outdir, "Figure_4_Nomogram.pdf"), width = 16, height = 14, onefile = TRUE)
par(mar = c(2, 1, 3, 1))
plot(nom, xfrac = 0.45, cex.axis = 0.95, cex.var = 1.0,
     col.grid = c("#2166AC", "#B2182B"),
     points.label = "Points", total.points.label = "Total Points",
     fun.side = c(1,1,1,1,3,3,3,3),
     lmgp = 0.3)
title(main = "Figure 4. Nomogram for Predicting Hospital Mortality in AKI Patients",
      cex.main = 1.4, font.main = 2)
dev.off()

png(file.path(outdir, "Figure_4_Nomogram.png"), width = 4800, height = 4200, res = 300)
par(mar = c(2, 1, 3, 1))
plot(nom, xfrac = 0.45, cex.axis = 0.95, cex.var = 1.0,
     col.grid = c("#2166AC", "#B2182B"),
     points.label = "Points", total.points.label = "Total Points",
     fun.side = c(1,1,1,1,3,3,3,3),
     lmgp = 0.3)
title(main = "Figure 4. Nomogram for Predicting Hospital Mortality in AKI Patients",
      cex.main = 1.4, font.main = 2)
dev.off()
cat("  Saved Figure 4\n")

# ============================================================
# FIGURE 5: ROC + LASSO Feature Importance
# UNIFIED model definitions matching Python compute_all_stats_v4.py
# ============================================================
cat("\n=== Figure 5: ROC + Feature Importance ===\n")

# Model A: age, sex, comorbidities, SOFA, AKI stage (NO CKD)
# Model B: A + CKD main effect
# Model C: B + CKD x Stage interactions
# Model D: Replace total SOFA with 6 individual components

models <- list(
  "Model A\n(Base)" = glm(died ~ age + male + htn + dm + hf + sofa_total + stage2 + stage3,
                           data = aki, family = binomial),
  "Model B\n(+CKD)" = glm(died ~ age + male + htn + dm + hf + sofa_total + stage2 + stage3 + ckd,
                           data = aki, family = binomial),
  "Model C\n(+Interact)" = glm(died ~ age + male + htn + dm + hf + sofa_total + stage2 + stage3 + ckd + ckd_s2 + ckd_s3,
                                data = aki, family = binomial),
  "Model D\n(+SOFA comp)" = glm(died ~ age + male + htn + dm + hf + sofa_resp + sofa_coag + sofa_liver + sofa_cv + sofa_cns + sofa_renal + stage2 + stage3 + ckd + ckd_s2 + ckd_s3,
                                 data = aki, family = binomial)
)

# Compute ROC curves
roc_data <- list()
auc_values <- numeric()
for (nm in names(models)) {
  pred <- predict(models[[nm]], type = "response")
  roc_obj <- roc(aki$died, pred, quiet = TRUE)
  auc_values[nm] <- as.numeric(auc(roc_obj))
  roc_data[[nm]] <- data.table(
    model = nm,
    specificity = rev(roc_obj$specificities),
    sensitivity = rev(roc_obj$sensitivities),
    fpr = 1 - rev(roc_obj$specificities)
  )
  cat(sprintf("  %s: AUC = %.4f\n", gsub("\n", " ", nm), auc_values[nm]))
}
roc_dt <- rbindlist(roc_data)
roc_dt[, model := factor(model, levels = names(models))]

# Use verified AUC CIs from Python bootstrap
auc_ci_l <- c(0.709, 0.709, 0.710, 0.738)
auc_ci_u <- c(0.725, 0.726, 0.727, 0.755)
auc_labels <- paste0(gsub("\n", " ", names(models)),
                     "  AUC=", sprintf("%.3f", auc_values),
                     " (", sprintf("%.3f", auc_ci_l), "-", sprintf("%.3f", auc_ci_u), ")")
names(auc_labels) <- names(models)

p5a <- ggplot(roc_dt, aes(x = fpr, y = sensitivity, color = model)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("#999999","#FC8D62","#2166AC","#B2182B"), labels = auc_labels) +
  coord_equal() +
  labs(title = "A. ROC Curves: Four Nested Models", x = "1 - Specificity", y = "Sensitivity", color = "") +
  theme(legend.position = c(0.55, 0.25),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill = "white", color = "grey80"))

# Panel B: Feature importance (|z-value| from Model D)
md_summary <- summary(models[["Model D\n(+SOFA comp)"]])$coefficients
feat_imp <- data.table(
  variable = rownames(md_summary),
  z = md_summary[, "z value"]
)
feat_imp <- feat_imp[!variable %in% c("(Intercept)")]
feat_imp[, variable := c("Age","Male","HTN","DM","HF","SOFA Resp","SOFA Coag","SOFA Liver","SOFA CV","SOFA CNS","SOFA Renal","Stage 2","Stage 3","CKD","CKD x Stage2","CKD x Stage3")]
feat_imp[, abs_z := abs(z)]
feat_imp[, direction := ifelse(z > 0, "Harmful", "Protective")]
feat_imp <- feat_imp[order(abs_z)]

p5b <- ggplot(feat_imp, aes(x = reorder(variable, abs_z), y = abs_z, fill = direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Protective" = C_PROT, "Harmful" = C_HARM)) +
  labs(title = "B. Feature Importance (|z-value|, Model D)",
       x = "", y = "|z-value|", fill = "Direction") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 9))

p5 <- (p5a + p5b) + plot_layout(widths = c(1, 1))
ggsave(file.path(outdir, "Figure_5_ROC_Importance.pdf"), p5, width = 12, height = 6, device = cairo_pdf)
ggsave(file.path(outdir, "Figure_5_ROC_Importance.png"), p5, width = 12, height = 6, dpi = 300)
cat("  Saved Figure 5\n")

# ============================================================
# FIGURE 6: Calibration Curves (4 models, 5-fold CV)
# UNIFIED model definitions matching Python
# ============================================================
cat("\n=== Figure 6: Calibration Curves ===\n")

model_names <- c("Model A","Model B","Model C","Model D")
model_objs <- list(
  lrm(died ~ age + male + htn + dm + hf + sofa_total + stage2 + stage3, data = aki, x = TRUE, y = TRUE),
  lrm(died ~ age + male + htn + dm + hf + sofa_total + stage2 + stage3 + ckd, data = aki, x = TRUE, y = TRUE),
  lrm(died ~ age + male + htn + dm + hf + sofa_total + stage2 + stage3 + ckd + ckd_s2 + ckd_s3, data = aki, x = TRUE, y = TRUE),
  lrm(died ~ age + male + htn + dm + hf + sofa_resp + sofa_coag + sofa_liver + sofa_cv + sofa_cns + sofa_renal + stage2 + stage3 + ckd + ckd_s2 + ckd_s3, data = aki, x = TRUE, y = TRUE)
)

cal_data <- list()
for (i in 1:4) {
  cal <- calibrate(model_objs[[i]], method = "crossvalidation", B = 5, bw = FALSE, plots = FALSE)
  cal_data[[i]] <- data.table(
    model = model_names[i],
    predicted = cal[, "predy"],
    actual = cal[, "calibrated.orig"],
    calibrated = cal[, "calibrated.corrected"]
  )
}
cal_dt <- rbindlist(cal_data)
cal_dt[, model := factor(model, levels = model_names)]

# Verified Brier scores
brier_labels <- c("Model A  Brier=0.133", "Model B  Brier=0.133", "Model C  Brier=0.133", "Model D  Brier=0.129")

p6 <- ggplot(cal_dt, aes(x = predicted, y = actual)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(color = model), linewidth = 0.8) +
  geom_point(aes(color = model, shape = model), size = 2) +
  scale_color_manual(values = c("#999999","#FC8D62","#2166AC","#B2182B"), labels = brier_labels) +
  scale_shape_manual(values = c(1, 2, 15, 18)) +
  labs(
    title = "Figure 6. Calibration Curves (5-Fold Cross-Validation)",
    x = "Predicted Mortality Probability",
    y = "Observed Mortality Probability",
    color = "", shape = ""
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text = element_text(size = 9))

ggsave(file.path(outdir, "Figure_6_Calibration.pdf"), p6, width = 8, height = 7, device = cairo_pdf)
ggsave(file.path(outdir, "Figure_6_Calibration.png"), p6, width = 8, height = 7, dpi = 300)
cat("  Saved Figure 6\n")

# ============================================================
# FIGURE S1: Subgroup Forest Plot (Stage 3 patients)
# Uses cutoff <=5 for SOFA strata
# ============================================================
cat("\n=== Figure S1: Subgroup Forest Plot ===\n")

s3dat <- aki[max_aki_stage == 3]
subgroups <- list(
  "Age < 65" = s3dat[age < 65],
  "Age >= 65" = s3dat[age >= 65],
  "Male" = s3dat[male == 1],
  "Female" = s3dat[male == 0],
  "Hypertension" = s3dat[htn == 1],
  "No Hypertension" = s3dat[htn == 0],
  "Diabetes" = s3dat[dm == 1],
  "No Diabetes" = s3dat[dm == 0],
  "Heart Failure" = s3dat[hf == 1],
  "No Heart Failure" = s3dat[hf == 0],
  "SOFA Low (0-5)" = s3dat[sofa_total <= 5],
  "SOFA Medium (6-10)" = s3dat[sofa_total > 5 & sofa_total <= 10],
  "SOFA High (>10)" = s3dat[sofa_total > 10]
)

sg_data <- list()
for (nm in names(subgroups)) {
  sub <- subgroups[[nm]]
  if (nrow(sub) > 10 && length(unique(sub$ckd)) > 1 && min(table(sub$ckd)) > 3) {
    m <- glm(died ~ ckd, data = sub, family = binomial)
    c <- coef(m)["ckd"]; se <- summary(m)$coefficients["ckd", "Std. Error"]
    sg_data[[nm]] <- data.table(
      subgroup = nm, or = exp(c), lo = exp(c - 1.96*se), hi = exp(c + 1.96*se),
      p = summary(m)$coefficients["ckd", "Pr(>|z|)"], n = nrow(sub)
    )
  }
}
sgd <- rbindlist(sg_data)
sgd[, label := sprintf("%.2f (%.2f-%.2f)", or, lo, hi)]
sgd[, sig := ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))]
sgd[, subgroup := factor(subgroup, levels = rev(subgroup))]

pS1 <- ggplot(sgd, aes(x = or, y = subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, color = C_PROT, linewidth = 0.6) +
  geom_point(color = C_PROT, size = 3, shape = 18) +
  geom_text(aes(x = max(hi) + 0.12, label = paste0(label, "  P=", sig)),
            hjust = 0, size = 3, family = "sans") +
  geom_text(aes(x = 0.35, label = paste0("n=", n)), hjust = 0, size = 2.8, color = "grey40") +
  scale_x_log10(limits = c(0.3, 3.0), breaks = c(0.5, 1, 2)) +
  labs(
    title = "Figure S1. Subgroup Analysis: Stage 3 AKI - CKD Effect on Mortality",
    x = "Odds Ratio (95% CI), log scale",
    y = "Subgroup"
  ) +
  theme(plot.margin = margin(15, 100, 15, 10))

ggsave(file.path(outdir, "Figure_S1_Subgroup_Forest.pdf"), pS1, width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(outdir, "Figure_S1_Subgroup_Forest.png"), pS1, width = 10, height = 8, dpi = 300)
cat("  Saved Figure S1\n")

# ============================================================
# FIGURE S2: Kaplan-Meier Curves (by AKI stage)
# ============================================================
cat("\n=== Figure S2: Kaplan-Meier Curves ===\n")

km_plots <- list()
for (s in 1:3) {
  sub <- sv[max_aki_stage == s]
  fit <- survfit(Surv(survival_days, died) ~ group, data = sub)
  sd <- survdiff(Surv(survival_days, died) ~ group, data = sub)
  pval <- 1 - pchisq(sd$chisq, length(sd$n) - 1)

  km_plots[[s]] <- ggsurvplot(
    fit, data = sub,
    pval = TRUE, pval.size = 3.5,
    conf.int = TRUE, conf.int.alpha = 0.15,
    risk.table = TRUE, risk.table.col = "strata",
    risk.table.height = 0.25,
    risk.table.fontsize = 2.5,
    risk.table.y.text = FALSE,
    palette = c(C_AKI, C_CKD),
    legend.labs = c("Pure AKI", "AKI+CKD"),
    legend.title = "",
    title = paste0("Stage ", s, " AKI"),
    xlab = "Days", ylab = "Survival Probability",
    ggtheme = theme_pub,
    surv.plot.height = 0.75,
    break.x.by = 14
  )
}

# Arrange 3 panels
pS2_arr <- arrange_ggsurvplots(km_plots, nrow = 3, ncol = 1,
                                title = "Figure S2. Kaplan-Meier Survival Curves by AKI Stage",
                                title.size = 1.2,
                                font.title = 2)
ggsave(file.path(outdir, "Figure_S2_KM_Curves.pdf"), pS2_arr, width = 9, height = 16, device = cairo_pdf)
ggsave(file.path(outdir, "Figure_S2_KM_Curves.png"), pS2_arr, width = 9, height = 16, dpi = 300)
cat("  Saved Figure S2\n")

# ============================================================
# FIGURE S3: Aalen-Johansen Competing Risk
# ============================================================
cat("\n=== Figure S3: Competing Risk Curves ===\n")

cr_plots <- list()
for (s in 1:3) {
  sub <- sv[max_aki_stage == s]
  sub[, event := ifelse(died == 1, 1, 2)]
  sub[, group_f := factor(ifelse(ckd == 1, "AKI+CKD", "Pure AKI"))]

  pure <- sub[group_f == "Pure AKI"]
  ckdg <- sub[group_f == "AKI+CKD"]

  ci_pure <- cuminc(pure$survival_days, pure$event, cencode = 0)
  ci_ckd <- cuminc(ckdg$survival_days, ckdg$event, cencode = 0)

  pure_cif <- ci_pure[["1 1"]]
  ckd_cif <- ci_ckd[["1 1"]]

  cif_dt <- rbind(
    data.table(time = pure_cif$time, cif = pure_cif$est, group = "Pure AKI"),
    data.table(time = ckd_cif$time, cif = ckd_cif$est, group = "AKI+CKD")
  )

  p <- ggplot(cif_dt, aes(x = time, y = cif, color = group)) +
    geom_step(linewidth = 0.8) +
    scale_color_manual(values = c("Pure AKI" = C_AKI, "AKI+CKD" = C_CKD)) +
    labs(
      title = paste0("Stage ", s, " AKI"),
      x = "Days", y = "Cumulative Incidence of Death",
      color = ""
    ) +
    theme(legend.position = "bottom")

  cr_plots[[s]] <- p
}

pS3 <- (cr_plots[[1]] | cr_plots[[2]] | cr_plots[[3]]) +
  plot_annotation(
    title = "Figure S3. Aalen-Johansen Cumulative Incidence Functions (Competing Risk)",
    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
  )

ggsave(file.path(outdir, "Figure_S3_Competing_Risk.pdf"), pS3, width = 14, height = 5, device = cairo_pdf)
ggsave(file.path(outdir, "Figure_S3_Competing_Risk.png"), pS3, width = 14, height = 5, dpi = 300)
cat("  Saved Figure S3\n")

# ============================================================
# FIGURE S4: RCS — Age Effect on Mortality (rms)
# Reference age = 65
# ============================================================
cat("\n=== Figure S4: RCS Age Effect ===\n")

dd2 <- datadist(aki)
options(datadist = "dd2")

rcs_model <- lrm(died ~ rcs(age, 4) + male + htn + dm + hf + sofa_total +
                  aki_stage_num + ckd, data = aki, x = TRUE, y = TRUE)

age_range <- seq(min(aki$age), max(aki$age), length.out = 200)
pred_dt <- data.table(
  age = age_range,
  log_odds = Predict(rcs_model, age = age_range)$yhat,
  lo = Predict(rcs_model, age = age_range)$lower,
  hi = Predict(rcs_model, age = age_range)$upper
)

ref_age <- 65
ref_idx <- which.min(abs(age_range - ref_age))
ref_lo <- pred_dt$log_odds[ref_idx]
pred_dt[, or := exp(log_odds - ref_lo)]
pred_dt[, or_lo := exp(lo - ref_lo)]
pred_dt[, or_hi := exp(hi - ref_lo)]

knots <- quantile(aki$age, probs = c(0.05, 0.35, 0.65, 0.95))

pS4 <- ggplot(pred_dt, aes(x = age, y = or)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = or_lo, ymax = or_hi), fill = C_PROT, alpha = 0.15) +
  geom_line(color = C_PROT, linewidth = 1) +
  geom_vline(xintercept = knots, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  annotate("text", x = knots, y = max(pred_dt$or_hi) * 0.95,
           label = sprintf("%.0f", knots), size = 3, color = "grey40", vjust = -0.5) +
  scale_y_continuous(breaks = c(0.5, 1, 2, 4, 8, 16), trans = "log2") +
  coord_cartesian(ylim = c(0.5, 16)) +
  labs(
    title = "Figure S4. Restricted Cubic Splines: Age Effect on Mortality",
    subtitle = "Reference: age 65 | 4 knots at 5th/35th/65th/95th percentiles",
    x = "Age (years)",
    y = "Odds Ratio (95% CI), log2 scale"
  )

ggsave(file.path(outdir, "Figure_S4_RCS_Age.pdf"), pS4, width = 8, height = 6, device = cairo_pdf)
ggsave(file.path(outdir, "Figure_S4_RCS_Age.png"), pS4, width = 8, height = 6, dpi = 300)
cat("  Saved Figure S4\n")

# ============================================================
# Summary
# ============================================================
cat("\n\n=== ALL FIGURES COMPLETE (v4) ===\n")
cat("Output directory:", outdir, "\n")
files <- list.files(outdir, pattern = "\\.(pdf|png)$")
for (f in files) cat(" ", f, "\n")
