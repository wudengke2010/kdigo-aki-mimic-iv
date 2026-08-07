###################################################
# Figure 1: Study Flow Diagram (R ggplot2 version)
# Publication-quality flowchart with actual MIMIC-IV data
###################################################

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(data.table)
}))

outdir <- "submission_figures_R"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data to get exact counts ----
dat <- fread("aki_paradox_cohort_with_sofa.csv")
total_icu <- nrow(dat)
aki <- dat[max_aki_stage > 0]
no_aki <- dat[max_aki_stage == 0]
n_aki <- nrow(aki)
n_no_aki <- nrow(no_aki)

# Stage counts
s1 <- aki[max_aki_stage == 1]; s2 <- aki[max_aki_stage == 2]; s3 <- aki[max_aki_stage == 3]
s1_pure <- s1[ckd == 0]; s1_ckd <- s1[ckd == 1]
s2_pure <- s2[ckd == 0]; s2_ckd <- s2[ckd == 1]
s3_pure <- s3[ckd == 0]; s3_ckd <- s3[ckd == 1]

# Mortality
m_s3_pure <- s3_pure[, mean(hospital_expire_flag) * 100]
m_s3_ckd  <- s3_ckd[, mean(hospital_expire_flag) * 100]

cat(sprintf("Total ICU: %d | No AKI: %d | AKI: %d\n", total_icu, n_no_aki, n_aki))
cat(sprintf("S1: %d (Pure %d / CKD %d)\n", nrow(s1), nrow(s1_pure), nrow(s1_ckd)))
cat(sprintf("S2: %d (Pure %d / CKD %d)\n", nrow(s2), nrow(s2_pure), nrow(s2_ckd)))
cat(sprintf("S3: %d (Pure %d / CKD %d) | Mortality: Pure %.1f%% / CKD %.1f%%\n",
            nrow(s3), nrow(s3_pure), nrow(s3_ckd), m_s3_pure, m_s3_ckd))

# ---- Colors ----
C_MAIN   <- "#2166AC"
C_AKI    <- "#E08214"
C_STAGE1 <- "#66C2A5"
C_STAGE2 <- "#FC8D62"
C_STAGE3 <- "#8DA0CB"
C_PURE   <- "#FFF3E0"
C_CKD    <- "#EDE7F6"
C_HARM  <- "#B2182B"
C_HIGHLIGHT <- "#B2182B"
C_PROT   <- "#2166AC"
C_TITLE  <- "#1A1A2E"

# ---- Build flowchart with ggplot2 ----
# Canvas: 14 x 18 (width x height in units)
# Layout: top-down flow

# Helper function to create a box
make_box <- function(xmin, xmax, ymin, ymax, fill, color, label, sublabel = NULL,
                     text_size = 3.8, sub_size = 3.2, text_color = "black",
                     label_face = "bold", line_width = 0.8) {
  bx <- data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  list(
    geom_rect(data = bx, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = fill, color = color, linewidth = line_width, inherit.aes = FALSE),
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2 + ifelse(is.null(sublabel), 0, 0.3),
             label = label, size = text_size, fontface = label_face, color = text_color)
  )
}

# Arrow helper
make_arrow <- function(x, y1, y2, color = "grey30", linewidth = 0.6, linetype = "solid") {
  annotate("segment", x = x, xend = x, y = y1, yend = y2,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           color = color, linewidth = linewidth, linetype = linetype)
}

# L-shaped arrow (horizontal then vertical)
make_larrow <- function(x1, y1, x2, y2, color = "grey30", linewidth = 0.6) {
  list(
    annotate("segment", x = x1, xend = x2, y = y1, yend = y1,
             color = color, linewidth = linewidth),
    annotate("segment", x = x2, xend = x2, y = y1, yend = y2,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             color = color, linewidth = linewidth)
  )
}

# Build the plot
p <- ggplot() + xlim(0, 20) + ylim(0, 22) +
  coord_fixed(ratio = 1) +

  # ===== Title =====
  annotate("text", x = 10, y = 21.7, label = "Figure 1. Study Flow Diagram",
           size = 5.5, fontface = "bold", color = C_TITLE) +
  annotate("text", x = 10, y = 20.95, label = "MIMIC-IV v3.1 (2008-2022) - AKI Paradox Cohort",
           size = 3.5, fontface = "italic", color = "grey40") +

  # ===== Level 1: MIMIC-IV Total =====
  geom_rect(aes(xmin = 6.5, xmax = 13.5, ymin = 19.2, ymax = 20.3),
            fill = "#E3F2FD", color = C_MAIN, linewidth = 1.0) +
  annotate("text", x = 10, y = 20.0, label = "MIMIC-IV v3.1 Database",
           size = 4.2, fontface = "bold", color = C_MAIN) +
  annotate("text", x = 10, y = 19.5, label = sprintf("Total ICU Patients: %s", format(total_icu, big.mark = ",")),
           size = 3.8, color = "black") +

  # Arrow down
  make_arrow(10, 19.2, 18.5) +

  # ===== Level 2: Inclusion Criteria =====
  geom_rect(aes(xmin = 5.5, xmax = 14.5, ymin = 17.2, ymax = 18.5),
            fill = "#FFF8E1", color = "#F9A825", linewidth = 0.8) +
  annotate("text", x = 10, y = 18.15, label = "Inclusion Criteria",
           size = 3.8, fontface = "bold", color = "#E65100") +
  annotate("text", x = 10, y = 17.55,
           label = "Age >= 18 | First ICU admission | ICU LOS >= 24h\nExcluded: ESRD, missing creatinine/urine output",
           size = 3.0, color = "grey30", lineheight = 0.9) +

  # Arrow down
  make_arrow(10, 17.2, 16.5) +

  # ===== Level 3: KDIGO Assessment =====
  geom_rect(aes(xmin = 6, xmax = 14, ymin = 15.2, ymax = 16.5),
            fill = "#E8F5E9", color = "#43A047", linewidth = 0.8) +
  annotate("text", x = 10, y = 16.15, label = "KDIGO AKI Assessment",
           size = 3.8, fontface = "bold", color = "#2E7D32") +
  annotate("text", x = 10, y = 15.55,
           label = "Creatinine + Urine Output criteria\nMax stage from either criterion",
           size = 3.0, color = "grey30", lineheight = 0.9) +

  # Arrow down (splits into two)
  make_arrow(10, 15.2, 14.5) +
  annotate("segment", x = 5, xend = 15, y = 14.5, yend = 14.5, color = "grey30", linewidth = 0.6) +
  make_arrow(5, 14.5, 13.5) +
  make_arrow(15, 14.5, 13.5) +

  # ===== Level 4: No AKI / AKI split =====
  # No AKI box (left)
  geom_rect(aes(xmin = 1.5, xmax = 8.5, ymin = 12.2, ymax = 13.5),
            fill = "#F5F5F5", color = "#9E9E9E", linewidth = 0.8) +
  annotate("text", x = 5, y = 13.15, label = "No AKI",
           size = 3.8, fontface = "bold", color = "#616161") +
  annotate("text", x = 5, y = 12.55, label = sprintf("n = %s (62.1%%)", format(n_no_aki, big.mark = ",")),
           size = 3.5, color = "black") +
  annotate("text", x = 5, y = 12.2, label = "Excluded from AKI analysis", size = 2.8, color = "grey50") +

  # AKI box (right)
  geom_rect(aes(xmin = 11.5, xmax = 18.5, ymin = 12.2, ymax = 13.5),
            fill = "#FFF3E0", color = C_AKI, linewidth = 1.0) +
  annotate("text", x = 15, y = 13.15, label = "AKI Cohort",
           size = 4.0, fontface = "bold", color = C_AKI) +
  annotate("text", x = 15, y = 12.55, label = sprintf("n = %s (37.9%%)", format(n_aki, big.mark = ",")),
           size = 3.5, color = "black") +
  annotate("text", x = 15, y = 12.2, label = sprintf("CKD: %s (22.2%%)", format(5372, big.mark = ",")),
           size = 3.0, color = C_CKD) +

  # Arrow down from AKI to stage split
  make_arrow(15, 12.2, 11.3) +
  # Horizontal split line
  annotate("segment", x = 5, xend = 15, y = 11.3, yend = 11.3, color = "grey30", linewidth = 0.6) +
  make_arrow(5, 11.3, 10.3) +
  make_arrow(10, 11.3, 10.3) +
  make_arrow(15, 11.3, 10.3) +

  # ===== Level 5: Three AKI Stages =====
  # Stage 1 (left)
  geom_rect(aes(xmin = 2, xmax = 8, ymin = 9.0, ymax = 10.3),
            fill = "#E0F2F1", color = C_STAGE1, linewidth = 0.8) +
  annotate("text", x = 5, y = 10.05, label = "Stage 1 AKI",
           size = 3.5, fontface = "bold", color = C_STAGE1) +
  annotate("text", x = 5, y = 9.5, label = sprintf("n = %s", format(nrow(s1), big.mark = ",")),
           size = 3.2, color = "black") +
  annotate("text", x = 5, y = 9.1, label = "Mild AKI", size = 2.8, color = "grey50") +

  # Stage 2 (center)
  geom_rect(aes(xmin = 8.5, xmax = 11.5, ymin = 9.0, ymax = 10.3),
            fill = "#FBE9E7", color = C_STAGE2, linewidth = 0.8) +
  annotate("text", x = 10, y = 10.05, label = "Stage 2 AKI",
           size = 3.5, fontface = "bold", color = C_STAGE2) +
  annotate("text", x = 10, y = 9.5, label = sprintf("n = %s", format(nrow(s2), big.mark = ",")),
           size = 3.2, color = "black") +
  annotate("text", x = 10, y = 9.1, label = "Moderate AKI", size = 2.8, color = "grey50") +

  # Stage 3 (right)
  geom_rect(aes(xmin = 12, xmax = 18, ymin = 9.0, ymax = 10.3),
            fill = "#E8EAF6", color = C_STAGE3, linewidth = 1.0) +
  annotate("text", x = 15, y = 10.05, label = "Stage 3 AKI",
           size = 3.5, fontface = "bold", color = C_STAGE3) +
  annotate("text", x = 15, y = 9.5, label = sprintf("n = %s", format(nrow(s3), big.mark = ",")),
           size = 3.2, color = "black") +

  # Arrows from each stage down to Pure/CKD split
  make_arrow(5, 9.0, 8.2) +
  make_arrow(10, 9.0, 8.2) +
  make_arrow(15, 9.0, 8.2) +

  # ===== Level 6: Pure AKI / AKI+CKD for each stage =====
  # Stage 1: Pure / CKD
  geom_rect(aes(xmin = 2, xmax = 4.8, ymin = 6.8, ymax = 8.2),
            fill = C_PURE, color = C_STAGE1, linewidth = 0.6) +
  annotate("text", x = 3.4, y = 7.95, label = "Pure AKI", size = 3.0, fontface = "bold", color = "#555") +
  annotate("text", x = 3.4, y = 7.35, label = sprintf("n = %s", format(nrow(s1_pure), big.mark = ",")), size = 2.8) +
  annotate("text", x = 3.4, y = 6.95, label = sprintf("Mort: %.1f%%", s1_pure[, mean(hospital_expire_flag)*100]), size = 2.6, color = "grey40") +

  geom_rect(aes(xmin = 5.2, xmax = 8, ymin = 6.8, ymax = 8.2),
            fill = C_CKD, color = C_STAGE1, linewidth = 0.6) +
  annotate("text", x = 6.6, y = 7.95, label = "AKI+CKD", size = 3.0, fontface = "bold", color = C_CKD) +
  annotate("text", x = 6.6, y = 7.35, label = sprintf("n = %s", format(nrow(s1_ckd), big.mark = ",")), size = 2.8) +
  annotate("text", x = 6.6, y = 6.95, label = sprintf("Mort: %.1f%%", s1_ckd[, mean(hospital_expire_flag)*100]), size = 2.6, color = "grey40") +

  # Stage 2: Pure / CKD
  geom_rect(aes(xmin = 8.5, xmax = 10, ymin = 6.8, ymax = 8.2),
            fill = C_PURE, color = C_STAGE2, linewidth = 0.6) +
  annotate("text", x = 9.25, y = 7.95, label = "Pure", size = 2.8, fontface = "bold", color = "#555") +
  annotate("text", x = 9.25, y = 7.35, label = format(nrow(s2_pure), big.mark = ","), size = 2.6) +
  annotate("text", x = 9.25, y = 6.95, label = sprintf("%.1f%%", s2_pure[, mean(hospital_expire_flag)*100]), size = 2.4, color = "grey40") +

  geom_rect(aes(xmin = 10, xmax = 11.5, ymin = 6.8, ymax = 8.2),
            fill = C_CKD, color = C_STAGE2, linewidth = 0.6) +
  annotate("text", x = 10.75, y = 7.95, label = "CKD", size = 2.8, fontface = "bold", color = C_CKD) +
  annotate("text", x = 10.75, y = 7.35, label = format(nrow(s2_ckd), big.mark = ","), size = 2.6) +
  annotate("text", x = 10.75, y = 6.95, label = sprintf("%.1f%%", s2_ckd[, mean(hospital_expire_flag)*100]), size = 2.4, color = "grey40") +

  # Stage 3: Pure / CKD (HIGHLIGHTED - this is the paradox)
  geom_rect(aes(xmin = 12, xmax = 14.8, ymin = 6.8, ymax = 8.2),
            fill = "#FFEBEE", color = C_HIGHLIGHT, linewidth = 1.0) +
  annotate("text", x = 13.4, y = 7.95, label = "Pure AKI", size = 3.0, fontface = "bold", color = C_HARM) +
  annotate("text", x = 13.4, y = 7.35, label = sprintf("n = %s", format(nrow(s3_pure), big.mark = ",")), size = 2.8) +
  annotate("text", x = 13.4, y = 6.95, label = sprintf("Mort: %.1f%%", m_s3_pure), size = 2.8, fontface = "bold", color = C_HARM) +

  geom_rect(aes(xmin = 15.2, xmax = 18, ymin = 6.8, ymax = 8.2),
            fill = "#E8F5E9", color = C_PROT, linewidth = 1.0) +
  annotate("text", x = 16.6, y = 7.95, label = "AKI+CKD", size = 3.0, fontface = "bold", color = C_PROT) +
  annotate("text", x = 16.6, y = 7.35, label = sprintf("n = %s", format(nrow(s3_ckd), big.mark = ",")), size = 2.8) +
  annotate("text", x = 16.6, y = 6.95, label = sprintf("Mort: %.1f%%", m_s3_ckd), size = 2.8, fontface = "bold", color = C_PROT) +

  # Arrow from Stage 3 CKD down to key finding
  annotate("segment", x = 16.6, xend = 16.6, y = 6.8, yend = 5.8,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           color = C_HIGHLIGHT, linewidth = 0.8, linetype = "dashed") +

  # ===== Level 7: Key Finding =====
  geom_rect(aes(xmin = 10, xmax = 19, ymin = 4.0, ymax = 5.8),
            fill = "#FFF8E1", color = "#F9A825", linewidth = 1.2) +
  annotate("text", x = 14.5, y = 5.5, label = "KEY FINDING: The AKI Paradox",
           size = 4.0, fontface = "bold", color = C_HIGHLIGHT) +
  annotate("text", x = 14.5, y = 4.8,
           label = sprintf("Stage 3 AKI: CKD patients have LOWER mortality\nPure AKI: %.1f%%  vs  AKI+CKD: %.1f%%  (P < 0.0001)\nSOFA-adjusted OR: 0.643 (0.537-0.769) — paradox STRENGTHENS", m_s3_pure, m_s3_ckd),
           size = 3.2, color = "black", lineheight = 1.1) +

  # ===== Left side: Analysis box =====
  geom_rect(aes(xmin = 1, xmax = 8.5, ymin = 4.0, ymax = 5.8),
            fill = "#E3F2FD", color = C_MAIN, linewidth = 0.8) +
  annotate("text", x = 4.75, y = 5.5, label = "Statistical Analysis",
           size = 3.5, fontface = "bold", color = C_MAIN) +
  annotate("text", x = 4.75, y = 4.5,
           label = "- Multivariable logistic regression (4 nested models)\n- Propensity score matching (with/without SOFA)\n- CKD x AKI stage interaction testing\n- SOFA-stratified sensitivity analysis\n- Nomogram + LASSO + 5-fold CV",
           size = 2.8, color = "grey30", lineheight = 0.95, hjust = 0.5) +

  # ===== Bottom: CKD identification note =====
  geom_rect(aes(xmin = 3, xmax = 17, ymin = 2.0, ymax = 3.5),
            fill = "#F3E5F5", color = C_CKD, linewidth = 0.6) +
  annotate("text", x = 10, y = 3.15, label = "CKD Identification",
           size = 3.2, fontface = "bold", color = C_CKD) +
  annotate("text", x = 10, y = 2.4,
           label = "ICD-9-CM (585.x) and ICD-10-CM (N18.x) codes  |  CKD prevalence: 22.2% of AKI cohort\nPure AKI: n = 18,832 (77.8%)  |  AKI+CKD: n = 5,372 (22.2%)",
           size = 2.8, color = "grey30", lineheight = 0.9) +

  # ===== Bottom-most: SOFA note =====
  geom_rect(aes(xmin = 3, xmax = 17, ymin = 0.3, ymax = 1.7),
            fill = "#E0F7FA", color = "#00838F", linewidth = 0.6) +
  annotate("text", x = 10, y = 1.35, label = "SOFA Score (First 24h ICU)",
           size = 3.2, fontface = "bold", color = "#00838F") +
  annotate("text", x = 10, y = 0.65,
           label = "6 organ systems: Respiratory, Coagulation, Liver, Cardiovascular, CNS, Renal\nMean SOFA: 8.1 (Stage 1: 7.4, Stage 2: 8.3, Stage 3: 10.8)",
           size = 2.8, color = "grey30", lineheight = 0.9) +

  # Remove axes
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))

# ---- Save ----
ggsave(file.path(outdir, "Figure_1_Flow_Diagram.pdf"), p, width = 12, height = 16, device = cairo_pdf)
ggsave(file.path(outdir, "Figure_1_Flow_Diagram.png"), p, width = 12, height = 16, dpi = 300)

cat("\nFigure 1 saved to submission_figures_R/\n")
