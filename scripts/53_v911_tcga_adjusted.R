# v9.11 Analysis A: confounder-adjusted sensitivity for the 26 probes (TCGA-UCEC)
# Models: unadjusted vs age-adjusted vs age+histology-adjusted (limma on M-values;
# beta-scale OLS for effect sizes). Tumour-only model tests subtype association.
library(limma)

df <- read.delim("results/v911_tcga_adj_input.tsv", stringsAsFactors = FALSE)
probes <- setdiff(colnames(df), c("sampleID","group","age","histology_simple","subtype_simple","menopause_simple"))
stopifnot(length(probes) == 26)

beta <- as.matrix(df[, probes]); rownames(df) <- df$sampleID
M <- log2(beta / (1 - beta))
M[beta <= 0] <- NA; M[beta >= 1] <- NA

df$tissue3 <- ifelse(df$group == "normal", "normal",
              ifelse(df$histology_simple == "endometrioid", "tum_eec", "tum_other"))
df$tissue3 <- factor(df$tissue3, levels = c("normal","tum_eec","tum_other"))

run_limma <- function(design, contrast, Msub) {
  fit <- lmFit(t(Msub), design)
  cm <- makeContrasts(contrasts = contrast, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  data.frame(probe = colnames(Msub),
             logFC = fit2$coefficients[, 1],
             P = fit2$p.value[, 1],
             FDR = p.adjust(fit2$p.value[, 1], method = "BH"))
}

# --- 1) unadjusted tumour vs normal ---
grp <- factor(df$group, levels = c("normal","tumour"))
d0 <- model.matrix(~ grp)
r_unadj <- run_limma(d0, "grptumour", M)

# --- 2) age-adjusted (samples with age) ---
ok <- !is.na(df$age)
d1 <- model.matrix(~ grp[ok] + df$age[ok])
colnames(d1) <- c("int","tumour","age")
r_age <- run_limma(d1, "tumour", M[ok, ])

# --- 3) age + histology (3-level tissue factor) ---
d2 <- model.matrix(~ 0 + df$tissue3[ok] + df$age[ok])
colnames(d2) <- c("normal","tum_eec","tum_other","age")
fit <- lmFit(t(M[ok, ]), d2)
cm <- makeContrasts(eec = tum_eec - normal, other = tum_other - normal, levels = d2)
fit2 <- eBayes(contrasts.fit(fit, cm))
r_hist <- data.frame(probe = colnames(M),
                     logFC_eec = fit2$coefficients[, "eec"], FDR_eec = p.adjust(fit2$p.value[, "eec"], "BH"),
                     logFC_other = fit2$coefficients[, "other"], FDR_other = p.adjust(fit2$p.value[, "other"], "BH"))

# --- beta-scale effect sizes (unadjusted vs age+histology adjusted) ---
db_unadj <- sapply(probes, function(p) mean(beta[df$group=="tumour", p], na.rm=TRUE) - mean(beta[df$group=="normal", p], na.rm=TRUE))
db_adj <- sapply(probes, function(p) {
  d <- df[ok & !is.na(df$histology_simple) | (ok & df$group=="normal"), ]
  d <- df[ok, ]
  m <- lm(beta[ok, p] ~ tissue3 + age, data = d)
  mean(coef(m)[c("tissue3tum_eec","tissue3tum_other")], na.rm = TRUE)  # average adjusted dbeta vs normal
})

out <- data.frame(probe = probes, dbeta_unadj = db_unadj[probes], dbeta_adj = db_adj[probes])
out <- merge(out, r_unadj[, c("probe","logFC","FDR")], by = "probe")
out <- merge(out, r_age[, c("probe","logFC","FDR")], by = "probe", suffixes = c("_unadj","_age"))
out <- merge(out, r_hist, by = "probe")
out <- out[order(-out$dbeta_unadj), ]
write.csv(out, "results/v911_tcga_adjusted_models.csv", row.names = FALSE)

cat("n samples:", nrow(df), "| with age:", sum(ok), "\n")
cat("unadj FDR<0.05:", sum(out$FDR_unadj < 0.05), "/26 | age-adj:", sum(out$FDR_age < 0.05),
    "| eec-adj:", sum(out$FDR_eec < 0.05), "| other-adj:", sum(out$FDR_other < 0.05), "\n")
cat("Pearson r (unadj vs age-adj logFC):", cor(out$logFC_unadj, out$logFC_age), "\n")
cat("Pearson r (unadj vs eec-adj logFC):", cor(out$logFC_unadj, out$logFC_eec), "\n")
cat("dbeta unadj range:", range(out$dbeta_unadj), "| adj range:", range(out$dbeta_adj), "\n")
cat("max |dbeta change|:", max(abs(out$dbeta_unadj - out$dbeta_adj)), "\n")

# --- 4) tumour-only: subtype/histology/age association ---
tum <- df$group == "tumour" & !is.na(df$subtype_simple) & !is.na(df$age)
sub_f <- factor(df$subtype_simple[tum])
d3 <- model.matrix(~ sub_f + df$age[tum])
fit3 <- eBayes(lmFit(t(M[tum, ]), d3))
subtype_p <- sapply(probes, function(p) {
  f <- topTable(fit3, coef = grep("sub_f", colnames(d3)), number = 1)$F[1]
  pf(f, 3, fit3$df.residual[1], lower.tail = FALSE)
})
tumour_only <- data.frame(probe = probes, subtype_Ftest_P = subtype_p[probes],
                          subtype_FDR = p.adjust(subtype_p[probes], "BH"))
write.csv(tumour_only, "results/v911_tcga_tumour_only_subtype.csv", row.names = FALSE)
cat("tumour-only subtype F-test: min FDR =", min(tumour_only$subtype_FDR),
    "| n FDR<0.05:", sum(tumour_only$subtype_FDR < 0.05), "\n")

# --- 5) endometrioid-only, age-adjusted ---
eec_ok <- ok & (df$group == "normal" | df$histology_simple == "endometrioid")
grp2 <- factor(ifelse(df$group[eec_ok] == "normal", "normal", "tumour"), levels = c("normal","tumour"))
d4 <- model.matrix(~ grp2 + df$age[eec_ok])
colnames(d4) <- c("int","tumour","age")
r_eec <- run_limma(d4, "tumour", M[eec_ok, ])
r_eec$dbeta <- sapply(r_eec$probe, function(p)
  mean(beta[eec_ok & df$group=="tumour", p], na.rm=TRUE) - mean(beta[eec_ok & df$group=="normal", p], na.rm=TRUE))
write.csv(r_eec, "results/v911_tcga_eec_only_adjusted.csv", row.names = FALSE)
cat("EEC-only age-adj: FDR<0.05:", sum(r_eec$FDR < 0.05), "/26 | dbeta range:", range(r_eec$dbeta), "\n")
