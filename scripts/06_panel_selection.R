# 06_panel_selection.R — 面板候选 v3：相关性去冗余 + 贪心增量覆盖 + 扩增子 SNP 注释
.libPaths(c("C:/Users/ld/Rlibs", .libPaths()))
suppressMessages({library(data.table); library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})
base <- "C:/Users/ld/WorkBuddy/2026-08-19-15-50-19"

tb <- fread(file.path(base, "results/panel_probes_tcga_beta.tsv"))
probes <- tb[[1]]; mat <- as.matrix(tb[, -1]); rownames(mat) <- probes
tumor <- grep("-01$", colnames(mat)); normal <- grep("-11$", colnames(mat))

elite <- fread(file.path(base, "results/elite_probes_v2.csv"))
clin <- fread(file.path(base, "data/UCEC_clinicalMatrix.tsv"))
meno <- setNames(clin$menopause_status, substr(clin$sampleID, 1, 15))
is_pre <- startsWith(meno[colnames(mat)], "Pre"); is_pre[is.na(is_pre)] <- FALSE
pre_idx <- which(colnames(mat) %in% colnames(mat)[tumor] & is_pre)
pre_tumor_cols <- intersect(tumor, which(is_pre))
cat("tumors:", length(tumor), "| Pre tumors:", length(pre_tumor_cols), "\n")

# 阳性矩阵（β>0.3）
pos <- mat[, tumor] > 0.3
pre_pos <- mat[, pre_tumor_cols, drop = FALSE] > 0.3

# 相关性冗余分组（癌样本 β 值，r>0.8）
cormat <- cor(t(mat[, tumor]), use = "pairwise.complete.obs")
grp <- list(); assigned <- c()
for (p in probes) {
  if (p %in% assigned) next
  members <- probes[cormat[p, probes] > 0.8]
  grp[[length(grp) + 1]] <- members
  assigned <- c(assigned, members)
}
redundant <- sapply(grp, function(g) if (length(g) > 1) paste(g, collapse = ";") else "")
cat("redundancy groups (>1 probe, r>0.8):", sum(nchar(redundant) > 0), "\n")
for (g in redundant[nchar(redundant) > 0]) cat("  group:", g, "\n")

# 贪心增量覆盖（以绝经前癌为主目标，全癌为次目标）
# 候选池：仅背景通过的精英探针（ne_pass==TRUE，含 CDO1）；原 5 位点中背景不过者不参与选择
score_map <- setNames(elite$score, elite$probe)
sc <- function(p) { v <- score_map[p]; ifelse(is.na(v), 0, v) }
pool <- intersect(elite[elite$ne_pass == TRUE]$probe, probes)
cat("greedy pool (ne_pass):", length(pool), "\n")
covered_pre <- rep(FALSE, length(pre_tumor_cols))
covered_all <- rep(FALSE, length(tumor))
panel <- c()
remaining <- pool
cat("\n--- greedy panel ---\n")
for (step in 1:12) {
  best <- NULL; best_gain <- -1; best_all_gain <- 0
  for (p in remaining) {
    new_cov <- pre_pos[p, ] & !covered_pre
    gain <- sum(new_cov, na.rm = TRUE)
    all_gain <- sum(pos[p, ] & !covered_all, na.rm = TRUE)
    if (gain > best_gain || (gain == best_gain && all_gain > best_all_gain) ||
        (gain == best_gain && all_gain == best_all_gain && !is.null(best) && sc(p) > sc(best))) {
      best <- p; best_gain <- gain; best_all_gain <- all_gain
    }
  }
  if (is.null(best) || (best_gain == 0 && best_all_gain == 0)) break
  panel <- c(panel, best)
  covered_pre <- covered_pre | (pre_pos[best, ] %in% TRUE)
  covered_all <- covered_all | (pos[best, ] %in% TRUE)
  remaining <- setdiff(remaining, best)
  cat(sprintf("%2d. %s  +Pre %d (累计 %.1f%%)  +All %d (累计 %.1f%%)\n",
              step, best, best_gain, 100 * mean(covered_pre), best_all_gain, 100 * mean(covered_all)))
  if (mean(covered_pre) >= 0.999 && step >= 6) { cat("Pre 覆盖已满，停\n"); break }
}

# SNP/扩增子注释（按可用列名自适应）
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
snp_col <- intersect(c("Probe_SNPs", "Probe_SNPs_10"), colnames(anno))
res <- data.frame(probe = panel,
                  gene = anno[panel, "UCSC_RefGene_Name"],
                  chr = anno[panel, "chr"], pos = anno[panel, "pos"],
                  island = anno[panel, "Relation_to_Island"],
                  pre_pos03 = round(rowMeans(pre_pos[panel, , drop = FALSE], na.rm = TRUE), 3),
                  all_pos03 = round(rowMeans(pos[panel, , drop = FALSE], na.rm = TRUE), 3),
                  stringsAsFactors = FALSE)
for (cc in snp_col) res[[cc]] <- anno[panel, cc]
fwrite(res, file.path(base, "results/panel_candidates_v3.csv"))
cat("\nsaved panel_candidates_v3.csv\n")
print(res)
