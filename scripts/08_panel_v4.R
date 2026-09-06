# 08_panel_v4.R — 排雷通过探针的贪心面板（v4 终版）
.libPaths(c("C:/Users/ld/Rlibs", .libPaths()))
suppressMessages({library(data.table); library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})
base <- "C:/Users/ld/WorkBuddy/2026-08-19-15-50-19"

tb <- fread(file.path(base, "results/panel_probes_tcga_beta_v4.tsv"))
probes <- tb[[1]]; mat <- as.matrix(tb[, -1]); rownames(mat) <- probes
tumor <- grep("-01$", colnames(mat))
elite <- fread(file.path(base, "results/elite_probes_v2.csv"))
score_map <- setNames(elite$score, elite$probe)
sc <- function(p) { v <- score_map[p]; ifelse(is.na(v), 0, v) }
clin <- fread(file.path(base, "data/UCEC_clinicalMatrix.tsv"))
meno <- setNames(clin$menopause_status, substr(clin$sampleID, 1, 15))
is_pre <- startsWith(meno[colnames(mat)], "Pre"); is_pre[is.na(is_pre)] <- FALSE
pre_tumor_cols <- intersect(tumor, which(is_pre))

pos <- mat[, tumor] > 0.3
pre_pos <- mat[, pre_tumor_cols, drop = FALSE] > 0.3

# 冗余
cormat <- cor(t(mat[, tumor]), use = "pairwise.complete.obs")
for (i in seq_len(nrow(cormat))) for (j in seq_len(ncol(cormat)))
  if (i < j && cormat[i, j] > 0.8) cat("redundant pair:", rownames(cormat)[i], "&", colnames(cormat)[j], "r=", round(cormat[i, j], 3), "\n")

covered_pre <- rep(FALSE, length(pre_tumor_cols)); covered_all <- rep(FALSE, length(tumor))
panel <- c(); remaining <- probes
cat("\n--- greedy (v4, triage-passed pool) ---\n")
for (step in 1:12) {
  best <- NULL; best_gain <- -1; best_all_gain <- 0
  for (p in remaining) {
    gain <- sum(pre_pos[p, ] & !covered_pre, na.rm = TRUE)
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
  if (mean(covered_all) >= 0.999 && step >= 4) { cat("覆盖已满\n"); break }
}

# 备份：按评分取未入选的前 4 个
backup <- setdiff(probes, panel)
backup <- backup[order(-sapply(backup, sc))][1:min(4, length(backup))]
final <- c(panel, backup)
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
tri <- fread(file.path(base, "results/triage_closedloop_v4.csv"))
res <- data.frame(role = c(rep("core", length(panel)), rep("backup", length(backup))),
                  probe = final,
                  gene = anno[final, "UCSC_RefGene_Name"],
                  chr = anno[final, "chr"],
                  island = anno[final, "Relation_to_Island"],
                  pre_pos03 = round(rowMeans(pre_pos[final, , drop = FALSE], na.rm = TRUE), 3),
                  all_pos03 = round(rowMeans(pos[final, , drop = FALSE], na.rm = TRUE), 3),
                  elite_score = round(sapply(final, sc), 3))
tri_sub <- tri[match(final, tri$probe), .(cycle_mean_abs_d, cycle_frac_gt01, cervix_p95, blood_max_p95, blood_worst_type)]
res <- cbind(res, tri_sub)
fwrite(res, file.path(base, "results/panel_candidates_v4.csv"))
cat("\n=== 推荐面板 v4 ===\n"); print(res)
