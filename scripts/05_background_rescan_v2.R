# 05_background_rescan_v2.R — 用 GSE73949 真实育龄正常内膜背景重排名 + 精英单探针轨
# 输出: DMR_candidates_ranked_v2.csv, elite_probes_v2.csv
.libPaths(c("C:/Users/ld/Rlibs", .libPaths()))
suppressMessages(library(data.table))
base <- "C:/Users/ld/WorkBuddy/2026-08-19-15-50-19/results"

ne <- readRDS(file.path(base, "GSE73949_betas.rds"))  # 17 例育龄正常内膜
ne_mean <- rowMeans(ne, na.rm = TRUE)
ne_p95 <- apply(ne, 1, quantile, 0.95, na.rm = TRUE)

cand <- fread(file.path(base, "DMR_candidates_ranked_v1.csv"))
gse_dmp <- fread(file.path(base, "GSE67116_dmp.csv"))

# ---- 轨 1: 候选 DMR 背景回扫 ----
ne_stats <- function(probes_) {
  pr <- unlist(strsplit(probes_, ";"))
  pr <- intersect(pr, names(ne_mean))
  if (!length(pr)) return(list(NA, NA))
  list(mean(ne_mean[pr]), mean(ne_p95[pr]))
}
st <- lapply(cand$probes, ne_stats)
cand$ne_bg_mean <- sapply(st, `[[`, 1)
cand$ne_bg_p95 <- sapply(st, `[[`, 2)

# 更新评分：良性低背景改用真实育龄内膜 P95（阈值 0.15），其余不变
cand$score_v2 <- with(cand, 0.25 * pmax(0, 1 - ne_bg_p95 / 0.15) +
                        0.25 * pre_pos03 + 0.20 * tumor_pos03 + 0.15 * 1 +
                        0.05 * ifelse(width >= 70 & width <= 500, 1, 0.5))
cand <- cand[order(-score_v2)]
cand$ne_pass <- cand$ne_bg_p95 <= 0.15
fwrite(cand, file.path(base, "DMR_candidates_ranked_v2.csv"))
cat("=== Track 1: 50 DMR 背景回扫 ===\n")
cat("ne_bg_p95<=0.15 通过:", sum(cand$ne_pass, na.rm = TRUE), "/", nrow(cand), "\n")
print(head(cand[, .(chr, start, end, n_cpgs, mean_dbeta = round(mean_dbeta, 3),
                    ne_bg_p95 = round(ne_bg_p95, 3), pre_pos03 = round(pre_pos03, 3),
                    genes, score_v2 = round(score_v2, 3))], 15))

# ---- 轨 2: 精英单探针（区域稀释救助） ----
tcga_dmp <- fread(file.path(base, "TCGA_dmp.csv"))
gse_db <- setNames(gse_dmp$delta_beta, gse_dmp$probe)
elite <- tcga_dmp[adj.P.Val < 0.05 & delta_beta >= 0.50 & case_pos03 >= 0.90 &
                  pre_pos03 >= 0.85 & ctrl_p95 <= 0.40]
elite$gse_dbeta <- gse_db[elite$probe]
elite <- elite[!is.na(gse_dbeta) & gse_dbeta > 0]
pr <- intersect(elite$probe, names(ne_p95))
elite <- elite[probe %in% pr]
elite$ne_bg_mean <- ne_mean[elite$probe]
elite$ne_bg_p95 <- ne_p95[elite$probe]
elite$ne_pass <- elite$ne_bg_p95 <= 0.15
elite$score <- with(elite, 0.25 * pmax(0, 1 - ne_bg_p95 / 0.15) +
                      0.25 * pre_pos03 + 0.20 * case_pos03 + 0.15 * 1 + 0.05 * 1)
elite <- elite[order(-score)]
five <- c("cg05016408", "cg01268824", "cg18675097", "cg23180938", "cg01580681")
elite$is_existing_5 <- elite$probe %in% five
fwrite(elite[, .(probe, chr, pos, gene, island, delta_beta = round(delta_beta, 3),
                 ctrl_p95 = round(ctrl_p95, 3), ne_bg_mean = round(ne_bg_mean, 3),
                 ne_bg_p95 = round(ne_bg_p95, 3), case_pos03 = round(case_pos03, 3),
                 pre_pos03 = round(pre_pos03, 3), gse_dbeta = round(gse_dbeta, 3),
                 ne_pass, is_existing_5, score = round(score, 3))],
       file.path(base, "elite_probes_v2.csv"))
cat("\n=== Track 2: 精英单探针 ===\n")
cat("入围探针:", nrow(elite), "| 其中 ne_bg_p95<=0.15:", sum(elite$ne_pass), "\n")
print(head(elite[, .(probe, gene, delta_beta = round(delta_beta, 3),
                     ne_bg_p95 = round(ne_bg_p95, 3), pre_pos03 = round(pre_pos03, 3),
                     ne_pass, is_existing_5, score = round(score, 3))], 20))
