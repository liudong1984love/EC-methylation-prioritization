# 04_rank_candidates.R — 跨队列一致性 + 入围过滤 + 综合评分排名（v2）
# 一致性定义（方法学说明）：GSE67116 增生组仅 8 例，正式 DMR 检出效能低，
# 故采用探针级效应量方向一致性（区域内 GSE67116 原发癌 vs 增生 mean Δβ > 0），
# 而非要求该区域在 GSE67116 中也被 formally 判为 DMR。
.libPaths(c("C:/Users/ld/Rlibs", .libPaths()))
suppressMessages(library(data.table))
base <- "C:/Users/ld/WorkBuddy/2026-08-19-15-50-19/results"

tcga_dmr <- fread(file.path(base, "TCGA_dmr.csv"))
tcga_dmp <- fread(file.path(base, "TCGA_dmp.csv"))
gse_dmp <- fread(file.path(base, "GSE67116_dmp.csv"))

# 只看高甲基化方向（qMSP 筛查标志物）
tcga_dmr <- tcga_dmr[mean_dbeta > 0]
cat("TCGA hyper DMRs:", nrow(tcga_dmr), "\n")

# 逐区域统计：用 GRanges findOverlaps 一次性完成（避免逐区域全表扫描）
suppressMessages(library(GenomicRanges))
gse_by_probe <- setNames(gse_dmp$delta_beta, gse_dmp$probe)
gse_hyper_p95 <- setNames(gse_dmp$ctrl_p95, gse_dmp$probe)  # 增生组 P95

g_dmp <- GRanges(tcga_dmp$chr, IRanges(tcga_dmp$pos, tcga_dmp$pos))
g_dmr <- GRanges(tcga_dmr$chr, IRanges(tcga_dmr$start, tcga_dmr$end))
ov <- findOverlaps(g_dmr, g_dmp)
dt <- data.table(region = queryHits(ov), tcga_dmp[subjectHits(ov), .(ctrl_p95, case_pos03, pre_pos03, gene, probe)])
dt$gse_db <- gse_by_probe[dt$probe]
dt$gse_hp95 <- gse_hyper_p95[dt$probe]
agg <- dt[, .(
  normal_p95_mean = mean(ctrl_p95, na.rm = TRUE),
  tumor_pos03 = mean(case_pos03, na.rm = TRUE),
  pre_pos03 = mean(pre_pos03, na.rm = TRUE),
  genes = paste(unique(na.omit(unlist(strsplit(gene, ";")))), collapse = ";"),
  gse_dbeta = mean(gse_db, na.rm = TRUE),
  gse_hyper_p95 = mean(gse_hp95, na.rm = TRUE)
), by = region]
setkey(agg, region)
tcga_dmr$region <- seq_len(nrow(tcga_dmr))
tcga_dmr <- merge(tcga_dmr, agg, by = "region", all.x = TRUE, sort = FALSE)

five <- c("cg05016408", "cg01268824", "cg18675097", "cg23180938", "cg01580681")
tcga_dmr$covers_5cpg <- sapply(strsplit(tcga_dmr$probes, ";"), function(p) paste(intersect(p, five), collapse = ";"))

# 诊断：各过滤条件通过数
d <- tcga_dmr
cat("filter diagnostics:\n")
cat("  gse 方向一致 (>0):", sum(d$gse_dbeta > 0, na.rm = TRUE), "\n")
cat("  mean_dbeta>=0.30:", sum(d$mean_dbeta >= 0.30), "\n")
cat("  normal_p95_mean<=0.15:", sum(d$normal_p95_mean <= 0.15, na.rm = TRUE), "\n")
cat("  tumor_pos03>=0.80:", sum(d$tumor_pos03 >= 0.80, na.rm = TRUE), "\n")
cat("  width<=2000:", sum(d$width <= 2000), "\n")

# 双轨入围（方案变更记录 v1.1）：
# 硬指标 Δβ>=0.30 与阳性率>=0.80 在低背景区域存在统计互斥（均值~0.32 的区域不可能 80% 样本>0.3），
# 改为：A 轨（mean_dbeta>=0.30）或 B 轨（区域阳性率>=0.50），其余硬条件不变；精细排序交给综合评分。
cand <- d[(mean_dbeta >= 0.30 | tumor_pos03 >= 0.50) & normal_p95_mean <= 0.15 &
          gse_dbeta > 0 & n_cpgs >= 3 & width <= 2000 & !is.na(gse_dbeta)]
cat("candidates after filters:", nrow(cand), "\n")

if (nrow(cand) > 0) {
  cand$score <- with(cand, 0.25 * pmax(0, 1 - normal_p95_mean / 0.15) +
                       0.25 * pre_pos03 + 0.20 * tumor_pos03 + 0.15 * 1 +
                       0.05 * ifelse(width >= 70 & width <= 500, 1, 0.5))
  cand <- cand[order(-score)]
  fwrite(cand, file.path(base, "DMR_candidates_ranked_v1.csv"))
  cat("saved ranked:", nrow(cand), "\n")
  print(head(cand[, .(chr, start, end, n_cpgs, width, mean_dbeta = round(mean_dbeta, 3),
                      normal_p95_mean = round(normal_p95_mean, 3), tumor_pos03 = round(tumor_pos03, 3),
                      pre_pos03 = round(pre_pos03, 3), gse_dbeta = round(gse_dbeta, 3),
                      genes, covers_5cpg, score = round(score, 3))], 25))
}
