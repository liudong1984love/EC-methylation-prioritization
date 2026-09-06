# 41_bumphunter_leads.R — cross-check BOLL/ZSCAN12 regions with bumphunter (established DMR method)
# Input: results/v97_bumphunter_chr26_M.tsv (chr, pos, M-values per sample), results/v97_bumphunter_pheno.tsv
suppressMessages({library(bumphunter); library(GenomicRanges)})
BASE <- "C:/Users/ld/WorkBuddy/2026-08-19-15-50-19"
M <- read.delim(file.path(BASE, "results/v97_bumphunter_chr26_M.tsv"), check.names = FALSE)
phe <- read.delim(file.path(BASE, "results/v97_bumphunter_pheno.tsv"))
# first col = probe id (unnamed), then chr, pos, then samples
colnames(M)[1] <- "probe"
stopifnot(all(colnames(M)[-(1:3)] == phe$sample))
mat <- as.matrix(M[, -(1:3)])
design <- model.matrix(~group, data = phe)

for (chr_ in c("chr2", "chr6")) {
  idx <- which(M$chr == chr_)
  cat("== running bumphunter on", chr_, "with", length(idx), "probes\n")
  cl <- clusterMaker(M$chr[idx], M$pos[idx], maxGap = 500)
  set.seed(42)
  bh <- bumphunter(mat[idx, ], design, chr = M$chr[idx], pos = M$pos[idx],
                   cluster = cl, B = 500, smooth = FALSE, pickCutoff = TRUE, verbose = FALSE)
  tab <- bh$table
  tab <- tab[order(tab$p.value), ]
  write.csv(tab, file.path(BASE, sprintf("results/v97_bumphunter_%s.csv", chr_)), row.names = FALSE)
  # report regions near lead loci
  if (chr_ == "chr2") { lo <- 197.5e6; hi <- 198.05e6; lab <- "BOLL" }
  else { lo <- 28.35e6; hi <- 28.45e6; lab <- "ZSCAN12" }
  hit <- tab[tab$start < hi & tab$end > lo, ]
  cat(sprintf("  %s: %d DMR(s) overlapping locus window\n", lab, nrow(hit)))
  if (nrow(hit)) print(hit[, c("chr", "start", "end", "value", "area", "cluster", "p.value", "fwer", "clusterL")])
}
cat("DONE\n")
