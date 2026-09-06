# 03_dmr.R — 统一 DMP + DMR 区域合并流程（limma + comb-p 思路，确定性、无外部依赖）
# 用法: Rscript 03_dmr.R GSE67116 | TCGA
.libPaths(c("C:/Users/ld/Rlibs", .libPaths()))
suppressMessages({library(limma); library(data.table);
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})

args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]
base <- "C:/Users/ld/WorkBuddy/2026-08-19-15-50-19"
outdir <- file.path(base, "results")

anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

# ---- DMR 区域合并函数 ----
find_dmr <- function(dmp, fdr = 0.05, maxgap = 500, mincpg = 3) {
  sig <- dmp[!is.na(dmp$adj.P.Val) & dmp$adj.P.Val < fdr & !is.na(dmp$chr) & !is.na(dmp$pos), ]
  sig$dir <- sign(sig$delta_beta)
  sig <- sig[sig$dir != 0, ]
  sig <- sig[order(sig$chr, sig$pos), ]
  regions <- list(); cur <- NULL
  for (i in seq_len(nrow(sig))) {
    p <- sig[i, ]
    if (is.null(cur) || p$chr != cur$chr || p$dir != cur$dir || (p$pos - cur$end) > maxgap) {
      if (!is.null(cur)) regions[[length(regions) + 1]] <- cur
      cur <- list(chr = p$chr, start = p$pos, end = p$pos, dir = p$dir, idx = i)
    } else {
      cur$end <- p$pos; cur$idx <- c(cur$idx, i)
    }
  }
  if (!is.null(cur)) regions[[length(regions) + 1]] <- cur
  res <- do.call(rbind, lapply(regions, function(r) {
    if (length(r$idx) < mincpg) return(NULL)
    sub <- sig[r$idx, ]
    z <- sum(sign(sub$delta_beta) * qnorm(sub$P.Value / 2, lower.tail = FALSE)) / sqrt(nrow(sub))
    data.frame(chr = r$chr, start = min(sub$pos), end = max(sub$pos),
               width = max(sub$pos) - min(sub$pos), n_cpgs = nrow(sub),
               mean_dbeta = mean(sub$delta_beta),
               stouffer_p = 2 * pnorm(abs(z), lower.tail = FALSE),
               probes = paste(sub$probe, collapse = ";"), stringsAsFactors = FALSE)
  }))
  if (!is.null(res)) res <- res[order(res$stouffer_p), ]
  res
}

add_probe_stats <- function(dmp, betas, case_idx, ctrl_idx) {
  dmp$beta_case <- rowMeans(betas[, case_idx, drop = FALSE], na.rm = TRUE)
  dmp$beta_ctrl <- rowMeans(betas[, ctrl_idx, drop = FALSE], na.rm = TRUE)
  dmp$delta_beta <- dmp$beta_case - dmp$beta_ctrl
  ctrl_sorted <- apply(betas[, ctrl_idx, drop = FALSE], 1, sort, na.last = TRUE)
  dmp$ctrl_p95 <- apply(betas[, ctrl_idx, drop = FALSE], 1, function(x) quantile(x, 0.95, na.rm = TRUE))
  dmp$case_pos03 <- rowMeans(betas[, case_idx, drop = FALSE] > 0.3, na.rm = TRUE)
  dmp
}

if (ds == "GSE67116") {
  betas <- readRDS(file.path(outdir, "GSE67116_betas.rds"))
  M <- readRDS(file.path(outdir, "GSE67116_M.rds"))
  sdrf <- read.delim(file.path(base, "data/GSE67116/E-GEOD-67116.sdrf.txt"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  sdrf$GSM <- sub(" .*", "", sdrf[["Source Name"]])
  sdrf <- sdrf[!duplicated(sdrf$GSM), ]
  tt <- setNames(sdrf[["Characteristics [tumor.type]"]], sdrf$GSM)
  gsm <- sub("_.*$", "", colnames(M))
  group <- ifelse(tt[gsm] == "Primary", "Primary",
           ifelse(tt[gsm] == "hyperplasia", "Hyperplasia", "Other"))
  keep <- group %in% c("Primary", "Hyperplasia")
  M <- M[, keep]; betas <- betas[, keep]; group <- factor(group[keep], levels = c("Hyperplasia", "Primary"))
  cat("contrast:", paste(table(group), names(table(group)), sep = "=", collapse = ", "), "\n")
  design <- model.matrix(~group)
  fit <- eBayes(lmFit(M, design))
  dmp <- topTable(fit, coef = "groupPrimary", number = Inf, sort.by = "none")
  dmp$probe <- rownames(dmp)
  dmp <- add_probe_stats(dmp, betas, group == "Primary", group == "Hyperplasia")
  common <- intersect(dmp$probe, rownames(anno))
  dmp <- dmp[dmp$probe %in% common, ]
  dmp$chr <- anno[dmp$probe, "chr"]; dmp$pos <- anno[dmp$probe, "pos"]
  dmp$gene <- anno[dmp$probe, "UCSC_RefGene_Name"]; dmp$island <- anno[dmp$probe, "Relation_to_Island"]
  fwrite(dmp, file.path(outdir, "GSE67116_dmp.csv"))
  dmr <- find_dmr(dmp)
  fwrite(dmr, file.path(outdir, "GSE67116_dmr.csv"))
  cat("DMRs:", nrow(dmr), "| with |mean dBeta|>=0.3:", sum(abs(dmr$mean_dbeta) >= 0.3), "\n")
}

if (ds == "TCGA") {
  cat("reading Xena beta matrix...\n")
  tb <- fread(file.path(base, "data/TCGA_UCEC_methylation450_beta.tsv"))
  probes <- tb[[1]]; tb[[1]] <- NULL
  betas <- as.matrix(tb); rownames(betas) <- probes; rm(tb, probes); gc()
  tumor <- grep("-01$", colnames(betas)); normal <- grep("-11$", colnames(betas))
  cat("tumor:", length(tumor), "normal:", length(normal), "\n")
  betas <- betas[, c(normal, tumor)]
  grp <- factor(rep(c("Normal", "Tumor"), c(length(normal), length(tumor))), levels = c("Normal", "Tumor"))
  # 绝经前亚组
  clin <- fread(file.path(base, "data/UCEC_clinicalMatrix.tsv"))
  meno <- setNames(clin$menopause_status, substr(clin$sampleID, 1, 15))
  is_pre <- startsWith(meno[colnames(betas)], "Pre"); is_pre[is.na(is_pre)] <- FALSE
  pre_idx <- which(grp == "Tumor" & is_pre)
  cat("Pre tumors:", length(pre_idx), "\n")
  M <- log2(betas / (1 - betas + 1e-4) + 1e-4)
  design <- model.matrix(~grp)
  fit <- eBayes(lmFit(M, design))
  dmp <- topTable(fit, coef = "grpTumor", number = Inf, sort.by = "none")
  dmp$probe <- rownames(dmp)
  dmp <- add_probe_stats(dmp, betas, grp == "Tumor", grp == "Normal")
  dmp$pre_pos03 <- rowMeans(betas[, pre_idx, drop = FALSE] > 0.3, na.rm = TRUE)
  common <- intersect(dmp$probe, rownames(anno))
  dmp <- dmp[dmp$probe %in% common, ]
  dmp$chr <- anno[dmp$probe, "chr"]; dmp$pos <- anno[dmp$probe, "pos"]
  dmp$gene <- anno[dmp$probe, "UCSC_RefGene_Name"]; dmp$island <- anno[dmp$probe, "Relation_to_Island"]
  fwrite(dmp, file.path(outdir, "TCGA_dmp.csv"))
  dmr <- find_dmr(dmp)
  fwrite(dmr, file.path(outdir, "TCGA_dmr.csv"))
  cat("DMRs:", nrow(dmr), "| with |mean dBeta|>=0.3:", sum(abs(dmr$mean_dbeta) >= 0.3), "\n")
}
