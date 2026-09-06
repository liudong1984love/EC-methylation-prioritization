# 01_minfi_preprocess.R — minfi 预处理：IDAT -> noob -> beta/M/detP + QC
# 用法: Rscript 01_minfi_preprocess.R <GSE> 
# 输入: data/<GSE>/*.idat + data/<GSE>/E-GEOD-*.sdrf.txt
# 输出: results/<GSE>_betas.rds, <GSE>_M.rds, <GSE>_detP.rds, <GSE>_pheno.csv, <GSE>_qc.csv

args <- commandArgs(trailingOnly = TRUE)
gse <- args[1]
.libPaths(c("C:/Users/ld/Rlibs", .libPaths()))
suppressMessages({library(minfi); library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})

base <- "C:/Users/ld/WorkBuddy/2026-08-19-15-50-19"
datadir <- file.path(base, "data", gse)
outdir <- file.path(base, "results")

# ---- 1. 读 IDAT ----
rgSet <- read.metharray.exp(datadir, recursive = TRUE, force = TRUE)
cat("read", ncol(rgSet), "samples x", nrow(rgSet), "probes\n")

# ---- 2. 读 SDRF 元数据 ----
sdrf_file <- list.files(datadir, pattern = "sdrf\\.txt$", full.names = TRUE)[1]
sdrf <- read.delim(sdrf_file, stringsAsFactors = FALSE, check.names = FALSE)
gsc <- grep("^Source Name$", colnames(sdrf))
sdrf$GSM <- sub(" .*", "", sdrf[[gsc]])
dedup <- !duplicated(sdrf$GSM)
sdrf <- sdrf[dedup, ]
type_col <- grep("sample type|condition|disease", colnames(sdrf), value = TRUE, ignore.case = TRUE)[1]
sentrix_col <- grep("sentrix id", colnames(sdrf), value = TRUE, ignore.case = TRUE)[1]
phase_col <- grep("menstrual phase|cycle phase", colnames(sdrf), value = TRUE, ignore.case = TRUE)[1]
patient_col <- grep("patient id", colnames(sdrf), value = TRUE, ignore.case = TRUE)[1]
pheno <- data.frame(GSM = sdrf$GSM,
                    sample_type = sdrf[[type_col]],
                    sentrix = if (!is.na(sentrix_col)) sdrf[[sentrix_col]] else NA,
                    menstrual_phase = if (!is.na(phase_col)) sdrf[[phase_col]] else NA,
                    patient_id = if (!is.na(patient_col)) sdrf[[patient_col]] else NA,
                    stringsAsFactors = FALSE)
# idat 列名形如 GSM1639323_8963303018_R06C02
pheno$idat_name <- sapply(pheno$GSM, function(g) {
  hit <- grep(paste0("^", g, "_"), colnames(rgSet), value = TRUE)
  if (length(hit) == 1) hit else NA
})
pheno <- pheno[!is.na(pheno$idat_name), ]
rgSet <- rgSet[, pheno$idat_name]
cat("matched", nrow(pheno), "samples with metadata\n")
print(table(pheno$sample_type))

# ---- 3. QC: 检测 p 值 ----
detP <- detectionP(rgSet)
qc <- data.frame(
  sample = pheno$GSM,
  sample_type = pheno$sample_type,
  mean_detP = colMeans(detP),
  frac_failed_0.01 = colMeans(detP > 0.01)
)
write.csv(qc, file.path(outdir, paste0(gse, "_qc.csv")), row.names = FALSE)
keep_samples <- qc$frac_failed_0.01 < 0.10
cat("samples passing QC:", sum(keep_samples), "/", length(keep_samples), "\n")
if (any(!keep_samples)) {
  cat("dropped:", paste(qc$sample[!keep_samples], collapse = ", "), "\n")
  rgSet <- rgSet[, keep_samples]; detP <- detP[, keep_samples]; pheno <- pheno[keep_samples, ]
}

# ---- 4. noob 归一化 ----
mset <- preprocessNoob(rgSet, dyeMethod = "single")
cat("noob done\n")

# ---- 5. 探针过滤 ----
betas <- getBeta(mset)
fail <- rowMeans(detP > 0.01) > 0.10
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sex <- rownames(betas) %in% rownames(anno)[anno$chr %in% c("chrX", "chrY")]
gset <- mapToGenome(mset)
mset_snp <- dropLociWithSnps(gset, snps = c("SBE", "CpG"), maf = 0.05)
snp <- !rownames(betas) %in% rownames(getBeta(mset_snp))
drop <- fail | sex | snp
cat("probes dropped: failP", sum(fail), "| sex", sum(sex), "| snp", sum(snp), "| total", sum(drop), "\n")

betas_f <- betas[!drop, ]
detP_f <- detP[!drop, ]
M_f <- log2(betas_f / (1 - betas_f + 1e-6) + 1e-6)

# ---- 6. 保存 ----
saveRDS(betas_f, file.path(outdir, paste0(gse, "_betas.rds")))
saveRDS(M_f, file.path(outdir, paste0(gse, "_M.rds")))
saveRDS(detP_f, file.path(outdir, paste0(gse, "_detP.rds")))
write.csv(pheno, file.path(outdir, paste0(gse, "_pheno.csv")), row.names = FALSE)
cat("saved:", nrow(betas_f), "probes x", ncol(betas_f), "samples\n")
