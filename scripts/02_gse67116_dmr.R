# 02_gse67116_dmr.R — GSE67116: 原发癌(33) vs 增生(8) 的 DMP + DMRcate DMR
.libPaths(c("C:/Users/ld/Rlibs", .libPaths()))
suppressMessages({library(minfi); library(limma); library(DMRcate);
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19); library(data.table)})

base <- "C:/Users/ld/WorkBuddy/2026-08-19-15-50-19"
outdir <- file.path(base, "results")

betas <- readRDS(file.path(outdir, "GSE67116_betas.rds"))
M <- readRDS(file.path(outdir, "GSE67116_M.rds"))

# 元数据: tumor.type 从 sdrf 取
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
cat("contrast samples:", paste(table(group), names(table(group)), sep="=", collapse=", "), "\n")

# limma DMP
design <- model.matrix(~group)
fit <- eBayes(lmFit(M, design))
dmp <- topTable(fit, coef = "groupPrimary", number = Inf, sort.by = "none")
dmp$probe <- rownames(dmp)
dmp$beta_primary <- rowMeans(betas[, group == "Primary"])
dmp$beta_hyper <- rowMeans(betas[, group == "Hyperplasia"])
dmp$delta_beta <- dmp$beta_primary - dmp$beta_hyper
fwrite(dmp, file.path(outdir, "GSE67116_dmp_primary_vs_hyperplasia.csv"))
cat("DMP table saved; FDR<0.05 & |dBeta|>0.2:", sum(dmp$adj.P.Val < 0.05 & abs(dmp$delta_beta) > 0.2), "\n")

# DMRcate
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
myanno <- cpg.annotate("array", M, what = "M", arraytype = "450K",
                       analysis.type = "differential", design = design,
                       coef = "groupPrimary", fdr = 0.05)
dmrs <- dmrcate(myanno, lambda = 1000, C = 2)
res <- extractRanges(dmrs, genome = "hg19")
res_df <- as.data.frame(res)
fwrite(res_df, file.path(outdir, "GSE67116_dmr_primary_vs_hyperplasia.csv"))
cat("DMRcate DMRs:", nrow(res_df), "\n")
if (nrow(res_df)) print(head(res_df[, c("seqnames","start","end","no.cpgs","meandiff","minfdr")], 10))
