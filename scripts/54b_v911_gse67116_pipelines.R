# v9.11 Analysis C part 2: GSE67116 tumour-signal rank consistency across pipelines
suppressMessages(library(minfi))
panel23 <- read.csv("results/wetlab_panel_23_primer_annotation.csv")$probe
subs <- c("cg23164203","cg20275132","cg25666433")
probes26 <- c(panel23, subs)

pheno <- read.csv("results/GSE67116_pheno.csv", stringsAsFactors = FALSE)
tum_idat <- pheno$idat_name[pheno$sample_type == "Tumor"]
hyp_idat <- pheno$idat_name[pheno$sample_type == "hyperplasia"]
cat("tumours:", length(tum_idat), "| hyperplasia:", length(hyp_idat), "\n")

rg <- read.metharray.exp("data/GSE67116", force = TRUE)
cat("RGset cols:", ncol(rg), "\n")
# match colnames to idat basenames
cn <- colnames(rg)
tum_idx <- which(cn %in% tum_idat | sub(".*_(GSM[0-9]+)$", "\\1", cn) %in% pheno$GSM[pheno$sample_type=="Tumor"])
# robust fallback: match by GSM prefix
gsm_of_col <- sub("^(GSM[0-9]+)_.*$", "\\1", cn)
tum_idx <- which(gsm_of_col %in% pheno$GSM[pheno$sample_type == "Tumor"])
hyp_idx <- which(gsm_of_col %in% pheno$GSM[pheno$sample_type == "hyperplasia"])
cat("matched tumour cols:", length(tum_idx), "| hyperplasia cols:", length(hyp_idx), "\n")

get_betas <- function(rg, method) {
  if (method == "raw")      mset <- preprocessRaw(rg)
  if (method == "noob")     mset <- preprocessNoob(rg)
  if (method == "quantile") mset <- suppressWarnings(preprocessQuantile(rg))
  getBeta(mset)
}

res <- list()
for (m in c("raw","noob","quantile")) {
  b <- get_betas(rg, m)
  b26 <- b[rownames(b) %in% probes26, , drop = FALSE]
  tum_mean <- rowMeans(b26[, tum_idx, drop = FALSE], na.rm = TRUE)
  hyp_mean <- rowMeans(b26[, hyp_idx, drop = FALSE], na.rm = TRUE)
  res[[m]] <- data.frame(probe = rownames(b26), tum = tum_mean, hyp = hyp_mean, dbeta = tum_mean - hyp_mean)
  cat(m, "probes:", nrow(b26), "\n")
}
df <- data.frame(probe = res$noob$probe,
                 tum_raw = res$raw$tum, tum_noob = res$noob$tum, tum_quantile = res$quantile$tum,
                 dbeta_raw = res$raw$dbeta, dbeta_noob = res$noob$dbeta, dbeta_quantile = res$quantile$dbeta)
write.csv(df, "results/v911_pipeline_tumour_gse67116.csv", row.names = FALSE)
cat("tumour mean Spearman: raw-noob", cor(df$tum_raw, df$tum_noob, method="spearman"),
    "| noob-quantile", cor(df$tum_noob, df$tum_quantile, method="spearman"),
    "| raw-quantile", cor(df$tum_raw, df$tum_quantile, method="spearman"), "\n")
cat("dbeta Spearman: raw-noob", cor(df$dbeta_raw, df$dbeta_noob, method="spearman"),
    "| noob-quantile", cor(df$dbeta_noob, df$dbeta_quantile, method="spearman"),
    "| raw-quantile", cor(df$dbeta_raw, df$dbeta_quantile, method="spearman"), "\n")
cat("dbeta range: raw", range(df$dbeta_raw), " noob", range(df$dbeta_noob), " quantile", range(df$dbeta_quantile), "\n")
