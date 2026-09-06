# v9.11: recompute ALL GSE73949 (healthy endometrium, n=17) background stats from the
# declared pipeline (minfi noob from IDAT, results/GSE73949_betas.rds), with bootstrap CIs,
# and check the funnel gate count (elite probes passing P95<=0.15) under noob.
b <- readRDS("results/GSE73949_betas.rds")
p95 <- function(x) as.numeric(quantile(x, 0.95, na.rm = TRUE, type = 7))

panel23 <- read.csv("results/wetlab_panel_23_primer_annotation.csv")$probe
subs <- c("cg23164203","cg20275132","cg25666433")
internal5 <- c("cg05016408","cg01268824","cg18675097","cg23180938","cg01580681")
kit9 <- c("cg23180938","cg13495205","cg11835068","cg15484532","cg03502002","cg04534765","cg10390058","cg25390440","cg04453971")
all_probes <- unique(c(panel23, subs, internal5, kit9))

set.seed(20260906)
B <- 10000
rows <- lapply(all_probes, function(p) {
  if (!p %in% rownames(b)) return(NULL)
  x <- as.numeric(b[p, ]); x <- x[!is.na(x)]
  boot <- replicate(B, p95(sample(x, replace = TRUE)))
  data.frame(probe = p, n = length(x),
             median = median(x), p90 = as.numeric(quantile(x, 0.9, type = 7)),
             p95 = p95(x), max = max(x), prop_gt_0.15 = mean(x > 0.15),
             boot_lo = as.numeric(quantile(boot, 0.025)), boot_hi = as.numeric(quantile(boot, 0.975)))
})
df <- do.call(rbind, rows)
df[, -1] <- round(df[, -1], 4)
write.csv(df, "results/v911_gse73949_noob_stats.csv", row.names = FALSE)
print(df[df$probe %in% c("cg24589459","cg07495363","cg27577527","cg23180938","cg23164203","cg20275132","cg25666433"), ])

# funnel check: elite probes passing healthy-endometrium P95<=0.15 under noob
elite <- read.csv("results/elite_probes_v2.csv")
ecol <- colnames(elite)[grepl("probe|cg", colnames(elite), ignore.case = TRUE)][1]
eprobes <- elite[[ecol]]
cat("elite probes:", length(eprobes), " column:", ecol, "\n")
present <- eprobes[eprobes %in% rownames(b)]
cat("present in GSE73949 rds:", length(present), "\n")
gate <- apply(b[present, , drop = FALSE], 1, p95) <= 0.15
cat("passing P95<=0.15 under noob:", sum(gate), "\n")
