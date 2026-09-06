# v9.9 -> v9.10: sixth-round reviewer revisions (12 points)
from docx import Document

SRC = "Manuscript_EC_methylation_panel_v9.9.docx"
DST = "Manuscript_EC_methylation_panel_v9.10.docx"
d = Document(SRC)

log = []

def set_text(p, text):
    for r in list(p.runs):
        r._element.getparent().remove(r._element)
    p.add_run(text)

def sub(old, new, tag, count_expected=None):
    n = 0
    for p in d.paragraphs:
        if old in p.text:
            set_text(p, p.text.replace(old, new))
            n += 1
    log.append(("OK  " if n else "MISS") + f" [{tag}] x{n}")

# --- P1: empirical P<0.002 -> P<=0.002 with +1 correction ---
sub("chromosome-wise empirical P<0.002; not a genome-wide FWER",
    "chromosome-wise empirical P\u22640.002 (k=0 of B=500 permutations, reported with the +1 finite-permutation correction (k+1)/(B+1)=1/501\u22480.002); not a genome-wide FWER",
    "p1-limitations")
sub("(chromosome-wise empirical P<0.002, the minimum reportable at B=500)",
    "(chromosome-wise empirical P\u22640.002 after the +1 finite-permutation correction)",
    "p1-results45")
sub("empirical P<0.002", "empirical P\u22640.002", "p1-residual")

# --- P2: Figure 5 legend — three lead probes undefined ---
sub("step indicators: pass/fail status of the three lead probes.",
    "step indicators: pass/fail status of the two anchor probes and the BOLL partner probe (cg24589459, cg27577527 and cg07495363, respectively).",
    "p2-fig5")

# --- P3: Table 7 legend thresholds overlap ---
sub("Table 7. Sampling-compartment suitability (background P95 \u03b2; \u2713\u22640.15, \u25b30.15\u20130.30, \u2717>0.30).",
    "Table 7. Sampling-compartment suitability (background P95 \u03b2; \u2713\u22640.15, \u25b3>0.15 and \u22640.30, \u2717>0.30). *Borderline: point-estimate pass whose bootstrap 95% confidence interval crosses the 0.15 gate.",
    "p3-table7-legend")

# --- P4: cycle sentence ---
sub("excluding cycle synchronisation as an explanation for the low benign background of BOLL.",
    "arguing against a substantial early- to mid-secretory phase effect on BOLL background; proliferative-phase, menstrual and anovulatory states were not assayed.",
    "p4-cycle")

# --- P9: SNP-free scoping ---
sub("all SNP-free and 22/23 embedded in multi-CpG DMRs",
    "none of them carrying annotated common SNPs at the CpG or single-base-extension sites under the Zhou et al. annotation [28], and 22/23 embedded in multi-CpG DMRs",
    "p9-snpfree")

# --- P8: platform-independent ---
sub("the qualitative conclusion is platform-independent, whereas the raw two- to three-fold contrast was inflated by processing.",
    "the qualitative difference persisted in this same-platform, same-reported-normalisation comparison, whereas the raw two- to three-fold contrast was inflated by processing.",
    "p8-platform")

# --- P6: face validity ---
sub("Convergence with independent wet-lab programmes is the strongest internal evidence of validity:",
    "Re-identification of established markers supports the face validity of the prioritisation framework:",
    "p6-facevalidity")

# --- P5: BOLL heterogeneity ---
sub("most likely reflects genuine regional methylation heterogeneity rather than probe artefact, reinforcing the need for region-level assay design.",
    "is consistent with regional methylation heterogeneity\u2014although probe-specific technical effects (probe chemistry, unannotated sequence variants or residual cross-reactivity) cannot be excluded\u2014reinforcing the need for region-level assay design.",
    "p5-heterogeneity")

# --- P7: GSE223817 gate computation in Methods ---
sub("The default elimination threshold was background P95 \u03b2>0.15 under a union rule requiring every available cohort to pass",
    "For GSE223817, the union-rule gate was computed on the 347 control samples; the 637 endometriosis samples were evaluated descriptively, with per-group P95 values reported in Supplementary Data. The default elimination threshold was background P95 \u03b2>0.15 under a union rule requiring every available cohort to pass",
    "p7-methods-gse223817")

# --- P10: reference 28 DOI/PMID ---
sub("Nucleic Acids Res. 2017;45(4):e22.",
    "Nucleic Acids Res. 2017;45(4):e22. doi:10.1093/nar/gkw967. PMID:27924034.",
    "p10-ref28")

# --- P11: code availability at submission ---
sub("and will be deposited in a public repository upon acceptance.",
    "and will be deposited in a public repository with a versioned DOI, sessionInfo output and environment lock files at the time of submission.",
    "p11-code")

d.save(DST)
print("\n".join(log))
print("saved", DST)
