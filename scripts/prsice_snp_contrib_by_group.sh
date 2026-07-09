#!/usr/bin/env bash
set -euo pipefail

# PRSice per-SNP contribution by group (mean beta*dosage)
# Uses PLINK2 --freq --within (no allele forcing). R maps PLINK allele freq -> effect allele.
#
# REQUIREMENTS: plink2, Rscript
#
# INPUTS:
#  --summary  PRSice .summary
#  --allscore PRSice .all_score (must include SNP, A1 (effect allele), BETA, P)
#  --pfile or --bfile  target genotypes (one only)
#  --MYGROUPS  MYGROUPS file: FID IID GROUP (no header)
#  --out    output prefix
#  [--pt]   optional override p-value threshold
#
# Example:
#  ./prsice_snp_contrib_by_group.sh --summary prsice.summary --allscore prsice.all_score \
#    --bfile target --MYGROUPS MYGROUPS.txt --out drivers

SUMMARY=""
ALLSCORE=""
PFILE=""
BFILE=""
MYGROUPS=""
OUT="prsice_drivers"
PT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --summary) SUMMARY="$2"; shift 2;;
    --allscore) ALLSCORE="$2"; shift 2;;
    --pfile) PFILE="$2"; shift 2;;
    --bfile) BFILE="$2"; shift 2;;
    --MYGROUPS) MYGROUPS="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --pt) PT="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

if [[ -z "$SUMMARY" || -z "$ALLSCORE" || -z "$MYGROUPS" ]]; then
  echo "Missing required args. Need --summary, --allscore, --MYGROUPS, and one of --pfile/--bfile." >&2
  exit 1
fi
if [[ -z "$PFILE" && -z "$BFILE" ]]; then
  echo "Provide either --pfile or --bfile." >&2
  exit 1
fi
if [[ -n "$PFILE" && -n "$BFILE" ]]; then
  echo "Provide only one of --pfile or --bfile, not both." >&2
  exit 1
fi

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT
# WORKDIR=$(mktemp -d -p . "prsice_debug_$(date +%Y%m%d_%H%M%S)_XXXX")
# echo "WORKDIR = $WORKDIR"


echo $MYGROUPS

# 1) Extract SNP list, effect allele (A1), BETA and PT from PRSice outputs
Rscript --vanilla - "$SUMMARY" "$ALLSCORE" "$WORKDIR" "$PT" <<'RSCRIPT'
args <- commandArgs(trailingOnly=TRUE)
summary_path <- args[1]; allscore_path <- args[2]; workdir <- args[3]; pt_override <- args[4]
pt_override <- if (nchar(pt_override) == 0) NA else as.numeric(pt_override)

read_flex <- function(path) {
  if (requireNamespace("data.table", quietly=TRUE)) data.table::fread(path, data.table=FALSE)
  else {
    x <- try(read.table(path, header=TRUE, stringsAsFactors=FALSE), silent=TRUE)
    if (!inherits(x, "try-error")) return(x)
    x <- try(read.delim(path, header=TRUE, stringsAsFactors=FALSE), silent=TRUE)
    if (!inherits(x, "try-error")) return(x)
    x <- try(read.csv(path, header=TRUE, stringsAsFactors=FALSE), silent=TRUE)
    if (!inherits(x, "try-error")) return(x)
    stop("Could not parse file: ", path)
  }
}

summ <- read_flex(summary_path)
alls <- read_flex(allscore_path)

# detect threshold column
thr_col <- intersect(names(summ), c("Threshold","P_T","PT","P","thresh","threshold"))
if (length(thr_col)==0) stop("Couldn't find threshold column in summary. Columns: ", paste(names(summ), collapse=", "))
thr_col <- thr_col[1]

pt <- pt_override
if (is.na(pt)) {
  if ("R2" %in% names(summ)) pt <- summ[[thr_col]][which.max(summ[["R2"]])]
  else if ("P" %in% names(summ)) pt <- summ[[thr_col]][which.min(summ[["P"]])]
  else pt <- summ[[thr_col]][1]
}
if (!is.finite(pt)) stop("PT parsed as non-finite.")
message("Using PT = ", pt)

# find columns in allscore
snp_col <- intersect(names(alls), c("SNP","rsID","rsid","ID","MarkerName","marker"))
a1_col  <- intersect(names(alls), c("A1","Allele1","Effect_Allele","EA","ALLELE1"))
beta_col<- intersect(names(alls), c("BETA","Beta","beta","Effect","OR","or"))
p_col   <- intersect(names(alls), c("P","p","PVAL","Pvalue","p.value","P-value"))

if (length(snp_col)==0) stop("Couldn't find SNP column in all_score. Columns: ", paste(names(alls), collapse=", "))
if (length(a1_col)==0)  stop("Couldn't find effect allele (A1) column in all_score. Columns: ", paste(names(alls), collapse=", "))
if (length(beta_col)==0) stop("Couldn't find beta column in all_score. Columns: ", paste(names(alls), collapse=", "))

snp_col <- snp_col[1]; a1_col <- a1_col[1]; beta_col <- beta_col[1]
p_col <- if (length(p_col)>0) p_col[1] else NA_character_

if (!is.na(p_col)) keep <- is.finite(alls[[p_col]]) & (alls[[p_col]] <= pt) else keep <- rep(TRUE, nrow(alls))
sub <- alls[keep, , drop=FALSE]
sub <- sub[is.finite(sub[[beta_col]]) & !is.na(sub[[snp_col]]) & !is.na(sub[[a1_col]]), , drop=FALSE]
if (nrow(sub) == 0) stop("After filtering, zero SNPs remain. Check PT and all_score p-value column.")

# write snplist and weights (SNP, A1, BETA)
snplist_path <- file.path(workdir, "snplist.txt"); writeLines(unique(sub[[snp_col]]), snplist_path)
weights_path <- file.path(workdir, "weights.txt")
w <- unique(sub[, c(snp_col, a1_col, beta_col)])
colnames(w) <- c("SNP","A1","BETA")
write.table(w, weights_path, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

writeLines(as.character(pt), file.path(workdir, "PT.txt"))
message("Wrote snplist (", nrow(w), " SNPs) and weights to ", workdir)
RSCRIPT

SNPLIST="$WORKDIR/snplist.txt"
WEIGHTS="$WORKDIR/weights.txt"
PTFILE="$WORKDIR/PT.txt"

# 2) Run plink2 to get allele frequencies per group WITHOUT --within.
# Create keep files for the first two groups listed in groups.txt.

# Identify the first two distinct group labels (column 3) in groups file
#echo "DEBUG: running plink2 for group $gA"
#read -r gA gB < <(awk '{print $3}' "$MYGROUPS" | awk '!seen[$0]++' | head -n 2 | tr '\n' ' ')
#echo "DEBUG: running plink2 for group $gA"

echo "DEBUG: detecting groups from $MYGROUPS"

[[ -f "$MYGROUPS" ]] || { echo "ERROR: MYGROUPS file not found: $MYGROUPS" >&2; exit 1; }
[[ -s "$MYGROUPS" ]] || { echo "ERROR: MYGROUPS file is empty: $MYGROUPS" >&2; exit 1; }

# Avoid silent exit under set -e by capturing output explicitly
set +e
groups_found=$(awk '{print $3}' "$MYGROUPS" | awk '!seen[$0]++' | head -n 2)
rc=$?
set -e
if [[ $rc -ne 0 ]]; then
  echo "ERROR: failed while extracting group labels from $MYGROUPS" >&2
  exit 1
fi

gA=$(echo "$groups_found" | sed -n '1p')
gB=$(echo "$groups_found" | sed -n '2p')

echo "DEBUG: raw groups_found:"
printf '%s\n' "$groups_found" | nl -ba

if [[ -z "${gA:-}" || -z "${gB:-}" ]]; then
  echo "ERROR: Could not detect two groups from $MYGROUPS (need >=2 distinct labels in column 3)." >&2
  exit 1
fi

echo "Detected groups: $gA vs $gB"







if [[ -z "${gA:-}" || -z "${gB:-}" ]]; then
  echo "ERROR: Could not detect two groups from $MYGROUPS (expected 3rd column to be group label)." >&2
  exit 1
fi

echo "Detected groups: $gA vs $gB"

KEEP_A="$WORKDIR/keep_${gA}.txt"
KEEP_B="$WORKDIR/keep_${gB}.txt"

# keep files are 2 columns: FID IID
awk -v G="$gA" '$3==G{print $1, $2}' "$MYGROUPS" > "$KEEP_A"
awk -v G="$gB" '$3==G{print $1, $2}' "$MYGROUPS" > "$KEEP_B"

PLINK_OUT_A="$WORKDIR/freq_${gA}"
PLINK_OUT_B="$WORKDIR/freq_${gB}"

if [[ -n "$PFILE" ]]; then
  plink2 --pfile "$PFILE" \
    --extract "$SNPLIST" \
    --keep "$KEEP_A" \
    --freq \
    --out "$PLINK_OUT_A"

  plink2 --pfile "$PFILE" \
    --extract "$SNPLIST" \
    --keep "$KEEP_B" \
    --freq \
    --out "$PLINK_OUT_B"
else
  plink2 --bfile "$BFILE" \
    --extract "$SNPLIST" \
    --keep "$KEEP_A" \
    --freq \
    --out "$PLINK_OUT_A"

  plink2 --bfile "$BFILE" \
    --extract "$SNPLIST" \
    --keep "$KEEP_B" \
    --freq \
    --out "$PLINK_OUT_B"
fi

# 3) Summarize into per-SNP mean(beta*dosage) per group + delta
Rscript --vanilla - "$PLINK_OUT_A.afreq" "$PLINK_OUT_B.afreq" "$WEIGHTS" "$OUT" "$gA" "$gB" "$WORKDIR/PT.txt" <<'RSCRIPT'
args <- commandArgs(trailingOnly=TRUE)
afA_path <- args[1]
afB_path <- args[2]
weights_path <- args[3]
out_prefix <- args[4]
gA <- args[5]
gB <- args[6]
pt_path <- args[7]

read_flex <- function(path) {
  if (requireNamespace("data.table", quietly=TRUE)) data.table::fread(path, data.table=FALSE)
  else read.table(path, header=TRUE, stringsAsFactors=FALSE, comment.char="", fill=TRUE)
}

afA <- read_flex(afA_path)
afB <- read_flex(afB_path)
w   <- read_flex(weights_path)
pt  <- suppressWarnings(as.numeric(readLines(pt_path, warn=FALSE)[1]))

# Expect PLINK2 .afreq columns like:
# #CHROM ID REF ALT PROVISIONAL_REF? ALT_FREQS OBS_CT
# We'll only rely on: ID, REF, ALT, ALT_FREQS
id_col <- if ("ID" %in% names(afA)) "ID" else stop("No ID column in afreq A.")
if (!all(c("REF","ALT","ALT_FREQS") %in% names(afA))) {
  stop("afreq A missing REF/ALT/ALT_FREQS. Columns: ", paste(names(afA), collapse=", "))
}
if (!all(c("REF","ALT","ALT_FREQS") %in% names(afB))) {
  stop("afreq B missing REF/ALT/ALT_FREQS. Columns: ", paste(names(afB), collapse=", "))
}

# weights must have SNP, A1(effect allele), BETA
if (!all(c("SNP","A1","BETA") %in% names(w))) {
  stop("weights.txt must contain columns: SNP, A1, BETA. Columns: ", paste(names(w), collapse=", "))
}

# Map effect allele frequency from REF/ALT/ALT_FREQS
# ALT may be comma-separated (multiallelic); ALT_FREQS similarly.
effect_freq_from_row <- function(ref, alt, alt_freqs, eff) {
  if (is.na(ref) || is.na(alt) || is.na(alt_freqs) || is.na(eff)) return(NA_real_)
  ref <- as.character(ref); alt <- as.character(alt); eff <- as.character(eff)
  alt_alleles <- strsplit(alt, ",", fixed=TRUE)[[1]]
  freqs <- as.numeric(strsplit(as.character(alt_freqs), ",", fixed=TRUE)[[1]])
  if (length(freqs) != length(alt_alleles)) return(NA_real_)

  if (eff == ref) {
    # ref freq = 1 - sum(alt freqs)
    return(1.0 - sum(freqs))
  }
  hit <- which(alt_alleles == eff)
  if (length(hit) == 1) return(freqs[hit])
  return(NA_real_)
}

map_effect_freq <- function(af, w) {
  m <- merge(w, af, by.x="SNP", by.y="ID", all.x=TRUE)
  m$EFFECT_FREQ <- mapply(effect_freq_from_row, m$REF, m$ALT, m$ALT_FREQS, m$A1)
  return(m[, c("SNP","A1","BETA","EFFECT_FREQ")])
}

mA <- map_effect_freq(afA, w)
mB <- map_effect_freq(afB, w)

colnames(mA)[colnames(mA)=="EFFECT_FREQ"] <- paste0("FREQ_", gA)
colnames(mB)[colnames(mB)=="EFFECT_FREQ"] <- paste0("FREQ_", gB)

m <- merge(mA, mB[, c("SNP", paste0("FREQ_", gB))], by="SNP", all.x=TRUE)

# Compute mean dosage and contribution
m[[paste0("MEAN_DOSAGE_", gA)]] <- 2 * m[[paste0("FREQ_", gA)]]
m[[paste0("MEAN_DOSAGE_", gB)]] <- 2 * m[[paste0("FREQ_", gB)]]

m[[paste0("MEAN_CONTRIB_", gA)]] <- m$BETA * m[[paste0("MEAN_DOSAGE_", gA)]]
m[[paste0("MEAN_CONTRIB_", gB)]] <- m$BETA * m[[paste0("MEAN_DOSAGE_", gB)]]

m$DELTA_CONTRIB <- m[[paste0("MEAN_CONTRIB_", gA)]] - m[[paste0("MEAN_CONTRIB_", gB)]]
m$DELTA_AF      <- m[[paste0("FREQ_", gA)]] - m[[paste0("FREQ_", gB)]]
m$ABS_DELTA     <- abs(m$DELTA_CONTRIB)

m <- m[order(-m$ABS_DELTA), ]

out_csv <- paste0(out_prefix, ".per_snp_contrib.csv")
top_csv <- paste0(out_prefix, ".top50_drivers.csv")
write.csv(m, out_csv, row.names=FALSE)
write.csv(head(m, 50), top_csv, row.names=FALSE)

cat("PT used:", pt, "\n")
cat("Groups:", gA, "vs", gB, "\n")
cat("Wrote:", out_csv, "\n")
cat("Wrote:", top_csv, "\n")

# Optional: warn if many NAs (allele mismatches)
naA <- mean(is.na(m[[paste0("FREQ_", gA)]]))
naB <- mean(is.na(m[[paste0("FREQ_", gB)]]))
if (naA > 0.05 || naB > 0.05) {
  cat("WARNING: High NA rate in mapped frequencies. Check SNP IDs and effect allele coding.\n")
  cat("NA rate ", gA, ": ", naA, "\n", sep="")
  cat("NA rate ", gB, ": ", naB, "\n", sep="")
}
RSCRIPT


echo "Done. Main output: ${OUT}.per_snp_contrib.csv"
echo "Top drivers: ${OUT}.top50_drivers.csv"
