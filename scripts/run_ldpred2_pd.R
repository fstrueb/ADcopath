#!/usr/bin/env Rscript
# run_ldpred2_pd.R
# Genome-wide shrinkage PRS (LDpred2-auto) for the PD GWAS, as a modern-method
# benchmark for the PRSice-2 clumping+threshold PD PRS (reviewer request).
#
# Pipeline:
#   1. Download the canonical bigsnpr European LD reference (HapMap3 + LD blocks).
#   2. Build df_beta from the PD GWAS base (BETA/P -> beta_se; n_eff from Nalls 2019).
#      Match GWAS <-> reference by rsID (base is hg38, reference is hg19; rsID match
#      sidesteps the build difference).
#   3. Assemble the genome-wide sparse correlation (SFBM) from the precomputed blocks.
#   4. LDpred2-auto -> genome-wide posterior effect sizes (averaged over good chains).
#   5. Score both target cohorts with plink2 --score (variant IDs set to chr:pos hg38,
#      so weights anchored on the base's hg38 positions match the targets).
#
# Output: output/ldpred2_pd_prs.tsv     (cohort, IID, ldpred2_prs)  -> read by Figure_1.Rmd
#         output/ldpred2_pd_betas.tsv   (per-SNP posterior weights, diagnostics)
#
# Run from project root:
#   cd /home/fstruebi/projects/ADcopath_final && Rscript scripts/run_ldpred2_pd.R

suppressPackageStartupMessages({
    library(bigsnpr)
    library(data.table)
    library(tidyverse)
})

options(timeout = 7200)                # LD reference zip is ~7.7 GB; default 300s is too short

adcopath  <- '/home/fstruebi/projects/ADcopath'
ref_dir   <- 'temp/ldpred2_ref'        # large LD reference (keep out of git)
work_dir  <- 'temp/ldpred2_work'       # SFBM + plink2 scratch
plink2    <- '/opt/plink2_binary/plink2'
NCORES    <- bigstatsr::nb_cores()
dir.create('output', showWarnings = FALSE)
dir.create(ref_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

# Nalls et al. 2019 (excluding 23andMe) effective sample size:
# n_eff = 4 / (1/Ncase + 1/Ncontrol), Ncase = 33674, Ncontrol = 449056.
N_EFF <- 4 / (1 / 33674 + 1 / 449056)

# ---------------------------------------------------------------------------
# 1. Download canonical European LD reference (HapMap3 + LD blocks)
# ---------------------------------------------------------------------------
message('Downloading LD reference (HapMap3 + blocks) if needed...')
map_rds <- runonce::download_file(
    'https://ndownloader.figshare.com/files/36360900',
    dir = ref_dir, fname = 'map_hm3_with_blocks.rds')
zip_f <- runonce::download_file(
    'https://ndownloader.figshare.com/files/36363087',
    dir = ref_dir, fname = 'ldref_with_blocks.zip')
if (!file.exists(file.path(ref_dir, 'ldref', 'LD_with_blocks_chr1.rds')) &&
    length(list.files(ref_dir, pattern = 'LD_with_blocks_chr1.rds', recursive = TRUE)) == 0)
    unzip(zip_f, exdir = ref_dir)
ld_chr1 <- list.files(ref_dir, pattern = 'LD_with_blocks_chr1.rds',
                      recursive = TRUE, full.names = TRUE)[1]
ld_path <- dirname(ld_chr1)
message('LD matrices in: ', ld_path)

map_ldref <- readRDS(map_rds)
map_ldref$ref_id <- seq_len(nrow(map_ldref))   # global row index into the reference

# ---------------------------------------------------------------------------
# 2. GWAS sumstats -> df_beta, matched to the reference by rsID
# ---------------------------------------------------------------------------
message('Reading PD GWAS base and matching to reference by rsID...')
base <- fread(file.path(adcopath, 'PRSice/input/PD_hg38_forPRSice.tsv'),
              select = c('SNP','A1','A2','BETA','P','CHR','BP')) %>%
    as_tibble() %>%
    filter(is.finite(BETA), is.finite(P), P > 0, P <= 1) %>%
    mutate(beta_se = abs(BETA) / qnorm(P / 2, lower.tail = FALSE)) %>%
    filter(is.finite(beta_se), beta_se > 0)

# match by rsID, then align effect allele to the reference a1 (flip beta if needed),
# drop allele mismatches and strand-ambiguous (A/T, C/G) SNPs.
flip <- c(A = 'T', T = 'A', C = 'G', G = 'C')
df_beta <- base %>%
    inner_join(map_ldref %>% dplyr::select(ref_id, r_chr = chr, r_a0 = a0, r_a1 = a1,
                                           rsid, ld),
               by = c('SNP' = 'rsid')) %>%
    mutate(
        A1 = toupper(A1), A2 = toupper(A2),
        ambiguous = (A2 == flip[A1]),
        orient = dplyr::case_when(
            A1 == r_a1 & A2 == r_a0 ~  1L,   # already aligned to reference a1
            A1 == r_a0 & A2 == r_a1 ~ -1L,   # swapped -> flip beta
            TRUE                    ~ NA_integer_)
    ) %>%
    filter(!ambiguous, !is.na(orient)) %>%
    mutate(
        beta    = BETA * orient,             # effect now refers to reference a1 (= r_a1)
        n_eff   = N_EFF,
        eff_allele = r_a1,                   # nucleotide the posterior weight refers to
        chr     = r_chr
    ) %>%
    distinct(ref_id, .keep_all = TRUE) %>%   # one GWAS SNP per reference position
    arrange(ref_id)                          # reference (chr, pos) order == corr order

message('Matched ', nrow(df_beta), ' SNPs to the LD reference.')

# ---------------------------------------------------------------------------
# 3. Genome-wide sparse correlation (SFBM) from precomputed blocks
# ---------------------------------------------------------------------------
message('Assembling genome-wide SFBM correlation...')
tmp_sfbm <- file.path(work_dir, paste0('corr_', Sys.getpid()))
corr <- NULL
for (chr in 1:22) {
    ind.chr  <- which(df_beta$chr == chr)
    if (length(ind.chr) == 0) next
    ind.glob <- df_beta$ref_id[ind.chr]                       # rows in map_ldref
    ind.loc  <- match(ind.glob, which(map_ldref$chr == chr))  # rows in this chr's matrix
    corr_chr <- readRDS(file.path(ld_path, paste0('LD_with_blocks_chr', chr, '.rds')))[ind.loc, ind.loc]
    if (is.null(corr)) {
        corr <- as_SFBM(corr_chr, tmp_sfbm, compact = TRUE)
    } else {
        corr$add_columns(corr_chr, nrow(corr))
    }
}
stopifnot(nrow(corr) == nrow(df_beta))

# ---------------------------------------------------------------------------
# 4. LDpred2-auto
# ---------------------------------------------------------------------------
message('LDSC h2 init...')
ldsc   <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                 chi2 = (beta / beta_se)^2,
                                 sample_size = n_eff, blocks = NULL))
h2_est <- ldsc[['h2']]
message(sprintf('LDSC h2 = %.4f', h2_est))

message('Running LDpred2-auto (', NCORES, ' cores)...')
set.seed(1)
multi_auto <- snp_ldpred2_auto(
    corr, df_beta, h2_init = h2_est,
    vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
    allow_jump_sign = FALSE, shrink_corr = 0.95, ncores = NCORES)

# keep well-behaved chains (stable posterior scale), average their effect sizes
rg    <- sapply(multi_auto, function(a)
    if (is.null(a$corr_est) || anyNA(a$corr_est)) NA_real_ else diff(range(a$corr_est)))
keep  <- which(rg > (0.95 * quantile(rg, 0.95, na.rm = TRUE)))
message('Kept ', length(keep), '/', length(multi_auto), ' LDpred2-auto chains.')
beta_auto <- rowMeans(sapply(multi_auto[keep], function(a) a$beta_est))
p_est  <- mean(sapply(multi_auto[keep], function(a) a$p_est))
message(sprintf('LDpred2-auto: p_est = %.2e, h2 = %.4f, %d SNPs',
                p_est, h2_est, length(beta_auto)))

df_beta$beta_auto <- beta_auto
betas_out <- df_beta %>%
    dplyr::transmute(rsid = SNP, chr = CHR, pos_hg38 = BP,
                     effect_allele = eff_allele, beta_auto)
write_tsv(betas_out, 'output/ldpred2_pd_betas.tsv')

# ---------------------------------------------------------------------------
# 5. Score both cohorts with plink2 (variant IDs -> chr:pos hg38)
# ---------------------------------------------------------------------------
weights_f <- file.path(work_dir, 'ldpred2_weights.txt')
betas_out %>%
    dplyr::transmute(ID = paste0(chr, ':', pos_hg38), A1 = effect_allele, BETA = beta_auto) %>%
    write_tsv(weights_f)

score_cohort <- function(pfile, out_prefix) {
    args <- c('--pfile', pfile,
              '--set-all-var-ids', '@:#',
              '--rm-dup', 'force-first',
              '--score', weights_f, '1', '2', '3', 'header', 'cols=+scoresums',
              '--out', out_prefix, '--silent')
    status <- system2(plink2, args)
    message('plink2 --score exit ', status, ' for ', basename(out_prefix))
    fread(paste0(out_prefix, '.sscore'))
}

message('Scoring SyNergy...')
syn_sc <- score_cohort('/neptune/ADWGS_hg38/plink2/synergy_GATKjoint_hg38',
                       file.path(work_dir, 'synergy_PD_ldpred2'))
message('Scoring AMP-AD...')
amp_sc <- score_cohort('/neptune/AMPAD_data/genomicVariants/plink2/AMPAD_filtered_hg38',
                       file.path(work_dir, 'ampad_PD_ldpred2'))

iid_col   <- function(d) names(d)[grepl('IID', names(d))][1]
score_col <- function(d) {
    cn <- names(d)
    cn[match(TRUE, cn %in% c('SCORE1_AVG','SCORE1_SUM'))]
}
tidy_sc <- function(d, cohort) {
    tibble(cohort = cohort,
           IID = d[[iid_col(d)]],
           ldpred2_prs = d[[score_col(d)]])
}
prs_out <- bind_rows(tidy_sc(syn_sc, 'SyNergy'), tidy_sc(amp_sc, 'AMP-AD'))
write_tsv(prs_out, 'output/ldpred2_pd_prs.tsv')

message('Done. Wrote output/ldpred2_pd_prs.tsv (', nrow(prs_out), ' individuals) ',
        'and output/ldpred2_pd_betas.tsv (', nrow(betas_out), ' SNPs).')
