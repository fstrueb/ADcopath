#!/usr/bin/env Rscript
# run_ld_stouffer.R
# Build mydrivers with LD-corrected Stouffer Z, parallelised over blocks.
#
# Run order:
#   1. Rscript scripts/run_ld_stouffer.R          # writes SNP lists + uncorrected mydrivers
#   2. (terminal) bash: loop plink2 --r-phased square over temp/ld_blocks/*.snps
#   3. Rscript scripts/run_ld_stouffer.R          # re-runs with LD correction applied
#
# Output: output/mydrivers_ld.rds
#
# Run from project root: cd /home/fstruebi/projects/ADcopath_final && Rscript scripts/run_ld_stouffer.R

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(broom)
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
    library(Matrix)
    library(parallel)
})

source('scripts/helper_functions.R')

adcopath <- '/home/fstruebi/projects/ADcopath'
ld_dir   <- 'temp/ld_blocks'
out_file <- 'output/mydrivers_ld.rds'
dir.create('output',   showWarnings = FALSE)
dir.create(ld_dir,     showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load metadata and group sizes
# ---------------------------------------------------------------------------
message('Loading metadata...')
synergy_metad <- read_tsv(file.path(adcopath, 'resources_paper/synergy_cohort_metad.tsv'),
                           show_col_types = FALSE)
ampad_metad   <- read_tsv(file.path(adcopath, 'resources_paper/AMPAD_cohort_metad.tsv'),
                           show_col_types = FALSE)

n_syn_AD    <- sum(synergy_metad$group == 'AD')
n_syn_ADLBP <- sum(synergy_metad$group == 'AD+LBP')
n_amp_AD    <- sum(ampad_metad$group == 'AD')
n_amp_ADLBP <- sum(ampad_metad$group == 'AD+LBP')

# ---------------------------------------------------------------------------
# 2. Load and annotate driver SNPs
# ---------------------------------------------------------------------------
message('Loading PD GWAS coords...')
pd_gwas_coords <- data.table::fread(
    file.path(adcopath, 'PRSice/input/PD_hg38_forPRSice.tsv'),
    select = c('SNP', 'CHR', 'BP'))

message('Annotating SyNergy drivers...')
drivers_synergy_joint <- read_csv(
    file.path(adcopath, 'PRSice/drivers/synergy_joint_PD.per_snp_contrib.csv'),
    show_col_types = FALSE) %>%
    left_join(annotate_rsids_band_nearest_gene(.$SNP, build = 'hg38',
                                               gwas_coords = pd_gwas_coords),
              by = c('SNP' = 'rsid'))

message('Annotating AMP-AD drivers...')
drivers_ampad <- read_csv(
    file.path(adcopath, 'PRSice/drivers/ampad_PD.per_snp_contrib.csv'),
    show_col_types = FALSE) %>%
    left_join(annotate_rsids_band_nearest_gene(.$SNP, build = 'hg38',
                                               gwas_coords = pd_gwas_coords),
              by = c('SNP' = 'rsid'))

# ---------------------------------------------------------------------------
# 3. Block-level aggregation (independence-assumed Stouffer Z)
# ---------------------------------------------------------------------------
message('Computing block scores...')
drivers_synergy_scored <- calc_contribution_scores(drivers_synergy_joint,
                                                   n_A = n_syn_AD,
                                                   n_B = n_syn_ADLBP)
drivers_ampad_scored   <- calc_contribution_scores(drivers_ampad,
                                                   n_A = n_amp_AD,
                                                   n_B = n_amp_ADLBP)

mydrivers <- drivers_synergy_scored %>%
    left_join(drivers_ampad_scored, by = 'block_id', suffix = c('.SYN', '.AMP')) %>%
    mutate(
        block_id = paste0('chr', chr.SYN, ':', start.SYN, '-', end.SYN),
        mean_z   = rowMeans(cbind(stouffer_z.SYN, stouffer_z.AMP), na.rm = TRUE),
        meta_z   = mean_z * sqrt(2),
        meta_p   = 2 * pnorm(-abs(meta_z)),
        meta_q   = p.adjust(meta_p, method = 'BH')
    ) %>%
    filter(block_id != 'chr17:45383525-50162864')

# ---------------------------------------------------------------------------
# 4. Write per-block SNP lists for plink2 LD computation
# ---------------------------------------------------------------------------
message('Writing SNP lists to ', ld_dir, '...')
pd_gwas_full <- data.table::fread(
    file.path(adcopath, 'PRSice/input/PD_hg38_forPRSice.tsv')
) %>%
    as_tibble() %>%
    filter(is.finite(BETA), is.finite(P), P > 0, P <= 1) %>%
    mutate(se  = abs(BETA) / qnorm(P / 2, lower.tail = FALSE),
           CHR = as.character(CHR)) %>%
    filter(is.finite(se), se > 0)

driver_blocks <- mydrivers %>%
    filter(!is.na(stouffer_z.SYN) | !is.na(stouffer_z.AMP)) %>%
    pull(block_id)

for (blk in driver_blocks) {
    coords  <- strsplit(sub('^chr', '', blk), '[:-]')[[1]]
    chr_b   <- coords[1]; start_b <- as.integer(coords[2]); end_b <- as.integer(coords[3])
    snps    <- pd_gwas_full %>%
        filter(CHR == chr_b, BP >= start_b, BP <= end_b) %>%
        distinct(SNP) %>% pull(SNP)
    slug    <- gsub(':', '_', gsub('^chr', '', blk))
    writeLines(snps, file.path(ld_dir, paste0(slug, '.snps')))
}
message('SNP lists written. Run plink2 loop before proceeding:')
message('  PFILE=/neptune/ADWGS_hg38/plink2/synergy_GATKjoint_hg38')
message('  for f in temp/ld_blocks/*.snps; do')
message('    plink2 --pfile $PFILE --extract $f --r-phased square --out "${f%.snps}" --silent')
message('  done')

# ---------------------------------------------------------------------------
# 5. LD-corrected Stouffer Z (runs only if LD matrices are present)
# ---------------------------------------------------------------------------
ld_files <- list.files(ld_dir, pattern = '\\.phased\\.vcor1$', full.names = FALSE)

if (length(ld_files) == 0L) {
    message('No LD matrices found — saving mydrivers with independence-assumed Stouffer Z.')
    message('Re-run this script after computing LD matrices with plink2.')
    saveRDS(mydrivers, out_file)
    message('Saved: ', out_file)
    quit(save = 'no', status = 0L)
}

message(length(ld_files), ' LD matrix files found. Applying LD correction...')

# Prepare SNP-level data with per-SNP Z-scores and block string IDs
prepare_drv <- function(driver_annotated, n_A, n_B) {
    as_tibble(driver_annotated) %>%
        filter(!is.na(block_id), !is.na(chr), !is.na(start), !is.na(end),
               !is.na(FREQ_0), !is.na(FREQ_1), !is.na(BETA), !is.na(SNP)) %>%
        mutate(
            freq_pool = (FREQ_0 * 2 * n_A + FREQ_1 * 2 * n_B) / (2 * (n_A + n_B)),
            se_i      = sqrt(freq_pool * (1 - freq_pool) * (1/n_A + 1/n_B)),
            z_i       = ifelse(se_i > 0, (FREQ_0 - FREQ_1) / se_i, 0),
            w         = abs(BETA),
            block_str = paste0('chr', chr, ':', start, '-', end)
        )
}

# Per-block LD-corrected Z (called via mclapply)
# Key design: subset to driver SNPs BEFORE nearPD — avoids O(n^3) on 16 k-SNP blocks.
# Uses fread(select=) to read only the driver-SNP columns from disk (~1 MB vs 2 GB).
compute_block_z_ld <- function(blk, drv, ld_dir) {
    d      <- drv[drv$block_str == blk, ]
    sz_ind <- sum(d$w * d$z_i) / sqrt(sum(d$w^2))

    slug   <- gsub(':', '_', gsub('^chr', '', blk))
    mat_f  <- file.path(ld_dir, paste0(slug, '.phased.vcor1'))
    vars_f <- file.path(ld_dir, paste0(slug, '.phased.vcor1.vars'))

    if (!file.exists(mat_f) || !file.exists(vars_f))
        return(data.frame(block_id = blk, stouffer_z = sz_ind, ld_corrected = FALSE))

    ld_ids   <- readLines(vars_f)
    keep     <- intersect(d$SNP, ld_ids)

    if (length(keep) < 2L)
        return(data.frame(block_id = blk, stouffer_z = sz_ind, ld_corrected = FALSE))

    keep_idx <- which(ld_ids %in% keep)

    # Read only the driver-SNP columns — O(n * m) not O(n^2) for the read,
    # then take only driver rows. n = all block SNPs, m = driver SNPs (<<n).
    sub_cols <- as.matrix(data.table::fread(mat_f, header = FALSE,
                                             select = keep_idx))
    R_k <- sub_cols[keep_idx, ]
    rm(sub_cols)

    rownames(R_k) <- colnames(R_k) <- ld_ids[keep_idx]
    R_k[!is.finite(R_k)] <- 0
    diag(R_k) <- 1

    # nearPD only on the small driver submatrix — trivially fast for m << n
    R_k <- as.matrix(Matrix::nearPD(R_k, keepDiag = TRUE)$mat)

    d_k      <- d[d$SNP %in% keep, ]
    d_k      <- d_k[order(match(d_k$SNP, rownames(R_k))), ]
    w_k      <- d_k$w
    denom_ld <- sqrt(as.numeric(t(w_k) %*% R_k[d_k$SNP, d_k$SNP] %*% w_k))
    sz_ld    <- if (denom_ld > 0) sum(w_k * d_k$z_i) / denom_ld else NA_real_

    data.frame(block_id = blk, stouffer_z = sz_ld, ld_corrected = TRUE)
}

drv_syn <- prepare_drv(annotate_ld_locus(drivers_synergy_joint), n_syn_AD,    n_syn_ADLBP)
drv_amp <- prepare_drv(annotate_ld_locus(drivers_ampad),         n_amp_AD,    n_amp_ADLBP)

n_cores <- max(1L, detectCores() - 1L)
message('Parallelising over ', n_cores, ' cores')

message('SyNergy: ', length(unique(drv_syn$block_str)), ' blocks')
ld_syn <- bind_rows(mclapply(unique(drv_syn$block_str), compute_block_z_ld,
                              drv = drv_syn, ld_dir = ld_dir, mc.cores = n_cores))

message('AMP-AD: ', length(unique(drv_amp$block_str)), ' blocks')
ld_amp <- bind_rows(mclapply(unique(drv_amp$block_str), compute_block_z_ld,
                              drv = drv_amp, ld_dir = ld_dir, mc.cores = n_cores))

n_corrected_syn <- sum(ld_syn$ld_corrected, na.rm = TRUE)
n_corrected_amp <- sum(ld_amp$ld_corrected, na.rm = TRUE)
message('LD-corrected blocks — SyNergy: ', n_corrected_syn, '  AMP-AD: ', n_corrected_amp)

# Update mydrivers in-place
mydrivers <- mydrivers %>%
    left_join(ld_syn %>% dplyr::select(block_id, sz_syn = stouffer_z), by = 'block_id') %>%
    left_join(ld_amp %>% dplyr::select(block_id, sz_amp = stouffer_z), by = 'block_id') %>%
    mutate(
        stouffer_z.SYN = coalesce(sz_syn, stouffer_z.SYN),
        stouffer_z.AMP = coalesce(sz_amp, stouffer_z.AMP),
        stouffer_p.SYN = 2 * pnorm(-abs(stouffer_z.SYN)),
        stouffer_p.AMP = 2 * pnorm(-abs(stouffer_z.AMP)),
        stouffer_q.SYN = p.adjust(stouffer_p.SYN, method = 'BH'),
        stouffer_q.AMP = p.adjust(stouffer_p.AMP, method = 'BH'),
        mean_z = rowMeans(cbind(stouffer_z.SYN, stouffer_z.AMP), na.rm = TRUE),
        meta_z = mean_z * sqrt(2),
        meta_p = 2 * pnorm(-abs(meta_z)),
        meta_q = p.adjust(meta_p, method = 'BH')
    ) %>%
    dplyr::select(-sz_syn, -sz_amp)

saveRDS(mydrivers, out_file)
message('Saved: ', out_file)
