#!/usr/bin/env Rscript
# run_loo_loci_pd.R
# Leave-one-locus-out (LOO) sensitivity for the PD PRS, generalising the SNCA
# leave-one-out to every LD locus carrying PRS weight (reviewer supplement).
#
# A PRS is additive (PRS = sum over loci), so LOO needs no per-locus PRSice re-run:
# we score each deCODE LD block's partial sum once (one multi-column plink2 --score
# per cohort) and subtract from the full score. The SNCA block's LOO effect should
# reproduce the independent --x-range noSNCA result (built-in validation).
#
# Loci = the same deCODE EUR LD blocks used for SNCA and Fig. 1d.
# PRS SNP set = PRSice's own clumped+thresholded SNPs (the .snp file, P < threshold).
#
# Output: output/loo_pd_effects.tsv  (cohort, locus, block, n_snps, estimate/CI/p)
# Run from project root:
#   cd /home/fstruebi/projects/ADcopath_final && Rscript scripts/run_loo_loci_pd.R

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(GenomicRanges)
    library(broom)
})
source('scripts/helper_functions.R')

adcopath <- '/home/fstruebi/projects/ADcopath'
plink2   <- '/opt/plink2_binary/plink2'
work     <- 'temp/loo_work'
dir.create(work,     showWarnings = FALSE, recursive = TRUE)
dir.create('output', showWarnings = FALSE)

# --- metadata / covariates (PRS IID conventions match the .all_score joins) ---
synergy_metad <- read_tsv(file.path(adcopath, 'resources_paper/synergy_cohort_metad.tsv'),
                          show_col_types = FALSE) %>% mutate(IID = paste0(IID, '_', IID))
ampad_metad   <- read_tsv(file.path(adcopath, 'resources_paper/AMPAD_cohort_metad.tsv'),
                          show_col_types = FALSE)
covars <- c('age', 'sex', paste0('PC', 1:10))

# --- PD GWAS base (effect allele + beta) and deCODE LD blocks ---
base <- fread(file.path(adcopath, 'PRSice/input/PD_hg38_forPRSice.tsv'),
              select = c('SNP','A1','A2','BETA','P','CHR','BP'))
pd_gwas_coords <- as.data.frame(base[, .(SNP, CHR, BP, P, BETA, A1)])

blocks <- fread('/earth/public_data/ldetect/LDblocks_GRCh38/data/deCODE_EUR_LD_blocks.bed',
                header = TRUE)
setnames(blocks, c('chr','start','end'))
blocks[, `:=`(chr = sub('^chr', '', chr), block_id = .I)]
blocks_gr <- GRanges(blocks$chr, IRanges(blocks$start + 1L, blocks$end), block_id = blocks$block_id)

std_effect <- function(x, covars) {
    x$PRS_z <- as.numeric(scale(x$PRS))
    broom::tidy(lm(reformulate(c('group', covars), 'PRS_z'), data = x, na.action = na.exclude),
                conf.int = TRUE) %>%
        dplyr::filter(term == 'groupAD+LBP') %>%
        dplyr::select(estimate, conf.low, conf.high, p.value)
}

run_cohort <- function(snp_file, threshold, target, metad, cohort) {
    message('=== ', cohort, ' (P < ', threshold, ') ===')
    snp <- fread(file.path(adcopath, snp_file))[P < threshold, .(SNP, CHR, BP, P)]
    prs <- merge(snp, base[, .(SNP, A1, BETA)], by = 'SNP')

    # assign each PRS SNP to a deCODE block; SNPs in gaps become singleton loci
    g  <- GRanges(as.character(prs$CHR), IRanges(prs$BP, prs$BP))
    ov <- findOverlaps(g, blocks_gr)
    prs[, block := NA_integer_]
    prs[queryHits(ov), block := blocks_gr$block_id[subjectHits(ov)]]
    na_idx <- which(is.na(prs$block))                    # SNPs in gaps -> singleton loci
    if (length(na_idx)) prs$block[na_idx] <- -seq_along(na_idx)
    prs[, grp := match(block, unique(block))]
    prs[, col := sprintf('S%03d', grp)]
    prs[, ID  := paste0(CHR, ':', BP)]

    # per-locus lead SNP + deCODE LD-block coordinates (matching the Fig. 1d labels);
    # SNPs in inter-block gaps keep their own position as the locus span.
    blk_coords <- blocks[, .(blk = block_id, b_start = start, b_end = end)]
    leads <- prs[order(P), .(lead = SNP[1], n_snps = .N, chr = CHR[1],
                             snp_lo = min(BP), snp_hi = max(BP), blk = block[1]),
                 by = .(grp, col)]
    leads <- merge(leads, blk_coords, by = 'blk', all.x = TRUE)
    leads[, `:=`(start = data.table::fifelse(is.na(b_start), snp_lo, b_start),
                 end   = data.table::fifelse(is.na(b_end),   snp_hi, b_end))]
    ann <- unique(as.data.table(annotate_rsids_band_nearest_gene(
        leads$lead, build = 'hg38', gwas_coords = pd_gwas_coords)), by = 'rsid')
    leads <- merge(leads, ann[, .(lead = rsid, gene = closest_gene, cytoband)],
                   by = 'lead', all.x = TRUE)

    # weight matrix (one column per locus) -> plink2 multi-column --score (sums)
    wide  <- dcast(prs, ID + A1 ~ col, value.var = 'BETA', fill = 0)
    setcolorder(wide, c('ID', 'A1', sort(setdiff(names(wide), c('ID','A1')))))
    wfile <- file.path(work, paste0(cohort, '_weights.txt'))
    fwrite(wide, wfile, sep = '\t')
    out   <- file.path(work, paste0(cohort, '_loo'))
    last  <- ncol(wide)
    system2(plink2, c('--pfile', target, '--set-all-var-ids', '@:#',
                      '--rm-dup', 'force-first',
                      '--score', wfile, '1', '2', 'header', 'cols=+scoresums',
                      '--score-col-nums', paste0('3-', last),
                      '--out', out, '--silent'))

    sc        <- fread(paste0(out, '.sscore'))
    iid_col   <- names(sc)[grepl('IID', names(sc))][1]
    sum_cols  <- grep('^SCORE[0-9]+_SUM$', names(sc), value = TRUE)   # per-locus sums, column order
    stopifnot(length(sum_cols) == nrow(leads))
    M         <- as.matrix(sc[, ..sum_cols])
    colnames(M) <- sort(unique(prs$col))                      # S001.. == grp order
    full_vec  <- rowSums(M)

    meta <- tibble(IID = sc[[iid_col]]) %>% left_join(metad, by = 'IID')

    # full reconstructed PRS effect (sanity check vs PRSice) + per-locus LOO
    eff_full <- std_effect(meta %>% mutate(PRS = full_vec), covars) %>%
        mutate(locus = '(full PRS)', n_snps = nrow(prs), block = NA_character_, .before = 1)

    eff_loo <- purrr::map_dfr(seq_len(nrow(leads)), function(i) {
        cj  <- leads$col[i]
        std_effect(meta %>% mutate(PRS = full_vec - M[, cj]), covars) %>%
            mutate(locus  = coalesce(leads$gene[i], leads$lead[i]),
                   n_snps = leads$n_snps[i],
                   block  = paste0('chr', leads$chr[i], ':', leads$start[i], '-', leads$end[i]),
                   .before = 1)
    })

    res <- bind_rows(eff_full, eff_loo) %>% mutate(cohort = cohort, .before = 1)
    snca <- res %>% dplyr::filter(locus == 'SNCA')
    if (nrow(snca))
        message(sprintf('  SNCA LOO beta = %.3f (compare to --x-range noSNCA); full = %.3f',
                        snca$estimate[1], eff_full$estimate))
    res
}

loo <- bind_rows(
    run_cohort('PRSice/fixed/synergy_GATKjoint_PD_fixed.snp', 5e-8,
               '/neptune/ADWGS_hg38/plink2/synergy_GATKjoint_hg38', synergy_metad, 'SyNergy'),
    run_cohort('PRSice/fixed/AMPAD_PD_fixed.snp', 1e-4,
               '/neptune/AMPAD_data/genomicVariants/plink2/AMPAD_filtered_hg38', ampad_metad, 'AMP-AD')
)
write_tsv(loo, 'output/loo_pd_effects.tsv')
message('Wrote output/loo_pd_effects.tsv (', nrow(loo), ' rows).')
