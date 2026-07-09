library(tidyverse)
library(data.table)
library(broom)
library(biomaRt)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)
library(rtracklayer)
library(CMplot)
library(GenomeInfoDb)

AD_color <- ggsci::pal_igv('alternating')(2)[1]
ASYN_color <- ggsci::pal_igv('alternating')(2)[2]

p_to_label <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    # paste0("p=", formatC(p, format = "f", digits = 2))
    'ns'
}

plot_violin_with_lm_p <- function(df,
                                  y = "PRS_adj",
                                  group = "group",
                                  covars = c("age","sex", paste0("PC", 1:10)),
                                  group_term = NULL,
                                  title) {
    
    # pick the coefficient name for group (handles factor coding)
    # default: first non-intercept term that starts with "group"
    if (is.null(group_term)) {
        # ensure group is factor with desired order
        df[[group]] <- factor(df[[group]], levels = c("AD", "AD+LBP"))
    }
    
    fml <- as.formula(paste(y, "~", group, "+", paste(covars, collapse = " + ")))
    fit <- lm(fml, data = df)
    
    tt <- tidy(fit)
    
    if (is.null(group_term)) {
        group_term <- tt$term[grepl(paste0("^", group), tt$term)][1]
    }
    
    p <- tt$p.value[match(group_term, tt$term)]
    lab <- p_to_label(p)
    
    # fixed annotation height (keeps it stable & compact)
    yvals <- df[[y]]
    y_top <- max(yvals, na.rm = TRUE)
    rng   <- diff(range(yvals, na.rm = TRUE))
    if (!is.finite(rng) || rng == 0) rng <- 1
    y_bracket <- y_top + 0.08 * rng
    y_text    <- y_top + 0.11 * rng
    
    # x positions for the 2 groups
    x1 <- 1
    x2 <- 2
    
    ggplot(df, aes(x = .data[[group]], y = .data[[y]], fill = .data[[group]])) +
        geom_violin(width = 0.85, trim = TRUE, linewidth = 0.25) +
        geom_boxplot(
            width = 0.25,
            outlier.shape = NA,
            staplewidth = 0,
            coef = 0
        ) +
        geom_point(position = position_jitter(width = 0.03, height = 0),
                   alpha = 0.3, size = 0.8) +
        scale_fill_manual(values = c(AD_color, ASYN_color), guide = FALSE) +
        annotate("segment", x = x1, xend = x1, y = y_top, yend = y_bracket, linewidth = 0.25) +
        annotate("segment", x = x1, xend = x2, y = y_bracket, yend = y_bracket, linewidth = 0.25) +
        annotate("segment", x = x2, xend = x2, y = y_top, yend = y_bracket, linewidth = 0.25) +
        annotate("text", x = (x1 + x2)/2, y = y_text, label = lab, vjust = 0, size = 3) +
        labs(title = basename(title)) +
        coord_cartesian(clip = "off") +
        theme_light(base_size = 14) +
        theme(
            axis.title = element_blank(),
            plot.margin = margin(2, 2, 10, 2)  # extra top margin so label isn't clipped
        )
}

annotate_rsids_band_nearest_gene <- function(rsids, build = c("hg38","hg19"),
                                             gwas_coords = NULL,
                                             cytoband_file = 'resources/cytoBand_hg38.txt.gz') {
    build <- match.arg(build)

    # --- 1) rsID -> chr + pos ---
    # Prefer a pre-loaded GWAS coords table (data.frame with columns SNP, CHR, BP)
    # over a live biomaRt query. The PD GWAS base file covers all driver SNPs.
    if (!is.null(gwas_coords)) {
        snp <- gwas_coords |>
            filter(SNP %in% unique(rsids)) |>
            transmute(rsid = SNP, chr = as.character(CHR), pos = as.integer(BP)) |>
            distinct() |>
            filter(chr %in% c(as.character(1:22), "X", "Y"))
        missing <- setdiff(unique(rsids), snp$rsid)
        if (length(missing) > 0)
            message(length(missing), " rsIDs not found in gwas_coords (e.g. ", missing[1], ")")
    } else {
        # Fall back to biomaRt when no local coords are provided
        if (build == "hg19") {
            snp_mart <- useEnsembl("snp", dataset = "hsapiens_snp", host = "grch37.ensembl.org")
        } else {
            snp_mart <- useEnsembl("snp", dataset = "hsapiens_snp")
        }
        snp <- getBM(
            attributes = c("refsnp_id", "chr_name", "chrom_start"),
            filters    = "snp_filter",
            values     = unique(rsids),
            mart       = snp_mart
        ) |>
            transmute(rsid = refsnp_id, chr = chr_name, pos = chrom_start) |>
            distinct() |>
            filter(chr %in% c(as.character(1:22), "X", "Y"))
    }

    if (nrow(snp) == 0) return(tibble())

    snp_gr <- GRanges(seqnames = snp$chr, ranges = IRanges(snp$pos, snp$pos), rsid = snp$rsid)

    # --- 2) cytoband via local cached file (UCSC hg38 format) ---
    # Falls back to live UCSC query only if the local file is absent.
    if (file.exists(cytoband_file)) {
        cyto <- read.table(gzfile(cytoband_file), sep = "\t",
                           col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
                           stringsAsFactors = FALSE) |>
            dplyr::transmute(
                chr   = sub("^chr", "", chrom),
                start = chromStart + 1L,
                end   = chromEnd,
                band  = name
            )
    } else {
        message("Local cytoband file not found; trying live UCSC query.")
        session <- rtracklayer::browserSession("UCSC")
        rtracklayer::genome(session) <- if (build == "hg19") "hg19" else "hg38"
        cyto <- rtracklayer::getTable(rtracklayer::ucscTableQuery(session, table = "cytoBand")) |>
            dplyr::transmute(chr = sub("^chr", "", chrom), start = chromStart + 1L,
                             end = chromEnd, band = name)
    }

    cyto_gr <- GRanges(seqnames = cyto$chr, ranges = IRanges(cyto$start, cyto$end), band = cyto$band)

    ov <- findOverlaps(snp_gr, cyto_gr, ignore.strand = TRUE)
    snp$cytoband <- NA_character_
    snp$cytoband[queryHits(ov)] <-
        paste0(as.character(seqnames(snp_gr))[queryHits(ov)], mcols(cyto_gr)$band[subjectHits(ov)])

    # --- 3) nearest gene via TxDb (no network required) ---
    suppressPackageStartupMessages(
        library(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly = TRUE)
    )
    suppressPackageStartupMessages(
        library(org.Hs.eg.db, quietly = TRUE)
    )

    txdb   <- TxDb.Hsapiens.UCSC.hg38.knownGene
    gene_ranges <- suppressMessages(genes(txdb))   # GRanges with entrez IDs
    # Map entrez -> gene symbol
    sym_map <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys  = as.character(gene_ranges$gene_id),
                                     columns = "SYMBOL",
                                     keytype = "ENTREZID")
    mcols(gene_ranges)$gene <- sym_map$SYMBOL[match(gene_ranges$gene_id, sym_map$ENTREZID)]

    # Strip "chr" prefix on snp_gr seqnames so they match TxDb's "chr1" format
    seqlevels(snp_gr) <- paste0("chr", seqlevels(snp_gr))
    gene_ranges <- keepStandardChromosomes(gene_ranges, pruning.mode = "coarse")

    idx <- nearest(snp_gr, gene_ranges, ignore.strand = TRUE)
    snp$closest_gene <- mcols(gene_ranges)$gene[idx]
    snp$distance_bp  <- mcols(distanceToNearest(snp_gr, gene_ranges[idx],
                                                 ignore.strand = TRUE))$distance

    gene_gr <- gene_ranges   # kept for the return below; seqnames carry "chr" prefix
    
    idx <- nearest(snp_gr, gene_gr, ignore.strand=TRUE)
    snp$closest_gene <- mcols(gene_gr)$gene[idx]
    snp$distance_bp  <- mcols(distanceToNearest(snp_gr, gene_gr[idx], ignore.strand=TRUE))$distance
    
    snp |> arrange(rsid)
}

plot_cm <- function(df) {
    plot_df <- df[, c("SNP", "chr", "pos", "ABS_DELTA")]
    colnames(plot_df) <- c("SNP", "Chromosome", "Position", "ABS_DELTA")
    present <- sort(unique(plot_df$Chromosome))
    missing <- setdiff(1:22, present)
    if (length(missing) > 0) {
        dummy <- data.frame(
            SNP = paste0("dummy_chr", missing),
            Chromosome = missing,
            Position = 1,
            ABS_DELTA = 0
        )
        plot_df <- rbind(plot_df, dummy)
    }
    
    CMplot(
        plot_df,
        file.output = FALSE,
        plot.type = "m",                 # Manhattan
        col = c('#F8766D', '#00BFC4'),
        cex = 0.6,
        ylab = "|Δ mean PRS contribution|",
        threshold = NULL,                # no significance line
        LOG10 = FALSE,                   # crucial
        chr.den.col = NULL,              # turn off density track
        chr.labels = c(1:22),
        chr.border = TRUE,
        ylim = c(0, 0.13)
    )
}

calc_contribution_scores <- function(driver_table, n_A = NULL, n_B = NULL) {
    d <- setDT(driver_table %>% filter(!is.na(pos)))
    # LD blocks BED: columns = chr start end block_id  (0-based BED start; end exclusive)
    b <- fread("/earth/public_data/ldetect/LDblocks_GRCh38/data/deCODE_EUR_LD_blocks.bed", header = TRUE)
    b[, block_id := .I]
    setnames(b, c("chr","start","end", "block_id"))
    b[, chr := as.character(gsub("^chr","", chr))]

    # Convert SNP BP to 0-based for BED overlap logic
    d[, pos0 := pos - 1L]
    d[, pos1 := pos - 1L]

    # data.table interval overlap
    setkey(b, chr, start, end)
    setkey(d, chr, pos0, pos0)

    ov <- foverlaps(
        d, b,
        by.x = c("chr","pos0","pos1"),
        by.y = c("chr","start","end"),
        type = "within",
        nomatch = 0L
    )

    x <- as.data.table(ov)
    x <- x[!is.na(DELTA_CONTRIB)]

    topk_abs_sum <- function(v, k=10L) {
        v <- abs(v)
        v <- v[is.finite(v)]
        if (!length(v)) return(NA_real_)
        sum(sort(v, decreasing = TRUE)[1:min(k, length(v))])
    }

    block_metrics <- x[, .(
        n_snps        = .N,
        sum_delta     = sum(DELTA_CONTRIB, na.rm = TRUE),
        l2            = sqrt(sum(DELTA_CONTRIB^2, na.rm = TRUE)),
        signed_l2     = sign(sum(DELTA_CONTRIB, na.rm = TRUE)) *
                        sqrt(sum(DELTA_CONTRIB^2, na.rm = TRUE)),
        max_abs       = max(abs(DELTA_CONTRIB), na.rm = TRUE),
        top10_abs_sum = topk_abs_sum(DELTA_CONTRIB, k = 5L),
        mean_abs      = mean(abs(DELTA_CONTRIB), na.rm = TRUE)
    ), by = .(block_id)]

    # BETA-weighted Stouffer Z per block.
    # Z_i = (FREQ_0 - FREQ_1) / SE_i, where SE_i is the Binomial SE under the
    # pooled-frequency null.  Block Z = sum(|BETA| * Z_i) / sqrt(sum(BETA^2)).
    # Positive Z => enriched in group 0 (AD); negative => enriched in group 1 (AD+LBP).
    # Requires n_A (group 0 size) and n_B (group 1 size).
    if (!is.null(n_A) && !is.null(n_B) &&
        all(c("FREQ_0", "FREQ_1", "BETA") %in% names(x))) {
        stouffer_metrics <- x[
            !is.na(FREQ_0) & !is.na(FREQ_1) & !is.na(BETA),
            {
                fp    <- (FREQ_0 * 2 * n_A + FREQ_1 * 2 * n_B) / (2L * (n_A + n_B))
                se_i  <- sqrt(fp * (1 - fp) * (1/n_A + 1/n_B))
                z_i   <- ifelse(se_i > 0, (FREQ_0 - FREQ_1) / se_i, 0)
                w     <- abs(BETA)
                denom <- sqrt(sum(w^2))
                sz    <- if (denom > 0) sum(w * z_i) / denom else NA_real_
                .(stouffer_z = sz, stouffer_p = 2 * pnorm(-abs(sz)))
            },
            by = .(block_id)
        ]
        block_metrics <- merge(block_metrics, stouffer_metrics, by = "block_id", all.x = TRUE)
        block_metrics[, stouffer_q := p.adjust(stouffer_p, method = "BH")]
    }

    coord_cols <- intersect(names(x), c("CHR","chr","start","end"))
    if (length(coord_cols) > 0) {
        coords <- x[, lapply(.SD, function(v) v[which(!is.na(v))[1]]), by=.(block_id), .SDcols=coord_cols]
        block_metrics <- merge(coords, block_metrics, by="block_id", all.y=TRUE)
    }
    block_metrics
}

map_tx <- function(x) {
    if (is.null(.mart_cache)) map_gene(character(0))   # initialise cache
    if (is.null(.mart_cache)) return(data.frame())
    getBM(attributes = c("ensembl_transcript_id", "external_gene_name", "external_transcript_name"),
          filters = "ensembl_transcript_id", values = x, mart = .mart_cache)
}

# mart / map_gene: lazy-initialised; only connects to Ensembl when called
.mart_cache <- NULL
map_gene <- function(x) {
    if (is.null(.mart_cache)) {
        .mart_cache <<- tryCatch(
            biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl"),
            error = function(e) {
                message("biomaRt unavailable (", conditionMessage(e), "). map_gene() will return NA gene names.")
                NULL
            }
        )
    }
    if (is.null(.mart_cache)) {
        return(data.frame(ensembl_gene_id = x, external_gene_name = NA_character_,
                          stringsAsFactors = FALSE))
    }
    getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
          filters = "ensembl_gene_id", values = x, mart = .mart_cache)
}

annotate_ld_locus <- function(driver_table,
                              ld_bed = "/earth/public_data/ldetect/LDblocks_GRCh38/data/deCODE_EUR_LD_blocks.bed") {
    
    # Convert to data.table and require valid coordinates
    d <- data.table::as.data.table(driver_table)
    d <- d[!is.na(chr) & !is.na(pos)]
    
    # Read LD blocks (BED: 0-based start, end exclusive)
    b <- data.table::fread(ld_bed, header = TRUE)
    
    # Standardize column names
    if ("block_id" %in% names(b)) {
        data.table::setnames(b, c("chr","start","end","block_id"))
    } else {
        data.table::setnames(b, names(b)[1:3], c("chr","start","end"))
        b[, block_id := .I]
    }
    
    # Harmonize chromosome format
    b[, chr := as.character(gsub("^chr","", chr))]
    d[, chr := as.character(gsub("^chr","", chr))]
    
    # Convert SNP position (1-based) → 0-based for BED overlap
    d[, pos0 := as.integer(pos) - 1L]
    d[, pos1 := as.integer(pos) - 1L]
    
    # Set keys for interval overlap
    data.table::setkey(b, chr, start, end)
    data.table::setkey(d, chr, pos0, pos1)
    
    # Overlap: keep all SNPs (even if no LD block match)
    annotated <- data.table::foverlaps(
        d, b,
        by.x = c("chr","pos0","pos1"),
        by.y = c("chr","start","end"),
        type = "within",
        nomatch = NA
    )
    
    # Drop temporary interval columns but keep everything else
    annotated[, c("pos0","pos1") := NULL]
    
    return(annotated[])
}
