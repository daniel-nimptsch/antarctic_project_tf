#
# Antarctic_phyloseq_heatmap_custom_functions.R
#
# This function is based on the 'plot_heatmap' function from the 'phyloseq' package,
# with modifications made by Daniel Nimptsch

plot_heatmap_ypos_right <- function(physeq, method = "NMDS", distance = "bray", sample.label = NULL,
                                    taxa.label = NULL, low = "#000033", high = "#66CCFF", na.value = "black",
                                    trans = scales::log_trans(4), max.label = 250, title = NULL, sample.order = NULL,
                                    taxa.order = NULL, first.sample = NULL, first.taxa = NULL,
                                    ...) {
    if (!is.null(taxa.order) & length(taxa.order) == 1) {
        rankcol <- which(rank_names(physeq) %in% taxa.order)
        taxmat <- as(tax_table(physeq)[, 1:rankcol], "matrix")
        taxa.order <- apply(taxmat, 1, paste, sep = "", collapse = "")
        names(taxa.order) <- taxa_names(physeq)
        taxa.order <- names(sort(taxa.order, na.last = TRUE))
    }
    if (!is.null(sample.order) & length(sample.order) == 1) {
        sample.order <- as.character(get_variable(physeq, sample.order))
        names(sample.order) <- sample_names(physeq)
        sample.order <- names(sort(sample.order, na.last = TRUE))
    }
    if (!is.null(method) & (is.null(taxa.order) | is.null(sample.order))) {
        junk <- capture.output(ps.ord <- ordinate(
            physeq, method,
            distance, ...
        ), file = NULL)
        if (is.null(sample.order)) {
            siteDF <- NULL
            trash1 <- try(
                {
                    siteDF <- scores(ps.ord,
                        choices = c(1, 2), display = "sites",
                        physeq = physeq
                    )
                },
                silent = TRUE
            )
            if (inherits(trash1, "try-error")) {
                warning(
                    "Attempt to access ordination coordinates for sample ordering failed.\n",
                    "Using default sample ordering."
                )
            }
            if (!is.null(siteDF)) {
                sample.order <- sample_names(physeq)[order(RadialTheta(siteDF))]
            }
        }
        if (is.null(taxa.order)) {
            specDF <- NULL
            trash2 <- try(
                {
                    specDF <- scores(ps.ord,
                        choices = c(1, 2), display = "species",
                        physeq = physeq
                    )
                },
                silent = TRUE
            )
            if (inherits(trash2, "try-error")) {
                warning(
                    "Attempt to access ordination coordinates for feature/species/taxa/OTU ordering failed.\n",
                    "Using default feature/species/taxa/OTU ordering."
                )
            }
            if (!is.null(specDF)) {
                taxa.order <- taxa_names(physeq)[order(RadialTheta(specDF))]
            }
        }
    }
    if (!is.null(first.sample)) {
        sample.order <- chunkReOrder(sample.order, first.sample)
    }
    if (!is.null(first.taxa)) {
        taxa.order <- chunkReOrder(taxa.order, first.taxa)
    }
    adf <- psmelt(physeq)
    adf$OTU <- as(adf$OTU, "character")
    adf$Sample <- as(adf$Sample, "character")
    if (!is.null(sample.order)) {
        adf$Sample <- factor(adf$Sample, levels = sample.order)
    } else {
        adf$Sample <- factor(adf$Sample)
    }
    if (!is.null(taxa.order)) {
        adf$OTU <- factor(adf$OTU, levels = taxa.order)
    } else {
        adf$OTU <- factor(adf$OTU)
    }
    p <- ggplot(adf, aes(x = Sample, y = OTU, fill = Abundance)) +
        geom_raster()
    if (nsamples(physeq) <= max.label) {
        p <- p + theme(axis.text.x = element_text(size = manytextsize(
            nsamples(physeq),
            4, 30, 12
        ), angle = -90, vjust = 0.5, hjust = 0))
    } else {
        p <- p + scale_x_discrete("Sample", labels = "")
    }
    if (ntaxa(physeq) <= max.label) {
        p <- p + theme(axis.text.y = element_text(size = manytextsize(
            ntaxa(physeq),
            4, 30, 12
        )))
    } else {
        p <- p + scale_y_discrete("OTU", labels = "", position = "right")
    }
    if (!is.null(sample.label) & nsamples(physeq) <= max.label) {
        labvec <- as(get_variable(physeq, sample.label), "character")
        names(labvec) <- sample_names(physeq)
        if (!is.null(sample.order)) {
            labvec <- labvec[sample.order]
        }
        labvec[is.na(labvec)] <- ""
        p <- p + scale_x_discrete(sample.label, labels = labvec)
    }
    if (!is.null(taxa.label) & ntaxa(physeq) <= max.label) {
        labvec <- as(tax_table(physeq)[, taxa.label], "character")
        names(labvec) <- taxa_names(physeq)
        if (!is.null(taxa.order)) {
            labvec <- labvec[taxa.order]
        }
        labvec[is.na(labvec)] <- ""
        p <- p + scale_y_discrete(taxa.label, labels = labvec, position = "right")
    }
    if (!is.null(trans)) {
        p <- p + scale_fill_gradient(
            low = low, high = high, trans = trans,
            na.value = na.value
        )
    } else {
        p <- p + scale_fill_gradient(low = low, high = high, na.value = na.value)
    }
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    return(p)
}

# Define an internal function for determining what the text-size should be
#' @keywords internal
manytextsize <- function(n, mins = 0.5, maxs = 4, B = 6, D = 100) {
    # empirically selected size-value calculator.
    s <- B * exp(-n / D)
    # enforce a floor.
    s <- ifelse(s > mins, s, mins)
    # enforce a max
    s <- ifelse(s < maxs, s, maxs)
    return(s)
}
