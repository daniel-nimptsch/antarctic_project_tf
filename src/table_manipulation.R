#
# table_manipulation.R
#
# Written by Daniel Nimptsch
#

project_output_path <- function(proejct_name) {
    if (!dir.exists(str_glue("output/{proejct_name}"))) {
        dir.create(str_glue("output/{proejct_name}"))
    }
    return(str_glue("output/{proejct_name}"))
}

add_sample_names <- function(tidy_table, sample_names) {
    tidy_table <- bind_cols(tidy_table, rep(NA, nrow(tidy_table)))
    colnames(tidy_table)[4] <- "type"
    for (i in 1:nrow(tidy_table)) {
        ind <- which(tidy_table$sample[i] == sample_names[, 1])
        tidy_table$type[i] <- sample_names[ind, 2]
    }
    levels <- as_vector(unique(sample_names[, 2]))
    tidy_table$type <- factor(unlist(tidy_table$type), levels = levels, ordered = T)
    return(tidy_table)
}

get_main_tax_table_phy <- function(phy) {
    main_tax_table <- as_tibble(as.data.frame(tax_table(phy)))
    main_tax_table <- bind_cols(taxa_names(phy), main_tax_table)
    colnames(main_tax_table)[1] <- "seq_ID"
    main_tax_table <- main_tax_table %>% arrange(seq_ID)
    return(main_tax_table)
}

get_otu_table_phy <- function(phy) {
    otu_table <- as_tibble(as.data.frame(otu_table(phy)))
    otu_table <- bind_cols(taxa_names(phy), otu_table)
    colnames(otu_table)[1] <- "asv_ID"
    otu_table <- otu_table %>% arrange(asv_ID)
    return(otu_table)
}

get_taxorder <- function(main_tax_table, taxgroup_name) {
    main_tax_table <- main_tax_table %>% arrange(path)
    col_tax_table <- which(colnames(main_tax_table) == taxgroup_name)
    taxorder <- as_vector(unique(main_tax_table[, col_tax_table]))
    return(taxorder)
}

get_colors <- function(length) {
    require(RColorBrewer)
    if (length <= 12) {
        col_vector <- brewer.pal(n = length, name = "Paired")
    } else {
        n <- length
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
        col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        col_vector <- col_vector[1:n]
    }
    return(col_vector)
}
