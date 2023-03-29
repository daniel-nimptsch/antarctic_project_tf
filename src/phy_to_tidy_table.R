#
# phy_to_tidy_table.R
#
# Written by Daniel Nimptsch
#

require(phyloseq)
require(tidyverse)
source("src/table_manipulation.R")

# returns tidy_table
# this is the main function for requesting tidy table out of phyloseq obj
# required is a phyloseq object with tax_table columns path and taxgroup_name
# if you want presence/absent counts spcify count = T
# relative abundance: percent = T
get_tidy_table_phy_object <- function(phy, taxgroup_name, count = FALSE, percent = FALSE) {
    # Functions

    # returns phyloseq obj
    # this function takes a phyloseq object and the taxgroup_name
    # which should be a string indicating the column of the tax_table to be utilized for
    # aglomarating taxa. All taxa beloging to the species of said column will be
    # merged with the tax_glom() function
    tax_glom_strict <- function(phy, taxgroup_name) {
        main_tax_table <- get_main_tax_table_phy(phy)
        # taxgroup column
        col_rank_maintaxtable <- which(colnames(main_tax_table) == taxgroup_name)
        # keep only the taxgroup column and the path column
        phy_tax_table <- main_tax_table[, c(col_rank_maintaxtable, which(colnames(main_tax_table) == "path"))]
        # new tax_table phyloseq object
        phy_tax_table <- tax_table(as.matrix(phy_tax_table))
        taxa_names(phy_tax_table) <- taxa_names(phy)
        tax_table(phy) <- phy_tax_table
        phy <- tax_glom(phy, taxrank = taxgroup_name)
    }

    # return tidy table
    # For a phyloseq object generate a tibble in the tidy column format
    phy_object_to_tidy_table <- function(phy, taxgroup_name) {
        # extract otu_table and tax_table from phyloseq obj
        otu_table <- as.data.frame(otu_table(phy))
        tax_table <- as.data.frame(tax_table(phy))
        # check if otu_table needs to be transposed
        if (dim(otu_table)[1] != dim(tax_table)[1]) {
            otu_table <- t(otu_table)
        }
        # generate a combined table
        tidy <- cbind(otu_table, tax_table)
        # remove the path column
        tidy <- tidy[, -which(colnames(tidy) == "path")]
        # tidy_table
        tidy <- pivot_longer(tidy, which(colnames(tidy) != taxgroup_name), names_to = "sample", values_to = "value")
        # remove entiries with value 0
        ind <- which(tidy$value == 0)
        if (!is_empty(ind)) {
            tidy <- tidy[-ind, ]
        }
        return(tidy)
    }

    # use decostand to transform otu_table from pyh obc to presence/absent
    if (count) {
        otu_table(phy) <- vegan::decostand(as.matrix(otu_table(phy)), "pa")
    }
    # glomerate taxa to the specified column from the tax_table
    phy <- tax_glom_strict(phy, taxgroup_name)
    # transform sample counts if relative abundance is desired
    if (percent) {
        phy <- transform_sample_counts(phy, function(x) x / sum(x) * 100)
    }
    # transform phy obj to tidy tibble with columns
    # taxgroup_name sample and value
    tidy_table <- phy_object_to_tidy_table(phy, taxgroup_name)
    return(tidy_table)
}

get_tidy_table_list <- function(phy, taxgroup_name, dataset_path) {
    if (file.exists(file.path(dataset_path, paste("tidy_table_list_", taxgroup_name, ".rds", sep = "")))) {
        tidy_table_list <- read_rds(file.path(dataset_path, paste("tidy_table_list_", taxgroup_name, ".rds", sep = "")))
    } else {
        tidy_table <- get_tidy_table_phy_object(phy, taxgroup_name = taxgroup_name)
        tidy_table_percent <- get_tidy_table_phy_object(phy, taxgroup_name, percent = T)
        tidy_table_absolut <- get_tidy_table_phy_object(phy, taxgroup_name, count = T)
        tidy_table_absolut_percent <- get_tidy_table_phy_object(phy, taxgroup_name, count = T, percent = T)

        tidy_table <- add_sample_names(tidy_table, sample_data(phy))
        tidy_table_percent <- add_sample_names(tidy_table_percent, sample_data(phy))
        tidy_table_absolut <- add_sample_names(tidy_table_absolut, sample_data(phy))

        phy_cat <- merge_samples(phy, "sample_categories")
        tidy_table_sample_category <- get_tidy_table_phy_object(phy_cat, taxgroup_name)
        tidy_table_percent_sample_category <- get_tidy_table_phy_object(phy_cat, taxgroup_name, percent = T)
        tidy_table_absolut_sample_category <- get_tidy_table_phy_object(phy_cat, taxgroup_name, count = T)
        tidy_table_absolut_sample_category_percent <- get_tidy_table_phy_object(phy_cat, taxgroup_name, count = T, percent = T)

        tidy_table_list <- list(
            tidy_table,
            tidy_table_percent,
            tidy_table_absolut,
            tidy_table_absolut_percent,
            tidy_table_sample_category,
            tidy_table_percent_sample_category,
            tidy_table_absolut_sample_category,
            tidy_table_absolut_sample_category_percent
        )

        names(tidy_table_list) <- c(
            "tidy_table",
            "tidy_table_percent",
            "tidy_table_absolut",
            "tidy_table_absolut_percent",
            "tidy_table_sample_category",
            "tidy_table_percent_sample_category",
            "tidy_table_absolut_sample_category",
            "tidy_table_absolut_sample_category_percent"
        )

        saveRDS(tidy_table_list, file.path(dataset_path, paste("tidy_table_list_", taxgroup_name, ".rds", sep = "")))
    }
    return(tidy_table_list)
}
