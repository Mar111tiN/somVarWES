
############# IMPORTS ######################
suppressMessages({
    library(glue)
    library(tidyverse, warn.conflicts = FALSE)
    library(gtools, warn.conflicts = FALSE)
    library(ASCAT, warn.conflicts = FALSE)
})


# load the helper functions
source(glue("{snakemake@scriptdir}/cnv_helpers.R"))


########## SNAKEMAKE GENERICS #####################
{
    i <- snakemake@input
    o <- snakemake@output
    w <- snakemake@wildcards
    p <- snakemake@params
    c.a <- snakemake@config[['CNV']][['ascat']]

    tumor <- w[['tumor']]
    normal <- w[['normal']]
    sample <- w[["sample"]]
    type <- w[["type"]]
}


########## PARAMS/CONFIG #########################
{
    gender <- c.a[['gender']]
    ascat.penalty <- c(25,75)
    gamma <- c.a[['gamma']]
    calculate.normal <- c.a[['calculate_normal']]
}


## get the file paths and names
{
    # create the sample directory
    sample.dir <- dirname(i[["sample_baf"]])
    # snakemake should take care of folder creation
    # dir.create(sample.dir, showWarnings = FALSE)
    output.dir <- dirname(o[["tumor"]])

    # get the sample names from wildcards
    sample.name <- glue(sample, "_", tumor, "-", normal)
    # get the full path to sample base name
    full.name <- glue(sample.name, "_", type)
    # e.g. "cnv/67/67_A-B_A"
    sample.file <- get.path(sample.dir, full.name)
}


# ------------------------------------------------------------------------------
# RUN ASCAT
# ------------------------------------------------------------------------------
ascat.output <- run.ascat2(
    sample=sample.name,   # 03_A-B
    type=type,       # A|B
    sample.dir=sample.dir,
    output.dir=output.dir
    )
