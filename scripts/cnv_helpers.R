# library(tidyverse)
# library(glue)

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS FOR ASCAT WORKFLOW
# ------------------------------------------------------------------------------
# 
# -----------------------------UTILS--------------------------------------------
get.path <- function(...) {
  # convenience function
  return(glue(..., .sep="/"))
}


# ---------------------LOAD FILES FOR ASCAT OBJECT ---------
ascat.load <- function(normal.base="", tumor.base="", type="tumor", input.dir=".", gender="XX") {
  # helper for loading the ascat file
  # make the filename base

  tumor.base <- get.path(input.dir, tumor.base)
  normal.base <- get.path(input.dir, normal.base)
  # sample.file <- glue(sample.base, "_", type)

  ascat.bc <- ascat.loadData(Tumor_LogR_file = glue(tumor.base, "_logr.tsv"),
                             Tumor_BAF_file = glue(tumor.base, "_baf.tsv"),
                             Germline_LogR_file = glue(normal.base, "_logr.tsv"),
                             Germline_BAF_file = glue(normal.base, "_baf.tsv"),
                             gender = gender) 
  return(ascat.bc)
}

# ------------------------------------------------------------------------------
# RUN ASCAT
# ------------------------------------------------------------------------------

run.ascat <- function(
    sample="",
    tumor="tumor",
    normal="normal",
    sample.dir=".",
    output.dir="."
    ) {
    sample.name <- glue(sample, "_", tumor, "-", normal)
    normal.base <- glue(sample.name, "_", normal)
    tumor.base <- glue(sample.name, "_", tumor)

    ascat.bc <- ascat.load(normal.base=normal.base, tumor.base=tumor.base, input.dir=sample.dir)
    
    ascat.plotRawData(ascat.bc, img.dir=output.dir)


    for (pen in ascat.penalty) {
        out.path <- glue(output.dir, "/", pen)
        dir.create(out.path, showWarnings = FALSE)
        ascat.bc <- ascat.aspcf(ascat.bc, penalty = pen, out.dir=out.path)
    
        ascat.plotSegmentedData(ascat.bc, img.dir=out.path)
    
        ascat.output <- ascat.runAscat(
            ascat.bc, 
            gamma = 1.0,    # actual log-diff between n=1 and n=2 --> 1 for NGS
            img.dir=out.path
        )
        # write the segment file from ascat output
        ascat.output$segments %>% 
        dplyr::rename(chrom=chr, start=startpos, end=endpos, cn_major=nMajor, cn_minor=nMinor) %>% 
        drop_na() %>% 
        select(-c(sample)) %>% 
        write_tsv(glue(out.path, "/", sample.name, "_", pen, "_ascat_segs.tsv"))
        rphase_sample_data <- data.frame(sample_id = sample.name, 
                                 segmentation = toString(glue(sample.file, "_ascat_segs.tsv.gz")), 
                                 snps = toString(glue(sample.file, "_snps.tsv.gz")), 
                                 purity = ascat.output$aberrantcellfraction[[1]],
                                 ploidy = ascat.output$ploidy[[1]]
                                 ) %>% 
        write_tsv(get.path(out.path, "rphase_sample_data.tsv"))
    }
}

run.ascat2 <- function(
    sample="",
    type="tumor",
    sample.dir=".",
    output.dir="."
    ) {
        germline.base <- glue(sample, "_", "germline")
        sample.base <- glue(sample, "_", type)

    ascat.bc <- ascat.load(normal.base=germline.base, tumor.base=sample.base, input.dir=sample.dir)
    
    ascat.plotRawData(ascat.bc, img.dir=output.dir)


    for (pen in ascat.penalty) {
        out.path <- glue(output.dir, "/", pen)
        dir.create(out.path, showWarnings = FALSE)
        ascat.bc <- ascat.aspcf(ascat.bc, penalty = pen, out.dir=out.path)
    
        ascat.plotSegmentedData(ascat.bc, img.dir=out.path)
    
        ascat.output <- ascat.runAscat(
            ascat.bc, 
            gamma = 1.0,    # actual log-diff between n=1 and n=2 --> 1 for NGS
            img.dir=out.path
        )
        # write the segment file from ascat output
        ascat.output$segments %>% 
        dplyr::rename(chrom=chr, start=startpos, end=endpos, cn_major=nMajor, cn_minor=nMinor) %>% 
        drop_na() %>% 
        select(-c(sample)) %>% 
        write_tsv(glue(output.dir, "/", sample.base, "_", pen, "_ascat_segs.tsv"))
    }
}