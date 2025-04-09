log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Save as calculate_tmb.R
library(maftools)

# Load MAF file
maf <- read.maf(snakemake@input$maf)

# Read bed file and cacluate capture size
twist_regions_file <- snakemake@input$twist_regions
data <- read.table(twist_regions_file, header = FALSE, sep = "\t")
twist_length <- sum(data[[3]] - data [[2]]) / 1000000

# Compute TMB (Assuming exome capture size ~38Mb)
tmb_result <- tmb(maf = maf, captureSize = twist_length)

# Save output
write.table(tmb_result, file=snakemake@output$tmb, sep="\t", quote=FALSE, row.names=FALSE)
