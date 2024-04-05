library(devtools)
library(sequenza)
devtools::install_github("aroneklund/copynumber")
#devtools::install_github("aroneklund/copynumber", ref= "2f54c08")
print(packageVersion("iotools"))

description_file <- system.file("DESCRIPTION", package = "copynumber")
description <- readLines(description_file)
package_version <- description[grep("^Version:", description)]
print(package_version)

## Function sequenza2PyClone from Sequenza 2.1.2 was added since this function is not available in 3.0.0. 
## Source code for the function is located at https://github.com/cran/sequenza/archive/refs/tags/2.1.2.zip (sequenza-2.1.2/R/next.R)
sequenza2PyClone <- function(mut.tab, seg.cn, sample.id, norm.cn = 2) {
    mut.tab <- cbind(mut.tab[,c('chromosome', 'position', 'good.reads','F', 'mutation')], CNt = NA, A = NA, B = NA)
    for (i in 1:nrow(seg.cn)) {
        pos.filt <- mut.tab$chromosome == seg.cn$chromosome[i] & mut.tab$position >= seg.cn$start.pos[i] & mut.tab$position <= seg.cn$end.pos[i]
        mut.tab[pos.filt, c("CNt", "A", "B")] <- seg.cn[i, c("CNt", "A", "B")]
    }
    id <- paste(sample.id, mut.tab$chromosome, mut.tab$position, sep = ':')
    var.counts  <- round(mut.tab$good.reads * mut.tab$F, 0)
    nor.counts  <- mut.tab$good.reads - var.counts
    pyclone.tsv <- data.frame(mutation_id = id, ref_counts = nor.counts, var_counts = var.counts,
			      normal_cn = norm.cn, minor_cn = mut.tab$B, major_cn = mut.tab$A,
			      variant_case = sample.id, variant_freq = mut.tab$F, genotype = mut.tab$mutation)
    na.exclude(pyclone.tsv)
}

run_sequenza <- function(arg_in, arg_out, arg_name) {
    options("scipen"=100, "digits"=4)
    
    ## Set the connection buffer to fit a complete line. Refer to https://github.com/tidyverse/vroom/issues/364 for the error received originally
    lines <- readLines(arg_in)
    buf_size <- sum(nchar(lines))
    print(buf_size)

    Sys.setenv(VROOM_CONNECTION_SIZE=as.integer(buf_size))
    ## Process seqz data, normalization, and segmentation
    test <- sequenza.extract(arg_in,assembly="hg38", chromosome.list=paste0("chr", c(1:22)))
    print(test)
    ## Perform parameter inference using the calculated B allele frequency and depth ratio of the obtained segments
    ##Run grid-search approach to estimate cellularity and ploidy
    CP <- sequenza.fit(test)
    ## Write files and plots using suggested or selected solution
    sequenza.results(sequenza.extract = test,cp.table = CP, sample.id = arg_name, out.dir=arg_out)
    mut.tab <- read.table(paste0(arg_out,"/",arg_name,"_mutations.txt"), header=TRUE)
    seg.res <- read.table(paste0(arg_out,"/",arg_name,"_segments.txt"),header=TRUE)
    ## fix chr list. Remark: do not use '$'
    chrs = names(test[['mutations']])
    mut.tab[['chromosome']] = factor(as.character(mut.tab[['chromosome']]), levels=chrs)
    seg.res[['chromosome']] = factor(as.character(seg.res[['chromosome']]), levels=chrs)
    
    tsv <- sequenza2PyClone(mut.tab, seg.res, arg_name, norm.cn = 2)
    tsv <- tsv[tsv[,'major_cn'] != 0 ,] # Required for Pyclone
    write.table(tsv,paste0(arg_out,"/",arg_name,'_pyclone.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
}

args <- commandArgs(trailingOnly=TRUE)
arg_in <- args[1]
arg_out <- args[2]
arg_name <- args[3]

run_sequenza(arg_in, arg_out, arg_name)
