library(devtools)
library(sequenza)
## TAR variable set; otherwise, "sh: 1: /bin/gtar:not found" error occurs during download of GitHub repo aroneklund/copynumber
Sys.setenv(TAR = '/bin/tar')
devtools::install_github("aroneklund/copynumber")

run_sequenza <- function(arg_in, arg_out, arg_name) {
    options("scipen"=100, "digits"=4) 
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
    
    tsv <- sequenza:::sequenza2PyClone(mut.tab, seg.res, arg_name, norm.cn=2)
    tsv <- tsv[tsv[,'major_cn'] != 0 ,] # Required for Pyclone
    write.table(tsv,paste0(arg_out,"/",arg_name,'_pyclone.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
}

args <- commandArgs(trailingOnly=TRUE)
arg_in <- args[1]
arg_out <- args[2]
arg_name <- args[3]

run_sequenza(arg_in, arg_out, arg_name)
