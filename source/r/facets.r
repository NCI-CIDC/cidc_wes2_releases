#!/usr/bin/env Rscript

library(pctGCdata)
library(facets)

run_facets <- function(arg_in, arg_out, arg_name) {
    ## Read in and format the snp-pileup file
    rcmat = readSnpMatrix(arg_in, perl.pileup = F)
    ## Generate lists to hold values
    filenames <- c()
    modes <- c()
    puritys <- c()
    ploidys <- c()
    dipLogRs <- c()
    purity.vals <- c()
    ploidy.vals <- c()
    dipLogR.vals <- c()
    cncf.vals <- list()

    ## Function to find the mode of a set of values
    getmode <- function(v) {
        uniqv <- unique(v)
	uniqv[which.max(tabulate(match(v, uniqv)))]
    }
  
    ## Performs process 10 iterations
    for(i in seq(10)) {
	## Set random generator seed for each iteration since preProcSample samples the SNPs to be used
	set.seed(i)
	## Processes a snp read count matrix and generates a segmentation tree
        xx <- preProcSample(rcmat)
        ## Processes the output from preProcSample for given cval and min.nhet
        oo <- procSample(xx, cval=150)
	## Uses genotype mixture model to estimate the cluster specific copy number and cellular fraction.
	## Uses estimates based on the cnlr.median and mafR as initial values for the EM iteration
        fit <- emcncf(oo, maxiter = 1000)
	## Assigns output values to lists
        purity.vals[i] <- signif(fit[['purity']], 3)
        ploidy.vals[i] <- fit[['ploidy']]
        dipLogR.vals[i] <- oo[['dipLogR']]
        cncf.vals[[i]] <- fit[['cncf']]
	## Plots copy number log-ratio, variant allele log-odds ratio as well as the copy number and cellular fraction fits
        plotSample(x=oo, emfit=fit)
    }

    ## Obtains the purity value that occurs most often in the iterations (mode) and its index
    mode <- getmode(purity.vals)
    index <- which(purity.vals == mode)[1]
    purity <- purity.vals[index]
    ploidy <- ploidy.vals[index]
    dipLogR <- dipLogR.vals[index]
    ## Writes the selected mode and its cncf data to the output file
    write.table(data.frame(cncf.vals[[index]]), paste0(arg_out,arg_name,'.cncf'), quote=F, sep="\t", row.names=F)

    filenames <- c(filenames, arg_in)
    modes <- c(modes, mode)
    puritys <- c(puritys, purity)
    ploidys <- c(ploidys, ploidy)
    dipLogRs <- c(dipLogRs, dipLogR)

    df <- data.frame(name=filenames, mode=modes, purity=puritys, ploidy=ploidys, dipLogR=dipLogRs)
    df.iter <- data.frame(name=filenames, purity=purity.vals, ploidy=ploidy.vals, dipLogR=dipLogR.vals)
    write.table(df,paste0(arg_out,arg_name,'_optimalpurityvalue.txt'), quote=F, sep="\t", row.names=F)
    write.table(df.iter, paste0(arg_out,arg_name,'_iterpurityvalues.txt') , quote=F, sep="\t", row.names=F)
    }

args <- commandArgs(trailingOnly=TRUE)
arg_in <- args[1]
arg_out <- args[2]
arg_name <- args[3]

run_facets(arg_in, arg_out, arg_name)
