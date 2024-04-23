# The copynumber module needs to be loaded first in an initial R session to update the Sequenza conda environment.
# If the user tries to start running Sequenza immediately after installation of the copynumber module in the same R session, 
# an error will occur indicating that the lazy-load database copynumber.rdb is corrupt (see below), and execution 
# will be halted. The solution is to simply restart the R session, which will occur in the subsequent rule sequenza 
# (Rscript sequenza.r) for each sample. The error seems to be related to the version of devtools that has to be 
# utilized in this Sequenza conda environment for the pipeline.
#
# Error example snippet:
# Processing chr1: Error in get0(oNam, envir = ns) :
#   lazy-load database '/media/scratch/wes2/wes2_output/.snakemake/conda/cccc14ba54cb95430a9fac0e84f9236f_/lib/R/library/copynumber/R/copynumber.rdb' is corrupt
#
# Solutions referencing to restart the R session:
# https://cran.r-project.org/web/packages/phecodemap/vignettes/trouble.html
# https://stackoverflow.com/questions/30424608/error-in-fetchkey-lazy-load-database
# https://github.com/lme4/lme4/issues/407
#
# This script serves as the intial R session to load the copynumber module to the Sequenza conda environment.
# After successful download, an empty output file is generated to serve as input for the subsequent rule sequenza.

library(devtools)
library(sequenza)
## TAR variable set; otherwise, "sh: 1: /bin/gtar:not found" error occurs during download of GitHub repo aroneklund/copynumber
Sys.setenv(TAR = '/bin/tar')
devtools::install_github("aroneklund/copynumber")

args <- commandArgs(trailingOnly = TRUE)
output_file = args[1]

file.create(output_file)
