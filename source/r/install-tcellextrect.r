library(devtools)
install.packages('/home/trivittge_nih_gov/TcellExTRECT/', repos=NULL, type ='source')
#install_github("McGranahanLab/TcellExTRECT", ref="b5e8b27")

args <- commandArgs(trailingOnly = TRUE)
output_file = args[1]

file.create(output_file)
