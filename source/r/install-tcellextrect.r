# This script installs TcellExTRECT in the conda environment for future use in rule tcellextrect.

args <- commandArgs(trailingOnly = TRUE)
repo_path = args[1]
output_file = args[2]

install.packages(repo_path, repos=NULL, type ='source')

file.create(output_file)
