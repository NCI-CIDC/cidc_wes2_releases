#############################################################################################################
# pipeline_template
#
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
#
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#
# Program:  Snakefile
# Version:  0.1
# Author:   Sami R. Cherikh, Travis L. Jensen
# Purpose:  Main snakefile for workflow template.
# Input:    Sample fastq files stored in a cloud location (google cloud, aws)
# Output:   'sample_metadata.csv','rseqc/bam_qc_parsed.tab', 'rseqc/bam_gc_parsed.tab'
#############################################################################################################

## Import modules
import shutil
import logging as _logging
import psutil
import os
from pathlib import Path
from box import Box
import yaml
import pandas as pd

wildcard_constraints:
    sample="[^_]+"


##############################
#       CONFIGURATION        #
##############################
# Specify YAML config file location
configfile: "config/config.yaml"

# Directories
# working output dir
PREDIR     = config["predir"]
# source dir for supporting scripts
SOURCEDIR  = config["srcdir"]
# analysis data results and reporting within working dir
DATADIR    = PREDIR+'/analysis/data'
REPDIR     = PREDIR+'/analysis/report'


# Use the source dir to import helper modules
try:
    sys.path.append(SOURCEDIR+'/python')
    import trimadapters
    import getfile
    import putfile
    import utils
    import time
except:
    print("####\nThe srcdir value in the config file has not been properly configured. \n \
           Please configure the config/config.yaml file and try again.\n####")

#added back in for to_log and to_benchmark functions
include: "./rules/common.smk"


## create file accessor
paths = create_path_accessor()

## read in reference genome locations file
reference_df = pd.read_table(config["reference"],
	     sep=",",
	     comment = "#")

## read in sample metadata file
sample_metadata_df = pd.read_table(config["sample_metadata"],
		   sep=",",
		   keep_default_na=False,
		   comment = "#")

pairings_df = pd.read_table(config["pairings"],
	     sep=",",
	     comment = "#").set_index("run_name", drop = False)

tumor_only_df = pairings_df.query("type == 'TO'")
tumor_normal_df = pairings_df.query("type == 'TN'")

# Reference genome gcloud URI locations
GENOME_FA_URI = reference_df.loc[reference_df["ref_file_name"]=="genome_fa", "google_bucket_URI"].item()
GENOME_GTF_URI = reference_df.loc[reference_df["ref_file_name"]=="genome_gtf", "google_bucket_URI"].item()
GENOME_BWA_URI = reference_df.loc[reference_df["ref_file_name"]=="genome_bwa_index", "google_bucket_URI"].item()
GENOME_BLACKLIST_URI = reference_df.loc[reference_df["ref_file_name"]=="genome_blacklist", "google_bucket_URI"].item()
GENOME_DHS_URI = reference_df.loc[reference_df["ref_file_name"]=="genome_dhs", "google_bucket_URI"].item()
GENOME_MILLS_URI =reference_df.loc[reference_df["ref_file_name"]=="mills_ref", "google_bucket_URI"].item()
GENOME_G1000_URI =reference_df.loc[reference_df["ref_file_name"]=="g1000_ref", "google_bucket_URI"].item()
GENOME_MILLS_INDEX_URI =reference_df.loc[reference_df["ref_file_name"]=="mills_index", "google_bucket_URI"].item()
GENOME_G1000_INDEX_URI =reference_df.loc[reference_df["ref_file_name"]=="g1000_index", "google_bucket_URI"].item()
GENOME_DBSNP_URI = reference_df.loc[reference_df["ref_file_name"]=="dbsnp", "google_bucket_URI"].item()
GENOME_DBSNP_INDEX_URI = reference_df.loc[reference_df["ref_file_name"]=="dbsnp_index", "google_bucket_URI"].item()
HLA_BED_URI = reference_df.loc[reference_df["ref_file_name"]=="hla_bed", "google_bucket_URI"].item()

# Sample info
## List of samples to process
SAMID = utils.toList(sample_metadata_df['samid'])
## List of input files
FASTQ_1 = utils.toList(sample_metadata_df['fastq_file_1'])
FASTQ_2 = utils.toList(sample_metadata_df['fastq_file_2'])
BAM = utils.toList(sample_metadata_df['bam_file'])
## Adapter sequences
FP_ADAPTERS   = [x.strip() for x in utils.toList(sample_metadata_df['fivep_adapter_seq'])]
TP_ADAPTERS   = [x.strip() for x in utils.toList(sample_metadata_df['threep_adapter_seq'])]


# Set workflow working (output) dir
workdir: PREDIR



##############################
#       SET GLOBAL VARS      #
##############################
# Workflow info
## number of cores dedicated to run
NCORES  = int(config["ncores"])
## initial sub folders
SUBDIRS  = 'benchmark log info progress genome annot input analysis analysis/data analysis/report'

## Set single or paired end
if (FASTQ_2[0] != '' or BAM[0] != ''):
    ENDS  = ['1','2']
else:
    ENDS  = ['1']
## Determine whether adapters should be trimmed or not
TRIM_FP = sum([x == ''  for x in FP_ADAPTERS]) == 0
TRIM_TP = sum([x == ''  for x in TP_ADAPTERS]) == 0
TRIM_ADAPTERS_OUTPUT = '.fastq.gz' if (TRIM_FP or TRIM_TP) else '.skipped'
## Set whether using FASTQ or BAM input
if (FASTQ_1 == [''] * len(SAMID)):
    BAM_INPUT = True
else:
    BAM_INPUT = False

# Preprocessing options
## Quality trimming option
QUAL_CUTOFF = config["quality_trim"]


# Cloud options retrieving files and archiving results
CLOUD  = config["cloud_prog"]
ARCHIVE = config["archive_bucket"]
DOARCHIVE = len(ARCHIVE) != 0



# Set up logging to log file
LOG_FILE  = config["log_file"]
_logging.basicConfig(level=_logging.INFO,
                    format='[cidc_wes2] %(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    handlers=[_logging.FileHandler(LOG_FILE)])



################################
#     DEFINE TARGET OUTPUT     #
################################
OUTPUT = [
          expand(paths.rseqc.bamqc_txt, sample=SAMID),
          expand(paths.rseqc.bamgc_txt, sample=SAMID),
          expand(paths.fastqc.targz, sample=SAMID),
          expand(paths.bqsr.report, sample = SAMID),
#          expand(paths.xhla.report, sample = SAMID),
]



#########################################
#    Define any onstart or onsuccess    #
#########################################
onsuccess:
    ## Copy sample_metadata.csv to the PREDIR
    shell('cp '+SOURCEDIR+'/../'+config["sample_metadata"]+' '+PREDIR)

    ## Merge sample rseqc results into single result files
    merged_results = utils.mergeRSEQC(SOURCEDIR)

    ## Copy some results to analysis data dir
    [shutil.copy2(x, DATADIR) for x in merged_results]

    ## knit rmarkdown html report
    #shell('Rscript --vanilla '+SOURCEDIR+'/r/run-report.r '+SOURCEDIR+'/r/cidc_atac-report-slidy.Rmd '+PREDIR+' '+DATADIR+'/../report')
    #merged_results.append('analysis/report/cidc_atac-report-slidy.html')

    ## Upload main results if needed
    if DOARCHIVE:
        [putfile.upload(file=x, destination=ARCHIVE, prog=CLOUD) for x in merged_results]

    shell("echo 'Pipeline complete!'")



################################
#   PASS OUTPUT TO all RULE    #
################################
rule all:
    input:
        OUTPUT

#just for testing ingestion of reference files
rule references:
    input:
        [
        paths.genome.mills,
        ]

################################
#        PIPELINE RULES        #
################################
include: "./rules/initialization.smk"
include: "./rules/ingest.smk"
include: "./rules/mapping.smk"
include: "./rules/realignment_recalibration_gatk.smk"
include: "./rules/xHLA.smk"