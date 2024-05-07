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

### read in reference genome locations file
ref_df = load_reference_data()

## read in sample metadata file
sample_metadata_df = pd.read_table(config["sample_metadata"],
		   sep=",",
		   keep_default_na=False,
		   comment = "#")

pairings_df = pd.read_table(config["pairings"],
	     sep=",",
	     comment = "#").set_index("run", drop = False)

tumor_only_df = pairings_df.query("type == 'TO'")
tumor_normal_df = pairings_df.query("type == 'TN'")
print(pairings_df)
print(tumor_normal_df)


#todo: wrap these in a function in rules/common.smk
GENOME_FA_URI = grab_ref_URI(ref_df, "genome_fa")
GENOME_GTF_URI = grab_ref_URI(ref_df,"genome_gtf")
GENOME_BWA_URI = grab_ref_URI(ref_df,"genome_bwa_index")

#reference SNVs for BQSR
GENOME_MILLS_URI = grab_ref_URI(ref_df,"mills_ref")
GENOME_G1000_URI = grab_ref_URI(ref_df,"g1000_ref")
GENOME_MILLS_INDEX_URI = grab_ref_URI(ref_df,"mills_index")
GENOME_G1000_INDEX_URI = grab_ref_URI(ref_df,"g1000_index")
GENOME_DBSNP_URI = grab_ref_URI(ref_df,"dbsnp")
GENOME_DBSNP_INDEX_URI = grab_ref_URI(ref_df,"dbsnp_index")

#reference data for xHLA
HLA_BED_URI = grab_ref_URI(ref_df,"hla_bed")
HLA_TSV_URI = grab_ref_URI(ref_df,"hla_tsv")
HLA_FNA_URI = grab_ref_URI(ref_df,"hla_fna")
HLA_DMND_URI = grab_ref_URI(ref_df,"hla_dmnd")
HLA_FAA_URI = grab_ref_URI(ref_df,"hla_faa")
HLA_FREQ_URI = grab_ref_URI(ref_df,"hla_freq")
HLA_SHIFT_URI = grab_ref_URI(ref_df,"hla_shift")
GENOME_COVERAGE_TARGETS = grab_ref_URI(ref_df,"coverage_targets")

#reference data for HLAHD
HLAHD_DICT_URI = grab_ref_URI(ref_df,"hlahd_dict")
HLAHD_FREQ_URI = grab_ref_URI(ref_df,"hlahd_freq")
HLAHD_SPLIT_URI = grab_ref_URI(ref_df,"hlahd_split")

# Reference data for MSIsensor2
MSISENSOR2_MODELS_URI = grab_ref_URI(ref_df,"msisensor2_ref")

# Reference data for HaplotypeCaller (Germline module)
GERMLINE_DBSNP_URI = grab_ref_URI(ref_df,"germline_dbsnp")
GERMLINE_INDEX_URI = grab_ref_URI(ref_df,"germline_index")

# CIMAC Center variants BED used in the Germline module
if config["cimac"] == "mocha":
    TARGETS_BED_URI = grab_ref_URI(ref_df,"targets_mocha")
elif config["cimac"] == "mda":
    TARGETS_BED_URI = grab_ref_URI(ref_df,"targets_mda")
else:
    TARGETS_BED_URI = grab_ref_URI(ref_df,"targets_broad")

# Reference data for use in FACETS (Purity and Copy Number modules)
FACETS_VCF_URI = grab_ref_URI(ref_df,"facets_vcf") 
FACETS_TBI_URI = grab_ref_URI(ref_df,"facets_tbi") 

# Reference data for use in Sequenza (Clonality and Copy Number modules)
SEQUENZA_WIG_URI = grab_ref_URI(ref_df,"sequenza_wig")

# Reference data for use in TcellExTRECT module
TCELLEXTRECT_BED_URI = grab_ref_URI(ref_df,"tcellextrect_bed")

# Reference data for Mutect2
AF_VCF_URI = grab_ref_URI(ref_df,"af_vcf")
AF_INDEX_URI = grab_ref_URI(ref_df,"af_index")

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
## List of the run names for each pairing or entry (tumor normal or tumor only)
RUN = utils.toList(pairings_df['run'])
## List of run names for each tumor normal pairing
TN = utils.toList(tumor_normal_df['run'])

# Set workflow working (output) dir
workdir: PREDIR



##############################
#       SET GLOBAL VARS      #
##############################
# Workflow info
## number of cores dedicated to run
NCORES  = int(config["ncores"])
## initial sub folders
SUBDIRS  = 'benchmark log info progress genome annot input analysis analysis/data analysis/report hlahd_references resources/vep/cache'

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
          expand(paths.bqsr.report, sample=SAMID),
          expand(paths.optitype.tsv, sample=SAMID),
          expand(paths.xhla.report, sample = SAMID),
          expand(paths.coverage.depth, sample=SAMID),
          expand(paths.coverage.bw, sample=SAMID),
          expand(paths.cnv.csv, sample=SAMID),
      	  paths.hlahd_references.dict_done,
          expand(paths.hlahd.done, sample= SAMID),
          expand(paths.msisensor2.output, sample=RUN),
          expand(paths.germline.tbi, sample=SAMID),
          expand(paths.germline.txt, sample=TN),
          expand(paths.germline.pdf, sample=TN),
          expand(paths.facets.opt, sample=TN), 
          expand(paths.sequenza.segments, sample=TN),
          expand(paths.pyclone6.summary, sample=TN),
          expand(paths.copynumber.seq_fac, sample=TN),
          expand(paths.tcellextrect.pdf, sample=RUN),
          expand(paths.mutect2.pon, sample=TN),
          expand(paths.mutect2.somatic_calls_vcf, sample=TN),
          expand(paths.mutect2.filtered_somatic_calls_vcf, sample=TN)
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
include: "./rules/bqsr.smk"
include: "./rules/xHLA.smk"
include: "./rules/coverage.smk"
include: "./rules/cnv.smk"
include: "./rules/hlahd.smk"
include: "./rules/optitype.smk"
include: "./rules/msisensor2.smk"
include: "./rules/germline.smk"
include: "./rules/facets.smk"
include: "./rules/vep_plugins.smk"
include: "./rules/vep_cache.smk"
include: "./rules/sequenza.smk"
include: "./rules/pyclone6.smk"
include: "./rules/copynumber.smk"
include: "./rules/tcellextrect.smk"
include: "./rules/mutect2.smk"
