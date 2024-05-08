## Set up directory structure based on dirs supplied in config
## Ignores non-zero exit status returned when any directories already exist
rule directory_setup:
    output:
        'progress/dirs.done'
    params:
        subdirs=SUBDIRS
    threads:1
    shell:
        '''
          mkdir {params.subdirs} -p 2> /dev/null
          touch {output}
        '''


## Download reference genome fa and supporting gtf annotation from ncbi ftp site
rule retrieve_reference_genome:
    input:
        rules.directory_setup.output
    output:
        fa=paths.genome.fa,
        gtf=paths.annot.gtf
    benchmark:
        'benchmark/retrieve_reference_genome.tab'
    log:
        'log/retrieve_reference_genome.log'
    params:
        fa_uri=GENOME_FA_URI,
        gtf_uri=GENOME_GTF_URI
    priority: 1000
    threads: 1
    shell:
        '''
          echo " gsutil cp {params.fa_uri} {output.fa}" | tee {log}
          gsutil cp {params.fa_uri} {output.fa} 2>> {log}

          echo "gsutil cp {params.gtf_uri} {output.gtf}" | tee {log}
          gsutil cp {params.gtf_uri} {output.gtf} 2>> {log}
        '''

rule retrieve_hlahd_reference_data:
    output:
        dict_done = paths.hlahd_references.dict_done,
        freq_done = paths.hlahd_references.freq_done,
	split = paths.hlahd_references.split,
    params:
        dest_folder = Path(paths.hlahd_references.dict_done).parent,
        HLAHD_DICT_URI = HLAHD_DICT_URI,
        HLAHD_FREQ_URI = HLAHD_FREQ_URI,
        HLAHD_SPLIT_URI = HLAHD_SPLIT_URI,
        dict_filename = Path(grab_field(ref_df, "hlahd_dict","google_bucket_URI")).name,
        freq_filename = Path(grab_field(ref_df, "hlahd_freq","google_bucket_URI")).name
    benchmark: "benchmark/retrieve_hlahd_refs.tab"
    log: "log/retrieve_hlahd_refs.log"
    shell:
       "echo {params.dict_filename};"
       "gsutil cp {params.HLAHD_DICT_URI} hlahd_references && "
       "gsutil cp {params.HLAHD_FREQ_URI} hlahd_references && "
       "gsutil cp {params.HLAHD_SPLIT_URI} hlahd_references && "
       "echo 'done copying files\n';"
       "tar zxvf {params.dest_folder}/{params.dict_filename} --directory {params.dest_folder} &&"
       "tar zxvf {params.dest_folder}/{params.freq_filename} --directory {params.dest_folder}  && "
       "touch {output.dict_done} && "
       "touch {output.freq_done} "


rule retrieve_ref_indel_sets:
    output:
        mills = paths.genome.mills,
        g1000 = paths.genome.g1000,
        dbsnp = paths.genome.dbsnp,
        mills_index = paths.genome.mills_index,
        g1000_index = paths.genome.g1000_index,
	dbsnp_index = paths.genome.dbsnp_index
    params:
        mills_gcp_uri = GENOME_MILLS_URI,
        g1000_gcp_uri = GENOME_G1000_URI,
        mills_index_gcp_url = GENOME_MILLS_INDEX_URI,
        g1000_index_gcp_uri = GENOME_G1000_INDEX_URI,
        dbsnp = GENOME_DBSNP_URI,
        dbsnp_index = GENOME_DBSNP_INDEX_URI,
	output_path = Path(PREDIR) / Path(paths.genome.dbsnp).parent
    shell:
        '''
        echo "Downloading gold standard indel files G1000 and Mills, along with their respective indices..." | tee {log}
        gsutil -m cp {params} 
        '''


rule retrieve_xHLA_ref_Data:
    output:
        hla_tsv = paths.genome.hla_tsv,
        hla_bed = paths.genome.hla_bed,
        hla_fna = paths.genome.hla_fna,
        hla_shift = paths.genome.hla_shift,
        hla_dmnd = paths.genome.hla_dmnd,
        hla_faa = paths.genome.hla_faa,
        hla_freq = paths.genome.hla_freq	
    params:
        hla_bed = HLA_BED_URI,
        hla_tsv = HLA_TSV_URI,
        hla_fna = HLA_FNA_URI,
        hla_shift = HLA_SHIFT_URI,
        hla_dmnd = HLA_DMND_URI,
        hla_faa = HLA_FAA_URI,
        hla_freq = HLA_FREQ_URI,	
	output_path = Path(PREDIR) / Path(paths.genome.hla_fna).parent
    shell:
        '''
        echo "Downloading reference data for xHLA..." | tee {log}
        gsutil -m cp {params} 
        '''

## Download built bwa_index files for the specified genome
## If using different genome, need to edit rule to build using 'bwa index'
rule build_bwa_index:
    input:
        rules.directory_setup.output,
        rules.retrieve_reference_genome.output.fa
    output:
        'progress/bwa_index_built.done'
    benchmark:
        'benchmark/build_bwa_index.tab'
    log:
        'log/build_bwa_index.log'
    conda:
        "../envs/bwa.yaml"
    params:
        bwa_uri=GENOME_BWA_URI,

    priority: 1000
    threads: 1
    shell:
        '''
          echo "Downloading bwa_index files from ncbi ftp associated with genome for mapping reads to GRCh38 or hg38..." | tee {log}
          gsutil cp {params.bwa_uri}/* genome &&
          touch {output}
	  
          ## export rule env details
          conda env export --no-builds > info/bwa.info
        '''

## Get genome chrom sizes for bw generation
rule genome_size:
    input:
        genome_fa=rules.retrieve_reference_genome.output.fa
    output:
        size=paths.genome.size,
        fai=paths.genome.fai
    benchmark:
        'benchmark/genome_size.tab'
    log:
        'log/genome_size.log'
    conda:
        SOURCEDIR+"/../envs/samtools.yaml"
    threads: 1
    shell:
        '''
          ## get genome chrom size
          echo "samtools faidx {input.genome_fa}" | tee {log}
          samtools faidx {input.genome_fa} 2>> {log}
          echo "cut -f1,2 {input.genome_fa}.fai > {output.size}" | tee -a {log}
          cut -f1,2 {input.genome_fa}.fai > {output.size} 2>> {log}
          
          ## export rule env details
          conda env export --no-builds > info/samtools.info
        '''

## Makes a Picard reference dictionary for use by picard re-aligner
rule picard_dictionary:
    input:
        paths.genome.fa
    output:
        paths.genome.picard_dict
    conda:
        "../envs/picard.yaml"
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output}"


## Retrieve coverage target regions
rule retrieve_coverage_targets_bed:
    output:
        paths.genome.coverage_targets
    benchmark:
        'benchmark/retrieve_coverage_targets_bed.tab'
    params:
        file_uri=GENOME_COVERAGE_TARGETS
    threads: 1
    shell:
       '''
          gsutil cp {params.file_uri} {output}
        '''

## Retrieve MSIsensor2 models:
rule retrieve_msisensor2_models:
    output:
        models=directory(paths.genome.models),
        tch='progress/msisensor2_models_downloaded.done'
    benchmark:
        'benchmark/retrieve_msisensor2_models.tab'
    log:
        'log/retrieve_msisensor2_models.log'
    params:
        uri=MSISENSOR2_MODELS_URI
    shell:
        '''
          echo "gsutil -m cp -R {params.uri} genome && touch {output.tch}" | tee {log}
          gsutil -m cp -R {params.uri} genome && touch {output.tch} 2>> {log}
        '''

## Set OptiType config.init to not delete intermediate bam files produced by RazerS3
rule optitype_config:
    output:
        tch='progress/optitype_config.done'
    benchmark:
        'benchmark/retrieve_coverage_targets_bed.tab'
    log:
        'log/optitype_config.log' 
    conda:
        "../envs/optitype.yaml"
    shell:
        '''
          script_path=$(conda info --envs | grep -E '\*' | awk '{{print $NF}}')

          echo "sed -i 's/deletebam=true/deletebam=false/' ${{script_path}}/bin/config.ini && touch {output.tch}" | tee {log}
          sed -i 's/deletebam=true/deletebam=false/' ${{script_path}}/bin/config.ini && touch {output.tch} 2>> {log}
        '''

## Retrieve dbSNP VCF from the GATK resource bundle for use in the Germline module (HaplotypeCaller)
rule retrieve_dbsnp_vcf:
    output:
        vcf=paths.genome.dbsnp_vcf,
        idx=paths.genome.dbsnp_vcf_idx
    benchmark:
        'benchmark/retrieve_dbsnp_vcf.tab'
    log:
        'log/retrieve_dbsnp_vcf.log'
    params:
        vcf_uri=GERMLINE_DBSNP_URI,
        idx_uri=GERMLINE_INDEX_URI
    shell:
        '''
          echo "gsutil cp {params.vcf_uri} {params.idx_uri} genome" | tee {log}
          gsutil cp {params.vcf_uri} {params.idx_uri} genome 2>> {log}
        '''

## Retrieve the specific CIMAC Center's targets BED for use in the Germline module
rule retrieve_targets_bed:
    output:
        bed=paths.genome.mocha if config["cimac"]=='mocha' else paths.genome.mda if config["cimac"]=='mda' else paths.genome.broad
    benchmark:
        'benchmark/retrieve_targets_bed.tab'
    log:
        'log/retrieve_targets_bed.log'
    params:
        bed_uri=TARGETS_BED_URI
    shell:
        '''
          echo "gsutil cp {params.bed_uri} {output.bed}" | tee {log}
          gsutil cp {params.bed_uri} {output.bed} 2>> {log}
        '''

## Retrieve the SNP VCF and its TBI for use with FACETS (Purity and Copy Number modules)
rule retrieve_facets_vcf:
    output:
        vcf=paths.genome.facets_vcf,
        tbi=paths.genome.facets_tbi
    benchmark:
        'benchmark/retrieve_facets_vcf.tab'
    log:
        'log/retrieve_facets_vcf.log'
    params:
        vcf_uri=FACETS_VCF_URI,
        tbi_uri=FACETS_TBI_URI
    shell:
        '''
          echo "gsutil -m cp {params.vcf_uri} {params.tbi_uri} genome" | tee {log}
          gsutil -m cp {params.vcf_uri} {params.tbi_uri} genome 2>> {log}
        '''

## Retrieve the GC WIG for use with Sequenza (Clonality and Copy Number modules)
rule retrieve_sequenza_wig:
    output:
        wig=paths.genome.wig
    benchmark:
        'benchmark/retrieve_sequenza_wig.tab'
    log:
        'log/retrieve_sequenza_wig.log'
    params:
        wig_uri=SEQUENZA_WIG_URI
    shell:
        '''
          echo "gsutil cp {params.wig_uri} {output.wig}" | tee {log}
          gsutil cp {params.wig_uri} {output.wig} 2>> {log}
        '''

## Install copynumber module in the Sequenza conda environment
rule install_copynumber:
    output:
        done='progress/install_copynumber.done'
    benchmark:
        'benchmark/install_copynumber.tab'
    log:
        'log/install_copynumber.log'
    conda:
        "../envs/sequenza.yaml"
    params:
        r=Path(SOURCEDIR) / "r" / "install-copynumber.r",
        done=PREDIR+"/progress/install_copynumber.done"
    shell:
        '''
          echo "Rscript {params.r} {params.done}" | tee {log}
          Rscript {params.r} {params.done} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/sequenza.info
        '''

## Retrieve the BED for use in the TcellExTRECT module
rule retrieve_tcellextrect_bed:
    output:
        bed=paths.genome.tcellextrect_bed
    benchmark:
        'benchmark/retrieve_tcellextrect_bed.tab'
    log:
        'log/retrieve_tcellextrect_bed.log'
    params:
        bed_uri=TCELLEXTRECT_BED_URI
    shell:
        '''
          echo "gsutil cp {params.bed_uri} {output.bed}" | tee {log}
          gsutil cp {params.bed_uri} {output.bed} 2>> {log}
        '''

## Clone the TcellExTRECT repo, reset the repo to specified commit, and install TcellExTRECT in its respective conda environment
rule install_tcellextrect:
    output:
        done='progress/install_tcellextrect.done'
    benchmark:
        'benchmark/install_tcellextrect.tab'
    log:
        'log/install_tcellextrect.log'
    conda:
        "../envs/tcellextrect.yaml"
    params:
        repo=PREDIR+"/repo/TcellExTRECT",
        commit=config["commit"],
        predir=PREDIR,
        r=Path(SOURCEDIR) / "r/install-tcellextrect.r",
        done=PREDIR+"/progress/install_tcellextrect.done"
    shell:
        '''
          echo "git clone https://github.com/McGranahanLab/TcellExTRECT.git {params.repo} \
          && cd {params.repo} \
          && git reset {params.commit} --hard \
          && cd {params.predir} \
          && Rscript {params.r} {params.repo} {params.done}" | tee {log}

          git clone https://github.com/McGranahanLab/TcellExTRECT.git {params.repo} \
          && cd {params.repo} \
          && git reset {params.commit} --hard \
          && cd {params.predir} \
          && Rscript {params.r} {params.repo} {params.done} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/tcellextrect.info
        '''

## Retrieve the allele frequency and gnomad annotations for use with mutect2
rule retrieve_mutect2_ref:
    output:
        wig=paths.annot.af_vcf
    benchmark:
        'benchmark/retrieve_sequenza_wig.tab'
    log:
        'log/retrieve_sequenza_wig.log'
    params:
       vcf_uri = AF_VCF_URI,
       idx_uri = AF_INDEX_URI       
    shell:
        '''
          echo "gsutil cp {params} annot " | tee {log}
          gsutil cp {params} annot  2>> {log}
        '''

## This is used in variant calling on tumor-only samples
rule retrieve_1kg_pon_file:
    output:
        wig=paths.annot.kg_pon
    benchmark:
        'benchmark/retrieve_kg_pon.tab'
    log:
        'log/retrieve_kg_pon.log'
    params:
       pon_uri = KG_PON_UR
    shell:
        '''
          echo "gsutil cp {params} annot " | tee {log}
          gsutil cp {params} annot  2>> {log}
        '''