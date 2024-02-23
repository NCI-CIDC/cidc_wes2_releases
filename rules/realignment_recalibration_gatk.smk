
##########################################
###basecall quality score recalibration###
##########################################
rule Base_recalibration_precal_gatk:
    """base recalibration for realigned files"""
    input:
        bam=paths.bam.filtered_bam,
        fa = paths.genome.fa,
        dict = paths.genome.picard_dict,
        dbsnp= paths.genome.dbsnp,
        mills= paths.genome.mills,
        g1000= paths.genome.g1000,	
    output:
        prerecaltable=paths.bqsr.prerecaltable,
        bam=paths.bqsr.recal_bam,
        bai=paths.bqsr.recal_index
    message:
        " PRE BASE RECALIBRATION: base recalibration for  realigned files"
    params:
        sentieon_path=config['sentieon_path'],
    threads: 4 #_realigner_threads
    group: "recalibration"
    conda: "../envs/gatk.yaml"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_precal_sentieon.txt"
    shell:
         """
         gatk BaseRecalibrator \
           --input {input.bam} \
           --known-sites {input.dbsnp} \
           --known-sites {input.mills} \
           --known-sites {input.g1000} \
           --output {output.prerecaltable} \
           --reference {input.fa} \
           --add-output-sam-program-record \
           --create-output-bam-index
        """

rule Base_recalibration_postcal_gatk:
    """post recalibration for realigned files"""
    input:
        bam=paths.bqsr.recal_bam,
        bai=paths.bqsr.recal_bam,
        prerecaltable=paths.bqsr.prerecaltable,
        index = paths.genome.fa,
        dbsnp= paths.genome.dbsnp,
        mills= paths.genome.mills,
        g1000= paths.genome.g1000        
    output:
        postrecaltable=paths.bqsr.postrecaltable
    message:
        "POST BASE RECALIBRATION: post base recalibration for  realigned files"
    params:
        sentieon_path=config['sentieon_path'],
    threads: 16 #_realigner_threads
    group: "recalibration"
    conda: "../envs/gatk.yaml"    
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_postcal_sentieon.txt"
    shell:
        """
        gatk ApplyBQSR \
        --bqsr-recal-file {input.precaltable} \
        --input {input.bam} \
        --output {output.postrecaltable} \
        --exclude-interval
        """

rule Base_recalibration_gatk:
    """ recalibration for realigned files"""
    input:
        postrecaltable= paths.bqsr.postrecaltable,
        prerecaltable=paths.bqsr.prerecaltable,
        index = paths.genome.fa,
    output:
        recalfile=paths.bqsr.recal_table,
        report=paths.bqsr.report,	
    message:
        "DIFF BASE RECALIBRATION: Difference in pre and post processing of realigned files"
    params:
        sentieon_path=config['sentieon_path'],
    threads: 1 #_realigner_threads
    group: "recalibration"
    conda: "../envs/gatk.yaml"    
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_sentieon.txt"
    shell:
        """
        gatk AnalyzeCovariates \
        --after-report-file {input.postrecaltable} \
        --before-report-file {input.prerecaltable} \
        --plots-report-file {output.report} \
        --intermediate-csv-file {output.recalfile}
        """