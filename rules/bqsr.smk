
##########################################
###basecall quality score recalibration###
##########################################
rule Base_recalibration_round1_gatk:
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


rule applyBQSR_gatk_round1:
    """post recalibration for realigned files"""
    input:
        bam=paths.bam.filtered_bam,
        precaltable=paths.bqsr.prerecaltable,
        ref_fa= paths.genome.fa,
    output:
        bam=paths.bqsr.recal_bam,
        bai=paths.bqsr.recal_index
    message:
        "POST BASE RECALIBRATION: post base recalibration for  realigned files"
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
        --output {output.bam} \
        --emit-original-quals \
        --reference {input.ref_fa} \
        --add-output-sam-program-record \
        && samtools index {output.bam}
        """


rule Base_recalibration_round2_gatk:
    """base recalibration for realigned files"""
    input:
        bam=paths.bqsr.recal_bam,
        fa = paths.genome.fa,
        dict = paths.genome.picard_dict,
        dbsnp= paths.genome.dbsnp,
        mills= paths.genome.mills,
        g1000= paths.genome.g1000,	
    output:
        postrecaltable=paths.bqsr.postrecaltable,
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
           --output {output.postrecaltable} \
           --reference {input.fa} \
           --add-output-sam-program-record \
           --create-output-bam-index
        """

rule applyBQSR_gatk_round2:
    """post recalibration for realigned files"""
    input:
        bam=paths.bqsr.recal_bam,
        precaltable=paths.bqsr.prerecaltable,
        ref_fa= paths.genome.fa,
    output:
        bam=paths.bqsr.recal_bam_round2,
        bai=paths.bqsr.recal_index_round2
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
        --output {output.bam} \
        --emit-original-quals \
        --reference {input.ref_fa} \
        --add-output-sam-program-record \
        && samtools index {output.bam}
        """


rule applyBQSR_gatk_round2:
    """post recalibration for realigned files"""
    input:
        bam=paths.bqsr.recal_bam,
        precaltable=paths.bqsr.postrecaltable,
        ref_fa= paths.genome.fa,
    output:
        bam=paths.bqsr.recal_round2_bam,
        bai=paths.bqsr.recal_round2_index
    message:
        "POST BASE RECALIBRATION: post base recalibration for  realigned files"
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
        --output {output.bam} \
        --emit-original-quals \
        --reference {input.ref_fa} \
        --add-output-sam-program-record \
        && samtools index {output.bam}
        """


rule Analyze_covariates_gatk:
    """ recalibration for realigned files"""
    input:
        round2_table = paths.bqsr.postrecaltable,
        round1_table = paths.bqsr.prerecaltable,
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
    conda: "../envs/bqsr_analyze_covariates.yml"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_sentieon.txt"
    shell:
        """
        gatk AnalyzeCovariates \
        --after-report-file {input.round2_table} \
        --before-report-file {input.round1_table} \
        --plots-report-file {output.report} \
        --intermediate-csv-file {output.recalfile}
        """


rule extract_chr6:
    """Extract chr6"""
    input:
        in_bam = paths.bqsr.recal_bam_round2,
    output:
        bam = paths.bam.filtered_chr6_bam,
        bai = paths.bam.filtered_chr6_bam_index
    threads: 16
    group: "optitype"
    conda: "../envs/samtools.yaml"
    log:
        "log/{sample}_extract_chr6.log"
    benchmark:
        "benchmark/{sample}_extract_chr6.tab"
    shell:
        """samtools view -@ {threads} {input.in_bam} chr6 -b -o {output.bam} &&
           samtools index -@ {threads} {output.bam}
        """

rule make_chr6_fastqs:
    input: bam = paths.bam.filtered_chr6_bam
    output:
        r1 = paths.bqsr.recal_chr6_fq_r1,
        r2 = paths.bqsr.recal_chr6_fq_r2
    params:
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n"