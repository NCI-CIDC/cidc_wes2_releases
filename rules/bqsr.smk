## Performs the first pass of the base quality score recalibration (BQSR). Generates a recalibration table based on various covariates
rule base_recalibrator:
    input:
        bam=rules.dedup_bam.output.dedup,
        bai=rules.dedup_bam.output.bai,
        fa=paths.genome.fa,
        dict=paths.genome.picard_dict,
        dbsnp=paths.genome.dbsnp,
        mills=paths.genome.mills,
        g1000=paths.genome.g1000,	
    output:
        prerecaltable=paths.bqsr.prerecaltable,
    benchmark:
        'benchmark/{sample}_base_recalibrator.tab'
    log:
        'log/{sample}_base_recalibrator.log'
    conda:
        "../envs/gatk.yaml"
    threads: max(1,min(4,NCORES))
    shell:
        '''
          echo "gatk BaseRecalibrator \
          --input {input.bam} \
          --known-sites {input.dbsnp} \
          --known-sites {input.mills} \
          --known-sites {input.g1000} \
          --output {output.prerecaltable} \
          --reference {input.fa} \
          --add-output-sam-program-record" | tee {log}

          gatk BaseRecalibrator \
          --input {input.bam} \
          --known-sites {input.dbsnp} \
          --known-sites {input.mills} \
          --known-sites {input.g1000} \
          --output {output.prerecaltable} \
          --reference {input.fa} \
          --add-output-sam-program-record 2>> {log}
        '''

## Performs the second pass of the BQSR process and recalibrates the base qualities of the input reads based on the recalibration table produced by the BaseRecalibrator tool and outputs a recalibrated BAM or CRAM file
rule apply_bqsr:
    input:
        bam=rules.dedup_bam.output.dedup,
        bai=rules.dedup_bam.output.bai,
        precaltable=rules.base_recalibrator.output.prerecaltable,
        ref_fa=paths.genome.fa,
    output:
        bam=paths.bqsr.recal_bam,
        bai=paths.bqsr.recal_index
    benchmark:
        'benchmark/{sample}_apply_bqsr.tab'
    log:
        'log/{sample}_apply_bqsr.log'
    conda:
        "../envs/gatk.yaml"    
    threads: max(1,min(16,NCORES))
    shell:
        '''
          echo "gatk ApplyBQSR \
          --bqsr-recal-file {input.precaltable} \
          --input {input.bam} \
          --output {output.bam} \
          --emit-original-quals \
          --reference {input.ref_fa} \
          --add-output-sam-program-record \
          --create-output-bam-index" | tee {log}

          gatk ApplyBQSR \
          --bqsr-recal-file {input.precaltable} \
          --input {input.bam} \
          --output {output.bam} \
          --emit-original-quals \
          --reference {input.ref_fa} \
          --add-output-sam-program-record \
          --create-output-bam-index 2>> {log}
        '''

## Generates a recalibration table on the recalibrated bam based on various covariates
rule base_recalibrator_post:
    input:
        bam=rules.apply_bqsr.output.bam,
        bai=rules.apply_bqsr.output.bai,
        fa=paths.genome.fa,
        dict=paths.genome.picard_dict,
        dbsnp=paths.genome.dbsnp,
        mills=paths.genome.mills,
        g1000=paths.genome.g1000
    output:
        postrecaltable=paths.bqsr.postrecaltable,
    benchmark:
        'benchmark/{sample}_base_recalibrator_post.tab'
    log:
        'log/{sample}_base_recalibrator_post.log'
    conda:
        "../envs/gatk.yaml"
    threads: max(1,min(4,NCORES))
    shell:
        '''
          echo "gatk BaseRecalibrator \
          --input {input.bam} \
          --known-sites {input.dbsnp} \
          --known-sites {input.mills} \
          --known-sites {input.g1000} \
          --output {output.postrecaltable} \
          --reference {input.fa} \
          --add-output-sam-program-record" | tee {log}

          gatk BaseRecalibrator \
          --input {input.bam} \
          --known-sites {input.dbsnp} \
          --known-sites {input.mills} \
          --known-sites {input.g1000} \
          --output {output.postrecaltable} \
          --reference {input.fa} \
          --add-output-sam-program-record 2>> {log}
        '''

## Evaluate and compare base quality score recalibration tables. Generates plots to assess the quality of a recalibration run as part of the BQSR procedure
rule analyze_covariates:
    input:
        prerecal=rules.base_recalibrator.output.prerecaltable,
        postrecal=rules.base_recalibrator_post.output.postrecaltable,
        index=paths.genome.fa
    output:
        recalfile=paths.bqsr.recal_table,
        report=paths.bqsr.report
    benchmark:
        'benchmark/{sample}_analyze_covariates.tab'
    log:
        'log/{sample}_analyze_covariates.log'
    conda:
        "../envs/analyze_covariates_gatk.yaml"
    shell:
        '''
          echo "gatk AnalyzeCovariates \
          --after-report-file {input.postrecal} \
          --before-report-file {input.prerecal} \
          --plots-report-file {output.report} \
          --intermediate-csv-file {output.recalfile}" | tee {log}

          gatk AnalyzeCovariates \
          --after-report-file {input.postrecal} \
          --before-report-file {input.prerecal} \
          --plots-report-file {output.report} \
          --intermediate-csv-file {output.recalfile} 2>> {log}
        '''

## Extract reads from chromosome 6 for HLA typing
rule extract_chr6:
    input:
        bam=rules.apply_bqsr.output.bam,
        bai=rules.apply_bqsr.output.bai
    output:
        bam=paths.bam.filtered_chr6_bam,
        bai=paths.bam.filtered_chr6_bam_index
    benchmark:
        'benchmark/{sample}_extract_chr6.tab'
    log:
        'log/{sample}_extract_chr6.log'
    conda:
        "../envs/samtools.yaml"
    threads: max(1,min(16,NCORES))
    shell:
        '''
          echo "samtools view -@ {threads} {input.bam} chr6 -b -o {output.bam} && samtools index -@ {threads} {output.bam}" | tee {log}
          samtools view -@ {threads} {input.bam} chr6 -b -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}
        '''

## Convert chromosome 6 bam to paired-end fastqs for HLA typing input
rule make_chr6_fastqs:
    input:
        bam=rules.extract_chr6.output.bam,
        bai=rules.extract_chr6.output.bai
    output:
        r1=paths.bqsr.recal_chr6_fq_r1,
        r2=paths.bqsr.recal_chr6_fq_r2
    benchmark:
        'benchmark/{sample}_make_chr6_fastqs.tab'
    log:
        'log/{sample}_make_chr6_fastqs.log'
    conda:
        "../envs/samtools.yaml"
    threads: max(1,min(16,NCORES))
    shell:
        '''
          echo "samtools fastq -@ {threads} {input.bam} -1 {output.r1} -2 {output.r2} -n" | tee {log}
          samtools fastq -@ {threads} {input.bam} -1 {output.r1} -2 {output.r2} -n 2>> {log}
        '''
