## Map reads to the reference genome using BWA and output sorted bam
rule run_bwa:
    input:
        tch=rules.build_bwa_index.output,
        fa=rules.qualityfilter.output
    output:
        paths.bam.bam
    benchmark:
        'benchmark/{sample}_run_bwa.tab'
    log:
        'log/{sample}_run_bwa.log'
    conda:
        SOURCEDIR+"/../envs/bwa.yaml"
    params:
        sample='{sample}',
        indexseq=paths.genome.fa,
        in_fa_str=expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])[0] + ' ' + expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])[2] if len(ENDS) == 2 else expand(paths.rqual_filter.qfilter_fastq_single, read=ENDS)[0],
        read_group= lambda wildcards: "@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA" % (wildcards.sample, wildcards.sample)
    priority: 4
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "bwa mem -t {threads} -R "{params.read_group}" {params.indexseq} {params.in_fa_str} | samtools view -@ {threads} -Sbh | samtools sort -@ {threads} > {output}" | tee {log}
          bwa mem -t {threads} -R \"{params.read_group}\" {params.indexseq} {params.in_fa_str} | samtools view -@ {threads} -Sbh | samtools sort -@ {threads} > {output} 2>> {log}
        '''

## Index BAM
rule index_bam:
    input:
        bam=rules.run_bwa.output
    output:
        paths.bam.index
    benchmark:
        'benchmark/{sample}_index_bam.tab'
    log:
        'log/{sample}_index_bam.log'
    conda:
        SOURCEDIR+"/../envs/samtools.yaml"
    params:
        sample='{sample}'
    priority: 5
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "samtools index -@ {threads} {input.bam}" > {log}
          samtools index -@ {threads} {input.bam} 2>> {log}
        '''

## Mark and remove duplicates in the sorted bam using a Spark implementation of Picard from GATK
rule dedup_bam:
    input:
        bam=rules.run_bwa.output,
        bai=rules.index_bam.output
    output:
        dedup=paths.bam.dedup_bam,
        bai=paths.bam.dedup_bai,
        metrics=paths.bam.dedup_metrics
    benchmark:
        'benchmark/{sample}_dedup_bam.tab'
    log:
        'log/{sample}_dedup_bam.log'
    conda:
        "../envs/gatk.yaml"
    params:
        tmp=Path(PREDIR) / "bam"
    resources:
        mem_mb=config["mem"]
    threads: max(1,min(16,NCORES))
    shell:
        '''
          echo "gatk MarkDuplicatesSpark -I {input.bam} -O {output.dedup} -M {output.metrics} --remove-all-duplicates true --tmp-dir {params.tmp} --conf 'spark.executor.cores={threads}'" | tee {log}
          gatk MarkDuplicatesSpark -I {input.bam} -O {output.dedup} -M {output.metrics} --remove-all-duplicates true --tmp-dir {params.tmp} --conf 'spark.executor.cores={threads}' 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/gatk.info
        '''

## Run FASTQC
rule fastqc:
    input:
        bam=rules.dedup_bam.output.dedup,
        idx=rules.dedup_bam.output.bai
    output:
        paths.fastqc.targz
    benchmark:
        'benchmark/{sample}_fastqc.tab'
    log:
        'log/{sample}_fastqc.log'
    conda:
        "../envs/fastqc.yaml"
    params:
        sample='{sample}',
        fq_base='fastqc/{sample}_dedup_fastqc',
        fq_zip='fastqc/{sample}_dedup_fastqc.zip',
        fq_html='fastqc/{sample}_dedup_fastqc.html'
    priority: 1
    threads: 1
    shell:
        '''
          echo "fastqc {input.bam} -q -o fastqc" > {log}
          fastqc {input.bam} -q -o fastqc 2>> {log}

          ## Unzip, remove zipped results, HTML duplicate, and tarball results
          unzip -qq {params.fq_zip} -d {params.fq_base} && tar -zcf {output} {params.fq_base} && rm -r {params.fq_zip} {params.fq_html} {params.fq_base}

          ## Export rule env details
          conda env export --no-builds > info/fastqc.info
        '''

## Run RSEQC bam_stat.py
rule bam_qc:
    input:
        bam=rules.dedup_bam.output.dedup,
        idx=rules.dedup_bam.output.bai
    output:
        paths.rseqc.bamqc_txt
    benchmark:
        'benchmark/{sample}_bam_qc.tab'
    log:
        'log/{sample}_bam_qc.log'
    conda:
        "../envs/rseqc.yaml"
    params:
        sample='{sample}'
    priority: 1
    threads: 1
    shell:
        '''
          echo "bam_stat.py -i {input.bam} > {output}" | tee {log}
          bam_stat.py -i {input.bam} > {output} 2>> {log}
        '''

## Run RSEQC read_gc.py
rule bam_gc:
    input:
        bam=rules.dedup_bam.output.dedup,
        idx=rules.dedup_bam.output.bai
    output:
        r=paths.rseqc.bamgc_r,
        txt=paths.rseqc.bamgc_txt
    benchmark:
        'benchmark/{sample}_bam_gc.tab'
    log:
        'log/{sample}_bam_gc.log'
    conda:
        "../envs/rseqc.yaml"
    params:
        sample='{sample}'
    priority: 1
    threads: 1
    shell:
      '''
        echo "read_GC.py -i {input.bam} -o rseqc/{params.sample}" | tee {log}
        read_GC.py -i {input.bam} -o rseqc/{params.sample} 2>> {log}

        ## R script to get txt output info
        echo "out=as.vector(summary(gc));dta = data.frame('{params.sample}',out[1],out[2],out[3],out[4],out[5],out[6]);write.table(dta,file='{output.txt}',sep="\t",row.names=F,col.names=F,quote=F);" >> {output.r}
        sed -i "s/pdf/png/g" {output.r} 
        sed -i 's/main=""/main="{params.sample}"/g' {output.r} 
        Rscript --vanilla --quiet {output.r}

        ## Export rule env details
        conda env export --no-builds > info/rseqc.info
      '''
