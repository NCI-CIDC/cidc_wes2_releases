## Generates the output files for rule bam2seqz and the input for rule mrg_seq_chr. Returns the list of files as one string
def seqz_chr(wildcards, file_type="seqz"):
    sample = wildcards
    chrs=["chr%s" % str(i) for i in range(1,23)]

    if file_type == "seqz":
        paths = " ".join([PREDIR+"/sequenza/%s/%s.seqz_%s.txt.gz" % (sample, sample, chr) for chr in chrs])
    elif file_type == "tbi":     
        paths = " ".join([PREDIR+"/sequenza/%s/%s.seqz_%s.txt.gz.tbi" % (sample, sample, chr) for chr in chrs])

    print(paths) 
    return paths

## Process a paired set of BAM/pileup files (tumor and matching normal) and GC-content genome-wide
## information to extract the common positions with A and B alleles frequencies
rule bam2seqz:
    input:
        fa=rules.retrieve_reference_genome.output.fa,
        wig=rules.retrieve_sequenza_wig.output.wig,
        tumor=rules.apply_bqsr.output.bam,
        tumor_bai=rules.apply_bqsr.output.bai,
        normal=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam",
        normal_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bai"
    output:
        tch=paths.sequenza.tch
    benchmark:
        'benchmark/{sample}_bam2seqz.tab'
    log:
        'log/{sample}_bam2seqz.log'
    conda:
        "../envs/sequenza_utils.yaml"
    threads: max(1,min(18,NCORES))
    params:
        seqz=paths.sequenza.seqz
    shell:
        '''
          echo "sequenza-utils bam2seqz -n {input.normal} -t {input.tumor} -gc {input.wig} -F {input.fa} -o {params.seqz} --parallel {threads} -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 && touch {output.tch}" | tee {log}
          sequenza-utils bam2seqz -n {input.normal} -t {input.tumor} -gc {input.wig} -F {input.fa} -o {params.seqz} --parallel {threads} -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 && touch {output.tch} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/sequenza_utils.info
        '''

## Merge all seqz outputs from each chromosome into one seqz file
rule merge_seqz_chr:
    input:
        tch=rules.bam2seqz.output.tch
    output:
        seqz_merged=paths.sequenza.seqz
    benchmark:
        'benchmark/{sample}_merge_chr.tab'
    log:
        'log/{sample}_merge_chr.log'
    conda:
        "../envs/sequenza_utils.yaml"
    threads: max(1,min(8,NCORES))
    params:
        seqz=lambda wildcards: seqz_chr(wildcards.sample, file_type="seqz")
    shell:
        '''
          echo "zcat {params.seqz} | gawk '{{if (NR!=1 && \$1 != \"chromosome\") {{print $0}}}}' | bgzip -@ {threads} > {output.seqz_merged}" | tee {log}
          zcat {params.seqz} | gawk '{{if (NR!=1 && $1 != "chromosome") {{print $0}}}}' | bgzip -@ {threads} > {output.seqz_merged} 2>> {log}
        '''

## Perform the binning of the seqz file to reduce file size and memory requirement for the analysis
rule seqz_binning:
    input:
        seqz=rules.merge_seqz_chr.output.seqz_merged
    output:
        bin50=paths.sequenza.bin50
    benchmark:
        'benchmark/{sample}_seqz_binning.tab'
    log:
        'log/{sample}_seqz_binning.log'
    conda:
        "../envs/sequenza_utils.yaml"
    shell:
        '''
          echo "sequenza-utils seqz_binning --seqz {input.seqz} --window 50 -o {output.bin50}" | tee {log}
          sequenza-utils seqz_binning --seqz {input.seqz} --window 50 -o {output.bin50} 2>> {log}
        '''

## Add the header to the final seqz file
rule seqz_header:
    input:
        bin50=rules.seqz_binning.output.bin50
    output:
        final=paths.sequenza.final
    benchmark:
        'benchmark/{sample}_seqz_header.tab'
    log:
        'log/{sample}_seqz_header.log'
    conda:
        "../envs/sequenza_utils.yaml"
    shell:
        '''
          echo "gunzip -c {input.bin50} | (echo "chromosome\tposition\tbase.ref\tdepth.normal\tdepth.tumor\tdepth.ratio\tAf\tBf\tzygosity.normal\tGC.percent\tgood.reads\tAB.normal\tAB.tumor\ttumor.strand"; cat) | gzip > {output.final}" | tee {log}
          gunzip -c {input.bin50} | (echo "chromosome\tposition\tbase.ref\tdepth.normal\tdepth.tumor\tdepth.ratio\tAf\tBf\tzygosity.normal\tGC.percent\tgood.reads\tAB.normal\tAB.tumor\ttumor.strand"; cat) | gzip > {output.final} 2>> {log}
        '''

## Runs Sequenza to extract the relevant information from the raw seqz file, fit the sequenza model to infer cellularity and ploidy,
## and apply the inferred parameters to estimate the copy number profile
rule sequenza:
    input:
        seqz=rules.seqz_header.output.final
    output:
        tsv=paths.sequenza.tsv,
        pdf=paths.sequenza.pdf,
        segments=paths.sequenza.segments,
        CP_contours=paths.sequenza.CP_contours,
        alt=paths.sequenza.alt,
        chr=paths.sequenza.chr
    benchmark:
        'benchmark/{sample}_sequenza.tab'
    log:
        'log/{sample}_sequenza.log'
    conda:
        "../envs/sequenza.yaml"
    params:
        r=Path(SOURCEDIR) / "r" / "sequenza.r",
        dir=PREDIR+"/sequenza/{sample}",
        sample="{sample}"
    shell:
        '''
          echo "Rscript {params.r} {input.seqz} {params.dir} {params.sample}" | tee {log}
          Rscript {params.r} {input.seqz} {params.dir} {params.sample} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/sequenza.info
        '''
