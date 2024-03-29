## Outputs for each SNP the counts of the reference nucleotide, alternative nucleotide, errors, and deletions. These counts can then be used in FACETS
rule snp_pileup:
    input:
        vcf=rules.retrieve_facets_vcf.output.vcf,
        tbi=rules.retrieve_facets_vcf.output.tbi,
        tumor=rules.apply_bqsr.output.bam,
        tumor_bai=rules.apply_bqsr.output.bai,
        normal=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam",
        normal_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bai"
    output:
        txt=paths.facets.txt
    benchmark:
        'benchmark/{sample}_snp_pileup.tab'
    log:
        'log/{sample}_snp_pileup.log'
    conda:
        "../envs/facets.yaml"
    shell:
        '''
          echo "snp-pileup -q15 -Q20 {input.vcf} {output.txt} {input.normal} {input.tumor}" | tee {log}
          snp-pileup -q15 -Q20 {input.vcf} {output.txt} {input.normal} {input.tumor} 2>> {log}
  
          ## Export rule env details
          conda env export --no-builds > info/facets.info
        '''

## Remove the "chr" prefix in column 1 of the results file from rule snp_pileup
rule snp_pileup_processing:
    input:
        txt=rules.snp_pileup.output.txt
    output:
        postprocessed=paths.facets.postprocessed
    benchmark:
        'benchmark/{sample}_snp_pileup_processing.tab'
    log:
        'log/{sample}_snp_pileup_processing.log'
    shell:
        '''
          echo "cat {input.txt} | sed 's/chr//g' > {output.postprocessed}" | tee {log}
          cat {input.txt} | sed 's/chr//g' > {output.postprocessed} 2>> {log}
        '''

## Run FACETS to estimate fraction and allele specific copy number from paired tumor/normal sequencing. As a result, also estimates purity and ploidy of a given tumor sample.
## Produces cncf output for futher analysis in the Copy Number module
rule facets:
    input:
        txt=rules.snp_pileup_processing.output.postprocessed
    output:
        cncf=paths.facets.cncf,
        opt=paths.facets.opt,
        iter=paths.facets.iter
    benchmark:
        'benchmark/{sample}_facets.tab'
    log:
        'log/{sample}_facets.log'
    conda:
        "../envs/facets.yaml"
    params:
        r=Path(SOURCEDIR) / "r" / "facets.r",
        dir=PREDIR+"/facets/"
    shell:
        '''
          echo " Rscript --vanilla {params.r} {input.txt} {params.dir} {wildcards.sample} " | tee {log}
          Rscript --vanilla {params.r} {input.txt} {params.dir} {wildcards.sample} 2>> {log}
        '''
