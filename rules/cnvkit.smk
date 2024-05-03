## Run CNVkit to call copy number variations for the tumor samples 
rule cnvkit:
    input:
        cnn=rules.retrieve_cnvkit_cnn.output.cnn,
        bam=rules.apply_bqsr.output.bam,
        bai=rules.apply_bqsr.output.bai,
    output:
        cns=paths.cnvkit.cns,
        call_cns=paths.cnvkit.call_cns,
        scatter=paths.cnvkit.scatter
    benchmark:
        'benchmark/{sample}_cnvkit.tab'
    log:
        'log/{sample}_cnvkit.log'
    conda:
        "../envs/cnvkit.yaml"
    params:
        output_dir=PREDIR+"/cnvkit/{sample}"
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "cnvkit.py batch {input.bam} -r {input.cnn} -p {threads} --scatter --diagram -d {params.output_dir}" | tee {log}
          cnvkit.py batch {input.bam} -r {input.cnn} -p {threads} --scatter --diagram -d {params.output_dir} 2>> {log}
  
          ## Export rule env details
          conda env export --no-builds > info/cnvkit.info
        '''

## Run FACETS to estimate fraction and allele specific copy number from paired tumor/normal sequencing. As a result, also estimates purity and ploidy of a given tumor sample.
## Produces cncf output for futher analysis in the Copy Number module
#rule facets:
#    input:
#        txt=rules.snp_pileup_processing.output.postprocessed
#    output:
#        cncf=paths.facets.cncf,
#        opt=paths.facets.opt,
#        iter=paths.facets.iter
#    benchmark:
#        'benchmark/{sample}_facets.tab'
#    log:
#        'log/{sample}_facets.log'
#    conda:
#        "../envs/facets.yaml"
#    params:
#        r=Path(SOURCEDIR) / "r" / "facets.r",
#        dir=PREDIR+"/facets/"
#    shell:
#        '''
#          echo " Rscript --vanilla {params.r} {input.txt} {params.dir} {wildcards.sample} " | tee {log}
#          Rscript --vanilla {params.r} {input.txt} {params.dir} {wildcards.sample} 2>> {log}
#        '''
