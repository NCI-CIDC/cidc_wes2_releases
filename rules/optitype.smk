## Runs OptiType to predict the HLA genotype
rule optitype:
    input:
        tch=rules.optitype_config.output.tch,
        fq1=rules.make_chr6_fastqs.output.r1,
        fq2=rules.make_chr6_fastqs.output.r2
    output:
        tsv=paths.optitype.tsv,
        pdf=paths.optitype.pdf
    benchmark:
        'benchmark/{sample}_optitype.tab'
    log:
        'log/{sample}_optitype.log'
    conda:
        "../envs/optitype.yaml"
    params:
        predir=PREDIR+"/optitype"
    threads: 16 
    shell:
        '''
          echo "OptiTypePipeline.py -i {input.fq1} {input.fq2} --dna -v -o {params.predir} -p {wildcards.sample}" | tee {log}
          OptiTypePipeline.py -i {input.fq1} {input.fq2} --dna -v -o {params.predir} -p {wildcards.sample} 2>> {log}
        '''
