rule calculate_tmb:
    input:
        maf=rules.vcf2maf.output.maf, 
        twist_regions=rules.retrieve_tcellextrect_bed.output.bed
    output:
        tmb=paths.tmb.tmb_Results,
    benchmark:
        'benchmark/{sample}_calculate_tmb.tab'
    log:
        'log/{sample}_calculate_tmb.log'
    conda:
        "../envs/calctmb.yaml"
    script:
        "../source/r/calculate_tmb.R"
