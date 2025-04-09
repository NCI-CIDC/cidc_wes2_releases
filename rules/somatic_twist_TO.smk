rule somatic_twist_TO:
    ## Takes output.vcf.gz and intersects it with the twist capture regions
    input:
        mutect2_TO=rules.filter_mutect2_TO.output.vcf,
        twist_regions=rules.retrieve_tcellextrect_bed.output.bed
    output:
        twist=paths.twist_TO.filtered_vcf
    benchmark:
        'benchmark/{sample}_somatic_twist_TO.txt'
    log:
        'log/{sample}_somatic_twist_TO.log'
    conda:
        "../envs/bcftools.yaml"
    shell:
        '''
          echo "bcftools view -R {input.twist_regions} {input.mutect2_TO} | bcftools sort | bcftools view -Ov > {output.twist}" | tee {log}
          bcftools view -R {input.twist_regions} {input.mutect2_TO} | bcftools sort | bcftools view -Ov > {output.twist}
        '''
