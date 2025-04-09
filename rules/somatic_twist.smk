rule somatic_twist:
    ## Takes output.vcf.gz and intersects it with the twist capture regions
    input:
        mutect2=rules.filter_mutect2.output.filtered_vcf,
        twist_regions=rules.retrieve_tcellextrect_bed.output.bed
    output:
        twist=paths.twist.filtered_vcf
    benchmark:
        'benchmark/{sample}_somatic_twist.txt'
    log:
        'log/{sample}_somatic_twist.log'
    conda:
        "../envs/bcftools.yaml"
    shell:
        '''
          echo "bcftools view -R {input.twist_regions} {input.mutect2} | bcftools sort | bcftools view -Ov > {output.twist}" | tee {log}
          bcftools view -R {input.twist_regions} {input.mutect2} | bcftools sort | bcftools view -Ov > {output.twist}
        '''
