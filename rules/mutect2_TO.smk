rule mutect2_TO:
    input:
        pon=rules.retrieve_1kg_pon_file.output.vcf,
        tbi=rules.retrieve_1kg_pon_file.output.tbi,
        ref=paths.genome.fa,
        tumor=paths.bqsr.recal_bam
    output:
        vcf=paths.mutect2_TO.vcf,
        tbi=paths.mutect2_TO.tbi
    benchmark:
        'benchmark/{sample}_mutect2_TO.tab'
    log:
        'log/{sample}_mutect2_TO.log'
    conda:
        "../envs/gatk.yaml"
    params:
        max_mnp_distance = config["mutect2"]["max_mnp_distance"],
        tumor_name='{sample}'
    shell:
        '''
          echo "gatk Mutect2 \
          -R {input.ref} \
          -I {input.tumor} -tumor {params.tumor_name} \
          -pon {input.pon} \
          -O {output.vcf} \
          --max-mnp-distance {params.max_mnp_distance}" | tee {log}

          gatk Mutect2 \
          -R {input.ref} \
          -I {input.tumor} -tumor {params.tumor_name} \
          -pon {input.pon} \
          -O {output.vcf} \
          --max-mnp-distance {params.max_mnp_distance} 2>> {log}
        '''

rule filter_mutect2_TO:
    input:
        vcf=paths.mutect2_TO.vcf,
        tbi=paths.mutect2_TO.tbi,
        ref=paths.genome.fa
    output:
        vcf=paths.mutect2_TO.filtered_vcf,
        tbi=paths.mutect2_TO.filtered_tbi
    benchmark:
        'benchmark/{sample}_filter_mutect2_TO.tab'
    log:
        'log/{sample}_filter_mutect2_TO.log'
    conda:
        "../envs/gatk.yaml"
    shell:
        '''
          echo "gatk FilterMutectCalls -V {input.vcf} -R {input.ref} -O {output.vcf}" | tee {log}
          gatk FilterMutectCalls -V {input.vcf} -R {input.ref} -O {output.vcf} 2>> {log}
        '''
