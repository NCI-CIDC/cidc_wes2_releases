rule mutect2_TO:
    input:
        pon =paths.annot.kg_pon,
        ref= paths.genome.fa,
        tumor=paths.bqsr.recal_bam,
    output:
        vcf=paths.mutect2_TO.vcf
    benchmark:
        'benchmark/{sample}_mutect2_TO.tab'
    log:
        'log/{sample}_mutect2_TO.log'
    conda:
        "../envs/gatk.yaml"
    params:
        max_mnp_distance = config["mutect2"]["max_mnp_distance"]
    shell:
        'gatk Mutect2 '
        '-R {input.ref}  '
        '-I {input.tumor} -tumor {params.tumor_name} '
        '-pon {input.pon}  '
        '-O {output.vcf}'
        '--max-mnp-distance {params.max_mnp_distance} '

rule filter_mutect2_TO:
    input:
       vcf = paths.mutect2_TO.vcf,
       ref = paths.genome.fa
    output:
       vcf = paths.mutect2_TO.filtered_vcf,
    conda:
       "../envs/gatk.yaml"
    log:
       "log/{sample}_mutect2_filter.log"
    benchmark:
       "benchmark/{sample}_mutect2_filter.tab"
    shell:
       "gatk FilterMutectCalls -V {input.vcf} -R {input.ref} -O {output}"
