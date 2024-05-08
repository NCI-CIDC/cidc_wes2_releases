#this rule will create vcf calls using mutect2 for all bams as written,
#  but should ideally only run on normal bams, as it feeds in to the creation
#  of a panel of normal file.
rule create_vcf_of_normals:
    input:
        ref = paths.genome.fa,
        bam = lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam",
#        bam = paths.bqsr.recal_bam
    output:
        pon = paths.mutect2.normal_vcf,
        #pon =lambda wildcards: Path(PREDIR) / "mutect2" / f"{tumor_normal_df.at[wildcards.sample,'normal']}.pon.vcf.gz", #
    params:
        max_mnp_distance = config["mutect2"]["max_mnp_distance"]
    threads: 8
    conda: "../envs/gatk.yaml"
    log: "log/{sample}_vcf_of_normals.log"
    benchmark: "benchmark/{sample}_vcf_of_normals.tab"
    shell: "gatk Mutect2 -R {input.ref} -I {input.bam} -O {output.pon} --max-mnp-distance {params.max_mnp_distance}"


rule copynumber_create_pon:
    input:
        vcf = paths.mutect2.normal_vcf
    output:
        ponfile=paths.mutect2.pon
    threads: 8
    conda: "../envs/gatk.yaml"
    benchmark:
        "benchmark/{sample}_create_pon.tab"
    log:
        "log/{sample}_create_pon.log"
    shell:
         "gatk CreateSomaticPanelOfNormals --variant {input.vcf} --output {output.ponfile}"

#        """{params.index1}/sentieon driver  -t {threads} -r {params.index} -i {input.normal_recalibratedbam} --algo CNV  --target {input.targetbed} --target_padding 0 --create_pon {output.ponfile}"""
	


rule mutect2:
    input:
        pon =paths.mutect2.normal_vcf,
        ref= paths.genome.fa,
        tumor=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bam",
        tumor_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bam",
        normal=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam",
        normal_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bai"
    output:
        vcf=paths.mutect2.somatic_calls_vcf
    benchmark:
        'benchmark/{sample}_mutect2.tab'
    log:
        'log/{sample}_mutect2.log'
    conda:
        "../envs/gatk.yaml"
    params:
        tumor_name = lambda wildcards: f"{tumor_normal_df.at[wildcards.sample,'tumor']}",
        normal_name = lambda wildcards: f"{tumor_normal_df.at[wildcards.sample,'normal']}"
    shell:
        'gatk Mutect2 '
        '-R {input.ref}  '
        '-I {input.normal} -normal {params.normal_name}  '
        '-I {input.tumor} -tumor {params.tumor_name} '
        '-pon {input.pon}  '
        '-O {output.vcf}'

#this is a very simple filtering. An example of a more detailed filter
#  which generates contamination tables for tumor and normal samples
#  can be found here:
#
#      https://www.biostars.org/p/9591924/
#
#  we can implement these additional steps if the simple filter option
#  seems insufficient
#
rule filter_mutect2:
    input:
       vcf = paths.mutect2.normal_vcf,
       ref = paths.genome.fa
    output:
       vcf = paths.mutect2.filtered_somatic_calls_vcf,
    conda:
       "../envs/gatk.yaml"
#    params:
    log:
       "log/{sample}_mutect2_filter.log"
    benchmark:
       "benchmark/{sample}_mutect2_filter.tab"
    shell:
       "gatk FilterMutectCalls -V {input.vcf} -R {input.ref} -O {output}"
