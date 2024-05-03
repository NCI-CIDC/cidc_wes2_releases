#this rule will create vcf calls using mutect2 for all bams as written,
#  but should ideally only run on normal bams, as it feeds in to the creation
#  of a panel of normal file.
rule create_vcf_of_normals:
    input:
        ref = paths.genome.fa,
        bam = paths.bqsr.recal_bam
    output:
        pon = paths.mutect2.normal_vcf,
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
         "gatk CreateSomaticPanelOfNormals -vcfs {input.vcf} -O {output.ponfile}"

#        """{params.index1}/sentieon driver  -t {threads} -r {params.index} -i {input.normal_recalibratedbam} --algo CNV  --target {input.targetbed} --target_padding 0 --create_pon {output.ponfile}"""
	


rule mutect2:
    input:
        vcf=paths.annot.af_vcf,
        tbi=paths.annot.af_index,
        tumor=rules.apply_bqsr.output.bam,
        tumor_bai=rules.apply_bqsr.output.bai,
        normal=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam",
        normal_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bai"
    output:
        txt=paths.mutect2.calls_vcf
    benchmark:
        'benchmark/{sample}_snp_pileup.tab'
    log:
        'log/{sample}_snp_pileup.log'
    conda:
        "../envs/gatk.yaml"
    shell:
        '''
          echo "snp-pileup -q15 -Q20 {input.vcf} {output.txt} {input.normal} {input.tumor}" | tee {log}
          snp-pileup -q15 -Q20 {input.vcf} {output.txt} {input.normal} {input.tumor} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/facets.info
        '''