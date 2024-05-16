## This rule will create vcf calls using mutect2 for all bams as written,
## but should ideally only run on normal bams, as it feeds in to the creation
## of a panel of normal file.
rule create_vcf_of_normals:
    input:
        ref=paths.genome.fa,
        bam=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam"
    output:
        vcf=paths.mutect2.normal_vcf,
        tbi=paths.mutect2.normal_tbi
    params:
        max_mnp_distance = config["mutect2"]["max_mnp_distance"]
    benchmark:
        'benchmark/{sample}_vcf_of_normals.tab'
    log:
        'log/{sample}_vcf_of_normals.log'
    conda:
        "../envs/gatk.yaml"
    threads: max(1,min(8,NCORES))
    params:
        max_mnp_distance=config["mutect2"]["max_mnp_distance"]
    shell:
        '''
          echo "Mutect2 -R {input.ref} -I {input.bam} -O {output.vcf} --max-mnp-distance {params.max_mnp_distance}" | tee {log}
          gatk Mutect2 -R {input.ref} -I {input.bam} -O {output.vcf} --max-mnp-distance {params.max_mnp_distance} 2>> {log}
        '''

rule copynumber_create_pon:
    input:
        vcf=paths.mutect2.normal_vcf,
        tbi=paths.mutect2.normal_tbi
    output:
        pon=paths.mutect2.pon,
        tbi=paths.mutect2.pon_idx
    benchmark:
        'benchmark/{sample}_create_pon.tab'
    log:
        'log/{sample}_create_pon.log'
    conda:
        "../envs/gatk.yaml"
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "gatk CreateSomaticPanelOfNormals --variant {input.vcf} --output {output.pon}" | tee {log}
          gatk CreateSomaticPanelOfNormals --variant {input.vcf} --output {output.pon} 2>> {log}
        '''

rule mutect2:
    input:
        pon=paths.mutect2.pon,
        pon_idx=paths.mutect2.pon_idx,
        ref=paths.genome.fa,
        tumor=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bam",
        tumor_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bam",
        normal=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam",
        normal_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bai"
    output:
        vcf=paths.mutect2.somatic_vcf,
        tbi=paths.mutect2.somatic_tbi
    benchmark:
        'benchmark/{sample}_mutect2.tab'
    log:
        'log/{sample}_mutect2.log'
    conda:
        "../envs/gatk.yaml"
    params:
        tumor_name=lambda wildcards: f"{tumor_normal_df.at[wildcards.sample,'tumor']}",
        normal_name=lambda wildcards: f"{tumor_normal_df.at[wildcards.sample,'normal']}"
    shell:
        '''
          echo "gatk Mutect2 \
          -R {input.ref} \
          -I {input.normal} -normal {params.normal_name} \
          -I {input.tumor} -tumor {params.tumor_name} \
          -pon {input.pon} \
          -O {output.vcf}" | tee {log}

          gatk Mutect2 \
          -R {input.ref} \
          -I {input.normal} -normal {params.normal_name} \
          -I {input.tumor} -tumor {params.tumor_name} \
          -pon {input.pon} \
          -O {output.vcf} 2>> {log}
        '''

##  This is a very simple filtering. An example of a more detailed filter
##  which generates contamination tables for tumor and normal samples
##  can be found here:
##
##  https://www.biostars.org/p/9591924/
##
##  We can implement these additional steps if the simple filter option
##  seems insufficient.
rule filter_mutect2:
    input:
        vcf=paths.mutect2.somatic_vcf,
        tbi=paths.mutect2.somatic_tbi,
        ref=paths.genome.fa
    output:
        filtered_vcf=paths.mutect2.filtered_somatic_vcf,
        filtered_tbi=paths.mutect2.filtered_somatic_tbi
    benchmark:
        'benchmark/{sample}_filter_mutect2.tab'
    log:
        'log/{sample}_filter_mutect2.log'
    conda:
        "../envs/gatk.yaml"
    shell:
        '''
          echo "gatk FilterMutectCalls -V {input.vcf} -R {input.ref} -O {output.filtered_vcf}" | tee {log}
          gatk FilterMutectCalls -V {input.vcf} -R {input.ref} -O {output.filtered_vcf} 2>> {log}
        '''
