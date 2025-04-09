rule mutect2:
    input:
        pon=rules.retrieve_1kg_pon_file.output.vcf,
        pon_idx=rules.retrieve_1kg_pon_file.output.tbi,
        ref=paths.genome.fa,
        fai=paths.genome.fai,
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
rule filter_mutect2:
    input:
        vcf=paths.mutect2.somatic_vcf,
        tbi=paths.mutect2.somatic_tbi,
        ref=paths.genome.fa,
        fai=paths.genome.fai
    output:
        unfiltered_vcf=paths.mutect2.unfiltered_somatic_vcf,
        unfiltered_tbi=paths.mutect2.unfiltered_somatic_tbi,
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
          echo "gatk FilterMutectCalls -V {input.vcf} -R {input.ref} -O {output.unfiltered_vcf}" | tee {log}
          gatk FilterMutectCalls -V {input.vcf} -R {input.ref} -O {output.unfiltered_vcf} 2>> {log}

          echo "vcftools --gzvcf {output.unfiltered_vcf} --remove-filtered-all --recode --stdout > {output.filtered_vcf}.tmp" | tee -a {log}
          vcftools --gzvcf {output.unfiltered_vcf} --remove-filtered-all --recode --stdout > {output.filtered_vcf}.tmp 2>> {log}

          # Recompress the temporary VCF file with bgzip
          echo "bgzip -c {output.filtered_vcf}.tmp > {output.filtered_vcf}" | tee -a {log}
          bgzip -c {output.filtered_vcf}.tmp > {output.filtered_vcf} 2>> {log}

          # Index the compressed VCF using tabix
          echo "tabix -f -p vcf {output.filtered_vcf}" | tee -a {log}
          tabix -f -p vcf {output.filtered_vcf} 2>> {log}
        '''
