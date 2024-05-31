rule run_manta:
    input:
        fa = paths.genome.fa,
        fai = paths.genome.fai,
        tumor=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bam",
        tumor_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bai",
        normal=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam",
        normal_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bai"    
    output:
        paths.manta.vcf
    params:
        predir = PREDIR,
    log: "log/{sample}_manta.log"
    benchmark: "benchmark/{sample}_manta.tab"
    conda: "../envs/manta.yaml"
    threads: 55
    shell:
        """
configManta.py \
--normalBam {input.normal} \
--tumorBam {input.tumor} \
--referenceFasta {input.fa} \
--runDir {params.predir}/manta_temp/
"""

rule config_strelka:
    input:
        fa = paths.genome.fa,
        fai = paths.genome.fai,
        tumor=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bam",
        tumor_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bai",
        normal=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam",
        normal_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bai",
        candidates = paths.manta.vcf
    output:
        vcf = paths.strelka.vcf
    params:
        sample = "{sample}",
	predir = PREDIR,
	sourcedir = SOURCEDIR
    log:
        "log/{sample}_strelka.log"
    benchmark:
        "benchmark/{sample}_strelka.tab"
    threads: 55
    conda:
        "../envs/strelka.yaml"
    shell: """

configureStrelkaSomaticWorkflow.py \
--normalBam {input.normal} \
--tumorBam {input.tumor} \
--referenceFasta {input.fa} \
--indelCandidates {input.candidates} \
--runDir {params.predir}/strelka_temp
"""
