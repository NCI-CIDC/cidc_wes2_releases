##this rule does not function. Illumina has not updated manta so that it can use python3
#  and python3 is not supported any more with conda. Implementing this would be a mess
#  and is considered optional anyway.
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
#        candidates = paths.manta.vcf
    output:
        vcf = paths.strelka.vcf,
	vcf_idx = paths.strelka.tbi,
	runstats_xml = paths.strelka.runstats_xml,
        runstats_tsv = paths.strelka.runstats_tsv,
    params:
        sample = "{sample}",
	predir = PREDIR,
	sourcedir = SOURCEDIR,
    log:
        "log/{sample}_strelka.log"
    benchmark:
        "benchmark/{sample}_strelka.tab"
    threads: 55
    conda:
        "../envs/strelka.yaml"
    shell:
        "rm -rf {params.predir}/strelka/{params.sample}/ ; " #strelka complains if this folder exists already
        "configureStrelkaSomaticWorkflow.py "
        "--normalBam {input.normal} "
        "--tumorBam {input.tumor} "
        "--referenceFasta {input.fa} "
       #"--indelCandidates {input.candidates} " #this optional line relies on manta, which cannot be run bc of a python2 dependency
        "--runDir {params.predir}/strelka/{params.sample}/ && "
        "{params.predir}/strelka/{params.sample}/runWorkflow.py -m local -j {threads} && "
        "cp {params.predir}/strelka/{params.sample}/results/stats/runStats.xml {output.runstats_xml}; "
        "cp {params.predir}/strelka/{params.sample}/results/stats/runStats.tsv {output.runstats_tsv}; "
        "bcftools concat -a {params.predir}/strelka/{params.sample}/results/variants/*vcf.gz -o {output.vcf} --write-index=tbi "

