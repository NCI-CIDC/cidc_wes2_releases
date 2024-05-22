rule setup_lancet:
    input:
        paths.conda.setup_lancet
    output:
        paths.git.lancet
    conda:
        "../envs/lancet_build.yaml"
    params:
        lancet_url = config["lancet"] ,
	path = PREDIR
    shell:
        """
        cd {params.path}/git
        gh repo clone {params.lancet_url}
        cd lancet/
        sudo apt-get update && sudo apt-get install libssl-dev libcurl3-dev liblzma-dev libbz2-dev
        make
        touch {output}
        """


rule lancet:
    input:
        genome_size = paths.genome.size,
        fa = paths.genome.fa,
	make_check = papths.git.setup_lancet,
        tumor=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bam",
        tumor_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bam",
        normal=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam",
        normal_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bai"
    output:
        vcf = paths.lancet.vcf
    params:
        sample = "{sample}",
	predir = PREDIR,
    conda: "../envs/lancet.yaml"
    shell:
        "for chrom in cut -f 1 {input.genome.size} | paste -sd ' '; "
        "do lancet --tumor {input.tumor} --normal {input.normal} --ref {input.fa} --reg $chrom --num-threads 8 | bgzip > {params.predir}/lancet/{params.sample}/${chrom}.vcf.gz; "
	"done"
	"bcftools concat {params.predir}/lancet/{params.sample}/*vcf.gz -o {output}"


    