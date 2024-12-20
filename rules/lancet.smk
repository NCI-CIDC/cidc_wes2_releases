rule run_lancet:
    input:
        genome_size = paths.genome.size,
        fa = paths.genome.fa,
        fai = paths.genome.fai,
        tumor=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bam",
        tumor_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'tumor']}_recalibrated.bai",
        normal=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bam",
        normal_bai=lambda wildcards: Path(PREDIR) / "bqsr" / f"{tumor_normal_df.at[wildcards.sample,'normal']}_recalibrated.bai"
    output:
        vcf = paths.lancet.vcf
    params:
        sample = "{sample}",
	predir = PREDIR,
	sourcedir = SOURCEDIR
    log:
        "log/{sample}_lancet.log"
    benchmark:
        "benchmark/{sample}_lancet.tab"
    threads: 55
    conda:
        "../envs/lancet.yaml"
    shell: """
counter=0

while read -r line; do
  chrom=$(echo "$line" | cut -f1)
  echo "{params.sourcedir}/bin/lancet  --num-threads {threads} --tumor {input.tumor} --normal {input.normal} --ref {input.fa} --reg "$chrom" | bgzip | bcftools reheader --fai {input.fai} > {params.predir}/lancet/{params.sample}_"$counter".vcf.gz"| tee {log}

  {params.sourcedir}/bin/lancet  --num-threads {threads} --tumor {input.tumor} --normal {input.normal} --ref {input.fa} --reg "$chrom" | bgzip | bcftools reheader --fai {input.fai} > {params.predir}/lancet/{params.sample}_"$counter".vcf.gz 2>> {log}
  counter=$((counter + 1))
done < {input.genome_size} &&

#MP 2024-05-31
#this line threw an error bc lancet once generated a completely empty vcf (ie sans header)
# an error here would prevent the rm cleanup operation, so would preserve any evidence
# of this bug for future troubleshooting.
bcftools concat {params.predir}/lancet/{params.sample}_*vcf.gz -o {output} 2>> {log}

rm {params.predir}/lancet/{params.sample}_*.vcf.gz #clean up single chromosome files
"""


    