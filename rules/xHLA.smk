rule optitype_extract_chr6:
    """Extract chr6"""
    input:
        in_bam = paths.bam.filtered_bam,
    output:
        bam = paths.bam.filtered_chr6_bam,
        bai = paths.bam.filtered_chr6_bam_index
    threads: 16
    group: "optitype"
    conda: "../envs/samtools.yaml"
    log:
        "log/{sample}_extract_chr6.log"
    benchmark:
        "benchmark/{sample}_extract_chr6.tab"
    shell:
        """samtools view -@ {threads} {input.in_bam} chr6 -b -o {output.bam} &&
           samtools index -@ {threads} {output.bam}
        """

rule xhla:
    """calculate hlatyping by xhla"""
    input:
        bam = paths.bam.filtered_chr6_bam,
        bai = paths.bam.filtered_chr6_bam_index,
	hla_bed = paths.genome.hla_bed,
	hla_tsv = paths.genome.hla_tsv,
	hla_fna = paths.genome.hla_fna,	
	
    output:
        paths.xhla.report
    threads: 4 #_xhla_threads
    group: "xhla"
    params:
        sample = "{sample}",
        output_dir=Path(paths.xhla.report).parent,
        srcdir=SOURCEDIR,
        ref_data = Path(paths.genome.hla_fna).parent
    conda: "../envs/xHLA.yaml"
    benchmark:
        "benchmarks/{sample}_xhla.tab"
    log:
        "log/{sample}_getfile.log"
    shell:
        """
        {params.srcdir}/xHLA/run.py --sample_id {params.sample} --input_bam_path {input.bam} --output_path  {params.output_dir} --exec_path "{params.srcdir}/xHLA/" --temp_path {params.output_dir} --ref_data {params.ref_data}
        """