
rule xhla:
    """calculate hlatyping by xhla"""
    input:
        bam = paths.bam.filtered_chr6_bam,
        bai = paths.bam.filtered_chr6_bam_index,
	hla_bed = paths.genome.hla_bed,
	hla_tsv = paths.genome.hla_tsv,
	hla_fna = paths.genome.hla_fna,	
	hla_faa = paths.genome.hla_faa,
	hla_shift = paths.genome.hla_shift,
#	hla_freq = paths.genome.hla_freq,	
    output:
        paths.xhla.report
    threads: 4 #_xhla_threads
    group: "xhla"
    params:
        sample = "{sample}",
        output_dir= Path(PREDIR) / "xhla",
        srcdir=SOURCEDIR,
        ref_data = Path(PREDIR) / "genome"
    conda: "../envs/xHLA.yaml"
    benchmark:
        "benchmark/{sample}_xhla.tab"
    log:
        "log/{sample}_xHLA.log"
    shell:
        """
	echo  {params.srcdir}/xHLA/run.py --sample_id {params.sample} --input_bam_path {input.bam} --output_path  {params.output_dir} --exec_path "{params.srcdir}/xHLA/" --temp_path {params.output_dir}-{params.sample} --ref_data {params.ref_data} > {log};

        {params.srcdir}/xHLA/run.py --sample_id {params.sample} --input_bam_path {input.bam} --output_path  {params.output_dir} --exec_path "{params.srcdir}/xHLA/" --temp_path {params.output_dir}-{params.sample} --ref_data {params.ref_data}

        """
