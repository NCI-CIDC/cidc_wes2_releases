rule hlahd:
    """calculate hlatyping by hla-hd"""
    input:
        chr6fastqfile1 = paths.bqsr.recal_chr6_fq_r1,
        chr6fastqfile2 = paths.bqsr.recal_chr6_fq_r2,
	split_file = paths.hlahd_references.split,
	dict_file = paths.hlahd_references.dict_done,
	freq_file = paths.hlahd_references.freq_done
    output:
        paths.hlahd.report
    threads: 4
    group: "hlahd"
    params:
        sample = "{sample}",
        PREDIR = PREDIR,
        hlahd_path = config["hlahd"]["path"],
        output_dir= Path(paths.hlahd.report).parent,
        min_length = config["hlahd"]["min_length"],
        cut_percent = config["hlahd"]["cut_percent"],
        freq_data = Path(paths.hlahd_references.split).parent / "freq_data",
        split_file = Path(paths.hlahd_references.split),
        dictionary = Path(paths.hlahd_references.split).parent / "dictionary"
    log: "log/{sample}_run_hlahd.log"
    conda:
        "../envs/hlahd.yaml"
    benchmark:
        "benchmark/{sample}_run_hlahd.tab"
    shell:
        """{params.PREDIR}/{params.hlahd_path}/hlahd.sh -m {params.min_length} -c {params.cut_percent} -t {threads} -f {params.freq_data} {input.chr6fastqfile1} {input.chr6fastqfile2} {params.split_file} {params.dictionary} {params.sample} {params.output_dir} 2> {log}"""

