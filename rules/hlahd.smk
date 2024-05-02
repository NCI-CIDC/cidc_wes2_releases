rule hlahd:
    """calculate hlatyping by hla-hd"""
    input:
        chr6fastqfile1 = rules.make_chr6_fastqs.output.r1,
        chr6fastqfile2 = rules.make_chr6_fastqs.output.r2,
	split_file = paths.hlahd_references.split,
	dict_file = paths.hlahd_references.dict_done,
	freq_file = paths.hlahd_references.freq_done
    output:
        touch(paths.hlahd.done)
    threads: 4
    group: "hlahd"
    params:
        sample = "{sample}",
        SOURCEDIR = SOURCEDIR,
        hlahd_path = config["hlahd"]["path"],
        output_dir= Path(paths.hlahd.done).parent,
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
        "export PATH={params.SOURCEDIR}/{params.hlahd_path}:$PATH; "
        "hlahd.sh -m {params.min_length} -c {params.cut_percent} -t {threads} -f {params.freq_data} {input.chr6fastqfile1} {input.chr6fastqfile2} {params.split_file} {params.dictionary} {params.sample} {params.output_dir} 2> {log}"


