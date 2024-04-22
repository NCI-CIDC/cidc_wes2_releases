rule annotate_variants:
    input:
        calls="variants.bcf",  # .vcf, .vcf.gz or .bcf
        cache=directory(Path(PREDIR) / "resources/vep/cache") 
        plugins=directory(Path(PREDIR) / "resources/vep/pugins")
        # optionally add reference genome fasta
        # fasta="genome.fasta",
        # fai="genome.fasta.fai", # fasta index
        # gff="annotation.gff",
        # csi="annotation.gff.csi", # tabix index
        # add mandatory aux-files required by some plugins if not present in the VEP plugin directory specified above.
        # aux files must be defined as following: "<plugin> = /path/to/file" where plugin must be in lowercase
        # revel = path/to/revel_scores.tsv.gz
    output:
        calls="variants.annotated.bcf",  # .vcf, .vcf.gz or .bcf
        stats="variants.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool"],
        extra="--everything",  # optional: extra arguments
    log:
        "logs/vep/annotate.log",
    threads: 4
    wrapper:
        "v3.7.0/bio/vep/annotate"


from_wes1 = """
rule neoantigen_vep_annotate:
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.twist.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.output.twist.neoantigen.vep.vcf"
    params:
        index=config['genome_fasta'],
        vep_data=config['vep_data'],
        vep_plugins=config['vep_plugins'],
	vcf_bin_path="%s/bin/" % config['vcf_root'],

        #normal = lambda wildcards: config['runs'][wildcards.run][0],
        #tumor = lambda wildcards: config['runs'][wildcards.run][1],
    group: "neoantigen"
    benchmark:
        "benchmarks/neoantigen/{run}/{run}_{caller}.neoantigen_vep_annotate.txt"
    shell:
        "{params.vcf_bin_path}vep --input_file {input} --output_file {output} --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta {params.index} --offline --cache --dir_cache {params.vep_data} --plugin Frameshift --plugin Wildtype --dir_plugins {params.vep_plugins} --pick --transcript_version" #MP 2023-02-01: this should not be able to run without a conda env unless vep is installed in the user's environment I think. elsewhere, vep is brought in by a conda env, but for some reason not here.
"""