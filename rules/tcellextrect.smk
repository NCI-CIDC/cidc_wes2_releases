## Calculates T cell fractions from WES data from hg38 aligned genomes using TcellExTRECT
rule tcellextrect:
    input:
        done=rules.install_tcellextrect.output.done,
        bed=rules.retrieve_tcellextrect_bed.output.bed,
        dedup=rules.dedup_bam.output.dedup,
        bai=rules.dedup_bam.output.bai
    output:
        txt=paths.tcellextrect.txt,
        pdf=paths.tcellextrect.pdf
    benchmark:
        'benchmark/{sample}_tcellextrect.tab'
    log:
        'log/{sample}_tcellextrect.log'
    params:
        r=Path(SOURCEDIR) / "r" / "tcellextrect.r",
        output_dir=PREDIR+"/tcellextrect",
        sample="{sample}"
    conda:
        "../envs/tcellextrect.yaml"
    shell:
        '''
          echo "Rscript --vanilla {params.r} {input.dedup} {input.bed} {params.output_dir} {params.sample}" | tee {log}
          Rscript --vanilla {params.r} {input.dedup} {input.bed} {params.output_dir} {params.sample} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/tcellextrect_rule.info
        '''
