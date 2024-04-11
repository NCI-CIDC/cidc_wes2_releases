## Adapts the PyClone output genereated by the Sequenza function sequenza2PyClone for use with PyClone-VI 
rule sequenza2pyclone6:
    input:
        tsv=rules.sequenza.output.tsv
    output:
        tsv=paths.pyclone6.tsv
    benchmark:
        'benchmark/{sample}_sequenza2pyclone6.tab'
    log:
        'log/{sample}_sequenza2pyclone6.log'
    params:
        py=Path(SOURCEDIR) / "python" / "sequenza2pyclone6.py",
        sample="{sample}"
    shell:
        '''
          echo "{params.py} -f {input.tsv} -n {params.sample} -o {output.tsv} 2" | tee {log}
          {params.py} -f {input.tsv} -n {params.sample} -o {output.tsv} 2>> {log}
        '''

## Performs the inference step with PyClone-VI
rule pyclone6_fit:
    input:
        tsv=rules.sequenza2pyclone6.output.tsv
    output:
        h5=paths.pyclone6.h5
    benchmark:
        'benchmark/{sample}_pyclone6_fit.tab'
    log:
        'log/{sample}_pyclone6_fit.log'
    conda:
        "../envs/pyclone6.yaml"
    params:
        num_clusters=40,
        density="beta-binomial",
        num_restarts=10
    shell:
        '''
          echo "pyclone-vi fit -i {input.tsv} -o {output.h5} -c {params.num_clusters} -d {params.density} -r {params.num_restarts}" | tee {log}
          pyclone-vi fit -i {input.tsv} -o {output.h5} -c {params.num_clusters} -d {params.density} -r {params.num_restarts} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/pyclone6.info
        '''

## Selects the best solution found by the PyClone-VI fit command and post-process the results
rule pyclone6_write_results_file:
    input:
        h5=rules.pyclone6_fit.output.h5
    output:
        results=paths.pyclone6.results
    benchmark:
        'benchmark/{sample}_pyclone6_write_results_file.tab'
    log:
        'log/{sample}_pyclone6_write_results_file.log'
    conda:
        "../envs/pyclone6.yaml"
    shell:
        '''
          echo "pyclone-vi write-results-file -i {input.h5} -o {output.results}" | tee {log}
          pyclone-vi write-results-file -i {input.h5} -o {output.results} 2>> {log}
        '''

## Summarizes the results from PyClone-VI
rule pyclone6_summarize_results:
    input:
        results=rules.pyclone6_write_results_file.output.results
    output:
        summary=paths.pyclone6.summary
    benchmark:
        'benchmark/{sample}_pyclone6_summarize_results.tab'
    log:
        'log/{sample}_pyclone6_summarize_results.log'
    shell:
        '''
          echo "cut -f 3-5 {input.results} | uniq > {output.summary}" | tee {log}
          cut -f 3-5 {input.results} | uniq > {output.summary} 2>> {log}
        '''
