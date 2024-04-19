## Use cutoffs to call regions of gain/loss for the {sample}_segments.txt from rule sequenza
rule callGainLoss_sequenza:
    input:
        segments=rules.sequenza.output.segments
    output:
        bed=paths.copynumber.sequenza
    benchmark:
        'benchmark/{sample}_callGainLoss_sequenza.tab'
    log:
        'log/{sample}_callGainLoss_sequenza.log'
    params:
        py=Path(SOURCEDIR) / "python" / "copynumber_callGainLoss.py"
    shell:
        '''
          echo "{params.py} -f {input.segments} -o {output.bed}" | tee {log}
          {params.py} -f {input.segments} -o {output.bed} 2>> {log}
        '''

## Use cutoffs to call regions of gain/loss for the {sample}.cncf from rule facets
rule callGainLoss_facets:
    input:
        cncf=rules.facets.output.cncf
    output:
        bed=paths.copynumber.facets
    benchmark:
        'benchmark/{sample}_callGainLoss_facets.tab'
    log:
        'log/{sample}_callGainLoss_facets.log'
    params:
        py=Path(SOURCEDIR) / "python" / "copynumber_callGainLoss.py"
    shell:
        '''
          echo "{params.py} -f {input.cncf} -o {output.bed}" | tee {log}
          {params.py} -f {input.cncf} -o {output.bed} 2>> {log}
        '''

## Generates the consensus calls between Sequenza and FACETS
rule consensus_sequenza_facets:
    input:
        sequenza=rules.callGainLoss_sequenza.output.bed,
        facets=rules.callGainLoss_facets.output.bed
    output:
        bed=paths.copynumber.seq_fac
    benchmark:
        'benchmark/{sample}_consensus_sequenza_facets.tab'
    log:
        'log/{sample}_consensus_sequenza_facets.log'
    conda:
        "../envs/bedtools.yaml"
    params:
        awk_cmd="awk -v OFS=\'\\t\' \'{if ($4 == $5) {print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$6}}\'"
    shell:
        '''
          echo "intersectBed -a {input.sequenza} -b {input.facets} | intersectBed -a stdin -b {input.facets} -filenames -wo | cut -f 1-3,5,10,11 | {{params.awk_cmd}} > {output.bed}" | tee {log}
          intersectBed -a {input.sequenza} -b {input.facets} | intersectBed -a stdin -b {input.facets} -filenames -wo | cut -f 1-3,5,10,11 | {params.awk_cmd} > {output.bed} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/bedtools.info
        '''

## Summarizes the results from PyClone-VI
#rule pyclone6_summarize_results:
#    input:
#        results=rules.pyclone6_write_results_file.output.results
#    output:
#        summary=paths.pyclone6.summary
#    benchmark:
#        'benchmark/{sample}_pyclone6_summarize_results.tab'
#    log:
#        'log/{sample}_pyclone6_summarize_results.log'
#    shell:
#        '''
#          echo "cut -f 3-5 {input.results} | uniq > {output.summary}" | tee {log}
#          cut -f 3-5 {input.results} | uniq > {output.summary} 2>> {log}
#        '''
