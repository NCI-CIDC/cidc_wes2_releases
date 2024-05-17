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

rule callGainLoss_cnvkit:
    input:
        enhanced_cns=paths.cnvkit.enhanced_cns
    output:
        bed=paths.copynumber.cnvkit
    benchmark:
        'benchmark/{sample}_callGainLoss_cnvkit.tab'
    log:
        'log/{sample}_callGainLoss_cnvkit.log'
    params:
        py=Path(SOURCEDIR) / "python" / "copynumber_callGainLoss.py"
    shell:
        '''
          echo "{params.py} -f {input.enhanced_cns} -o {output.bed}" | tee {log}
          {params.py} -f {input.enhanced_cns} -o {output.bed} 2>> {log}
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

## Generates the consensus calls between CNVkit and Sequenza
rule consensus_cnvkit_sequenza:
    input:
        cnvkit=rules.callGainLoss_cnvkit.output.bed,
        sequenza=rules.callGainLoss_sequenza.output.bed
    output:
        bed=paths.copynumber.cnv_seq
    benchmark:
        'benchmark/{sample}_consensus_cnvkit_sequenza.tab'
    log:
        'log/{sample}_consensus_cnvkit_sequenza.log'
    conda:
        "../envs/bedtools.yaml"
    params:
        awk_cmd="awk -v OFS=\'\\t\' \'{if ($4 == $5) {print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$6}}\'"
    shell:
        '''
          echo "intersectBed -a {input.cnvkit} -b {input.sequenza} | intersectBed -a stdin -b {input.cnvkit} -filenames -wo | cut -f 1-3,5,10,11 | {{params.awk_cmd}} > {output.bed}" | tee {log}
          intersectBed -a {input.cnvkit} -b {input.sequenza} | intersectBed -a stdin -b {input.cnvkit} -filenames -wo | cut -f 1-3,5,10,11 | {params.awk_cmd} > {output.bed} 2>> {log}
        '''

## Generates the consensus calls between CNVkit and FACETS
rule consensus_cnvkit_facets:
    input:
        cnvkit=rules.callGainLoss_cnvkit.output.bed,
        facets=rules.callGainLoss_facets.output.bed
    output:
        bed=paths.copynumber.cnv_fac
    benchmark:
        'benchmark/{sample}_consensus_cnvkit_facets.tab'
    log:
        'log/{sample}_consensus_cnvkit_facets.log'
    conda:
        "../envs/bedtools.yaml"
    params:
        awk_cmd="awk -v OFS=\'\\t\' \'{if ($4 == $5) {print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$6}}\'"
    shell:
        '''
          echo "intersectBed -a {input.cnvkit} -b {input.facets} | intersectBed -a stdin -b {input.cnvkit} -filenames -wo | cut -f 1-3,5,10,11 | {{params.awk_cmd}} > {output.bed}" | tee {log}
          intersectBed -a {input.cnvkit} -b {input.facets} | intersectBed -a stdin -b {input.cnvkit} -filenames -wo | cut -f 1-3,5,10,11 | {params.awk_cmd} > {output.bed} 2>> {log}
        '''

## Combines all the consensus BEDs into a singular file and sorts the entries
rule consensus_all:
    input:
        seq_fac=rules.consensus_sequenza_facets.output.bed,
        cnv_seq=rules.consensus_cnvkit_sequenza.output.bed,
        cnv_fac=rules.consensus_cnvkit_facets.output.bed
    output:
        consensus=paths.copynumber.consensus
    benchmark:
        'benchmark/{sample}_consensus_all.tab'
    log:
        'log/{sample}_consensus_all.log'
    shell:
        '''
          echo "cat {input.seq_fac} {input.cnv_seq} {input.cnv_fac} | sort -V -k1,1 -k2,2n > {output.consensus}" | tee {log}
          cat {input.seq_fac} {input.cnv_seq} {input.cnv_fac} | sort -V -k1,1 -k2,2n > {output.consensus} 2>> {log}
        '''

## Merges the overlapping entries and splits the gain and loss entries into separate BEDs
rule consensus_all_merge:
    input:
        consensus=rules.consensus_all.output.consensus
    output:
        gain=paths.copynumber.merged_gain,
        loss=paths.copynumber.merged_loss
    benchmark:
        'benchmark/{sample}_consensus_all_merge.tab'
    log:
        'log/{sample}_consensus_all_merge.log'
    conda:
        "../envs/bedtools.yaml"
    params:
        sh=Path(SOURCEDIR) / "shell/copynumber-merge.sh",
        prefix=PREDIR+"/copynumber/{sample}"
    shell:
        '''
          echo "{params.sh} {input.consensus} {params.prefix}" | tee {log}
          {params.sh} {input.consensus} {params.prefix} 2>> {log}
        '''
