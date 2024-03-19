## Map reads to the reference genome using BWA and output sorted bam
rule exome_coverage:
    input:
        genome_size=rules.genome_size.output.size,
        bam=rules.dedup_bam.output.dedup,
        bai=rules.dedup_bam.output.bai,
        targets_bed=rules.retrieve_coverage_targets_bed.output
    output:
        depth=paths.coverage.depth,
        bdg=paths.coverage.bdg,
        bw=paths.coverage.bw
    benchmark:
        'benchmark/{sample}_exome_coverage.tab'
    log:
        'log/{sample}_exome_coverage.log'
    conda:
        SOURCEDIR+"/../envs/filter_bam.yaml"
    params:
        sample='{sample}',
        min_mapping_quality=20,
        min_base_quality=20
    priority: 4
    threads: 1
    shell:
        '''
          # Perform coverage analysis
          echo "samtools depth -b {input.targets_bed} -Q {params.min_mapping_quality} -q {params.min_base_quality} -a {input.bam} > {output.depth}" | tee {log}
          samtools depth -b {input.targets_bed} -Q {params.min_mapping_quality} -q {params.min_base_quality} -a {input.bam} > {output.depth} 2>> {log}
          
          # Generate coverage bedgraph file
          echo "bedtools genomecov -ibam {input.bam} -bg | sort -k1,1 -k2,2n > {output.bdg}" | tee -a {log}
          bedtools genomecov -ibam {input.bam} -bg | sort -k1,1 -k2,2n > {output.bdg} 2>> {log}
          
          # Generate coverage bw file for visualizations
          echo "bedGraphToBigWig {output.bdg} {input.genome_size} {output.bw}" | tee -a {log}
          bedGraphToBigWig {output.bdg} {input.genome_size} {output.bw} 2>> {log}
        '''
