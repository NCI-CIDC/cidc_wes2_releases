rule Indel_realigner_sentieon:
    """indel realigner for uniquely mapped reads"""
    input:
        bam=paths.bam.filtered_bam,
        bai=paths.bam.filtered_index
    output:
        realignbam=paths.bam.realigned_bam,
        realignbai=paths.bam.realigned_index
    message: "INDEL REALIGNER: indel realigner for mapped reads"
    params:
        index=paths.genome.fa,
        sentieon_path=config['sentieon_path'],
        mills=paths.genome.mills,
        g1000=paths.genome.g1000,
    group: "recalibration"
    threads: 8 #_realigner_threads
    benchmark:
        'benchmark/{sample}_Indel_realigner_sentieon.tab'
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.bam} --algo Realigner -k {params.mills} -k {params.g1000} {output.realignbam}"""
