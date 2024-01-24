rule Indel_realigner_sentieon:
    """indel realigner for uniquely mapped reads"""
    input:
        bam=paths.bam.filtered_bam,
        bai=paths.bam.filtered_index,
        mills=paths.genome.mills,
        g1000=paths.genome.g1000,
    output:
        realignbam=paths.bam.realigned_bam,
        realignbai=paths.bam.realigned_index
    message: "INDEL REALIGNER: indel realigner for mapped reads"
    params:
        index=paths.genome.fa,
        sentieon_path=config['sentieon_path'],
    group: "recalibration"
    threads: 8 #_realigner_threads
    benchmark:
        'benchmark/{sample}_Indel_realigner_sentieon.tab'
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.bam} --algo Realigner -k {input.mills} -k {input.g1000} {output.realignbam}"""

##########################################
###basecall quality score recalibration###
##########################################

rule Base_recalibration_precal_sentieon:
    """base recalibration for realigned files"""
    input:
        realignbam=paths.bam.realigned_bam,
        index = paths.genome.fa,
        dbsnp= paths.genome.dbsnp,
        mills= paths.genome.mills,
        g1000= paths.genome.g1000,	
    output:
        prerecaltable=paths.bqsr.prerecaltable,
        recalibratedbam=paths.bqsr.recal_bam,
        recalibratedbai=paths.bqsr.recal_index
    message:
        " PRE BASE RECALIBRATION: base recalibration for  realigned files"
    params:
        sentieon_path=config['sentieon_path'],
    threads: 4 #_realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_precal_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {input.index} -t {threads} -i {input.realignbam} --algo QualCal -k {input.dbsnp} -k {input.mills} -k {input.g1000}  {output.prerecaltable} --algo ReadWriter {output.recalibratedbam}"""

rule Base_recalibration_postcal_sentieon:
    """post recalibration for realigned files"""
    input:
        recalibratedbam=paths.bqsr.recal_bam,
        recalibrated_bai=paths.bqsr.recal_bam,
        prerecaltable=paths.bqsr.prerecaltable,
        index = paths.genome.fa,
        dbsnp= paths.genome.dbsnp,
        mills= paths.genome.mills,
        g1000= paths.genome.g1000        
    output:
        postrecaltable=paths.bqsr.postrecaltable
    message:
        "POST BASE RECALIBRATION: post base recalibration for  realigned files"
    params:
        sentieon_path=config['sentieon_path'],
    threads: 16 #_realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_postcal_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {input.index} -t {threads} -i {input.recalibratedbam} -q {input.prerecaltable} --algo QualCal -k {input.dbsnp} -k {input.mills} -k {input.g1000}  {output.postrecaltable}"""

rule Base_recalibration_sentieon:
    """ recalibration for realigned files"""
    input:
        postrecaltable= paths.bqsr.postrecaltable,
        prerecaltable=paths.bqsr.prerecaltable,
        index = paths.genome.fa,	
    output:
        recalfile=paths.bqsr.recal_table
    message:
        "DIFF BASE RECALIBRATION: Difference in pre and post processing of realigned files"
    params:
        sentieon_path=config['sentieon_path'],
    threads: 1 #_realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -t {threads} --algo QualCal --plot --before {input.prerecaltable} --after {input.postrecaltable} {output.recalfile}"""

