## Run CNV analysis with QDNAseq
rule cnv_analysis:
   input:
       bam=rules.apply_bqsr.output.bam,
       idx=rules.apply_bqsr.output.bai,
       gtf=rules.retrieve_reference_genome.output.gtf
   output:
       bed=paths.cnv.bed,
       igv=paths.cnv.igv,
       csv=paths.cnv.csv
   benchmark:
       'benchmark/{sample}_cnv_analysis.tab'
   log:
       'log/{sample}_cnv_analysis.log'
   conda:
       SOURCEDIR+"/../envs/cnv.yaml"
   params:
       sample='{sample}',
       srcdir=SOURCEDIR,
       predir=PREDIR
   priority: 1
   threads: 1
   shell:
       '''
         echo "Rscript --vanilla {params.srcdir}/r/init-qdnaseq-cnv-analysis.r {params.predir} {params.sample} {input.bam} {input.gtf}" | tee {log}
         Rscript --vanilla {params.srcdir}/r/init-qdnaseq-cnv-analysis.r {params.predir} {params.sample} {input.bam} {input.gtf} 2>> {log}
       '''
