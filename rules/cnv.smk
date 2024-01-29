## Run CNV analysis with QDNAseq
rule cnv_analysis:
   input:
       bam=rules.Base_recalibration_precal_sentieon.output.recalibratedbam,
       idx=rules.Base_recalibration_precal_sentieon.output.recalibratedbai
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
       predir=PREDIR,
       annot_gtf=paths.annot.gtf
   priority: 1
   threads: 1
   shell:
       '''
         echo "Rscript --vanilla {params.srcdir}/r/init-qdnaseq-cnv-analysis.r {params.predir} {params.sample} {input.bam} {params.annot_gtf}" | tee {log}
         Rscript --vanilla {params.srcdir}/r/init-qdnaseq-cnv-analysis.r {params.predir} {params.sample} {input.bam} {params.annot_gtf} 2>> {log}
       '''
