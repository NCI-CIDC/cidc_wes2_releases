## Calculates microsatellite instability (MSI) with MSIsensor2
##  CHECK IF YOU WANT TO ADD .TXT AT THE END OF THE MSISENSOR FILE
rule msisensor2:
    input:
        models=rules.retrieve_msisensor2_models.output.models,
        tch=rules.retrieve_msisensor2_models.output.tch,
        tumor=paths.bam.bam,
        normal=lambda wildcards: PREDIR+"/bam/{normal}.bam".format(normal=pairings_df.at[wildcards.sample, 'normal']) if pairings_df.at[wildcards.sample,'type'] == "TN" else []
    output:
        output=paths.msisensor2.output,
        dis=paths.msisensor2.dis,
        somatic=paths.msisensor2.somatic
    benchmark:
        'benchmark/{sample}_msisensor2.tab'
    log:
        'log/{sample}_msisensor2.log'
    conda:
        "../envs/msisensor2.yaml"
    threads: max(1,min(8,NCORES))
    params:
        prefix='msisensor2/{sample}_msisensor2',
        file=lambda wildcards, input: "-t %s" % input.tumor if pairings_df.at[wildcards.sample, 'type'] == "TO" else "-t %s -n %s" % (input.tumor, input.normal)
    shell:
        '''
           echo "msisensor2 msi -M {input.models} {params.file} -o {params.prefix} -b {threads}" | tee {log}
           msisensor2 msi -M {input.models} {params.file} -o {params.prefix} -b {threads} 2>> {log}

           ## Export rule env details
           conda env export --no-builds > info/msisensor2.info
        '''
