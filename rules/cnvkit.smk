## Extracts the purity value and generates the -m and --purity options in the CNVkit command
def purity_checker(sample):
    purity="" # Purity defaults to ''
    file=PREDIR+"/facets/"+sample+"_optimalpurityvalue.txt"
    if os.path.exists(file):
        df=pd.read_csv(file, na_filter=False, delimiter="\t")
        if df["purity"][0] != "NA": ## Valid purity value, set purity options in the command
            purity="-m clonal --purity %s" % df["purity"][0]
            print(purity)
    return purity # Either "-m clonal --purity VAL" or ""

## Run CNVkit to call copy number variations for the tumor samples 
rule cnvkit:
    input:
        cnn=rules.retrieve_cnvkit_cnn.output.cnn,
        bam=rules.apply_bqsr.output.bam,
        bai=rules.apply_bqsr.output.bai,
    output:
        cns=paths.cnvkit.cns,
        call_cns=paths.cnvkit.call_cns,
        scatter=paths.cnvkit.scatter
    benchmark:
        'benchmark/{sample}_cnvkit.tab'
    log:
        'log/{sample}_cnvkit.log'
    conda:
        "../envs/cnvkit.yaml"
    params:
        output_dir=PREDIR+"/cnvkit/{sample}"
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "cnvkit.py batch {input.bam} -r {input.cnn} -p {threads} --scatter --diagram -d {params.output_dir}" | tee {log}
          cnvkit.py batch {input.bam} -r {input.cnn} -p {threads} --scatter --diagram -d {params.output_dir} 2>> {log}
  
          ## Export rule env details
          conda env export --no-builds > info/cnvkit.info
        '''

## Add somatic SNP and purity information to CNVkit's refined call
## INPUT VCF AND TBI PATHS WILL NEED TO CHANGE ONCE SOMATIC OUTPUT IS AVAILABLE
rule cnvkit_enhance:
    input:
        cns=rules.cnvkit.output.call_cns,
        vcf=paths.cnvkit.vcf,
        tbi=paths.cnvkit.tbi,
        purity=lambda wildcards: Path(PREDIR) / "facets" / f"{pairings_df.at[wildcards.sample, 'tumor']}_optimalpurityvalue.txt" if pairings_df.at[wildcards.sample,'type'] == "TN" else []
    output:
        enhanced_cns=paths.cnvkit.enhanced_cns    
    benchmark:
        'benchmark/{sample}_cnvkit_enhance.tab'
    log:
        'log/{sample}_cnvkit_enhance.log'
    conda:
        "../envs/cnvkit.yaml"
    params:
        tumor="-i {sample}",
        normal=lambda wildcards: f"-n {pairings_df.at[wildcards.sample, 'normal']}" if pairings_df.at[wildcards.sample,'type'] == "TN" else "",
        purity=lambda wildcards: purity_checker(wildcards.sample)
    shell:
        '''
          echo "cnvkit.py call {input.cns} -y -v {input.vcf} {params.tumor} {params.normal} {params.purity} -o {output.enhanced_cns}" | tee {log}
          cnvkit.py call {input.cns} -y -v {input.vcf} {params.tumor} {params.normal} {params.purity} -o {output.enhanced_cns} 2>> {log}
        '''
