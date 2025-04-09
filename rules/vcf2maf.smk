rule vcf2maf:
    input:
        vcf=rules.somatic_twist.output.twist
    output:
        vcf_uncompressed=paths.vcf2maf.uncompress,
        maf=paths.vcf2maf.maf
    benchmark:
        'benchmark/{sample}_vcf2maf.tab'
    log:
        'log/{sample}_vcf2maf.log'
    params:
        indexseq=paths.genome.fa,
        vep_cache_dir="../wes2_output/resources/vep/cache",
        normal=lambda wildcards: {pairings_df.at[wildcards.sample,'normal']} if pairings_df.at[wildcards.sample,'type'] == "TN" else [],
        tumor=lambda wildcards: {pairings_df.at[wildcards.sample,'tumor']} if pairings_df.at[wildcards.sample,'type'] == "TN" else []
    conda:
        "../envs/vcf2maf.yaml"
    shell:
        '''
          echo "bcftools view {input.vcf} > {output.vcf_uncompressed}" | tee {log}
          bcftools view {input.vcf} > {output.vcf_uncompressed}
        
          echo "vcf2maf.pl --input-vcf {output.vcf_uncompressed} \
          --output-maf {output.maf} \
          --ref-fasta {params.indexseq} \
          --cache-version 111 \
          --vep-path $CONDA_PREFIX/share/ensembl-vep-111.0-0 \
          --ncbi-build GRCh38 \
          --species homo_sapiens \
          --tumor-id {params.tumor} \
          --normal-id {params.normal} \
          --vep-data {params.vep_cache_dir}" | tee {log}
          
          vcf2maf.pl --input-vcf {output.vcf_uncompressed} \
          --output-maf {output.maf} \
          --ref-fasta {params.indexseq} \
          --cache-version 111 \
          --vep-path $CONDA_PREFIX/share/ensembl-vep-111.0-0 \
          --ncbi-build GRCh38 \
          --species homo_sapiens \
          --tumor-id {params.tumor} \
          --normal-id {params.normal} \
          --vep-data {params.vep_cache_dir} 2>> {log}
        '''
