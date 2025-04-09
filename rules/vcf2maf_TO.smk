rule vcf2maf_TO:
    input:
        vcf=rules.somatic_twist_TO.output.twist
    output:
        vcf_uncompressed=paths.vcf2maf_TO.uncompress,
        maf=paths.vcf2maf_TO.maf
    benchmark:
        'benchmark/{sample}_vcf2maf_TO.tab'
    log:
        'log/{sample}_vcf2maf_TO.log'
    params:
        indexseq=paths.genome.fa,
        vep_cache_dir="../wes2_output/resources/vep/cache/",
        tumor=lambda wildcards: {pairings_df.at[wildcards.sample,'tumor']} if pairings_df.at[wildcards.sample,'type'] == "TO" else []
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
          --vep-data {params.vep_cache_dir}" | tee {log}

          vcf2maf.pl --input-vcf {output.vcf_uncompressed} \
          --output-maf {output.maf} \
          --ref-fasta {params.indexseq} \
          --cache-version 111 \
          --vep-path $CONDA_PREFIX/share/ensembl-vep-111.0-0 \
          --ncbi-build GRCh38 \
          --species homo_sapiens \
          --tumor-id {params.tumor} \
          --vep-data {params.vep_cache_dir} 2>> {log}
        '''
