# pcawg3-transcriptome-integration

Figure 10 is generated using [compare_rna_dna_ranks.R](compare_rna_dna_ranks.R)

Recurrence analysis is done using [recurrance_histogram.R](recurrance_histogram.R)

Master table generation and filtering is done using [make_master_table.R](make_master_table.R) followed by [filter_master_table.R](filter_master_table.R)

Metadata needed for filtering and analysis:

- To account for gene_length: [ensembl_gene_length_gc_path](https://github.com/ICGC-TCGA-PanCancer/pcawg3-transcriptome-integration/blob/master/metadata/mart_export_transcript_length_gc.txt)
    
- To translate gene names: [ensembl_annot_path](https://github.com/ICGC-TCGA-PanCancer/pcawg3-transcriptome-integration/blob/master/metadata/hgnc_ensembl_output2.tsv)

- Filter out the genes from the ranked list that don't pass permutation test: [filter_score_path](https://github.com/ICGC-TCGA-PanCancer/pcawg3-transcriptome-integration/blob/master/metadata/scores_to_filter.tsv)

- Cancer Census genes used in this analysis: [cancer_census_genes_path](https://github.com/ICGC-TCGA-PanCancer/pcawg3-transcriptome-integration/blob/master/metadata/gencodeV14.v7.pancan_subset.ensembleID.list)

- PCAWG driver genes used in this analysis: [driver_gene_path](https://github.com/ICGC-TCGA-PanCancer/pcawg3-transcriptome-integration/blob/master/metadata/cds_driver_genes_22_04_2017.tsv)
