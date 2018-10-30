source("../set_path_variables.R")
source("../process_table_functions.R")
source("../visualization_scripts.R")
require(reshape2)
require(magrittr)
require(tidyr)

#################################################
### FUNCTIONS FOR TRANSLATING SAMPLE & GENE NAMES
#################################################
make_RNA_cancer_type_set <- function(){

	in_table = fread(histo_meta_path, header=T, sep="\t")

	icgc_meta_table = data.frame(in_table, check.names=F)
	icgc_meta_table = icgc_meta_table[which(icgc_meta_table$wgs_white_black_gray == "Whitelist"),]
	icgc_meta_table = icgc_meta_table[which(!is.na(icgc_meta_table$histotype)),]

	#make cancer-specific labels
	#proj_codes = toChar(icgc_meta_table$project_code)
	#proj_codes = substr(proj_codes,1,nchar(proj_codes)-3)
	#icgc_meta_table$project_code = proj_codes

	#make the names list
	sample_type_list = data.frame(Cancer_Type=icgc_meta_table$histotype, Sample_ID=icgc_meta_table$aliquot_id, ICGC_DONOR_ID=icgc_meta_table$icgc_donor_id)

	return(sample_type_list)
}

#make_Analysis_cancer_type_set <- function(){
#
#	in_table = fread(rnaseq_meta_table_1_3, header=T, sep="\t")
#	icgc_meta_table = data.frame(in_table, check.names=F)
#	icgc_meta_table = icgc_meta_table[which(icgc_meta_table$wgs_white_black_gray == "Whitelist"),]
#
#
#	#make cancer-specific labels
#	proj_codes = toChar(icgc_meta_table$project_code)
#	proj_codes = substr(proj_codes,1,nchar(proj_codes)-3)
#	icgc_meta_table$project_code = proj_codes
#
#	sample_type_list = data.frame(Cancer_Type=icgc_meta_table$project_code, Analysis_ID=icgc_meta_table$analysis_id, ICGC_DONOR_ID=icgc_meta_table$icgc_donor_id)
#
#	# remove the replicates
#	#rep_ids = which(duplicated(sample_type_list$ICGC_DONOR_ID))
#	#sample_type_list = sample_type_list[-rep_ids,]
#
#
#	return(sample_type_list)
#}

make_WGS_cancer_type_set <- function(){

	icgc_meta_table = read.table(icgc_meta_path, header=T, sep="\t")
	icgc_meta_table = icgc_meta_table[which(icgc_meta_table$wgs_exclusion_white_gray == "Whitelist"), c("icgc_donor_id", "tumor_wgs_aliquot_id")]

	rna_table = fread(histo_meta_path, header=T, sep="\t")
	rna_table = data.frame(rna_table, check.names=F)
	rna_table = rna_table[which(rna_table$wgs_white_black_gray == "Whitelist"),]
	rna_table = rna_table[,c("histotype", "aliquot_id", "icgc_donor_id")]
	rna_table = rna_table[which(!is.na(rna_table$histotype)),]

	total_table = merge(icgc_meta_table, rna_table, by="icgc_donor_id")
	total_table = unique(total_table)

	#make cancer-specific labels
	#proj_codes = toChar(total_table$project_code)
	#proj_codes = substr(proj_codes,1,nchar(proj_codes)-3)
	#total_table$project_code = proj_codes

	#make the names list
	sample_type_list = data.frame(Cancer_Type=total_table$histotype, Sample_ID=total_table$tumor_wgs_aliquot_id, ICGC_DONOR_ID=total_table$icgc_donor_id)


	return(sample_type_list)


}

translate_NAMETYPE_TO_ENS <-function(genes_to_trans, org="hsapiens_gene_ensembl", name_type='hgnc_symbol'){

	# select mart and data set
	bm <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
	bm <- useDataset(org, mart=bm)

	# Get ensembl gene ids and hgnc_symbol
	# this is our translation table
	EG2HGNC_name <- getBM(mart=bm, attributes=c('ensembl_gene_id', name_type))
	EG2HGNC_name <- na.omit(EG2HGNC_name)


	genes_to_trans = merge(EG2HGNC_name, genes_to_trans, by=name_type, all.y=TRUE)
	genes_to_trans=unique(genes_to_trans)


	return(genes_to_trans)
}


translate_ENS_TO_HGNC <-function(genes_to_trans, org="hsapiens_gene_ensembl", name_type='ensembl_gene_id'){

	# select mart and data set
	bm <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
	bm <- useDataset(org, mart=bm)

	# Get ensembl gene ids and hgnc_symbol
	# this is our translation table
	EG2HGNC_name <- getBM(mart=bm, attributes=c('hgnc_symbol', name_type))
	EG2HGNC_name <- na.omit(EG2HGNC_name)


	genes_to_trans = merge(EG2HGNC_name, genes_to_trans, by=name_type, all.y=TRUE)
	genes_to_trans=unique(genes_to_trans)


	return(genes_to_trans)
}


#################################################
### FUNCTIONS FOR VALIDATING SAMPLES AND GENES
#################################################

get_valid_samples <- function(){

	in_table = fread(rnaseq_meta_table_1_3, header=T, sep="\t")
	in_table = data.frame(in_table, check.names=F)
	in_table = in_table[, c("project_code", "icgc_donor_id", "wgs_white_black_gray")]
	in_table = unique(in_table)
	icgc_samples_whitelist = in_table[which(in_table$wgs_white_black_gray == "Whitelist"), c("icgc_donor_id")]

	return(unique(icgc_samples_whitelist))

}

get_valid_genes <- function(){

	in_table = fread(prt_genes, header=F, sep="\t")
	in_table = data.frame(in_table, check.names=F)
	in_table[,1] = remove_period(in_table[,1])

	colnames(in_table) = c("Ensembl_Gene_ID", "HGNC_Gene_Symbol")

	return(in_table)

}

validate_table <- function(table_to_validate, doGenes=T){

	valid_samples = tolower(get_valid_samples())
	valid_genes = get_valid_genes()

	if(doGenes){
		gene_keep_idx = which(table_to_validate$Ensembl_Gene_ID %in% valid_genes$Ensembl_Gene_ID )
		table_to_validate = table_to_validate[gene_keep_idx,]
	}

	samp_keep_idx = which(tolower(toChar(table_to_validate$ICGC_DONOR_ID)) %in% valid_samples )
	table_to_validate = table_to_validate[samp_keep_idx,]

	return(list(clean_table=table_to_validate, num_genes=length(valid_genes$Ensembl_Gene_ID), num_samps=length(valid_samples) ))

}


#################################################
### TABLE PROCESSING
#################################################
read_format_table <- function(cancer_type_table, table_path, num_header, outfile){

	### This method is for processing the fusion tables, or other general tables that are binary and Aliquot IDs ##

	# PREREQUISITES #
	# gene id column must be named Ensembl_Gene_ID

	# read in table and melt it
	in_table = fread(table_path, header=T, sep="\t")
	in_table = data.frame(in_table, check.names=F)

	if(length(which(colnames(in_table) == "ensembl_gene_id")) > 0){
		colnames(in_table)[which(colnames(in_table) == "ensembl_gene_id")] = "Ensembl_Gene_ID"
	}

	if(length( "Ensembl_Gene_ID" %in% colnames(in_table)) == 0){
		print("gene id column must be EXACTLY named Ensembl_Gene_ID")
		return(0)
	}

	gene_id_col = which(colnames(in_table) == "Ensembl_Gene_ID")
	in_table = in_table[,c(gene_id_col, num_header:ncol(in_table))]
	melt_in_table = melt(in_table, id.vars="Ensembl_Gene_ID")
	colnames(melt_in_table) = c("Ensembl_Gene_ID", "Sample_ID", "isVal")
	melt_in_table$Sample_ID = gsub("[.]", "-", melt_in_table$Sample_ID)

	#merge with cancer type
	melt_in_table = data.table(melt_in_table)
	cancer_type_table = data.table(cancer_type_table)
	melt_in_table = merge(melt_in_table, cancer_type_table, by="Sample_ID")
	melt_in_table = data.frame(melt_in_table)

	#format and write
	long_in_table = melt_in_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", "isVal")]
	colnames(long_in_table) = c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", outfile)
	write.table(long_in_table, paste("master_tables/", outfile, "_full.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make gene_centric version
	gc_long_in_table = unique(long_in_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", outfile)])
	gc_long_in_table = dcast(data = gc_long_in_table,formula = Ensembl_Gene_ID~ICGC_DONOR_ID, fun.aggregate = sum, value.var = outfile)
	write.table(gc_long_in_table, paste("master_tables/", outfile, "_gene_centric.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make version with non-zero values
	long_in_table = long_in_table[which(long_in_table[,outfile] != 0),]
	write.table(long_in_table, paste("master_tables/", outfile, ".tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)
	return(long_in_table)
}

collapse_and_sum <- function(table_to_smash, col_smash_id){

	col_interest = which(! colnames(table_to_smash) %in% c(col_smash_id, "ICGC_DONOR_ID", "Ensembl_Gene_ID", "Cancer_Type", "hgnc_symbol", "gene_sample_pairs_curr") )

	DT <- as.data.table(table_to_smash)
	summed_table = DT[, lapply(.SD, sum), by = c(col_smash_id), .SDcols = colnames(table_to_smash)[col_interest]]

	summed_table = data.frame(summed_table)

	return(summed_table)

}


read_format_rna_edit_table <- function(cancer_type_table, rna_edit_table_path, num_header, outfile){

	### This method is for processing the fusion tables, or other general tables that are binary and Aliquot IDs ##

	# PREREQUISITES #
	# gene id column must be named Ensembl_Gene_ID

	# read in table and melt it
	in_table = fread(rna_edit_table_path, header=T, sep="\t")
	in_table = data.frame(in_table, check.names=F)

	# remove non-prt coding genes
	prt_genes_table = data.frame(fread(prt_genes, header=F, sep="\t"))
	prt_genes_table = remove_period(prt_genes_table[,1])
	in_table = in_table[which(in_table$ensembl_gene_id %in% prt_genes_table),]

	# remove events that occur in less than 3 samples
	in_table_sums = rowSums(in_table[,num_header:ncol(in_table)])
	in_table_idx = which(in_table_sums >= 3)
	in_table = in_table[in_table_idx,]


	if(length( "ensembl_gene_id" %in% colnames(in_table)) == 0){
		print("gene id column must be EXACTLY named Ensembl_Gene_ID")
		return(0)
	}

	# make sure the mutation is is classified functional
	in_table = in_table[which(!in_table$predicted_function %in% c("-","slient")),]
	in_table = in_table[which(in_table$Region %in% c("exonic")),]

	# make sure the mutation is NOT in already identified missense position
	missense_table = data.frame(fread(missense_table_path, header=F, sep="\t"), check.names=F)
	all_mut = paste(missense_table[,3], missense_table[,4], missense_table[,5], sep="_")
	in_table_pos = paste(in_table[,1], in_table[,2], in_table[,3], sep="_")
	in_table_remove = which(in_table_pos %in% all_mut)
	in_table[in_table_remove,num_header:ncol(in_table)] = 0

	# merge with cancer_type
	gene_id_col = colnames(in_table)[c(1:4,9)]
	in_table = in_table[,c(c(1:4,9), num_header:ncol(in_table))]
	melt_in_table = melt(in_table, id.vars=gene_id_col)
	colnames(melt_in_table) = c("chr", "start", "end", "Ensembl_Gene_ID", "genome_change", "Sample_ID", "isVal")
	melt_in_table$Sample_ID = gsub("[.]", "-", melt_in_table$Sample_ID)

	#merge with cancer type
	melt_in_table_dt = data.table(melt_in_table)
	cancer_type_table_dt = data.table(cancer_type_table)
	melt_in_table = merge(melt_in_table_dt, cancer_type_table_dt, by="Sample_ID")
	melt_in_table = data.frame(melt_in_table)


	# make sure its seen >3 times in a cancer type
	temp = data.table(melt_in_table)
	temp %<>% unite(pos, chr, start, genome_change, remove = FALSE)
	melt_in_table = data.frame(temp)

	collapse_ct = melt_in_table[,c("pos", "Cancer_Type", "isVal")]
	colnames(collapse_ct)[3] = c("rna_edit")
	collapse_ct = collapse_and_sum(collapse_ct, c("pos", "Cancer_Type"))
	pass_pos = collapse_ct$pos[which(collapse_ct$rna_edit > 2)]
	filtered_in_table = melt_in_table[melt_in_table$pos %in% pass_pos,]

	# now here filter out multiple genes occuring
	filtered_in_table = filtered_in_table[filtered_in_table$isVal == 1,]

	#format and write
	long_in_table = filtered_in_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", "isVal")]
	colnames(long_in_table) = c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", outfile)
	long_in_table = unique(long_in_table)
	write.table(long_in_table, paste("master_tables/", outfile, "_full.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make gene_centric version
	gc_long_in_table = unique(long_in_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", outfile)])
	gc_long_in_table = dcast(data = gc_long_in_table,formula = Ensembl_Gene_ID~ICGC_DONOR_ID, fun.aggregate = sum, value.var = outfile)
	write.table(gc_long_in_table, paste("master_tables/", outfile, "_gene_centric.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make version with non-zero values
	long_in_table = long_in_table[which(long_in_table[,outfile] != 0),]
	write.table(long_in_table, paste("master_tables/", outfile, ".tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)
	return(long_in_table)
}


read_format_splice_table <- function(cancer_type_set, splice_table_path, num_header, outfile, binarize=F){

	# read in table and melt it
	splice_table = fread(splice_table_path, header=T, sep="\t")
	splice_table = data.frame(splice_table, check.names=F)
	sample_names = colnames(splice_table)[num_header:ncol(splice_table)]
	melt_splice_table = melt(splice_table, id.vars="gene_id")
	#melt_splice_table = melt_splice_table[which(melt_splice_table$value != 0),]
	colnames(melt_splice_table) = c("Ensembl_Gene_ID", "Sample_ID", "isSplice")
	melt_splice_table$Sample_ID = gsub("[.]", "-", melt_splice_table$Sample_ID)

	#merge with cancer type
	melt_splice_table = data.table(melt_splice_table)
	cancer_type_set = data.table(cancer_type_set)
	melt_splice_table = merge(melt_splice_table, cancer_type_set, by="Sample_ID")
	melt_splice_table = data.frame(melt_splice_table)
	melt_splice_table_merge = melt_splice_table
	sum(melt_splice_table_merge[melt_splice_table_merge$ICGC_DONOR_ID == "DO52691", "isSplice"] > 4.5)

	# remove duplicates
	paste_table = paste(melt_splice_table[,1], melt_splice_table[,2], melt_splice_table[,4], melt_splice_table[,5], sep="_")
	duplicate_idx = duplicated(paste_table)
	melt_splice_table = melt_splice_table[which(!duplicate_idx),]
	undup_melt_splice_table = melt_splice_table
	sum(undup_melt_splice_table[undup_melt_splice_table$ICGC_DONOR_ID == "DO52691", "isSplice"] > 4.5)

	#format
	melt_splice_table = melt_splice_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", "isSplice")]
	melt_splice_table_orig = melt_splice_table
	head(melt_splice_table)

	if(binarize){
		top_5_percent = quantile(melt_splice_table$isSplice, 0.95)
		melt_splice_table[melt_splice_table_orig$isSplice > top_5_percent, c("isSplice")] = 1
		melt_splice_table[melt_splice_table_orig$isSplice <= top_5_percent, c("isSplice")] = 0
	}

	# validate tables
	#res = validate_table(melt_splice_table, doGenes=F)
	#alt_table_clean = res$clean_table
	sum(melt_splice_table[melt_splice_table$ICGC_DONOR_ID == "DO52691", "isSplice"] > 4.5)

	write.table(melt_splice_table, paste("master_tables/", outfile, "_full.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make gene_centric version
	gc_long_in_table = unique(melt_splice_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "isSplice")])
	gc_long_in_table = dcast(data = gc_long_in_table,formula = Ensembl_Gene_ID~ICGC_DONOR_ID, fun.aggregate = sum, value.var = "isSplice")
	write.table(gc_long_in_table, paste("master_tables/", outfile, "_gene_centric.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make version with non-zero values
	melt_splice_table = melt_splice_table[which(melt_splice_table$isSplice != 0),]
	sum(melt_splice_table[melt_splice_table$ICGC_DONOR_ID == "DO52691", "isSplice"] > 4.5)
	write.table(melt_splice_table, paste("master_tables/", outfile, ".tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)



}

read_format_de_table <- function(cancer_type_set, de_table_path, num_header, outfile, log2FC_cutoff, binarize=F){

	# read in table and melt it
	de_table = fread(de_table_path, header=T, sep="\t")
	de_table = data.frame(de_table, check.names=F)
	gene_id_col = which(colnames(de_table) == "Ensembl_Gene_ID")
	de_table = de_table[,c(gene_id_col, num_header:ncol(de_table))]
	de_table[is.na(de_table)] = 0

	melt_de_table = melt(de_table, id.vars="Ensembl_Gene_ID")
	#melt_de_table = melt_de_table[which(melt_de_table$value != 0),]
	colnames(melt_de_table) = c("Ensembl_Gene_ID", "Sample_ID", "log2FC")
	melt_de_table$Sample_ID = gsub("[.]", "-", melt_de_table$Sample_ID)

	melt_de_table$log2FC[which(melt_de_table$log2FC < -1000000)] = 0

	#merge with cancer type
	melt_de_table = merge(melt_de_table, cancer_type_set, by="Sample_ID")
	#melt_de_table = melt_de_table[!(melt_de_table$Cancer_Type == "BRCA" | melt_de_table$Cancer_Type == "STAD"), ]

	# write out the full table
	melt_de_table = unique(melt_de_table)
	bin_full_table = melt_de_table
	if(binarize){
		bin_full_table$log2FC = 0
		bin_full_table$log2FC[which(melt_de_table$log2FC > log2FC_cutoff)] = 1
		bin_full_table$log2FC[which(melt_de_table$log2FC < -1*log2FC_cutoff)] = 1
	}
	bin_full_table = bin_full_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", "log2FC")]

	colnames(bin_full_table)[4] = "isDE"
	write.table(bin_full_table, paste("master_tables/", outfile, "_full.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	bin_full_table = bin_full_table[which(bin_full_table$isDE != 0),]
	write.table(bin_full_table, paste("master_tables/", outfile, ".tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make gene_centric version
	gc_long_in_table = unique(bin_full_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "log2FC")])
	gc_long_in_table = dcast(data = gc_long_in_table,formula = Ensembl_Gene_ID~ICGC_DONOR_ID, fun.aggregate = sum, value.var = outfile)
	write.table(gc_long_in_table, paste("master_tables/", outfile, "_gene_centric.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)


	# filtering method
	melt_de_table_up = melt_de_table[which(melt_de_table$log2FC > log2FC_cutoff),]
	melt_de_table_down = melt_de_table[which(melt_de_table$log2FC < -1*log2FC_cutoff),]

	melt_de_table_up$log2FC = 1
	melt_de_table_down$log2FC = 1

	# format the table
	melt_de_table_up = melt_de_table_up[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", "log2FC")]
	melt_de_table_down = melt_de_table_down[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", "log2FC")]

	colnames(melt_de_table_up)[4] = "isDE_up"
	colnames(melt_de_table_down)[4] = "isDE_down"

	#format and write
	write.table(melt_de_table_up, paste("master_tables/", "de_up", ".tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)
	write.table(melt_de_table_down, paste("master_tables/", "de_down", ".tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)



}

read_format_eQTL_table <- function(cancer_type_set, table_path, num_header, outfile){

	eQTL_table = read.table(table_path, header=T, row.names=NULL, sep="\t")
	sample_names = colnames(eQTL_table)[num_header:ncol(eQTL_table)]
	gene_names = unique(remove_period(eQTL_table$Ensembl_Gene_ID))

	# combine the same genes
	reduced_eQTL_table = matrix(,ncol=length(9:ncol(eQTL_table)), nrow=length(gene_names))
	idx = 1
	for(curr_gene in gene_names){
		eQTL_table_gene = eQTL_table[which(remove_period(eQTL_table$Ensembl_Gene_ID)==curr_gene),9:ncol(eQTL_table)]
		eQTL_table_gene = colSums(eQTL_table_gene)
		reduced_eQTL_table[idx,] = eQTL_table_gene
		idx = idx + 1
	}
	colnames(reduced_eQTL_table) = colnames(eQTL_table)[9:ncol(eQTL_table)]
	rownames(reduced_eQTL_table) = gene_names


	# read in table and melt it
	melt_eQTL_table = melt(reduced_eQTL_table, id.vars="Ensembl_Gene_ID")
	melt_eQTL_table = melt_eQTL_table[which(melt_eQTL_table$value != 0),]
	colnames(melt_eQTL_table) = c("Ensembl_Gene_ID", "WGS_Aliquot_ID", "isVal")
	melt_eQTL_table$WGS_Aliquot_ID = gsub("[.]", "-", melt_eQTL_table$WGS_Aliquot_ID)

	#merge with cancer type
	melt_eQTL_table = merge(melt_eQTL_table, cancer_type_set, by="WGS_Aliquot_ID")

	#format and write
	long_eqtl_table = melt_eQTL_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", "isVal")]
	colnames(long_eqtl_table) = c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", outfile)
	write.table(long_eqtl_table, paste("master_tables/", outfile, ".tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)


}

read_format_variant_table <- function(cancer_type_set, variant_table_path, outfile){

	var_table = fread(variant_table_path, header=T, sep="\t")
	var_table = data.frame(var_table)
	colnames(var_table) = c("Cancer_Type", "Ensembl_Gene_ID", "ICGC_DONOR_ID")
	var_table = var_table[,-1]

	var_table = merge(cancer_type_set, var_table, by="ICGC_DONOR_ID")

	cancer_table = var_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type")]
	isVariant = rep(1, nrow(cancer_table))
	cancer_table = cbind(cancer_table, isVariant)
	colnames(cancer_table) = c("hgnc_symbol", "ICGC_DONOR_ID", "Cancer_Type", outfile)

	# translate gene ids to ensembl
	cancer_table_trans = translate_NAMETYPE_TO_ENS(cancer_table, org="hsapiens_gene_ensembl", name_type='hgnc_symbol')
	cancer_table_trans = na.omit(cancer_table_trans)
	colnames(cancer_table_trans)[2] = c("Ensembl_Gene_ID")
	cancer_table_trans_final = cancer_table_trans[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", outfile)]

	write.table(cancer_table_trans_final, paste("master_tables/", outfile, ".tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make gene_centric version
	gc_long_in_table = unique(cancer_table_trans_final[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", outfile)])
	gc_long_in_table = dcast(data = gc_long_in_table,formula = Ensembl_Gene_ID~ICGC_DONOR_ID, fun.aggregate = sum, value.var = outfile)
	write.table(gc_long_in_table, paste("master_tables/", outfile, "_gene_centric.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)


}


read_format_ase_table_roland <- function(cancer_type_set, ase_table_path_all, outfile){

	ase_table = fread(ase_table_path_all, header=T, sep="\t")
	ase_table = data.frame(ase_table, check.names=F)
	sample_names = toChar(ase_table$rna_seq_aliquot_id)
	gene_names = remove_period(ase_table$geneid)

	# filtering method
	values = as.integer(ase_table$padj_cc < 0.05 & (abs(ase_table$mean_allelic_imbalance-ase_table$mean_cn_ratio)>0.2))
	values = as.integer(ase_table$padj_cc < 0.05 )
	values = rep(1, nrow(ase_table))

	#pdf(paste("master_tables/", outfile, "_ratio.pdf", sep=""))
	#hist(abs(ase_table$mean_allelic_imbalance-ase_table$mean_cn_ratio))
	#dev.off()
	#pdf(paste("master_tables/", outfile, "_padj_cc.pdf", sep=""))
	#hist(ase_table$padj_cc)
	#dev.off()

	# filter out normal ASE genes
	ase_norm = data.frame(fread(ase_norm_table_path, header=T, sep="\t"), check.names=F)
	ase_norm$norm_ratio = ase_norm$padj_norm/ase_norm$nsamples

	# dont allow this filtering for now
	ase_norm = ase_norm[ase_norm$nsamples > 10 & ase_norm$norm_ratio > 0.5,"geneid"]

	long_table = data.frame(gene_names, sample_names, values)
	colnames(long_table) = c("Ensembl_Gene_ID", "Sample_ID", outfile)


	# plot distributions of ASE in genes we are removing
	pdf("/cluster/project/raetsch/lab/03/ICGC/genecentric_analysis/results/qc_files/qc_ase/hist_genes_removed.pdf")
	for(curr_gene in unique(ase_norm)){
		tab_interest = ase_table[which(ase_table$geneid == curr_gene),]
		hist(tab_interest$ase_phased, breaks=30, main=curr_gene)

		g = ggplot(tab_interest, aes(x=ase_phased, fill=geneid)) +
		    geom_density(alpha=.5, position="identity") +
		    ggtitle(curr_gene)
		print(g)

	}
	dev.off()



	# now filter
	long_table = long_table[which(!long_table$Ensembl_Gene_ID %in% ase_norm),]

	long_table = merge(cancer_type_set, long_table, by="Sample_ID")
	long_table = long_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", outfile)]
	write.table(long_table, paste("master_tables/", outfile, "_full.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make gene_centric version
	gc_long_in_table = unique(long_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", outfile)])
	gc_long_in_table = dcast(data = gc_long_in_table,formula = Ensembl_Gene_ID~ICGC_DONOR_ID, fun.aggregate = sum, value.var = outfile)
	write.table(gc_long_in_table, paste("master_tables/", outfile, "_gene_centric.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)


}



read_format_cn_table <- function(cancer_type_set, cn_table_path, outfile){

	#read
	cn_table = fread(cn_table_path, header=T, sep="\t")
	cn_table = data.frame(cn_table, check.names=F)
	colnames(cn_table)[1] = "Sample_ID"

	#merge
	cancer_type_set = data.table(cancer_type_set)
	cn_table = data.table(cn_table)
	cn_table = merge(cancer_type_set, cn_table, by=c("Sample_ID"))
	cn_table = data.frame(cn_table)

	#melt
	melt_cn_table = melt(cn_table, id.vars=c("Sample_ID", "Cancer_Type", "ICGC_DONOR_ID"))
	colnames(melt_cn_table) = c("Sample_ID", "Cancer_Type", "ICGC_DONOR_ID", "Ensembl_Gene_ID", "cn")
	melt_cn_table = melt_cn_table[which(melt_cn_table$cn < 1 |  melt_cn_table$cn > 6),]
	melt_cn_table$Ensembl_Gene_ID = remove_period(melt_cn_table$Ensembl_Gene_ID)

	melt_cn_table = melt_cn_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", "cn")]
	write.table(melt_cn_table, paste("master_tables/", outfile, "_full.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)


}


read_format_expr_outlier_table <- function(cancer_type_set, expr_outlier_table_path, outfile, binarize=F){

	expr_outlier_table = fread(expr_outlier_table_path, header=T, sep="\t")
	expr_outlier_table = data.frame(expr_outlier_table, check.names=F)

	# filtering method
	long_table = expr_outlier_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", "zscore")]
	if(binarize){
		values = rep(0, nrow(expr_outlier_table))
		top_percent = quantile(expr_outlier_table$zscore, 0.975, na.rm=T)
		bot_percent = quantile(expr_outlier_table$zscore, 0.025, na.rm=T)
		values[which(expr_outlier_table$zscore < bot_percent)] = 1
		values[which(expr_outlier_table$zscore > top_percent)] = 1

		long_table = data.frame(expr_outlier_table$Ensembl_Gene_ID, expr_outlier_table$ICGC_DONOR_ID, expr_outlier_table$Cancer_Type, values)
	}

	colnames(long_table) = c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "proj_id", outfile)
	long_table = long_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", outfile)]
	long_table = unique(long_table)

	# now we need to remove samples that come from many donor_ids
	# so we take the more extreme ones, by first sorting by zscore then removing the second duplicate
	long_table$abs_val = abs(long_table$expr_outlier)
	long_table = long_table[order(long_table$abs_val, decreasing=T),]
	dup_rows = which(duplicated(long_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID")]))
	long_table_no_dup = long_table
	if(length(dup_rows) > 0){
		long_table_no_dup = long_table[-dup_rows,]
	}

	cancer_type_set = cancer_type_set[,c("Cancer_Type", "ICGC_DONOR_ID")]
	cancer_type_set = unique(cancer_type_set)
	cancer_type_set_dt = data.table(cancer_type_set)
	long_table_dt = data.table(long_table_no_dup)
	long_table = merge(cancer_type_set_dt, long_table_dt, by="ICGC_DONOR_ID")

	long_table = data.frame(long_table)
	long_table = long_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", outfile)]


	write.table(long_table, paste("master_tables/", outfile, "_full.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	if(binarize){
		# make gene_centric version
		gc_long_in_table = unique(long_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", outfile)])
		gc_long_in_table = dcast(data = gc_long_in_table,formula = Ensembl_Gene_ID~ICGC_DONOR_ID, fun.aggregate = sum, value.var = outfile)
		write.table(gc_long_in_table, paste("master_tables/", outfile, "_gene_centric.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)
	}
}




read_format_ase_table <- function(cancer_type_set, ase_table_path, num_header, outfile){

	ase_table = fread(ase_table_path, header=T, sep="\t")
	ase_table = data.frame(ase_table, check.names=F)
	sample_names = toChar(ase_table$rna_seq_aliquot_id)
	gene_names = remove_period(ase_table$Ensembl_Gene_ID)

	gene_id_col = which(colnames(ase_table) == "Ensembl_Gene_ID")
	reduced_ase_table = ase_table[,c(gene_id_col, num_header:ncol(ase_table))]

	melt_ase_table = melt(reduced_ase_table, id.vars="Ensembl_Gene_ID")
	colnames(melt_ase_table) = c("Ensembl_Gene_ID", "Sample_ID", outfile)

	long_table = merge(cancer_type_set, melt_ase_table, by="Sample_ID")
	long_table = long_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "Cancer_Type", outfile)]

	write.table(long_table, paste("master_tables/", outfile, "_full.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make gene_centric version
	gc_long_in_table = unique(long_table[,c("Ensembl_Gene_ID", "ICGC_DONOR_ID", outfile)])
	gc_long_in_table = dcast(data = gc_long_in_table,formula = Ensembl_Gene_ID~ICGC_DONOR_ID, fun.aggregate = sum, value.var = outfile)
	write.table(gc_long_in_table, paste("master_tables/", outfile, "_gene_centric.tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

	# make version with non-zero values
	long_table = long_table[which(long_table[,c(outfile)] != 0),]
	write.table(long_table, paste("master_tables/", outfile, ".tsv", sep=""), sep="\t", row.name=F, col.names=T, quote=F)

}

######### RAW RUN ############
sum_all_tables <- function(list_of_alt_files){

	# get all gene_sample pairs
	alt_table = NA
	gene_sample_pairs = c()
	for(alt_file in list_of_alt_files){
		print(alt_file)
		curr_alt_table = data.frame(fread(alt_file, header=T), check.names=F)
		curr_alt_table$Ensembl_Gene_ID = remove_period(curr_alt_table$Ensembl_Gene_ID)

		print(unique(curr_alt_table$Cancer_Type))

		if(is.na(alt_table)){
			alt_table = curr_alt_table
		}else{
			alt_table = data.table(alt_table)
			curr_alt_table = data.table(curr_alt_table)
			alt_table = merge(curr_alt_table, alt_table, by=c("ICGC_DONOR_ID", "Ensembl_Gene_ID", "Cancer_Type"), all=T)
		}
		rm(curr_alt_table)
		gc()
	}



	# Unite
	alt_table %<>% unite(gene_sample_pairs_curr, Ensembl_Gene_ID, ICGC_DONOR_ID, remove = FALSE)

 	length(unique(alt_table$gene_sample_pairs_curr)) == length(alt_table$gene_sample_pairs_curr)

	write.table(alt_table, master_table_raw_path, sep="\t", row.name=F, col.names=T, quote=F)
	return(alt_table)


}


#################################################
### RUN
#################################################

run_make_master_table <- function(){

	setwd(paste(base_dir, "/results/", sep=""))

	RNA_sample_type_list = make_RNA_cancer_type_set()
	WGS_sample_type_list = make_WGS_cancer_type_set()


	if(FALSE){

		dom_trans_total = read_format_table(RNA_sample_type_list, dt_table_path, 6, "dom_trans")

		alt_prom = read_format_table(RNA_sample_type_list, alt_prom_table_path, num_header=9, "alt_prom")
		alt_pa = read_format_table(RNA_sample_type_list, apa_table_path, num_header=9, "alt_pa")


		read_format_ase_table_roland(RNA_sample_type_list, ase_table_path_all, "ase_all")

		read_format_cn_table(WGS_sample_type_list, cn_table_path, "cn")

		read_format_expr_outlier_table(RNA_sample_type_list, expr_outlier_table_path, "expr_outlier")

		fusion_total = read_format_table(RNA_sample_type_list, fusion_table_path, 6, "fusion")

		rna_edit = read_format_rna_edit_table(RNA_sample_type_list, rna_edit_table_path, num_header=10, "rna_edit")

		read_format_splice_table(RNA_sample_type_list, splice_allEvents_table_path, 1, "splice_allEvent")

		read_format_variant_table(WGS_sample_type_list, variant_table_path, "variants")
		rna_edit = read_format_rna_edit_table(RNA_sample_type_list, rna_edit_table_path, num_header=10, "rna_edit2")

	}


	fusion_file = "master_tables/fusion.tsv"
	eqtl_somatic_file = "master_tables/eqtl_somatic.tsv"
	ase_file = "master_tables/ase_all_full.tsv"
	de_file_up = "master_tables/de_up.tsv"
	de_file_down = "master_tables/de_down.tsv"
	de_file = "master_tables/de.tsv"
	rna_edit_file = "master_tables/rna_edit2.tsv"
	alt_prom_file = "master_tables/alt_prom.tsv"
	alt_pa_file = "master_tables/alt_pa.tsv"
	cn_file = "master_tables/cn_full.tsv"
	expr_outlier_file = "master_tables/expr_outlier_full.tsv"

	dom_trans = "master_tables/dom_trans.tsv"

	splice_file = "master_tables/splice.tsv"
	splice_allEvent_file = "master_tables/splice_allEvent.tsv"
	var_file = "master_tables/variants.tsv"
	#var_file = "/cluster/project/raetsch/lab/03/ICGC/starks_hackathon/gene_variant_info.tsv"


	#all_res_files = c(var_file, ase_file, cn_file, expr_outlier_file, fusion_file, rna_edit_file, alt_prom_file, splice_file)
	#all_res_files = c(var_file, ase_file, expr_outlier_file, fusion_file, rna_edit_file, alt_prom_file, splice_allEvent_file, cn_file)

	all_res_files = c(var_file, ase_file, expr_outlier_file, fusion_file, rna_edit_file, alt_prom_file, alt_pa_file, splice_allEvent_file)
	master_table_raw = sum_all_tables(all_res_files)

}

#run_make_master_table()
