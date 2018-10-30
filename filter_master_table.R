source("../set_path_variables.R")
source("../process_table_functions.R")
source("../visualization_scripts.R")
source("../make_master_table.R")
#source("../filter_master_table.R")
require(reshape2)
library(RColorBrewer)
library(pheatmap)


get_aus_paca_ov_donor_id <- function(){

	in_table = fread(histo_meta_path, header=T, sep="\t")
	icgc_meta_table = data.frame(in_table, check.names=F)
	icgc_meta_table = icgc_meta_table[which(icgc_meta_table$wgs_white_black_gray == "Whitelist"),]

	#make cancer-specific labels
	#paca_ov_meta_table = icgc_meta_table[which(icgc_meta_table$histotype %in% c("Panc-AdenoCA", "Ovary-AdenoCA")),]
	#paca_ov_meta_table = paca_ov_meta_table[which(!paca_ov_meta_table$project_code %in% c("OV-US")),]

	return(paca_ov_meta_table$icgc_donor_id)
}


binarize_table <- function(master_table_raw, col_to_binarize){

	for(curr_col in col_to_binarize){
		values = rep(0, nrow(master_table_raw))

		# take bottom and top 5%
		#top_percent = quantile(master_table_raw[,curr_col], 0.99, na.rm=T)
		#bot_percent = quantile(master_table_raw[,curr_col], 0.01, na.rm=T)

		if(curr_col == "isSplice"){
			values[which(abs(master_table_raw[,curr_col]) > 6)] = 1
			values[which(abs(master_table_raw[,curr_col]) < -6)] = 1

		}else if(curr_col == "expr_outlier"){

			values[which((master_table_raw[,curr_col]) > 3)] = 1
			values[which((master_table_raw[,curr_col]) < -3)] = 1
		}

		master_table_raw[,curr_col] = values
	}

	# replace all NA with zero
	master_table_raw[is.na(master_table_raw)] <- 0

	# remove rows with 0 sum across
	row_sums = rowSums(master_table_raw[,5:ncol(master_table_raw)])
	master_table_raw = master_table_raw[which(row_sums > 0),]

	return(master_table_raw)
}

translate_ENS_to_HGNC_fast <- function(genes_to_trans, org="hsapiens_gene_ensembl", name_type='ensembl_gene_id'){

        # select mart and data set
        #bm <- useMart("ENSEMBL_MART_ENSEMBL", host="useast.ensembl.org")
        #bm <- useDataset(org, mart=bm)

        # Get ensembl gene ids and hgnc_symbol
        # this is our translation table
        #EG2HGNC_name <- getBM(mart=bm, attributes=c('hgnc_symbol', name_type))
        #EG2HGNC_name <- na.omit(EG2HGNC_name)

				EG2HGNC_name = data.frame(fread(ensembl_annot_path))
				colnames(EG2HGNC_name) = c("ensembl_gene_id", "hgnc_symbol")
        genes_to_trans = merge(data.table(EG2HGNC_name), data.table(genes_to_trans), by=name_type, all.y=TRUE)
				genes_to_trans = data.frame(genes_to_trans)
        genes_to_trans=unique(genes_to_trans)


        return(genes_to_trans)
}


translate_validate_table <- function(master_table_bin){

	# validate tables
	res = validate_table(master_table_bin)
	master_table_bin_valid = res$clean_table
	num_genes = res$num_genes
	num_samps = res$num_samps

	# merge with hgnc_symbol
	colnames(master_table_bin_valid)[2] = "ensembl_gene_id"
	master_table_bin_valid_trans = translate_ENS_to_HGNC_fast(master_table_bin_valid, org="hsapiens_gene_ensembl", name_type='ensembl_gene_id')
	colnames(master_table_bin_valid_trans)[1] = "Ensembl_Gene_ID"

	idx_na = which(is.na(master_table_bin_valid_trans$hgnc_symbol))
	master_table_bin_valid_trans$hgnc_symbol[idx_na] = master_table_bin_valid_trans$Ensembl_Gene_ID[idx_na]

	idx_empty = which(master_table_bin_valid_trans$hgnc_symbol == "")
	master_table_bin_valid_trans$hgnc_symbol[idx_empty] = master_table_bin_valid_trans$Ensembl_Gene_ID[idx_empty]

	return(master_table_bin_valid_trans)

}


permute_samples <- function(master_table){

	all_genes = unique(master_table$Ensembl_Gene_ID)
	alt_interest = which(! colnames(master_table) %in% c("ICGC_DONOR_ID", "Ensembl_Gene_ID", "Cancer_Type", "hgnc_symbol", "gene_sample_pairs_curr") )
	shuffled_table = NA
	for(curr_alt in alt_interest){

		curr_tab = master_table[,c("ICGC_DONOR_ID", "Ensembl_Gene_ID", "hgnc_symbol", colnames(master_table)[curr_alt])]

		# make matrix
		tab1 = dcast(data = curr_tab,formula = Ensembl_Gene_ID+hgnc_symbol~ICGC_DONOR_ID, fun.aggregate = sum)

		# shuffle samples
		tab1_shuffle = apply(tab1[,3:ncol(tab1)], 2, sample)
		tab1_shuffle = data.frame(tab1[,c(1,2)], tab1_shuffle)
		colnames(tab1_shuffle) = colnames(tab1)

		# flatten
		tab1_flat = melt(data = tab1_shuffle,id =c("Ensembl_Gene_ID", "hgnc_symbol"))
		colnames(tab1_flat) = c("Ensembl_Gene_ID", "hgnc_symbol", "ICGC_DONOR_ID", colnames(master_table)[curr_alt])

		tab1_flat = data.table(tab1_flat)
		if(is.na(shuffled_table)){
			shuffled_table = tab1_flat
		}else{
			shuffled_table = data.table(shuffled_table)
			shuffled_table = merge(shuffled_table, tab1_flat, by=c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "hgnc_symbol"))
		}
	}

	perm_master_table = data.frame(shuffled_table)

	return(perm_master_table)
}


permute_samples_inner <- function(in_table, all_samp_ids){
	res = NA

	samp_missing = all_samp_ids[which(!all_samp_ids %in% in_table$ICGC_DONOR_ID)]

	padded_table = data.frame(samp_missing, unique(in_table$Ensembl_Gene_ID), NA, unique(in_table$hgnc_symbol), NA, 0, 0, 0, 0, 0, 0, 0, 0)
	colnames(padded_table) = colnames(in_table)

	full_table = rbind(in_table, padded_table)
	mini_table = full_table[,6:ncol(in_table)]

	res = apply(mini_table, 2, sample)
	res = cbind(full_table[,1:5], res)

	remove_col = rowSums(res[6:ncol(in_table)])
	res = res[which(remove_col > 0),]
	return( res )
}


final_process_master_table <- function(master_table){


	# merge with the CN
	cn_table = data.frame(fread(cn_processed_table_path))
	cn_table = cn_table[which(cn_table$cn < 1 | cn_table$cn > 6),]
	cn_table[is.na(cn_table)] = 0
	cn_table = cn_table[which(cn_table$Cancer_Type != 0),]
	cn_table = data.table(cn_table)
	cn_table = unique(cn_table)
	cn_table %<>% unite(gene_sample_pairs_curr, Ensembl_Gene_ID, ICGC_DONOR_ID, remove = FALSE)
	cn_table$cn = 1


	master_table = merge(data.table(master_table), data.table(cn_table), by=c("Ensembl_Gene_ID", "ICGC_DONOR_ID", "gene_sample_pairs_curr", "Cancer_Type"), all=T)
	master_table = data.frame(master_table)

	master_table = subset(master_table, select=-c(hgnc_symbol))

	# get hgnc_symbol
	colnames(master_table)[1] = "ensembl_gene_id"
	master_table = translate_ENS_to_HGNC_fast(master_table, org="hsapiens_gene_ensembl", name_type='ensembl_gene_id')
	colnames(master_table)[1] = "Ensembl_Gene_ID"

	master_table[is.na(master_table)] = 0
	master_table$hgnc_symbol[which(master_table$hgnc_symbol==0)] = 	master_table$Ensembl_Gene_ID[which(master_table$hgnc_symbol==0)]
	master_table$hgnc_symbol[which(master_table$hgnc_symbol=="")] = 	master_table$Ensembl_Gene_ID[which(master_table$hgnc_symbol=="")]


	# get only whitelist genes
	hla_genes = data.frame(fread(hla_genes_path, header=F))
	hla_genes = remove_period(toChar(hla_genes[,1]))

	fpkm_pass = data.frame(fread(fpkm_pass_path, header=F))
	fpkm_pass = remove_period(toChar(fpkm_pass[,1]))

	prt_genes_table = fread(prt_genes, header=F, sep="\t")
	prt_genes_table = data.frame(prt_genes_table, check.names=F)
	prt_genes_table[,1] = remove_period(prt_genes_table[,1])

	master_table = master_table[which(master_table$Ensembl_Gene_ID %in% fpkm_pass),]
	master_table = master_table[which(!master_table$Ensembl_Gene_ID %in% hla_genes),]
	master_table = master_table[which(master_table$Ensembl_Gene_ID %in% prt_genes_table[,1]),]
	master_table = master_table[which(!master_table$Cancer_Type %in% "ESAD"),]

	return(master_table)
}


collapse_and_sum_multi <- function(table_to_smash, col_smash_id){

	col_interest = which(! colnames(table_to_smash) %in% c("ICGC_DONOR_ID", "Ensembl_Gene_ID", "Cancer_Type", "hgnc_symbol", "gene_sample_pairs_curr") )

	DT <- as.data.table(table_to_smash)
	summed_table = DT[, lapply(.SD, sum), by = c("gene_sample_pairs_curr"), .SDcols = colnames(table_to_smash)[col_interest]]

	summed_table = data.frame(summed_table)

	id_table = unique(table_to_smash[,c("ICGC_DONOR_ID", "Ensembl_Gene_ID", "Cancer_Type", "hgnc_symbol", "gene_sample_pairs_curr")])
	merged_summed_table = merge(id_table, summed_table, by="gene_sample_pairs_curr")

	return(merged_summed_table)

}

merge_cancer_type <- function(master_table_CT){

	for(curr_samp in unique(master_table_CT$ICGC_DONOR_ID)){

		cancer_type_curr = unique(master_table_CT[master_table_CT$ICGC_DONOR_ID == curr_samp, "Cancer_Type"])
		cancer_type_curr = names(sort(sapply(cancer_type_curr, nchar), decreasing=T))[1]
		master_table_CT[master_table_CT$ICGC_DONOR_ID == curr_samp, "Cancer_Type"] = cancer_type_curr

	}

	return(master_table_CT)

}

make_pca_alt <- function(master_table, col_interest, outname){
	variants_tab = master_table[,c(col_interest)]
	variants_tab = unique(variants_tab)

	variants_matr = dcast(variants_tab, ICGC_DONOR_ID+Cancer_Type~Ensembl_Gene_ID)
	variants_matr[is.na(variants_matr)] = 0
	variants_matr_filtered = variants_matr[, 3:ncol(variants_matr)]
	variants_matr_filtered = variants_matr_filtered[ ,which(colSums(variants_matr_filtered) > 0) ]

	ir.pca <- prcomp(variants_matr_filtered,
			 center = T,
			 scale. = T,
			 rank.=2)

	pca_out = ir.pca$x
	pca_out = data.frame(pca_out[,1:3], Cancer_Type=variants_matr$Cancer_Type)
	 pdf(paste(base_dir, "/results/qc_files/qc_master_table/", col_interest[length(col_interest)], outname, ".pdf", sep=""))
	  print(ggplot(pca_out, aes(x=PC1, y=PC2, color=Cancer_Type)) + geom_point(shape=1))
	 dev.off()

}


collapse_and_sum <- function(table_to_smash, col_smash_id){

	col_interest = which(! colnames(table_to_smash) %in% c("ICGC_DONOR_ID", "Ensembl_Gene_ID", "Cancer_Type", "hgnc_symbol", "gene_sample_pairs_curr") )

	DT <- as.data.table(table_to_smash)
	summed_table = DT[, lapply(.SD, sum), by = c(col_smash_id), .SDcols = colnames(table_to_smash)[col_interest]]

	summed_table = data.frame(summed_table)

	return(summed_table)

}

plot_heatmap <- function(gene_ct_table, curr_alt, samp_table){
	curr_tab = gene_ct_table[,c("Ensembl_Gene_ID", "Cancer_Type", curr_alt)]
	val_matr = dcast(curr_tab, Ensembl_Gene_ID~Cancer_Type)
	val_matr[is.na(val_matr)] = 0
	row_val = rowSums(val_matr[,2:ncol(val_matr)])
	row_val_idx = which(row_val >= quantile(row_val, 0.5))
	val_matr = val_matr[row_val_idx,2:ncol(val_matr)]

	val_matr = sweep(val_matr, 2, samp_table, `/`)

	pdf(paste(base_dir, "/results/qc_files/qc_master_table/", curr_alt, "_heatmap_v.13.pdf", sep=""))
	gg = pheatmap(val_matr, main=paste(curr_alt, "num_genes ", nrow(val_matr), sep=""), scale="none", show_rownames=FALSE)
	print(gg)
	dev.off()

}

run_filter <- function(){
	# read in file
	master_table_raw = data.frame(fread(master_table_raw_path), check.names=F)
	col_order = c("ICGC_DONOR_ID", "Ensembl_Gene_ID", "Cancer_Type", "gene_sample_pairs_curr", "alt_prom", "alt_pa", "expr_outlier", "rna_edit", "variants", "ase_all", "fusion", "isSplice")
	master_table_raw = master_table_raw[,col_order]


	# binarize and validate
	col_to_binarize = c("isSplice", "expr_outlier")
	master_table_bin = binarize_table(master_table_raw, col_to_binarize)
	master_table = translate_validate_table(master_table_bin)

	# remove all duplicates and collapse if needed
	master_table = data.frame(unique(data.table(master_table)))
	master_table_collapse = collapse_and_sum_multi(master_table, "gene_sample_pairs_curr")

	# make sure binary after collapse
	master_table_numeric = master_table_collapse[,6:ncol(master_table_collapse)]
	master_table_numeric[master_table_numeric > 0] = 1
	master_table_collapse[,6:ncol(master_table_collapse)] = master_table_numeric




	sum(master_table_collapse[master_table_collapse$ICGC_DONOR_ID == "DO52691", "isSplice"]) == 6290
	master_table_collapse = data.frame(unique(data.table(master_table_collapse)))


	master_table_final = final_process_master_table(master_table_collapse)
	collapse_genes = collapse_and_sum(master_table_final, c("Ensembl_Gene_ID", "hgnc_symbol"))
	collapse_genes[collapse_genes$Ensembl_Gene_ID == "ENSG00000110321", ]

	# remove the imprinted genes
	imprint_table = data.frame(fread(imprint_path, header=F))
	colnames(imprint_table) = c("imprint_status", "Ensembl_Gene_ID")
	imprint_genes = imprint_table$Ensembl_Gene_ID[imprint_table$imprint_status == "I"]
	master_table_final$ase_all[which(master_table_final$Ensembl_Gene_ID %in% imprint_genes)] = 0

	#write it out
	write.table(master_table_final, master_table_path, sep="\t", row.name=F, col.names=T, quote=F)

	master_table_plot = master_table_final

	# now do QC
	make_pca_alt(master_table_plot, col_interest=c("ICGC_DONOR_ID", "Cancer_Type", "Ensembl_Gene_ID", "expr_outlier"), outname="_2_center_v.13")
	make_pca_alt(master_table_plot, col_interest=c("ICGC_DONOR_ID", "Cancer_Type", "Ensembl_Gene_ID", "isSplice"), outname="_2_center_v.13")
	make_pca_alt(master_table_plot, col_interest=c("ICGC_DONOR_ID", "Cancer_Type", "Ensembl_Gene_ID", "cn"), outname="_2_center_v.13")
	make_pca_alt(master_table_plot, col_interest=c("ICGC_DONOR_ID", "Cancer_Type", "Ensembl_Gene_ID", "alt_prom"), outname="_2_center_v.13")
	make_pca_alt(master_table_plot, col_interest=c("ICGC_DONOR_ID", "Cancer_Type", "Ensembl_Gene_ID", "alt_pa"), outname="_2_center_v.13")
	make_pca_alt(master_table_plot, col_interest=c("ICGC_DONOR_ID", "Cancer_Type", "Ensembl_Gene_ID", "rna_edit"), outname="_2_center_v.13")
	make_pca_alt(master_table_plot, col_interest=c("ICGC_DONOR_ID", "Cancer_Type", "Ensembl_Gene_ID", "fusion"), outname="_2_center_v.13")
	make_pca_alt(master_table_plot, col_interest=c("ICGC_DONOR_ID", "Cancer_Type", "Ensembl_Gene_ID", "ase_all"), outname="_2_center_v.13")
	make_pca_alt(master_table_plot, col_interest=c("ICGC_DONOR_ID", "Cancer_Type", "Ensembl_Gene_ID", "variants"), outname="_2_center_v.13")


	gene_ct_table = collapse_and_sum(master_table_plot, c("Ensembl_Gene_ID", "Cancer_Type"))

	samp_table = unique(master_table_plot[,c("ICGC_DONOR_ID", "Cancer_Type")])
	samp_table = table(samp_table$Cancer_Type)
	samp_table = samp_table[order(names(samp_table))]

	for(curr_alt in colnames(gene_ct_table)[3:12]){
		plot_heatmap(gene_ct_table, curr_alt, samp_table)
	}


	make_boxplot <- function(alt_table, outfile){

		plot_colors <- brewer.pal(ncol(alt_table),"Set3")
		pdf(outfile)
		boxplot((alt_table), las=2, main="Number of Samples with each alteration", col=plot_colors, ylim=c(0,700))
		dev.off()

	}

	collapse_genes = collapse_and_sum(master_table_plot, c("Ensembl_Gene_ID", "hgnc_symbol"))
	curr_alt_interest = c("alt_prom", "expr_outlier", "rna_edit", "variants", "ase_all", "fusion", "isSplice", "dom_trans")
	curr_alt_interest = c("alt_prom", "alt_pa", "expr_outlier", "rna_edit", "variants", "ase_all", "fusion", "isSplice", "cn")
	make_boxplot(collapse_genes[,curr_alt_interest], paste(base_dir, "/results/recurrance_analysis/alteration_boxplot_final_v.13.pdf", sep=""))

}
#run_filter()
