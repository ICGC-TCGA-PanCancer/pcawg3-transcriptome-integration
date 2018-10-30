source("../set_path_variables.R")
source("../process_table_functions.R")
source("../visualization_scripts.R")
require(reshape2)
require("RColorBrewer")



get_gene_lengths <- function(master_table){

	trans_len_bm = data.frame(fread(ensembl_gene_length_gc_path))
	colnames(trans_len_bm) = c('Ensembl_Gene_ID', "gc_content", 'transcript_length')
	trans_len_bm <- na.omit(trans_len_bm)

	#take the longest transcript
	trans_len_bm = trans_len_bm[order(trans_len_bm$transcript_length, decreasing=T),]
	trans_len_bm = trans_len_bm[!duplicated(trans_len_bm$Ensembl_Gene_ID), ]
	colnames(trans_len_bm)[1] = "Ensembl_Gene_ID"

	master_table_gl = join(master_table, trans_len_bm, by="Ensembl_Gene_ID", type="left")

	return(master_table_gl)

}



collapse_and_sum <- function(table_to_smash, col_smash_id){

	col_interest = which(! colnames(table_to_smash) %in% c("ICGC_DONOR_ID", "Ensembl_Gene_ID", "Cancer_Type", "hgnc_symbol", "gene_sample_pairs_curr") )

	DT <- as.data.table(table_to_smash)
	summed_table = DT[, lapply(.SD, sum), by = c(col_smash_id), .SDcols = colnames(table_to_smash)[col_interest]]

	summed_table = data.frame(summed_table)

	return(summed_table)

}


collapse_cancer_and_sum <- function(table_to_smash, col_smash_id){

	col_interest = which(! colnames(table_to_smash) %in% c("ICGC_DONOR_ID", "Ensembl_Gene_ID", "Cancer_Type", "hgnc_symbol", "gene_sample_pairs_curr") )

	DT <- as.data.table(table_to_smash)
	summed_table = DT[, lapply(.SD, sum), by = c(col_smash_id, "Cancer_Type"), .SDcols = colnames(table_to_smash)[col_interest]]

	summed_table = data.frame(summed_table)

	return(summed_table)

}

make_boxplot <- function(alt_table, outfile){

	plot_colors <- brewer.pal(ncol(alt_table),"Set3")
	pdf(outfile)
	boxplot((alt_table), las=2, main="Number of Samples with each alteration", col=plot_colors, ylim=c(0,700))
	dev.off()

}

make_counts_to_orders <- function(alt_table){

	for( curr_col in 1:ncol(alt_table)){

		# to normalize by the number of zero elements, we scale the rank based on the number of non-zero elements
		zero_idx = which(alt_table[,curr_col] == 0)
		alt_table[,curr_col] = rank(-as.numeric(alt_table[,curr_col]))
		alt_table[zero_idx,curr_col] = nrow(alt_table)

	}

	return(alt_table)
}


translate_ENS_TO_HGNC <-function(genes_to_trans, org="hsapiens_gene_ensembl", name_type='ensembl_gene_id'){
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

		genes_to_trans$hgnc_symbol[which(genes_to_trans$hgnc_symbol == "" | is.na(genes_to_trans$hgnc_symbol))] =
			genes_to_trans$ensembl_gene_id[which(genes_to_trans$hgnc_symbol == "" | is.na(genes_to_trans$hgnc_symbol))]

		return(genes_to_trans)
}

calc_cutoff_hypergeo <- function(total_genes, genes_pass, ranked_genes, scores_genes){

	# number of white balls in the urn
	pop_success = length(unique(genes_pass))

	# number of black balls in the urn
	total_space = unique(c(toChar(total_genes), genes_pass))
	print(head(total_space))
	print(head(genes_pass))

	pop_fail = sum(! toChar(total_space) %in% toChar(genes_pass))

	print("length of all genes")
	print(length(total_space))
	print("number of genes passing")
	print(pop_success)
	print("number of genes failing")
	print(pop_fail)


	# filter out the genes from the ranked list that don't pass permutation test
	percentile_pass = data.frame(fread(filter_score_path, header=T))
	score_pass = percentile_pass[percentile_pass$percentiles==5, "scores"]
	samp_size = sum(scores_genes <= score_pass)

	#samp_size = 1712

	samp_success = sum(ranked_genes[c(1:samp_size)] %in% genes_pass)
	print("CGC genes in enrich")
	print(ranked_genes[which(ranked_genes[c(1:samp_size)] %in% genes_pass)] )
	curr_padj = phyper(samp_success, pop_success, pop_fail, samp_size, lower.tail = FALSE)

	enrichment = (samp_success/samp_size) / (pop_success/length(total_genes))

	min_pval = c(enrichment, curr_padj, samp_success, samp_size)
	print(min_pval)
	print(c(samp_success, pop_success, pop_fail, samp_size))
	names(min_pval) = c("enrichment", "pval", "samp_success", "samp_size")

	min_pval = c(min_pval, ranked_genes[which(ranked_genes[c(1:samp_size)] %in% genes_pass)])

	return(min_pval)

}

calc_minimal_hypergeo <- function(total_genes, genes_pass, ranked_genes, dummy_var=T){

  max_genes = 1*length(unique(total_genes))

  # starting at 30, go 10 at a time until you hit max_genes
  tests = seq(10, max_genes, 10)

  pop_success = length(unique(genes_pass))
  pop_fail = length(total_genes) - pop_success

  all_pval = c()

  for(i in tests){
    samp_size = i
    samp_success = sum(ranked_genes[c(1:i)] %in% genes_pass)

    curr_padj = phyper(samp_success, pop_success, pop_fail, samp_size, lower.tail = FALSE)

    enrichment = (samp_success/samp_size) / (pop_success/length(total_genes))

    all_pval = rbind(all_pval, c(enrichment, curr_padj, samp_success, samp_size))
    colnames(all_pval) = c("enrichment", "pval", "samp_success", "samp_size")

  }
	all_pval = data.frame(all_pval)

	all_pval$padj = p.adjust(unlist(all_pval$pval), method="BH")

	min_pval = all_pval[which(all_pval$padj == min(all_pval$padj))[1], ]

	print(min_pval)
	print(length(tests))
	print(c(samp_success, pop_success, pop_fail, samp_size))

  return(min_pval)

}



calc_enrichment_cgc_driver <- function(collapse_genes, filename, write_file=T, run_minimal=F){

	cgc = data.frame(fread(cancer_census_genes_path, header=F), header=F)
	colnames(cgc) = "ensembl_gene_id"

	cgc_trans = translate_ENS_TO_HGNC(cgc)
	idx_cancer_genes = which(na.omit(collapse_genes$hgnc_symbol) %in% cgc_trans$hgnc_symbol)

	driver_genes = data.frame(fread(driver_gene_path), header=F)
	idx_driver_genes = which(na.omit(collapse_genes$hgnc_symbol) %in% driver_genes$gene)

	idx_cancer_driver_genes = unique(c(idx_cancer_genes, idx_driver_genes))
	all_cancer_driver_genes = unique(c(cgc_trans$hgnc_symbol, driver_genes$gene))

	if(run_minimal){
		min_pval = calc_minimal_hypergeo(na.omit(collapse_genes$hgnc_symbol), all_cancer_driver_genes, na.omit(collapse_genes$hgnc_symbol))
	}else{
		min_pval = calc_cutoff_hypergeo(na.omit(collapse_genes$hgnc_symbol), all_cancer_driver_genes, na.omit(collapse_genes$hgnc_symbol), collapse_genes[which(!is.na(collapse_genes$hgnc_symbol)),"sum_order"])
	}
	if(write_file){
		write.table(min_pval, filename, sep="\t", row.name=T, quote=F)
	}
	return(min_pval)

}


get_trans_table <- function(file_path){

	trans_tab = data.frame(fread(file_path))
	trans_tab = trans_tab$geneID[trans_tab$g_adj_pval < 0.1]
	trans_tab = data.frame(ensembl_gene_id = remove_period(trans_tab))
	trans_tab = translate_ENS_TO_HGNC(trans_tab)
	return(trans_tab$hgnc_symbol)

}

get_cis_table <- function(file_path){

	trans_tab = data.frame(fread(file_path))
	trans_tab = trans_tab$geneID[trans_tab$g_adj_pval < 0.05]
	trans_tab = data.frame(ensembl_gene_id = remove_period(trans_tab))
	trans_tab = translate_ENS_TO_HGNC(trans_tab)
	return(trans_tab$hgnc_symbol)

}

get_germ_table <- function(){

	trans_tab = data.frame(fread(germline_mut_path))
	trans_tab = data.frame(ensembl_gene_id = remove_period(trans_tab[,1]))
	trans_tab = translate_ENS_TO_HGNC(trans_tab)
	return(trans_tab$hgnc_symbol)

}


get_sv_table <- function(){

	trans_tab = data.frame(fread(sv_upreg_path))
	colnames(trans_tab) = c("ensembl_gene_id")
	trans_tab = translate_ENS_TO_HGNC(trans_tab)
	return(trans_tab$hgnc_symbol)

}


get_qtl_tables <- function(){

	pT = get_trans_table(promoter_trans_path)
	iT = get_trans_table(intron_trans_path)
	eT = get_trans_table(enhancer_trans_path)

	pC = get_cis_table(promoter_cis_path)
	iC = get_cis_table(intron_cis_path)
	eC = get_cis_table(enhancer_cis_path)
	germ = get_germ_table()
	sv = get_sv_table()

	return(list(pT=pT, iT=iT, eT=eT, pC=pC, iC=iC, eC=eC, g=germ, sv=sv))

}

calc_enrichment_qtl_inner <- function(hgnc_symbol_qtl, collapse_genes, filename, write_file){

	min_pval = calc_cutoff_hypergeo(na.omit(collapse_genes$hgnc_symbol), hgnc_symbol_qtl, na.omit(collapse_genes$hgnc_symbol), collapse_genes[which(!is.na(collapse_genes$hgnc_symbol)),"sum_order"])
	if(write_file){
		write.table(min_pval, filename, sep="\t", row.name=T, quote=F)
	}

}


calc_enrichment_qtl <- function(collapse_genes, filename, write_file=T){

	qtl_tables = get_qtl_tables()
	calc_enrichment_qtl_inner(qtl_tables$pT, collapse_genes, paste(filename, "_pT", sep=""), write_file)
	calc_enrichment_qtl_inner(qtl_tables$iT, collapse_genes, paste(filename, "_iT", sep=""), write_file)
	calc_enrichment_qtl_inner(qtl_tables$eT, collapse_genes, paste(filename, "_eT", sep=""), write_file)

	calc_enrichment_qtl_inner(qtl_tables$pC, collapse_genes, paste(filename, "_pC", sep=""), write_file)
	calc_enrichment_qtl_inner(qtl_tables$iC, collapse_genes, paste(filename, "_iC", sep=""), write_file)
	calc_enrichment_qtl_inner(qtl_tables$eC, collapse_genes, paste(filename, "_eC", sep=""), write_file)

	calc_enrichment_qtl_inner(qtl_tables$g, collapse_genes, paste(filename, "_g", sep=""), write_file)
	calc_enrichment_qtl_inner(qtl_tables$sv, collapse_genes, paste(filename, "_sv", sep=""), write_file)

	calc_enrichment_qtl_inner(unique(unlist(qtl_tables)), collapse_genes, paste(filename, "_TOT", sep=""), write_file)
}




calc_enrichment_cgc <- function(collapse_genes, filename, write_file=T){

	cgc = data.frame(fread(cancer_census_genes_path, header=F), header=F)
	colnames(cgc) = "ensembl_gene_id"

	cgc_trans = translate_ENS_TO_HGNC(cgc)

	idx_cancer_genes = which(na.omit(collapse_genes$hgnc_symbol) %in% cgc_trans$hgnc_symbol)

	min_pval = calc_cutoff_hypergeo(na.omit(collapse_genes$hgnc_symbol), cgc_trans$hgnc_symbol, na.omit(collapse_genes$hgnc_symbol), collapse_genes[which(!is.na(collapse_genes$hgnc_symbol)),"sum_order"])
	if(write_file){
		write.table(min_pval, filename, sep="\t", row.name=T, quote=F)
	}
	return(min_pval)


}


calc_enrichment_drivers <- function(collapse_genes, filename, write_file=T){

	driver_genes = data.frame(fread(driver_gene_path), header=F)

	idx_driver_genes = which(na.omit(collapse_genes$hgnc_symbol) %in% driver_genes$gene)

	min_pval = calc_cutoff_hypergeo(na.omit(collapse_genes$hgnc_symbol), driver_genes$gene, na.omit(collapse_genes$hgnc_symbol), collapse_genes[which(!is.na(collapse_genes$hgnc_symbol)),"sum_order"])
	if(write_file){
		write.table(min_pval, filename, sep="\t", row.name=T, quote=F)
	}
	return(min_pval)


}

write_driver_genes <-function(collapse_genes, filename){

	driver_genes = data.frame(fread(driver_gene_path), header=F)
	idx_driver_genes = which(na.omit(collapse_genes$hgnc_symbol) %in% driver_genes$gene)

	write.table(collapse_genes$hgnc_symbol[idx_driver_genes], filename, sep="\t", row.name=T, quote=F)


}


write_driver_cgc_genes <-function(collapse_genes, filename){

	driver_genes = data.frame(fread(driver_gene_path), header=F)
	idx_driver_genes = which(na.omit(collapse_genes$hgnc_symbol) %in% driver_genes$gene)

	cgc = data.frame(fread(cancer_census_genes_path, header=F), header=F)
	colnames(cgc) = "ensembl_gene_id"
	cgc_trans = translate_ENS_TO_HGNC(cgc)
	idx_cgc_genes = which(na.omit(collapse_genes$hgnc_symbol) %in% cgc_trans$hgnc_symbol)


	idx_cgc_driver_genes = unique(c(idx_cgc_genes, idx_driver_genes))
	idx_cgc_driver_genes = sort(idx_cgc_driver_genes)
	write.table(collapse_genes$hgnc_symbol[idx_cgc_driver_genes], filename, sep="\t", row.name=T, quote=F)

}

plot_alt_distr <- function(rank_matr, value_matr, col_interest, outdir){

	pallete = brewer.pal(length(col_interest),"Dark2")
	#pdf(paste(outdir, "/alteration_distr.pdf", sep=""))
	#barplot(table(rank_matr$alt_chosen[1:100])[col_interest], col=pallete)
	#dev.off()

	pdf(paste(outdir, "/rank_distr.pdf", sep=""))
	rank_matr_values = rank_matr[,col_interest]

	for(col_idx in 1:ncol(rank_matr_values)){
		if(col_idx==1){
			plot(density(rank_matr_values[,col_idx]), col=pallete[col_idx], main="Distr of # ranks")
		}else{
			lines(density(rank_matr_values[,col_idx]), col=pallete[col_idx])
		}
	}
	legend("topleft", legend = c(col_interest), col=pallete, title = "Alteration Types", pch=15)
	dev.off()

	pdf(paste(outdir, "/value_distr.pdf", sep=""))
	value_matr = value_matr[,col_interest]
	for(col_idx in 1:ncol(value_matr)){
		if(col_idx==1){
			plot(density(value_matr[,col_idx]), col=pallete[col_idx], main="Distr of # samples altered")
		}else{
			lines(density(value_matr[,col_idx]), col=pallete[col_idx])
		}
	}
	legend("topleft", legend = c(col_interest), col=pallete, title = "Alteration Types", pch=15)
	dev.off()

	pdf(paste(outdir, "/value_hist.pdf", sep=""))
	for(col_idx in 1:ncol(value_matr)){
		curr_val = value_matr[,col_idx]
		hist(curr_val, col=pallete[col_idx], main=paste(colnames(value_matr)[col_idx], "Distr of # samples altered"), breaks=length(curr_val)/100)
	}
	dev.off()


}

run_rank_method <- function(ordered_alts=NA, collapse_genes, curr_alt_interest, master_table, outname, order_decreasing){

	res_dir = paste(base_dir, "/results/recurrance_analysis/no_dom_trans/", outname, "/", sep="")
	dir.create(res_dir, showWarnings=F)

	collapse_genes_ord = collapse_genes[order(collapse_genes$sum_order, decreasing=order_decreasing),]
	write.table(collapse_genes_ord, paste(res_dir, outname, ".txt", sep=""), sep="\t", row.name=F, quote=F)

	if(!is.na(ordered_alts)){
		ordered_alts = ordered_alts[order(collapse_genes$sum_order, decreasing=order_decreasing),]
		write.table(ordered_alts, paste(res_dir, outname, "_rank.txt", sep=""), sep="\t", row.name=F, quote=F)

		if(length(curr_alt_interest) > 1){
			rank_table = ordered_alts[,curr_alt_interest]
			colnames(rank_table) = paste(colnames(rank_table), "_rank", sep="")
			full_table = cbind(collapse_genes_ord[,c(1,2)], collapse_genes_ord[,c("alt_chosen", "sum_order")], rank_table, collapse_genes_ord[,curr_alt_interest])
			write.table(full_table, paste(res_dir, outname, "_full.txt", sep=""), sep="\t", row.name=F, quote=F)

		}

	}



	#calc_enrichment_qtl(collapse_genes_ord, filename=paste(res_dir, outname, "_qtl_enrichment.txt", sep=""))

	print("---------------------------------")
	print(curr_alt_interest)
	print("cgc")
	calc_enrichment_cgc(collapse_genes_ord, filename=paste(res_dir, outname, "_cgc_enrichment.txt", sep=""))
	print("drivers")
	calc_enrichment_drivers(collapse_genes_ord, filename=paste(res_dir, outname, "_drivers_enrichment.txt", sep=""))
	print("cgc_drivers")
	calc_enrichment_cgc_driver(collapse_genes_ord, filename=paste(res_dir, outname, "_cgc_drivers_enrichment.txt", sep=""))
	print("---------------------------------")

	write_driver_genes(collapse_genes_ord, filename=paste(res_dir, outname, "_drivers.txt", sep=""))
	write_driver_cgc_genes(collapse_genes_ord, filename=paste(res_dir, outname, "_cgc_drivers.txt", sep=""))

	if(!is.na(ordered_alts)){
		plot_alt_distr(ordered_alts, collapse_genes_ord, curr_alt_interest, res_dir)
	}else{
		plot_alt_distr(collapse_genes_ord, collapse_genes_ord, curr_alt_interest, res_dir)
	}

	# write out relevent cancer types
	collapse_genes_cancer = collapse_cancer_and_sum(master_table, "hgnc_symbol")
	collapse_genes_cancer = collapse_genes_cancer[which(collapse_genes_cancer$hgnc_symbol %in% collapse_genes_ord$hgnc_symbol[1:35]),]
	collapse_genes_cancer = merge(collapse_genes_cancer, collapse_genes_ord[,c("hgnc_symbol", "alt_chosen", "sum_order")], by="hgnc_symbol")
	#colnames(collapse_genes_cancer)[13] = "rank"
	write.table(collapse_genes_cancer, paste(res_dir, outname, "_cancer_info.txt", sep=""), sep="\t", row.name=F, quote=F)



}

shuffle_rows_each_col <- function(curr_table){

	for(curr_col in 1:ncol(curr_table)){

		curr_table[,curr_col] = curr_table[sample(nrow(curr_table)), curr_col]
	}

	return(curr_table)
}

do_shuffle_qc <- function(table_to_run, num_shuffle, curr_alt_interest){

	all_pval_cgc = c()
	all_pval_driver = c()
	all_pval_cgc_driver = c()

	for(idx in 1:num_shuffle){
		table_to_run[,curr_alt_interest] = shuffle_rows_each_col(table_to_run[,curr_alt_interest])

		all_res = do_basic_qc(table_to_run, curr_alt_interest)
		all_pval_cgc = rbind(all_pval_cgc, all_res[[1]][1:4])

		all_pval_driver = rbind(all_pval_driver, all_res[[2]][1:4])

		all_pval_cgc_driver = rbind(all_pval_cgc_driver, all_res[[3]][1:4])
		print(idx)
		print(" ")
	}

	return(list(cgc_pval=all_pval_cgc, driver_pval=all_pval_driver, cgc_driver_pval=all_pval_cgc_driver))
}

do_basic_qc <- function(table_to_run, curr_alt_interest){

	collapse_genes = collapse_and_sum(table_to_run, "hgnc_symbol")

	# 1-> second smallest rank
	collapse_genes = collapse_and_sum(table_to_run, "hgnc_symbol")
	ordered_alts = make_counts_to_orders(collapse_genes[,curr_alt_interest])
	collapse_genes = data.frame(collapse_genes, alt_chosen=rep("NA", nrow(ordered_alts)), stringsAsFactors=FALSE)
	for( curr_row in 1:nrow(ordered_alts) ){

		collapse_genes$sum_order[curr_row] = unlist(sort(ordered_alts[curr_row,],partial=2)[2])
		collapse_genes$alt_chosen[curr_row] = unlist(colnames(ordered_alts)[order(ordered_alts[curr_row,])[2]])

	}
	ordered_alts = data.frame(hgnc_symbol=collapse_genes$hgnc_symbol, ordered_alts, alt_chosen=collapse_genes$alt_chosen)

	collapse_genes_ord = collapse_genes[order(collapse_genes$sum_order, decreasing=F),]

	if(!is.na(ordered_alts)){
		ordered_alts = ordered_alts[order(collapse_genes$sum_order, decreasing=F),]
	}

	min_pval_cgc = calc_enrichment_cgc(collapse_genes_ord, filename="NA", write_file=F)

	min_pval_driver = calc_enrichment_drivers(collapse_genes_ord, filename="NA", write_file=F)

	min_pval_cgc_driver = calc_enrichment_cgc_driver(collapse_genes_ord, filename="NA", write_file=F)

	return(list(cgc_pval=min_pval_cgc, driver_pval=min_pval_driver, cgc_driver_pval=min_pval_cgc_driver))


}

num_bubble_swap <-function(sorted_ranks_table, curr_alt_interest){

	searched_table = sorted_ranks_table[1:500,c(curr_alt_interest)]

	# calculate for each column
	col_sums = c()
	for(curr_col in 1:ncol(searched_table)){
		col_sum = 0
		# calculate the number of swaps needed
		for(outer_idx in 1:(nrow(searched_table)-1)){
			for(inner_idx in (outer_idx+1):nrow(searched_table)){
				if(searched_table[outer_idx,curr_col] < searched_table[inner_idx,curr_col]){
					col_sum = col_sum + 1
				}
			}
		}
		col_sums = c(col_sums, col_sum)
		print(curr_col)
	}
	return(col_sums)
}

do_qc <- function(){

	res_dir = paste(base_dir, "/results/qc_files/qc_on_ranking/", sep="")

	master_table = data.frame(fread(master_table_path), check.names=F)
	curr_alt_interest = c("alt_prom", "expr_outlier", "rna_edit", "variants", "ase_all", "fusion", "isSplice")

	# shuffle and run on entire table
	full_shuffle_res = do_shuffle_qc(master_table, num_shuffle=10, curr_alt_interest)
	write.table(full_shuffle_res, paste(res_dir, "shuffle_enrichment.tsv", sep=""), sep="\t", quote=F)

	# hold out cancer_type
	cancer_types = unique(master_table$Cancer_Type)
	all_cancer_res = c()
	all_genes_found_cancer = list()
	for(curr_ct in cancer_types){

		print(curr_ct)
		minus_cancer_type = master_table[which(!master_table$Cancer_Type == curr_ct),]
		curr_res = do_basic_qc(minus_cancer_type, curr_alt_interest)

		genes_found = toChar(curr_res[[3]][5:length(curr_res[[3]])])

		cancer_res = c(curr_res[[1]][1:4], curr_res[[2]][1:4], curr_res[[3]][1:4])
		names(cancer_res) = c(paste("cgc", names(cancer_res)[1:4], sep="_"), paste("driver", names(cancer_res)[1:4], sep="_"), paste("cgc_driver", names(cancer_res)[1:4], sep="_") )

		all_genes_found_cancer[curr_ct] = list(genes_found)
		all_cancer_res = rbind(all_cancer_res, cancer_res)
	}
	rownames(all_cancer_res) = cancer_types
	write.table(all_cancer_res, paste(res_dir, "cancer_enrichment.tsv", sep=""), sep="\t", quote=F)
	write(cancer_types, paste(res_dir, "cancer_genes_enrich.tsv", sep=""), append=FALSE)
	lapply(all_genes_found_cancer, write, paste(res_dir, "cancer_genes_enrich.tsv", sep=""), append=TRUE, ncolumns=1000)

	# hold out alt_type
	all_alt_res = c()
	all_genes_found_alt = list()
	for(curr_alt in curr_alt_interest){

		print(curr_alt)
		curr_alt_interest_minus = curr_alt_interest[which(!curr_alt_interest == curr_alt)]
		curr_res = do_basic_qc(master_table, curr_alt_interest_minus)

		genes_found = toChar(curr_res[[3]][5:length(curr_res[[3]])])

		alt_res = c(curr_res[[1]][1:4], curr_res[[2]][1:4], curr_res[[3]][1:4])
		names(alt_res) = c(paste("cgc", names(alt_res)[1:4], sep="_"), paste("driver", names(alt_res)[1:4], sep="_"), paste("cgc_driver", names(alt_res)[1:4], sep="_") )

		all_genes_found_alt[curr_alt] = list(genes_found)
		all_alt_res = rbind(all_alt_res, alt_res)
	}
	rownames(all_alt_res) = curr_alt_interest
	write.table(all_alt_res, paste(res_dir, "alt_enrichment.tsv", sep=""), sep="\t", quote=F)
	write(curr_alt_interest, paste(res_dir, "alt_genes_enrich.tsv", sep=""), append=FALSE)
	lapply(all_genes_found_alt, write, paste(res_dir, "alt_genes_enrich.tsv", sep=""), append=TRUE, ncolumns=1000)


}

qc_each_table <- function(master_table, curr_alt_interest){

	collapse_genes = collapse_and_sum(master_table, c("Ensembl_Gene_ID", "hgnc_symbol"))

	for(alt_event in curr_alt_interest){
		print(alt_event)
		# get the enrichment for each alteration type
		ordered_alts = make_counts_to_orders(collapse_genes[,colnames(collapse_genes)[3:ncol(collapse_genes)]])
		collapse_genes = data.frame(collapse_genes, alt_chosen=alt_event, stringsAsFactors=FALSE)
		collapse_genes$sum_order = ordered_alts[,alt_event]
		collapse_genes$alt_chosen = alt_event
		ordered_alts = data.frame(hgnc_symbol=collapse_genes$hgnc_symbol, ordered_alts, alt_chosen=collapse_genes$alt_chosen)

		run_rank_method(ordered_alts, collapse_genes, curr_alt_interest, master_table, paste("secondSmallest_full_v.13", alt_event), order_decreasing=F)

	}

}

correct_variants <- function(master_table_sum){

	master_table_gl = get_gene_lengths(master_table_sum)

	if(any(master_table_gl$Ensembl_Gene_ID != master_table_sum$Ensembl_Gene_ID)){
		print("ERROR NOT MATCHING!!!!!")
	}

	master_table_sum$variants = (master_table_sum$variants)/((master_table_gl$transcript_length+1)/1000)

	return(master_table_sum)

}

run_ranking_method <- function(){

	master_table = data.frame(fread(master_table_path), check.names=F)
	dim(master_table)

	# our alterations of interest
	curr_alt_interest = c("alt_prom", "expr_outlier", "rna_edit", "variants", "ase_all", "fusion", "isSplice", "cn")


	# get only whitelist genes
	hla_genes = data.frame(fread(hla_genes_path, header=F))
	hla_genes = remove_period(toChar(hla_genes[,1]))

	fpkm_pass = data.frame(fread(fpkm_pass_path, header=F))
	fpkm_pass = remove_period(toChar(fpkm_pass[,1]))


	master_table = master_table[which(master_table$Ensembl_Gene_ID %in% fpkm_pass),]
	master_table = master_table[which(!master_table$Ensembl_Gene_ID %in% hla_genes),]
	master_table = master_table[,colnames(master_table) != "alt_pa"]
	non_zero = which(!rowSums(master_table[,curr_alt_interest]) == 0)
	master_table = master_table[non_zero,]

	dim(master_table)


	#####   run on all genes  with var   #####

	collapse_genes = collapse_and_sum(master_table, c("Ensembl_Gene_ID", "hgnc_symbol"))
	make_boxplot(collapse_genes[,curr_alt_interest], paste(base_dir, "/results/recurrance_analysis/alteration_boxplot_final_v.13.pdf", sep=""))


	qc_each_table(master_table, curr_alt_interest)


	# RANK METHODS

	# 1-> second smallest rank
	collapse_genes = collapse_and_sum(master_table, c("Ensembl_Gene_ID", "hgnc_symbol"))
	collapse_genes = correct_variants(collapse_genes)

	ordered_alts = make_counts_to_orders(collapse_genes[,curr_alt_interest])
	collapse_genes = data.frame(collapse_genes, alt_chosen=rep("NA", nrow(ordered_alts)), stringsAsFactors=FALSE)
	for( curr_row in 1:nrow(ordered_alts) ){

		collapse_genes$sum_order[curr_row] = unlist(sort(ordered_alts[curr_row,],partial=2)[2])
		collapse_genes$alt_chosen[curr_row] = unlist(colnames(ordered_alts)[order(ordered_alts[curr_row,])[2]])

	}
	ordered_alts = data.frame(hgnc_symbol=collapse_genes$hgnc_symbol, ordered_alts, alt_chosen=collapse_genes$alt_chosen)

	run_rank_method(ordered_alts, collapse_genes, curr_alt_interest, master_table, "secondSmallest_full_v.13", order_decreasing=F)


	# 1-> second smallest rank only drivers
	driver_genes = data.frame(fread(driver_gene_path))$gene
	master_table_gene = master_table[master_table$hgnc_symbol %in% driver_genes, ]
	collapse_genes = collapse_and_sum(master_table_gene, c("Ensembl_Gene_ID", "hgnc_symbol"))
	collapse_genes = correct_variants(collapse_genes)
	ordered_alts = make_counts_to_orders(collapse_genes[,curr_alt_interest])
	collapse_genes = data.frame(collapse_genes, alt_chosen=rep("NA", nrow(ordered_alts)), stringsAsFactors=FALSE)
	for( curr_row in 1:nrow(ordered_alts) ){

		collapse_genes$sum_order[curr_row] = unlist(sort(ordered_alts[curr_row,],partial=2)[2])
		collapse_genes$alt_chosen[curr_row] = unlist(colnames(ordered_alts)[order(ordered_alts[curr_row,])[2]])

	}
	ordered_alts = data.frame(hgnc_symbol=collapse_genes$hgnc_symbol, ordered_alts, alt_chosen=collapse_genes$alt_chosen)

	run_rank_method(ordered_alts, collapse_genes, curr_alt_interest, master_table_gene, "secondSmallest_drivers_v.13", order_decreasing=F)




}

#run_ranking_method()
