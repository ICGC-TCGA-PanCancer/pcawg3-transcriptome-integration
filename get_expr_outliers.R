
source("../set_path_variables.R")
source("../make_master_table.R")

require(data.table)
require(RColorBrewer)
library("DESeq2")

get_hdf5_counts <- function(RNA_sample_type_list, pan_file){
    #this method is used after the reviews
    
    full_file <- h5file(pan_file, 'r')
    
    sample_ids = full_file["/row_header"][1:1188]
    gene_ids = remove_period(full_file["/col_header"][1:18898])
    full_resid = full_file["/peer35"][1:1188, 1:18898]
    full_resid = t(full_resid)
    colnames(full_resid) = sample_ids[,2]
    rownames(full_resid) = gene_ids
    
    # filter samples
    full_resid_whitelist  = full_resid[,which(colnames(full_resid) %in% RNA_sample_type_list$Sample_ID)]
    
    return(full_resid_whitelist)
    
    
}

read_icgc_data <- function(RNA_sample_type_list, icgc_counts_file){
    # this method was used before the reviews
    icgc_table = data.frame(fread(icgc_counts_file, sep="\t", header=T), check.names=F)
    
    gene_names = toChar(icgc_table$feature)
    icgc_table_whitelist  = icgc_table[,which(colnames(icgc_table) %in% RNA_sample_type_list$Sample_ID)]
    rownames(icgc_table_whitelist) = gene_names
    
    return(icgc_table_whitelist)
}

get_scale_values <- function(x){
    if(sd(x) < 1){
        return(rep(0, length(x)))
    }else{
        return(scale(log10(x+1)))
    }
}

zscore_icgc_table_per_cancer <- function(curr_cancer_table, cancer_type){
    
    curr_cancer_table_scale = apply(curr_cancer_table, 1, get_scale_values)
    
    curr_cancer_table_scale = t(curr_cancer_table_scale)
    colnames(curr_cancer_table_scale) = colnames(curr_cancer_table)
    rownames(curr_cancer_table_scale) = rownames(curr_cancer_table)
    
    curr_cancer_table_scale = data.frame(gene_id=remove_period(row.names(curr_cancer_table)), curr_cancer_table_scale, check.names=F)
    melt_zscore_table = melt(curr_cancer_table_scale, id.vars="gene_id")
    melt_zscore_table = data.frame(melt_zscore_table, Cancer_Type=cancer_type)
    colnames(melt_zscore_table) = c("Ensembl_Gene_ID", "Sample_ID", "zscore", "Cancer_Type")
    
    curr_cancer_table = data.frame(gene_id=remove_period(row.names(curr_cancer_table)), curr_cancer_table, check.names=F)
    melt_expr_table = melt(curr_cancer_table, id.vars="gene_id")
    melt_expr_table = data.frame(melt_expr_table, Cancer_Type=cancer_type)
    
    melt_zscore_table = data.frame(melt_zscore_table, fpkm_uq=melt_expr_table$value, check.names=F)
    
    return(melt_zscore_table)
    
}




### read in arguments
args = commandArgs(trailingOnly=TRUE)
#icgc_counts_uq_file = args[1] # something like "joint_fpkm_uq.tsv"
pan_file = args[1] # something like "Pan.phenotype.hdf5"
RESULTS_FOLDER = args[2] # something like "Pan.phenotype.hdf5"


setwd(RESULTS_FOLDER)

RNA_sample_type_list = make_RNA_cancer_type_set()
icgc_table  = get_hdf5_counts(RNA_sample_type_list, pan_file)
cancer_types = unique(RNA_sample_type_list$Cancer_Type)

icgc_table_mini = icgc_table

all_res = list()
for(cancer_type in cancer_types){
    
    print(cancer_type)
    
    samples_used = RNA_sample_type_list[RNA_sample_type_list$Cancer_Type == cancer_type, c("Sample_ID")]
    idx_samp = which(colnames(icgc_table_mini) %in% toChar(samples_used))
    cancer_samp = icgc_table_mini[,idx_samp]
    
    zscore_cancer_samp = zscore_icgc_table_per_cancer(cancer_samp, cancer_type)
    
    temp_RNA_sample_type_list = RNA_sample_type_list[RNA_sample_type_list$Cancer_Type == cancer_type, ]
    zscore_cancer_samp_merge = merge(data.table(zscore_cancer_samp), data.table(temp_RNA_sample_type_list), by=c("Sample_ID", "Cancer_Type"))
    
    zscore_cancer_samp_merge = data.frame(zscore_cancer_samp_merge)
    all_res[[cancer_type]] = zscore_cancer_samp_merge
}


all_res_df <- do.call("rbind", all_res)

write.table(all_res_df, expr_outlier_table_path, sep="\t", row.name=F, col.names=T, quote=F)

