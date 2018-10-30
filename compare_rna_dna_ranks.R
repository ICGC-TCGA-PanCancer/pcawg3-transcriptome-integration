# plot DNA vs RNA rank differences

## this is the code for figure 10 generation

source("../set_path_variables.R")
source("../process_table_functions.R")
source("../recurrance_histogram.R")
source("../visualization_scripts.R")
require(reshape2)
require("RColorBrewer")
library(GGally)
library(ggrepel)
library(gridExtra)
library(wesanderson)

cgc = data.frame(fread(cancer_census_genes_path, header=F), header=F)
colnames(cgc) = "ensembl_gene_id"

cgc_trans = translate_ENS_TO_HGNC(cgc)
cgc_trans = cgc_trans[,1:2]
colnames(cgc_trans)[1] = "Ensembl_Gene_ID"

get_cancer_genes <- function(collapse_genes){


  	driver_genes = data.frame(fread(driver_gene_path), header=F)

    collapse_genes$cancer_genes = "normal"
    collapse_genes$cancer_genes[collapse_genes$hgnc_symbol %in% cgc_trans$hgnc_symbol] = "cgc"

    collapse_genes$cancer_genes[collapse_genes$hgnc_symbol %in% driver_genes$gene] = "driver"

    genes_both = intersect(driver_genes$gene, cgc_trans$hgnc_symbol)

    collapse_genes$cancer_genes[collapse_genes$hgnc_symbol %in% genes_both] = "both"

    return(collapse_genes)
}


get_master_table <- function(curr_alt_interest){

  	master_table = data.frame(fread(master_table_path), check.names=F)

  	dim(master_table)

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


  	return(master_table)
}

get_ranks <- function(master_table, rna_interest, dna_interest){

  	#####   run on all genes  with var   #####

  	collapse_genes = collapse_and_sum(master_table, c("Ensembl_Gene_ID", "hgnc_symbol"))


  	# RANK METHODS

  	# 1-> second smallest rank
    curr_alt_interest = c(rna_interest, dna_interest)
  	collapse_genes = collapse_and_sum(master_table, c("Ensembl_Gene_ID", "hgnc_symbol"))
  	ordered_alts = make_counts_to_orders(collapse_genes[,c(rna_interest, dna_interest)])
  	collapse_genes = data.frame(collapse_genes, alt_chosen=rep("NA", nrow(ordered_alts)), stringsAsFactors=FALSE)
  	for( curr_row in 1:nrow(ordered_alts) ){

  		collapse_genes$rna_rank[curr_row] = unlist(sort(ordered_alts[curr_row,rna_interest],partial=2)[2])
  		collapse_genes$rna_alt_chosen[curr_row] = unlist(colnames(ordered_alts[,rna_interest])[order(ordered_alts[curr_row,rna_interest])[2]])

      collapse_genes$dna_rank[curr_row] = unlist(sort(ordered_alts[curr_row,dna_interest],partial=2)[2])
  		collapse_genes$dna_alt_chosen[curr_row] = unlist(colnames(ordered_alts[,dna_interest])[order(ordered_alts[curr_row,dna_interest])[2]])

      collapse_genes$rank[curr_row] = unlist(sort(ordered_alts[curr_row,curr_alt_interest],partial=2)[2])
  		collapse_genes$alt_chosen[curr_row] = unlist(colnames(ordered_alts[,curr_alt_interest])[order(ordered_alts[curr_row,curr_alt_interest])[2]])

  	}
    ordered_alts$rank = collapse_genes$rank
    ordered_alts$hgnc_symbol = collapse_genes$hgnc_symbol


    return(list(collapse_genes, ordered_alts))

}

plot_ranks <- function(plot_df, outdir, outname){

  plot_df$rna_contribution = (plot_df$dna_rank - plot_df$rank)
  #plot_df$rna_contribution[plot_df$rna_contribution < 0] = 0
  plot_df$dna_contribution = (plot_df$rna_rank - plot_df$rank)
  #plot_df$dna_contribution[plot_df$dna_contribution < 0] = 0

    pdf(paste(outdir, "/", outname, ".pdf", sep=""))

      # regular scatter plots
      gg = ggplot(plot_df, aes(16801-rna_rank, 16801-dna_rank, color=cancer_genes, alpha=rank)) +
        geom_point()+ geom_abline() +
        geom_text_repel(data=plot_df[plot_df$rank < 500 &
          plot_df$cancer_genes != "normal", ],
          aes(label=hgnc_symbol), force=5,
          segment.color = '#cccccc', max.iter=3e2)
      print(gg)
      gg = ggplot(plot_df, aes(16801-rank, 16801-dna_rank, color=cancer_genes, alpha=rank)) +
        geom_point()+ geom_abline() +
        geom_text_repel(data=plot_df[plot_df$rank < 500 &
          plot_df$cancer_genes != "normal", ],
          aes(label=hgnc_symbol), force=5,
          segment.color = '#cccccc', max.iter=3e2)
      print(gg)

      gg = ggplot(plot_df[plot_df$cancer_genes !="normal",], aes(rna_contribution, dna_contribution, color=cancer_genes)) +
        geom_point()+ geom_abline() +
        geom_text_repel(data=plot_df[plot_df$rank < 500 &
          plot_df$cancer_genes != "normal", ],
          aes(label=hgnc_symbol), force=5,
          segment.color = '#cccccc', max.iter=3e2)
      print(gg)


      gg = ggplot(plot_df[plot_df$rank < 1005 & plot_df$cancer_genes!="normal",], aes(rank, dna_rank, color=cancer_genes)) +
        geom_point()+ geom_abline() +
        geom_text_repel(data=plot_df[plot_df$rank < 500 &
          plot_df$cancer_genes != "normal", ],
          aes(label=hgnc_symbol), force=5,
          segment.color = '#cccccc', max.iter=3e2)
      print(gg)

      gg = ggplot(plot_df[plot_df$rank < 1005 & plot_df$cancer_genes!="normal",], aes(rank, rna_rank, color=cancer_genes)) +
        geom_point()+ geom_abline() +
        geom_text_repel(data=plot_df[plot_df$rank < 500 &
          plot_df$cancer_genes != "normal", ],
          aes(label=hgnc_symbol), force=5,
          segment.color = '#cccccc', max.iter=3e2)
      print(gg)

      melt_df = plot_df[order(plot_df$rank),]
      melt_df = melt_df[1:50,]
      melt_df$plot_rank = 1:nrow(melt_df)
      melt_df$hgnc_symbol <- factor(melt_df$hgnc_symbol, melt_df$hgnc_symbol[1:nrow(melt_df)])
      melt_plot = melt(melt_df[,c("plot_rank", "rna_rank", "dna_rank", "rank")], id="plot_rank")
      gg = ggplot(melt_plot, aes(plot_rank, value, color=variable)) +
        geom_line()+ geom_abline()
        print(gg)

    dev.off()

}


make_pie_charts <- function(outname, outdir, ordered_alts_melt, collapse_genes_rank, ngenes, rna_interest, dna_interest){

  curr_alt_interest = c(rna_interest, dna_interest)

  total_cgc = collapse_genes_rank$cancer_genes
  total_cgc = collapse_genes_rank$cancer_genes
  ordered_alts_melt = ordered_alts_melt[1:ngenes,]
  collapse_genes_rank = collapse_genes_rank[1:ngenes,]

  # now go through each row and get the scaled versions of top two alterations
  top_two = c()
  for( curr_row in 1:nrow(ordered_alts_melt) ){

    row_order = unlist(colnames(ordered_alts_melt[,curr_alt_interest])[order(ordered_alts_melt[curr_row,curr_alt_interest], decreasing=T)])
    top_two = c(top_two, row_order[1:2])
  }

  top_two_rna_dna = top_two
  top_two_rna_dna[top_two_rna_dna %in% rna_interest] = "rna"
  top_two_rna_dna[top_two_rna_dna %in% dna_interest] = "dna"
  top_two_df = data.frame(alt = top_two, rna_dna=top_two_rna_dna)
  top_two_df$alt = factor(top_two_df$alt, levels(top_two_df$alt)[c(1, 2, 3, 5, 6, 7, 8, 4, 9)])
  top_two_df$rna_dna = factor(top_two_df$rna_dna, levels(top_two_df$rna_dna)[c(2, 1)])

  print(table(top_two_df$rna_dna))

  pdf(paste(outdir, "/", outname, "_pie_charts_top2_alt.pdf", sep=""))
    # regular scatter plots
    gg1 = ggplot(top_two_df, aes(x="", fill=alt)) +
      geom_bar(width=1)  +
      theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
           axis.text.x=element_blank(), axis.text.y=element_blank(),
           axis.title.x=element_blank(), axis.title.y=element_blank()) +
           scale_fill_brewer(palette="Set1") +
      coord_polar(theta="y")
    print(gg1)
    gg2 = ggplot(top_two_df, aes(x="", fill=rna_dna)) +
      geom_bar(width=1)  +
      theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
           axis.text.x=element_blank(), axis.text.y=element_blank(),
           axis.title.x=element_blank(), axis.title.y=element_blank()) +
           scale_fill_manual(values=c("#009E73", "#F0E442")) +
      coord_polar(theta="y")
    print(gg2)

    gg = grid.arrange(gg1, gg2, ncol=2, nrow=1, widths=c(1,1), heights=c(1))
    print(gg)

    cgc_df = data.frame(genes="overall_genes", cancer_genes = total_cgc)
    cgc_df$cancer_genes = factor(cgc_df$cancer_genes,
      c("normal", "cgc", "driver", "both"))

    cgc_df_sig = data.frame(genes="sig_genes", cancer_genes = collapse_genes_rank$cancer_genes)
    cgc_df = rbind(cgc_df, cgc_df_sig)
    gg = ggplot(cgc_df, aes(x=genes, y=cancer_genes, fill=cancer_genes)) +
      geom_bar(stat="identity", position="fill", width=0.5) +
      theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
           axis.text.x=element_blank(), axis.text.y=element_blank(),
           axis.title.y=element_blank()) +
           scale_fill_brewer(palette="YlOrRd")
    print(gg)

  dev.off()

}

make_volcano_plot <- function(outname, outdir, collapse_genes_rank, ngenes, rna_interest, dna_interest){


  collapse_genes_rank = collapse_genes_rank[1:ngenes]

  collapse_genes_rank$pval = collapse_genes_rank$rank / 16801
  collapse_genes_rank$log2FC = log2((rowMeans(collapse_genes_rank[,rna_interest])+1) / (rowMeans(collapse_genes_rank[,dna_interest])+1))
  collapse_genes_rank$sig = factor("not_sig", levels =
                              c("abs(log2FC) < 1 & pval < 0.05",
                              "abs(log2FC) > 1 & pval > 0.05",
                              "not_sig",
                              "abs(log2FC) > 1 & pval < 0.05"))
  collapse_genes_rank$sig[abs(collapse_genes_rank$log2FC) > 1 & collapse_genes_rank$pval < 0.05] =  "abs(log2FC) > 1 & pval < 0.05"
  collapse_genes_rank$sig[abs(collapse_genes_rank$log2FC) < 1 & collapse_genes_rank$pval < 0.05] =  "abs(log2FC) < 1 & pval < 0.05"
  collapse_genes_rank$sig[abs(collapse_genes_rank$log2FC) > 1 & collapse_genes_rank$pval > 0.05] =  "abs(log2FC) > 1 & pval > 0.05"
  collapse_genes_rank$sig[abs(collapse_genes_rank$log2FC) < 1 & collapse_genes_rank$pval > 0.05] =  "not_sig"

  # Make a basic volcano plot
  gg = ggplot(collapse_genes_rank, aes(log2FC, -log10(pval), color=sig)) +
        geom_point() +
        geom_text_repel(data=collapse_genes_rank[collapse_genes_rank$cancer_genes != "normal" & collapse_genes_rank$rank < 320,], aes(label=hgnc_symbol))

  pdf(paste(outdir, "/", outname, "_volcano.pdf", sep=""))
    print(gg)
  dev.off()


}

make_curve_plot <- function(outname, outdir, ordered_alts_melt, collapse_genes_rank, ngenes, rna_interest, dna_interest){


    total_cgc = collapse_genes_rank$cancer_genes
    collapse_genes_rank = collapse_genes_rank[1:ngenes,]

    dna_sum = rowSums(collapse_genes_rank[,dna_interest])+1
    rna_sum = rowSums(collapse_genes_rank[,rna_interest])+1

    dna_prop = dna_sum/(dna_sum+rna_sum)
    rna_prop = rna_sum/(dna_sum+rna_sum)

    collapse_genes_rank_rna_dna = data.frame(rank=1:nrow(collapse_genes_rank), dna_prop,
      rna_prop, rna_sum, dna_sum, cancer_genes=collapse_genes_rank$cancer_genes, hgnc_symbol=collapse_genes_rank$hgnc_symbol,
      collapse_genes_rank[,rna_interest], collapse_genes_rank[,dna_interest])
    collapse_genes_rank_rna_dna = collapse_genes_rank_rna_dna[order(dna_prop, decreasing=T),]
    collapse_genes_rank_rna_dna$rank = 1:ngenes

    print(collapse_genes_rank_rna_dna[which(collapse_genes_rank_rna_dna$hgnc_symbol %in% c("TP53", "GAS7")),])

    collapse_genes_rank_rna_dna$cancer_text_annot = factor("", levels=c("", "*"))
    collapse_genes_rank_rna_dna$cancer_text_annot[collapse_genes_rank_rna_dna$cancer_genes == "normal"] = ""
    collapse_genes_rank_rna_dna$cancer_text_annot[collapse_genes_rank_rna_dna$cancer_genes != "normal"] = "*"

    pval_res_dir = "/cluster/project/raetsch/lab/03/ICGC/genecentric_analysis/results/qc_files/qc_master_table/perm_analysis/"
    p_vals = data.frame(fread(paste(pval_res_dir, "perm_res.tsv", sep="")))
    p_vals = p_vals[,c("hgnc_symbol", "pval")]
    collapse_genes_rank_rna_dna = merge(collapse_genes_rank_rna_dna, p_vals)
    collapse_genes_rank_rna_dna$pval = -1*log10(collapse_genes_rank_rna_dna$pval)

    ranked_prop_table = collapse_genes_rank_rna_dna[order(collapse_genes_rank_rna_dna$rank),]
    write.table(ranked_prop_table, paste(outdir, "/sig_ranked_table.tsv", sep=""), sep="\t", row.name=F, quote=F)

    collapse_genes_rank_rna_dna$fill_rna = "reg"
    collapse_genes_rank_rna_dna$fill_rna[which(collapse_genes_rank_rna_dna$hgnc_symbol %in% c("TP53", "CDK12", "GAS7"))] = "special"

    collapse_genes_rank_rna_dna$fill_dna = "reg"
    collapse_genes_rank_rna_dna$fill_dna[which(collapse_genes_rank_rna_dna$hgnc_symbol %in% c("TP53", "CDK12", "GAS7"))] = "special"

    color_names_rna = c("#009E73", "#000000")
    names(color_names_rna) = c("reg", "special")

    color_names_dna = c("#F0E442", "#000000")
    names(color_names_dna) = c("reg", "special")


    pdf(paste(outdir, "/", outname, "_barplot_sample_LR.pdf", sep=""))

      # regular scatter plots
      gg_rna = ggplot(collapse_genes_rank_rna_dna, aes(x=rank, y=rna_prop, fill=fill_rna)) +
        geom_bar(stat="identity")  +
        scale_fill_manual(values=color_names_rna) +
        theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
             axis.text.y=element_blank(),
             axis.title.x=element_blank(), axis.title.y=element_blank(),
             legend.position="none") + coord_flip() + scale_x_reverse()

      gg_dna = ggplot(collapse_genes_rank_rna_dna, aes(x=rank, y=dna_prop, fill=fill_dna)) +
        geom_bar(stat="identity") + coord_flip() + scale_x_reverse() +
        scale_y_reverse()  +
        scale_fill_manual(values=color_names_dna) +
        theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
             axis.text.y=element_blank(),
             axis.title.x=element_blank(), axis.title.y=element_blank(),
             legend.position="none")

      gg_pval = ggplot(collapse_genes_rank_rna_dna, aes(x=rank, y=0.2, fill=pval)) +
        geom_bar(stat="identity") + geom_text(aes(x=rank, y=0.1, label=cancer_text_annot, size=0.0001, vjust="top", hjust="top")) +
        coord_flip() + scale_x_reverse() +
        scale_y_reverse() +
        theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
             axis.text.y=element_blank(),
             axis.title.x=element_blank(), axis.title.y=element_blank(),
             legend.position="none") +
             scale_fill_gradientn(colours=brewer.pal(9,"YlOrRd"))

      gg = grid.arrange(gg_dna, gg_pval, gg_rna, ncol=3, nrow=1, widths=c(2, 1, 2), heights=c(1))
      print(gg)

      gg_pval = ggplot(collapse_genes_rank_rna_dna, aes(x=rank, y=0.2, fill=pval)) +
        geom_bar(stat="identity") +
             scale_fill_gradientn(colours=brewer.pal(9,"YlOrRd"))
    print(gg_pval)
    dev.off()

    collapse_genes_rank_rna_dna_cg = collapse_genes_rank_rna_dna[collapse_genes_rank_rna_dna$cancer_genes!="normal",]
    collapse_genes_rank_rna_dna_cg = collapse_genes_rank_rna_dna_cg[order(collapse_genes_rank_rna_dna_cg$dna_prop, decreasing=T),]
    collapse_genes_rank_rna_dna_cg$rank = 1:nrow(collapse_genes_rank_rna_dna_cg)

    # now melt to get the inside bars correct proportion...
    col_ids = setdiff(colnames(collapse_genes_rank_rna_dna_cg), c(dna_interest, rna_interest))
    collapse_genes_rank_rna_dna_cg = melt(collapse_genes_rank_rna_dna_cg, id.vars=col_ids)

    # change factor order
    collapse_genes_rank_rna_dna_cg$variable = factor(collapse_genes_rank_rna_dna_cg$variable,
        c("alt_prom", "ase_all", "expr_outlier", "fusion", "isSplice", "rna_edit", "cn", "variants"))


    ranked_prop_table = collapse_genes_rank_rna_dna_cg[order(collapse_genes_rank_rna_dna_cg$rank),]
    write.table(ranked_prop_table, paste(outdir, "/sig_cg_ranked_table.tsv", sep=""), sep="\t", row.name=F, quote=F)

    pdf(paste(outdir, "/", outname, "_barplot_sample_LR_cg.pdf", sep=""))

      # regular scatter plots
      gg_rna = ggplot(collapse_genes_rank_rna_dna_cg, aes(x=rank, y=rna_prop, fill=fill_rna)) +
        geom_bar(stat="identity")  +
        scale_fill_manual(values=color_names_rna)  +
        theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
             axis.text.y=element_blank(),
             axis.title.x=element_blank(), axis.title.y=element_blank(),
             legend.position="none")  +
             coord_flip() + scale_x_reverse()

      gg_dna = ggplot(collapse_genes_rank_rna_dna_cg, aes(x=rank, y=dna_prop, fill=fill_dna)) +
        geom_bar(stat="identity") + coord_flip() + scale_x_reverse() +
        scale_y_reverse()  +
        scale_fill_manual(values=color_names_dna) +
        theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
             axis.text.y=element_blank(),
             axis.title.x=element_blank(), axis.title.y=element_blank(),
             legend.position="none")

     gg_pval = ggplot(collapse_genes_rank_rna_dna_cg, aes(x=rank, y=value, fill=variable)) +
       geom_bar(stat="identity", position="fill") + #geom_text(aes(x=rank, y=0.1, label=cancer_text_annot, size=0.0001, vjust="top", hjust="top")) +
       coord_flip() + scale_x_reverse() +
       scale_y_reverse() +
       theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(), axis.title.y=element_blank(),
            legend.position="none") +
           scale_fill_brewer(palette="Set1")

      gg = grid.arrange(gg_dna, gg_pval, gg_rna, ncol=3, nrow=1, widths=c(2, 1, 2), heights=c(1))
      print(gg)

      gg_pval = ggplot(collapse_genes_rank_rna_dna_cg, aes(x=rank, y=0.2, fill=pval)) +
        geom_bar(stat="identity") +
             scale_fill_gradientn(colours=brewer.pal(9,"YlOrRd"))
    print(gg_pval)
    dev.off()


}


make_wing_plot <- function(outname, outdir, ordered_alts_melt, collapse_genes_rank, ngenes, rna_interest, dna_interest){

  curr_alt_interest = c(rna_interest, dna_interest)

  total_cgc = collapse_genes_rank$cancer_genes
  ordered_alts_melt = ordered_alts_melt[1:ngenes,]
  collapse_genes_rank = collapse_genes_rank[1:ngenes,]

  min_rank = min(ordered_alts_melt$rank)

  ordered_alts_melt$shifted_rank = ordered_alts_melt$rank - min_rank + 10

  # now go through each row and get the scaled versions of top two alterations
  for( curr_row in 1:nrow(ordered_alts_melt) ){

    top_values = unlist(sort(ordered_alts_melt[curr_row,curr_alt_interest],partial=2, decreasing=T)[1:2])
    top_values = (top_values/sum(top_values)) * ordered_alts_melt[curr_row, "shifted_rank"]
    row_order = unlist(colnames(ordered_alts_melt[,curr_alt_interest])[order(ordered_alts_melt[curr_row,curr_alt_interest], decreasing=T)])
    ordered_alts_melt[curr_row,row_order[1:2]] = top_values
    ordered_alts_melt[curr_row,row_order[3:length(curr_alt_interest)]] = 0

    # now go through the same thing and zero-out samples with zero contributions to rank
    top_values = unlist(sort(collapse_genes_rank[curr_row,curr_alt_interest],partial=2, decreasing=T)[1:2])
    top_values = (top_values/sum(top_values)) * ordered_alts_melt[curr_row, "shifted_rank"]
    collapse_genes_rank[curr_row,row_order[1:2]] = top_values
    collapse_genes_rank[curr_row,row_order[3:length(curr_alt_interest)]] = 0
  }

  dna_prop = rowSums(ordered_alts_melt[,dna_interest])
  rna_prop = rowSums(ordered_alts_melt[,rna_interest])
  ordered_alts_dna_rna = data.frame(rank=1:nrow(ordered_alts_melt), dna_prop, rna_prop, cancer_genes=collapse_genes_rank$cancer_genes)
  ordered_alts_dna_rna_melt = melt(ordered_alts_dna_rna, id=c("rank", "cancer_genes"))
  ordered_alts_dna_rna_melt$cancer_text_annot = factor("", levels=c("", "*"))
  ordered_alts_dna_rna_melt$cancer_text_annot[ordered_alts_dna_rna_melt$cancer_genes == "normal"] = ""
  ordered_alts_dna_rna_melt$cancer_text_annot[ordered_alts_dna_rna_melt$cancer_genes != "normal"] = "*"


  dna_prop = rowSums(collapse_genes_rank[,dna_interest])
  rna_prop = rowSums(collapse_genes_rank[,rna_interest])
  collapse_genes_rank_rna_dna = data.frame(rank=1:nrow(collapse_genes_rank), dna_prop, rna_prop, cancer_genes=collapse_genes_rank$cancer_genes)
  collapse_genes_rank_rna_dna = melt(collapse_genes_rank_rna_dna, id=c("rank", "cancer_genes"))
  collapse_genes_rank_rna_dna$cancer_text_annot = factor("", levels=c("", "*"))
  collapse_genes_rank_rna_dna$cancer_text_annot[collapse_genes_rank_rna_dna$cancer_genes == "normal"] = ""
  collapse_genes_rank_rna_dna$cancer_text_annot[collapse_genes_rank_rna_dna$cancer_genes != "normal"] = "*"


  ordered_alts_melt$rank = 1:nrow(ordered_alts_melt)
  ordered_alts_melt = melt(ordered_alts_melt[,c("rank", "cancer_genes", curr_alt_interest)], id=c("rank", "cancer_genes"))
  ordered_alts_melt$cancer_text_annot = factor("", levels=c("", "*"))
  ordered_alts_melt$cancer_text_annot[ordered_alts_melt$cancer_genes == "normal"] = ""
  ordered_alts_melt$cancer_text_annot[ordered_alts_melt$cancer_genes != "normal"] = "*"


  collapse_genes_rank$rank = 1:nrow(collapse_genes_rank)
  collapse_genes_rank_melt = melt(collapse_genes_rank[,c("rank", "cancer_genes", curr_alt_interest)], id=c("rank", "cancer_genes"))
  collapse_genes_rank_melt$cancer_text_annot = factor("", levels=c("", "*"))
  collapse_genes_rank_melt$cancer_text_annot[collapse_genes_rank_melt$cancer_genes == "normal"] = ""
  collapse_genes_rank_melt$cancer_text_annot[collapse_genes_rank_melt$cancer_genes != "normal"] = "*"

  pdf(paste(outdir, "/", outname, "_barplot_top2_alt.pdf", sep=""))
    # regular scatter plots
    gg_rank_norm = ggplot(ordered_alts_melt, aes(x=rank, y=value, fill=variable)) +
      geom_bar(stat="identity")  +
      theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
           axis.text.y=element_blank(),
           axis.title.x=element_blank(), axis.title.y=element_blank(),
           legend.position="none") + scale_fill_brewer(palette="Set1") +
      geom_text(aes(label=cancer_text_annot, y=-10)) + coord_flip() + scale_x_reverse()

    gg_samp_prop = ggplot(collapse_genes_rank_melt, aes(x=rank, y=value, fill=variable)) +
      geom_bar(stat="identity") + coord_flip() + scale_x_reverse() +
      scale_y_reverse() +
      theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
           axis.text.y=element_blank(),
           axis.title.x=element_blank(), axis.title.y=element_blank(),
           legend.position="none") + scale_fill_brewer(palette="Set1")

    gg = grid.arrange(gg_samp_prop, gg_rank_norm, ncol=2, nrow=1, widths=c(2, 2), heights=c(1))
    print(gg)


    gg_rna_dna_rank_normalized = ggplot(ordered_alts_dna_rna_melt, aes(x=rank, y=value, fill=variable)) +
      geom_bar(stat="identity") +
      theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
           axis.text.y=element_blank(),
           axis.title.x=element_blank(), axis.title.y=element_blank(),
           legend.position="none")  + scale_fill_manual(values=c("#009E73", "#F0E442")) +
      geom_text(aes(label=cancer_text_annot, y=-10)) + coord_flip() + scale_x_reverse()

    gg_rna_dna_samp_proportions = ggplot(collapse_genes_rank_rna_dna, aes(x=rank, y=value, fill=variable)) +
      geom_bar(stat="identity") + coord_flip() + scale_x_reverse() +
      scale_y_reverse() +
      theme(axis.ticks.y=element_blank(), panel.background=element_blank(),
           axis.text.y=element_blank(),
           axis.title.x=element_blank(), axis.title.y=element_blank(),
           legend.position="none") + scale_fill_manual(values=c("#009E73", "#F0E442"))

    gg = grid.arrange(gg_rna_dna_samp_proportions, gg_rna_dna_rank_normalized, ncol=2, nrow=1, widths=c(2, 2), heights=c(1))
    print(gg)

    gg = gg_rank_norm = ggplot(ordered_alts_melt, aes(x=rank, y=value, fill=variable)) +
      geom_bar(stat="identity") + scale_fill_brewer(palette="Set1")
    print(gg)

    gg = gg_rank_norm = ggplot(ordered_alts_dna_rna_melt, aes(x=rank, y=value, fill=variable)) +
      geom_bar(stat="identity") + scale_fill_manual(values=c("#009E73", "#F0E442"))
    print(gg)


  dev.off()

}

make_full_plot <- function(master_table, rna_interest, dna_interest, outdir, ngenes) {

  curr_alt_interest = c(rna_interest, dna_interest)

  rank_res = get_ranks(master_table, rna_interest, dna_interest)

  # rank and process collapse_genes
  collapse_genes_rank = rank_res[[1]]
  collapse_genes_rank = get_cancer_genes(collapse_genes_rank)
  collapse_genes_rank = collapse_genes_rank[order(collapse_genes_rank$rank, decreasing=F),]

  # rank and collapse ordered alts
  ordered_alts = rank_res[[2]]
  ordered_alts = get_cancer_genes(ordered_alts)
  ordered_alts_melt = ordered_alts
  ordered_alts_melt[,c("rank", curr_alt_interest)] = nrow(ordered_alts_melt)-ordered_alts_melt[,c("rank", curr_alt_interest)]
  ordered_alts_melt = ordered_alts_melt[order(ordered_alts_melt$rank, decreasing=T),]

  # make plots of basic stats
  outname="overall"
  make_pie_charts(outname, outdir, ordered_alts_melt, collapse_genes_rank, ngenes, rna_interest, dna_interest)


  # plot all significant CGC/DRIVER genes
  ordered_alts_cg = ordered_alts_melt[ordered_alts_melt$cancer_genes != "normal",]
  collapse_gene_cg = collapse_genes_rank[collapse_genes_rank$cancer_genes != "normal",]
  ngenes_cg = nrow(collapse_gene_cg[collapse_gene_cg$rank < ngenes,])
  outname = "significant_cg"
  make_pie_charts(outname, outdir, ordered_alts_cg, collapse_gene_cg, ngenes_cg,
    rna_interest, dna_interest)

  outname="all_cg"
  make_curve_plot(outname, outdir, ordered_alts_melt, collapse_genes_rank, ngenes, rna_interest, dna_interest)


}

run <- function(){

  rna_interest  = c("alt_prom", "expr_outlier", "rna_edit", "ase_all", "fusion", "isSplice")
  dna_interest  = c("variants", "cn")
  curr_alt_interest = c(rna_interest, dna_interest)

  # get number of genes
  outname = "secondSmallest_full_v.13"
  full_table_file = paste(base_dir, "/results/recurrance_analysis/no_dom_trans/", outname, "_full.txt", sep="")
  full_table = data.frame(fread(full_table_file))
  percentile_pass = data.frame(fread(filter_score_path, header=T))
  score_pass = percentile_pass[percentile_pass$percentiles==5, "scores"]
  ngenes = nrow(full_table[full_table$sum_order < score_pass,])


  master_table = get_master_table(curr_alt_interest)
  outdir = paste(base_dir, "/results/rna_dna_ranks/all_sample_rank", sep="")
  dir.create(outdir, showWarnings=F)
  make_full_plot(master_table, rna_interest, dna_interest, outdir, ngenes)


}


# plot within a cancer type
run()
