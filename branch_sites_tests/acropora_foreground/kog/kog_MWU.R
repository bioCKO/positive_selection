setwd("~/git_Repositories/positive_selection/branch_sites_tests/acropora_foreground/kog")
source('~/git_Repositories/positive_selection/scripts/kogMWU_functions.R')

### inname = two column table (same as GO) where first column is gene and second is a measure
### gene2kog = table of gene names and their associated KOGs
### alternative = "g" if the measure of interest gets bigger with interest (e.g., 0 and 1 for presence/absence), "l" if the measure of interest gets smaller with interest (e.g., untransformed pvalue), or "t" if the measure of interest is interesting a both ends (e.g., log-transformed pvalue)


#run kog
file.name = 'branch_site_LRT_results_for_KOGGmwu.csv'

M = kog.mwu(inname = file.name, gene2kog = "AdigNCBI_PRO_iso2kogClass.tab", Alternative = "g")
M
write.table(M, file="branch_sites_kog_enrichment_table.tsv", quote = F, row.names = F, sep = "\t")
################# Generate heatmap
# compiling a table of delta-ranks to compare these results:
kognames=M$term
kogtable=data.frame("M" = M$delta.rank)
row.names(kogtable)=kognames

#sort the dataframe
y = rownames(kogtable)
x = kogtable$M
z = data.frame("M" = x, 'names' = y)
z = z[with(z, order(M)),]
kogtable = data.frame('M' = z$M)
row.names(kogtable) = z$names



# making a heatmap. Un-remark and run the next two lines if you don't have packages pheatmap and RColorBrewer installed. 
# install.package("pheatmap")
# install.package("RColorBrewer")
library(pheatmap)
library(RColorBrewer)

# you may want to adust these colors to make sure 0 is in the yellow; this will depend on your scale of delta-ranks
color = colorRampPalette(rev(c(brewer.pal(n = 7, name = "RdYlBu"),"darkblue","darkblue")))(100)

quartz()
pheatmap(as.matrix(kogtable), clustering_distance_rows = "euclidean", cluster_rows = F, cluster_cols = F, color=color, treeheight_row=15, treeheight_col=20)

#write out the results
out.name = 'kogg_bs_acro_rcut100_out.txt'
write.table(M, out.name, row.names = F, quote = F, sep = "\t")
