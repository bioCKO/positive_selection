#LRT_for_branch_sites_models
#Groves Dixon
#8/28/15
#R script LRT_for_branch_site_models.R


#set working directory
setwd("~/git_Repositories/positive_selection/branch_sites_tests/acropora_foreground")
# load("~/git_Repositories/metaTranscriptomes/working_directory/MBD-seq_Image.R")
# source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")


#likelihood ratio test function
lrt = function(la, lo, df){
  G = 2*(la - lo)
  # print(G)
  p.values = pchisq(G, df, ncp = 0, lower.tail = F)
  return(p.values)
}

#read in the data from the null and alternative models for the branch-site test for positive selection
#see PAML manual for descriptions of these models
null = read.table("nullLikelihoods_branchSites.txt", header = T) #file generated using parse_codeml_branch_sites.py (see Positive_Selection_Walkthrough.txt)
alt = read.table("altLikelihoods_branchSites.txt", header = T)   #file generated using parse_codeml_branch_sites.py (see Positive_Selection_Walkthrough.txt)


#these should show the contig, the number of parameters for the model, and it's log likelihood
head(null); head(alt)



#run lrt between the models
contig = alt$contig
contigs.null = null$contig
sum(contig == contigs.null) == nrow(alt) #double-check the files were assembled correctly
la = alt$likelihood
lo = null$likelihood
p.values = lrt(la, lo, df=1) #note degrees of freedom is difference in number of parameters between models (in this case always 1)

#record the results of the test in a new dataframe
dat = data.frame(la, lo, p.values, contig, contigs.null)
head(dat)

#doublecheck the contigs match
sum(dat$contig == dat$contigs.null) == nrow(dat)


#get adjusted p values for multiple tests
dat$adj.p = p.adjust(dat$p.values, method = 'BH')

#get an idea of how many genes are significant
cuts = c(0.1, 0.05, 0.01, 0.001, 0.0001)
for (i in cuts){
	sub = dat[dat$adj.p < i,]
	unadjust = dat[dat$p.values < i,]
	print(paste(paste(paste(i, nrow(sub)), nrow(unadjust)), (nrow(sub)/nrow(dat))*100))
}

hist(dat$p.values)



#WRITE OUT THE RESULTS
#get gene names
lnames = load("~/git_Repositories/positive_selection/datasets/ProteinTable.Rdata")#assembled from here: https://www.ncbi.nlm.nih.gov/genome/proteins/10529?genome_assembly_id=263537
colnames(ptable) = c('genome_contig', 'locusName', 'contig', 'protein.name')
head(ptable)
head(dat)
out = merge(ptable, dat, by = 'contig', all.y=TRUE)
dim(out)
dim(dat)
out = out[order(out$p.values),]
out$geneName = out$protein.name
out$protein.name<-NULL
out$contigs.null<-NULL
head(out, n=30)
out[out$adj.p < 0.15,]#look at marginally significant genes
#write out text file of results
write.table(out, "branch_sites_LRT_results_acro.txt", quote = T, row.names = T, sep = "\t")
#save for other analyses
a.pdat=out
save(a.pdat, file='branch_sites_results.Rdata')



#export the data for GO and KOGG enrichment tests using MWU-tests
out2 = data.frame(out$contig, -log(out$p.values, 10))
colnames(out2) = c('gene', 'logp')
head(out2)
nrow(out2)
#write out for GO
#use this file as input for 
write.table(out2, 'branch_site_LRT_results_for_GOmwu.csv', quote = F, row.names = F, sep = ",")
ps = out$p.values
x = ps < 0.05
ps[x==TRUE]<-1
ps[x==FALSE]<-0
out3 = data.frame(out$contig, ps)
colnames(out3) = c('isogroup', 'sig')
write.table(out3, 'branch_site_LRT_results_for_Fisher_GOmwu.csv', quote = F, row.names = F, sep = ",")

#write out for KOGG
write.table(out2, 'branch_site_LRT_results_for_KOGGmwu.csv', quote = F, row.names = F, sep = ",")