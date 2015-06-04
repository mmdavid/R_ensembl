setwd("/Users/mmdavid/Documents/lifeboat/R_ensembl")

#Package load bioconductor http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library("biomaRt", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")

# Biomart citation
#Durinck S, Spellman P, Birney E and Huber W (2009). “Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.” Nature Protocols, 4, pp. 1184–1191.
#Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A and Huber W (2005). “BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis.” Bioinformatics, 21, pp. 3439–3440.

browseVignettes("biomaRt") #download documentation 

#stupid problem with evrison, here is the right address
install.packages("Matrix", repos = "http:http://cran.rstudio.com/bin/macosx/mavericks/contrib/3.1/Matrix_1.2-0.tgz", type="source")


#ok let's start to use it 
mart.hs <- useMart("ensembl","hsapiens_gene_ensembl")

#I'm using GRCh38.p2 => last verison
att<-listAttributes(mart.hs)
write.csv(att, file="att.csv")

grep("mmusculus_homolog_dn",listAttributes(mart.hs), ignore.case = T)

#ok that's for each gene gives me the name of variation to go grab in the next request
getBM(attributes = "variation_name", filters = "hgnc_symbol", values = "GABRA6", mart = mart.hs)

#that gives me al the snps for eahc gene
getBM(attributes = "ensembl", filters = "hgnc_symbol", values = "GABRA6", mart = mart.hs)

#ok now the satus of each snp
snpmart = useMart("snp", dataset="hsapiens_snp")

getBM(attributes = "clinical_significance", filters = "hgnc_symbol", values = "GABRA6", mart = snpmart)

#ok let's load mylist
load("mylist2.Rda")

mylistnoASD<-mylist2[-4]
unlist_noASD<-unlist(mylistnoASD)
#unique to ASD 
allASDgenes<-unlist(mylist2[4])
uniqueASD<-setdiff(allASDgenes,unlist_noASD)
nonuniqueASD<-intersect(allASDgenes,unlist_noASD)

#find the ds for all the genes
mice_human_dn=c()
for (i in 1:length(uniqueASD)){
a<-getBM(attributes = "mmusculus_homolog_dn", filters = "hgnc_symbol", values = uniqueASD[i], mart = mart.hs)
mice_human_dn<-c(mice_human_dn,a)
cat(i,"\t")
}
names(mice_human_dn)<-uniqueASD
unlist_mice_human_dn<-unlist(mice_human_dn)
unlist_mice_human_dn2<-as.numeric(unlist_mice_human_dn)
names(unlist_mice_human_dn2)<-names(unlist_mice_human_dn)

mice_human_ds=c()
for (i in 1:length(uniqueASD)){
  a<-getBM(attributes = "mmusculus_homolog_ds", filters = "hgnc_symbol", values = uniqueASD[i], mart = mart.hs)
  mice_human_ds<-c(mice_human_ds,a)
  cat(i,"\t")
}
names(mice_human_ds)<-uniqueASD
unlist_mice_human_ds<-unlist(mice_human_ds)
unlist_mice_human_ds2<-as.numeric(unlist_mice_human_ds)
names(unlist_mice_human_ds2)<-names(unlist_mice_human_ds)

mice_human_dn_ds<-unlist_mice_human_dn/unlist_mice_human_ds


nonuniqueASD_mice_human_dn=c()
for (i in 1:length(nonuniqueASD)){
  a<-getBM(attributes = "mmusculus_homolog_dn", filters = "hgnc_symbol", values = nonuniqueASD[i], mart = mart.hs)
  nonuniqueASD_mice_human_dn<-c(nonuniqueASD_mice_human_dn,a)
  cat(i,"\t")
}
names(nonuniqueASD_mice_human_dn)<-nonuniqueASD
unlist_nonuniqueASD_mice_human_dn<-unlist(nonuniqueASD_mice_human_dn)
unlist_nonuniqueASD_mice_human_dn<-as.numeric(unlist_nonuniqueASD_mice_human_dn)

nonuniqueASD_mice_human_ds=c()
for (i in 1:length(nonuniqueASD)){
  a<-getBM(attributes = "mmusculus_homolog_ds", filters = "hgnc_symbol", values = nonuniqueASD[i], mart = mart.hs)
  nonuniqueASD_mice_human_ds<-c(nonuniqueASD_mice_human_ds,a)
  cat(i,"\t")
}
names(nonuniqueASD_mice_human_ds)<-nonuniqueASD
unlist_nonuniqueASD_mice_human_ds<-unlist(nonuniqueASD_mice_human_ds)
unlist_nonuniqueASD_mice_human_ds<-as.numeric(unlist_nonuniqueASD_mice_human_ds)
#problem it seems ot have several ds or dn value per gene
#for example
#$RHOXF1
#[1] 1.1269 1.0964 0.9145 1.0531
#What is going on ?
getBM(attributes = "mmusculus_homolog_ds", filters = "hgnc_symbol", values = "RHOXF1", mart = mart.hs)
getBM(attributes = "hgnc_symbol", filters = "hgnc_symbol", values = "RHOXF1", mart = mart.hs)

#facing tons of problems with R studio and git syn
# I had to go to the config file in my folder: /Users/mmdavid/R_ensembl/.git #to see the hidden file/folder: ls -a
#then I change 
#OK so according to ensembl there is actually orthologs! So I'll take the avearge of the orhtologs 


nonuniqueASD_mice_human_dn_ds<-unlist_nonuniqueASD_mice_human_dn/unlist_nonuniqueASD_mice_human_ds



#now let's grab all the RS per genes



