setwd("/Users/mmdavid/Documents/lifeboat/R_ensembl")

--------------------------------------------------------------------------------------------------------
  #BIOMART API FOR R ANd ENSENBL DATABASE: extraction dsdn value for human and mice alignment 
  --------------------------------------------------------------------------------------------------------
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
#OK so according to ensembl there is actually orthologs! So I'll take the avearge of the orhtologs 
mice_human_dn_mean<-sapply(mice_human_dn, function(x) mean(x))
#remove NA
mice_human_dn_mean<-sapply(mice_human_dn_mean, function (x) x[!is.na(x)]) 
#stupid function put it back as a list, so just unlist again. 
unlist_mice_human_dn_mean<-unlist(mice_human_dn_mean)

mice_human_ds=c()
for (i in 1:length(uniqueASD)){
  a<-getBM(attributes = "mmusculus_homolog_ds", filters = "hgnc_symbol", values = uniqueASD[i], mart = mart.hs)
  mice_human_ds<-c(mice_human_ds,a)
  cat(i,"\t")
}
names(mice_human_ds)<-uniqueASD

--------------------------------------------------------------------------------------------------------
  #MICE DNDS CALCULATIONS
  --------------------------------------------------------------------------------------------------------
  #MEAN
  #sometimes several ds or dn value due to isomers or copy gene numbers, I took the mean, actually David indicated to take the lmax, so need to redo this later
  mice_human_ds_mean<-sapply(mice_human_ds, function(x) mean(x))
#remove NA
mice_human_ds_mean<-sapply(mice_human_ds_mean, function (x) x[!is.na(x)]) 
#stupid function put it back as a list, so just unlist again. 
unlist_mice_human_ds_mean<-unlist(mice_human_ds_mean)
#ok do the ratio here, I need ot make sure ht ename of the gene is the same 
uniqueASD_dn_ds=c()
for (i in 1:length(unlist_mice_human_dn_mean)){
  for (j in 1:length(unlist_mice_human_ds_mean)){
    if (names(unlist_mice_human_dn_mean[i]) == (names(unlist_mice_human_ds_mean[j]))) {
      a<-unlist_mice_human_dn_mean[i]/unlist_mice_human_ds_mean[j]
      names(a)<-names(unlist_mice_human_dn_mean[i])
      uniqueASD_dn_ds<-c(uniqueASD_dn_ds,a)
    }
  }
}

#same with genes involved in ASD and other comorbid to ASD
nonuniqueASD_mice_human_dn=c()
for (i in 1:length(nonuniqueASD)){
  a<-getBM(attributes = "mmusculus_homolog_dn", filters = "hgnc_symbol", values = nonuniqueASD[i], mart = mart.hs)
  nonuniqueASD_mice_human_dn<-c(nonuniqueASD_mice_human_dn,a)
  cat(i,"\t")
}
names(nonuniqueASD_mice_human_dn)<-nonuniqueASD
nonuniqueASD_mice_human_dn_mean<-sapply(nonuniqueASD_mice_human_dn, function(x) mean(x))
#remove NA #26 not calculated
nonuniqueASD_mice_human_dn_mean<-sapply(nonuniqueASD_mice_human_dn_mean, function (x) x[!is.na(x)]) 
#stupid function put it back as a list, so just unlist again. 
unlist_nonuniqueASD_mice_human_dn_mean<-unlist(nonuniqueASD_mice_human_dn_mean)

nonuniqueASD_mice_human_ds=c()
for (i in 1:length(nonuniqueASD)){
  a<-getBM(attributes = "mmusculus_homolog_ds", filters = "hgnc_symbol", values = nonuniqueASD[i], mart = mart.hs)
  nonuniqueASD_mice_human_ds<-c(nonuniqueASD_mice_human_ds,a)
  cat(i,"\t")
}
names(nonuniqueASD_mice_human_ds)<-nonuniqueASD
nonuniqueASD_mice_human_ds_mean<-sapply(nonuniqueASD_mice_human_ds, function(x) mean(x))
#remove NA #26 not calculated
nonuniqueASD_mice_human_ds_mean<-sapply(nonuniqueASD_mice_human_ds_mean, function (x) x[!is.na(x)]) 
#function put it back as a list, so just unlist again. 
unlist_nonuniqueASD_mice_human_ds_mean<-unlist(nonuniqueASD_mice_human_ds_mean)

#ok do the ratio here, I need ot make sure ht ename of the gene is the same 

nonuniqueASD_dn_ds=c()
for (i in 1:length(unlist_nonuniqueASD_mice_human_dn_mean)){
  for (j in 1:length(unlist_nonuniqueASD_mice_human_ds_mean)){
    if (names(unlist_nonuniqueASD_mice_human_dn_mean[i]) == (names(unlist_nonuniqueASD_mice_human_ds_mean[j]))) {
      a<-unlist_nonuniqueASD_mice_human_dn_mean[i]/unlist_nonuniqueASD_mice_human_ds_mean[j]
      names(a)<-names(unlist_nonuniqueASD_mice_human_dn_mean[i])
      nonuniqueASD_dn_ds<-c(nonuniqueASD_dn_ds,a)
    }
  }
}

#OK let's try to box plot that thing!!
boxplot(nonuniqueASD_dn_ds, uniqueASD_dn_ds)
plot(sort(nonuniqueASD_dn_ds))
points(sort(uniqueASD_dn_ds), col="red")
wilcox.test(nonuniqueASD_dn_ds,uniqueASD_dn_ds)

# Wilcoxon rank sum test with continuity correction
#data:  nonuniqueASD_dn_ds and uniqueASD_dn_ds
#W = 80975, p-value = 0.008392
#alternative hypothesis: true location shift is not equal to 0

#let's look at the density;
plot(nonuniqueASD_dn_ds)

#let's try to plot ggplot
install.packages("ggplot2")
library("ggplot2")

#seems like I need a serious dataframe here. Let's do that: 
nonuniqueASD_dn_ds
all_dn_ds<-append(nonuniqueASD_dn_ds,uniqueASD_dn_ds)
NONUNIQ<-rep("NONUNIQ",length(nonuniqueASD_dn_ds))
UNIQ<-rep("UNIQ",length(uniqueASD))
factorforall_dn_ds<-append(NONUNIQ,UNIQ)
all_dn_ds<-as.data.frame(all_dn_ds)
all_dn_ds$factorforall_dn_ds<-factorforall_dn_ds

#ggploting 
#usefulpage for ggplot http://www.ling.upenn.edu/~joseff/rstudy/summer2010_ggplot2_intro.html
#http://www.sthda.com/english/wiki/ggplot2-boxplot-easy-box-and-whisker-plots-maker-function
update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
ggplot(all_dn_ds, aes(factorforall_dn_ds, all_dn_ds))+
  geom_boxplot(aes(colour = factor(factorforall_dn_ds)), fill = c("green", "yellow"), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1)

ggplot(all_dn_ds, aes(factorforall_dn_ds, all_dn_ds))+
  geom_boxplot(aes(colour = factor(factorforall_dn_ds)), fill = c("green", "yellow"))

#with all the points: just not great 
ggplot(all_dn_ds, aes(factorforall_dn_ds, all_dn_ds))+
  geom_boxplot(aes(colour = factor(factorforall_dn_ds)), fill = c("green", "yellow"), color = c("chartreuse3","orange"))+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=0.2) 


#nicer boxplot
#ok let's try this density graph now

blou<-ggplot(all_dn_ds, aes(all_dn_ds, fill = factorforall_dn_ds)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow"))
ggplot(mpg, aes(displ, hwy))

