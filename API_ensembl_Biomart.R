setwd("/Users/mmdavid/Documents/lifeboat/R_ensembl")

--------------------------------------------------------------------------------------------------------------------------
#BIOMART API for R for Ensembl database
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


#same with non unique to ASD

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
#stupid function put it back as a list, so just unlist again. 
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

#useful page for ggplot
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

#usefulpage for ggplot http://www.ling.upenn.edu/~joseff/rstudy/summer2010_ggplot2_intro.html
#http://www.sthda.com/english/wiki/ggplot2-boxplot-easy-box-and-whisker-plots-maker-function


#OK so David sent me the value for several primate dn/ds and african american let's do this again
#I need ot grab all the names of Ensembl into my 
to_search<-attributes(mart.hs)
grep("ensembl", to_search)

to_search[1]
#let's try this attribute

getBM(attributes = "ensembl_gene_id", filters = "hgnc_symbol", values = "GABRA6", mart = mart.hs)

#that works! BIM! 
nonuniqueASD_emsblID=c()
for (i in 1:length(nonuniqueASD)){
  a<-getBM(attributes = "ensembl_gene_id", filters = "hgnc_symbol", values = nonuniqueASD[i], mart = mart.hs)
  nonuniqueASD_emsblID<-c(nonuniqueASD_emsblID,a)
  cat(i,"\t")
  
}

save(nonuniqueASD_emsblID,file="nonuniqueASD_emsblID.Rda")
nonuniqueASD_emsblID<-unlist(nonuniqueASD_emsblID)

uniqueASD_emsblID=c()
for (i in 1:length(uniqueASD)){
  a<-getBM(attributes = "ensembl_gene_id", filters = "hgnc_symbol", values = uniqueASD[i], mart = mart.hs)
  uniqueASD_emsblID<-c(uniqueASD_emsblID,a)
  cat(i,"\t")
}

save(uniqueASD_emsblID,file="uniqueASD_emsblID.Rda")

uniqueASD_emsblID<-unlist(uniqueASD_emsblID)
names(uniqueASD_emsblID)<-c()

#now import files from David
#first file: 
#primates_dnds_v69 is a file that contains the average dN/dS for nine primates. Better than the human since it is closer to what constraint is in human.
#pnps_file_1000Gafricanpops contains Ensembl Gene ID, PN, PS, DN and DS (order of the columns). PN and PS is from african pops from the 1000 Genomes project.

primates_dnds_v69<-read.table("primates_dnds_v69")
head(primates_dnds_v69)
primates_dnds_v69$V1<-as.character(primates_dnds_v69$V1)

primate_dnds_ASDonly<-subset(primates_dnds_v69,primates_dnds_v69$V1 %in% uniqueASD_emsblID)
primate_dnds_ASDcomord<-subset(primates_dnds_v69,primates_dnds_v69$V1 %in% nonuniqueASD_emsblID)

#I need ot od a df
ASDONLY<-rep("ASDONLY",dim(primate_dnds_ASDonly)[1])
ASDCOMOR<-rep("ASDCOMOR",dim(primate_dnds_ASDcomord)[1])

ASD_status<-c(ASDONLY,ASDCOMOR)
genes_names<-c(primate_dnds_ASDonly$V1,primate_dnds_ASDcomord$V1)
dnds_primate<-c(primate_dnds_ASDonly$V2,primate_dnds_ASDcomord$V2)

ASD_status<-as.data.frame(ASD_status)

ASD_status$genes_names<-genes_names

ASD_status$dnds_primate<-dnds_primate

primates<-ASD_status

colnames(ASD_status)<-c("ASD_status_onlyornot","genes_names","dnds_primate")

#super boxplot

update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
ggplot(ASD_status, aes(ASD_status_onlyornot, dnds_primate))+
  geom_boxplot(aes(colour = factor(ASD_status)), fill = c("green", "yellow"), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1)

#ggplot(all_dn_ds, aes(factorforall_dn_ds, all_dn_ds))+
 # geom_boxplot(aes(colour = factor(factorforall_dn_ds)), fill = c("green", "yellow"), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1)

#cool looking good
#test
wilcox.test( dnds_primate ~ ASD_status_onlyornot, data= primates)

#data:  dnds_primate by ASD_status_onlyornot
#W = 94548.5, p-value = 0.1574
#alternative hypothesis: true location shift is not equal to 0
#density graph


blou<-ggplot(primates, aes(dnds_primate, fill = ASD_status_onlyornot)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow"))

ggplot(mpg, aes(displ, hwy))

#now let's try with the 1000 genomes projects
#pnps_file_1000Gafricanpops contains Ensembl Gene ID, PN, PS, DN and DS (order of the columns). PN and PS is from african pops from the 1000 Genomes project.

pnps_file_1000Gafricanpops<-read.table("pnps_file_1000Gafricanpops")
colnames(pnps_file_1000Gafricanpops)<-c("Ensembl_Gene_ID", "PN", "PS", "DN", "DS")
pnps_file_1000Gafricanpops$PNPS<-(pnps_file_1000Gafricanpops$PN /pnps_file_1000Gafricanpops$PS)
pnps_file_1000Gafricanpops$DNDS<-(pnps_file_1000Gafricanpops$DN /pnps_file_1000Gafricanpops$DS)

#quick silly fix
names(nonuniqueASD_emsblID)<-c()
uniqueASD_emsblID

none<-rep("none",dim(pnps_file_1000Gafricanpops)[1])
pnps_file_1000Gafricanpops$ASD_status<-none

#fix fix fix
allgenes<-c(nonuniqueASD_emsblID,uniqueASD_emsblID)
pnps_file_1000Gafricanpops$Ensembl_Gene_ID<-as.character(pnps_file_1000Gafricanpops$Ensembl_Gene_ID)

pnps_file_1000GafricanpopsASD<-subset(pnps_file_1000Gafricanpops, pnps_file_1000Gafricanpops$Ensembl_Gene_ID %in% allgenes)

for (i in 1:dim(pnps_file_1000GafricanpopsASD)[1]){
  if (pnps_file_1000GafricanpopsASD$Ensembl_Gene_ID[i] %in% nonuniqueASD_emsblID) {pnps_file_1000GafricanpopsASD$ASD_status[i] <- "comorbdASD"} else {pnps_file_1000GafricanpopsASD$ASD_status[i] <- "ASDonly"}
  cat(i,"\t")
}

#ok let's do the dnds boxplot
pnps_file_1000GafricanpopsASD_save<-pnps_file_1000GafricanpopsASD
pnps_file_1000GafricanpopsASD[is.na(pnps_file_1000GafricanpopsASD)] <- 0
pnps_file_1000GafricanpopsASD1<- pnps_file_1000GafricanpopsASD

#remove the one with infinite
inf_pb<-pnps_file_1000GafricanpopsASD$DNDS == "Inf"
noinf_pnps_file_1000GafricanpopsASD<-pnps_file_1000GafricanpopsASD[!inf_pb,]

update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
ggplot(noinf_pnps_file_1000GafricanpopsASD, aes(ASD_status,DNDS))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)

ggplot(all_dn_ds, aes(factorforall_dn_ds, all_dn_ds))+
  geom_boxplot(aes(colour = factor(factorforall_dn_ds)), fill = c("green", "yellow"), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1)


blou<-ggplot(noinf_pnps_file_1000GafricanpopsASD, aes(DNDS, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow"))

#cool looking good
#test
wilcox.test( DNDS ~ ASD_status, data= noinf_pnps_file_1000GafricanpopsASD)

#Wilcoxon rank sum test with continuity correction
#data:  DNDS by ASD_status
#W = 79654, p-value = 0.9188
#alternative hypothesis: true location shift is not equal to 0

#ok how about the polymorphism ------------------------------------------------------------------

#ok let's do the dnds boxplot


#remove the one with infinite
inf_pb<-pnps_file_1000GafricanpopsASD$PNPS == "Inf"
psnoinf_pnps_file_1000GafricanpopsASD<-pnps_file_1000GafricanpopsASD[!inf_pb,]

update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
ggplot(psnoinf_pnps_file_1000GafricanpopsASD, aes(ASD_status,PNPS))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)


blou<-ggplot(psnoinf_pnps_file_1000GafricanpopsASD, aes(PNPS, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow"))

#cool looking good
#test
wilcox.test( PNPS ~ ASD_status, data= noinf_pnps_file_1000GafricanpopsASD)

#Wilcoxon rank sum test with continuity correction
#data:  PNPS by ASD_status
#W = 75719.5, p-value = 0.9077
#alternative hypothesis: true location shift is not equal to 0
#calculate the Mc Donald Kreitman  neutrality index (NI) 

inf_pb<-(pnps_file_1000GafricanpopsASD$PNPS == "Inf") 

noinf_pnps_file_1000GafricanpopsASD<-pnps_file_1000GafricanpopsASD[!inf_pb,]


inf_pb_dnds<-noinf_pnps_file_1000GafricanpopsASD$DNDS == "Inf"
noinf_dnda_noinf_pnps_file_1000GafricanpopsASD<-noinf_pnps_file_1000GafricanpopsASD[!inf_pb_dnds,]
head(noinf_dnda_noinf_pnps_file_1000GafricanpopsASD)


noinf_dnda_noinf_pnps_file_1000GafricanpopsASD$NI<-noinf_dnda_noinf_pnps_file_1000GafricanpopsASD$DNDS*noinf_dnda_noinf_pnps_file_1000GafricanpopsASD$PNPS
#test andplots
update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
ggplot(noinf_dnda_noinf_pnps_file_1000GafricanpopsASD, aes(ASD_status,NI))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)


blou<-ggplot(noinf_dnda_noinf_pnps_file_1000GafricanpopsASD, aes(NI, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow"))

#cool looking good
#test
wilcox.test( NI ~ ASD_status, data= noinf_dnda_noinf_pnps_file_1000GafricanpopsASD)

#data:  NI by ASD_status
#W = 77867.5, p-value = 0.338
#alternative hypothesis: true location shift is not equal to 0


# let's try again but adding 1 to NOT had the remove some genes 

pnps_file_1000GafricanpopsASD2<-subset(pnps_file_1000Gafricanpops, pnps_file_1000Gafricanpops$Ensembl_Gene_ID %in% allgenes)

for (i in 1:dim(pnps_file_1000GafricanpopsASD2)[1]){
  if (pnps_file_1000GafricanpopsASD2$Ensembl_Gene_ID[i] %in% nonuniqueASD_emsblID) {pnps_file_1000GafricanpopsASD2$ASD_status[i] <- "comorbdASD"} else {pnps_file_1000GafricanpopsASD2$ASD_status[i] <- "ASDonly"}
  cat(i,"\t")
}

#now add one to eveyry one

dim(pnps_file_1000GafricanpopsASD2)

pnps_file_1000GafricanpopsASD2$PS<-pnps_file_1000GafricanpopsASD2$PS + 1
pnps_file_1000GafricanpopsASD2$PN<-pnps_file_1000GafricanpopsASD2$PN + 1
pnps_file_1000GafricanpopsASD2$DS<-pnps_file_1000GafricanpopsASD2$DS + 1
pnps_file_1000GafricanpopsASD2$DN<-pnps_file_1000GafricanpopsASD2$DN + 1

pnps_file_1000GafricanpopsASD2$PNPS<-(pnps_file_1000GafricanpopsASD2$PN /pnps_file_1000GafricanpopsASD2$PS)
pnps_file_1000GafricanpopsASD2$DNDS<-(pnps_file_1000GafricanpopsASD2$DN /pnps_file_1000GafricanpopsASD2$DS)



update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
ggplot(pnps_file_1000GafricanpopsASD2, aes(ASD_status,DNDS))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)


blou<-ggplot(pnps_file_1000GafricanpopsASD2, aes(DNDS, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow"))

#cool looking good
#test
wilcox.test( DNDS ~ ASD_status, data= pnps_file_1000GafricanpopsASD2)
#data:  DNDS by ASD_status
#W = 89255.5, p-value = 0.2047
#alternative hypothesis: true location shift is not equal to 0
#Wilcoxon rank sum test with continuity correction


#data:  DNDS by ASD_status
#W = 89255.5, p-value = 0.2047
#alternative hypothesis: true location shift is not equal to 0

#ok how about the polymorphism ------------------------------------------------------------------

#ok let's do the dnds boxplot

update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
ggplot(pnps_file_1000GafricanpopsASD2, aes(ASD_status,PNPS))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)

blou<-ggplot(pnps_file_1000GafricanpopsASD2, aes(PNPS, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow"))

#cool looking good
#test
wilcox.test( PNPS ~ ASD_status, data= pnps_file_1000GafricanpopsASD2)


#Wilcoxon rank sum test with continuity correction

#data:  PNPS by ASD_status: YES IT'S SIGNIFICATIVE 
#W = 91933, p-value = 0.03767
#alternative hypothesis: true location shift is not equal to 0

#but is it becaus eI am adding +1? How about the test with inf and Na? Let's try, it's ranking, may work. 

pnps_file_1000GafricanpopsASD3<-subset(pnps_file_1000Gafricanpops, pnps_file_1000Gafricanpops$Ensembl_Gene_ID %in% allgenes)

for (i in 1:dim(pnps_file_1000GafricanpopsASD3)[1]){
  if (pnps_file_1000GafricanpopsASD3$Ensembl_Gene_ID[i] %in% nonuniqueASD_emsblID) {pnps_file_1000GafricanpopsASD3$ASD_status[i] <- "comorbdASD"} else {pnps_file_1000GafricanpopsASD3$ASD_status[i] <- "ASDonly"}
  cat(i,"\t")
}

#BTW: we found 864 / 946 genes with values in 1000 genomes project
wilcox.test( PNPS ~ ASD_status, data= pnps_file_1000GafricanpopsASD3)
#Wilcoxon rank sum test with continuity correction
#data:  PNPS by ASD_status
#W = 20661.5, p-value = 0.05914
#alternative hypothesis: true location shift is not equal to 0


#Apparently R remove the NaN and Inf on its own, wonderful

update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
p1<-ggplot(pnps_file_1000GafricanpopsASD3, aes(ASD_status,PNPS))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)

#trying to add dots
# if needed update_geom_defaults("point", list(colour = "black"))
#Making the effort for nice boxplots 

p1<-ggplot(pnps_file_1000GafricanpopsASD3, aes(ASD_status,PNPS))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 5, outlier.shape = 1, fill = c("green", "yellow"),)

p1 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 

#tried to change the scale but not great because it's infinite+ coord_cartesian(ylim=c(0, 50))



#usefull lnik http://docs.ggplot2.org/0.9.3.1/geom_point.html

blou<-ggplot(pnps_file_1000GafricanpopsASD3, aes(PNPS, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.)

blou + scale_fill_manual( values = c("green","yellow"))

#dnds all I can do now
wilcox.test( DNDS ~ ASD_status, data= pnps_file_1000GafricanpopsASD3)

#not significative 
#Wilcoxon rank sum test with continuity correction
#data:  DNDS by ASD_status
#W = 23366, p-value = 0.3547
#alternative hypothesis: true location shift is not equal to 0
