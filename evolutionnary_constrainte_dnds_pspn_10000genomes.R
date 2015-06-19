#usefull ggplot2 links:
#http://docs.ggplot2.org/0.9.3.1/geom_point.html

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
allgenes<-c(nonuniqueASD_emsblID,uniqueASD_emsblID)  #1208 genes
pnps_file_1000Gafricanpops$Ensembl_Gene_ID<-as.character(pnps_file_1000Gafricanpops$Ensembl_Gene_ID)

pnps_file_1000GafricanpopsASD<-subset(pnps_file_1000Gafricanpops, pnps_file_1000Gafricanpops$Ensembl_Gene_ID %in% allgenes)

#850 genes with DNDS or PNPS values

for (i in 1:dim(pnps_file_1000GafricanpopsASD)[1]){
  if (pnps_file_1000GafricanpopsASD$Ensembl_Gene_ID[i] %in% nonuniqueASD_emsblID) {pnps_file_1000GafricanpopsASD$ASD_status[i] <- "comorbdASD"} else {pnps_file_1000GafricanpopsASD$ASD_status[i] <- "ASDonly"}
  cat(i,"\t")
}

#ok let's do the dnds boxplot
#save a verison
pnps_file_1000GafricanpopsASD_save<-pnps_file_1000GafricanpopsASD
#one verison zero instead of na because it's zero divided by a number
pnps_file_1000GafricanpopsASD[is.na(pnps_file_1000GafricanpopsASD)] <- 0

------------------------------------------------------------------------------------------------------------------------------------------
  #DNDS
------------------------------------------------------------------------------------------------------------------------------------------
#one verison without any modification
pnps_file_1000GafricanpopsASD1<- pnps_file_1000GafricanpopsASD

#Lets try with no changes at all DNDS

update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
p2<-ggplot(pnps_file_1000GafricanpopsASD1, aes(ASD_status,DNDS))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)
p2 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 

#a few warning: Warning messages: due to NaN and Inf
# Removed 403 rows containing non-finite values (stat_boxplot).
# Removed 376 rows containing missing values (geom_point).

wilcox.test( DNDS ~ ASD_status, data=pnps_file_1000GafricanpopsASD1)

#Wilcoxon rank sum test with continuity correction $ not significatif but lots of Nan and Inf problem
#Wilcoxon rank sum test with continuity correction
#data:  DNDS by ASD_status
#W = 76363.5, p-value = 0.7743
#alternative hypothesis: true location shift is not equal to 0

#let's try removing all Inf 
#remove the one with infinite

inf_pb<-pnps_file_1000GafricanpopsASD$DNDS == "Inf"
noinf_pnps_file_1000GafricanpopsASD<-pnps_file_1000GafricanpopsASD[!inf_pb,]
#test
wilcox.test( DNDS ~ ASD_status, data= noinf_pnps_file_1000GafricanpopsASD)

#  Of course: Wilcoxon rank sum test with continuity correction
#data:  DNDS by ASD_status
#W = 76363.5, p-value = 0.7743
#alternative hypothesis: true location shift is not equal to 0

#let's remove zero and inf
no_zero_pnps_file_1000GafricanpopsASD2<-pnps_file_1000GafricanpopsASD_save

#one verison zero instead of na because it's zero divided by a number
inf_pb<-no_zero_pnps_file_1000GafricanpopsASD2$DNDS == "Inf"
no_zero_pnps_file_1000GafricanpopsASD2<-no_zero_pnps_file_1000GafricanpopsASD2[!inf_pb,]
zero<-is.na(no_zero_pnps_file_1000GafricanpopsASD2$DNDS)
no_zero_pnps_file_1000GafricanpopsASD2<-no_zero_pnps_file_1000GafricanpopsASD2[!zero,]

#test
wilcox.test( DNDS ~ ASD_status, data=no_zero_pnps_file_1000GafricanpopsASD2)

#still not significant 
#Wilcoxon rank sum test with continuity correction
#data:  DNDS by ASD_status
#W = 22827.5, p-value = 0.4562
#alternative hypothesis: true location shift is not equal to 0

#ok and now let's add 1 to avoid all the problems of zero in the ratios
#now add one to eveyry one

pnps_file_1000GafricanpopsASD2<-pnps_file_1000GafricanpopsASD_save
dim(pnps_file_1000GafricanpopsASD2)

pnps_file_1000GafricanpopsASD2$PS<-pnps_file_1000GafricanpopsASD2$PS + 1
pnps_file_1000GafricanpopsASD2$PN<-pnps_file_1000GafricanpopsASD2$PN + 1
pnps_file_1000GafricanpopsASD2$DS<-pnps_file_1000GafricanpopsASD2$DS + 1
pnps_file_1000GafricanpopsASD2$DN<-pnps_file_1000GafricanpopsASD2$DN + 1

pnps_file_1000GafricanpopsASD2$PNPS<-(pnps_file_1000GafricanpopsASD2$PN /pnps_file_1000GafricanpopsASD2$PS)
pnps_file_1000GafricanpopsASD2$DNDS<-(pnps_file_1000GafricanpopsASD2$DN /pnps_file_1000GafricanpopsASD2$DS)

#test
wilcox.test( DNDS ~ ASD_status, data= pnps_file_1000GafricanpopsASD2)

#plot of +1
update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
p2<-ggplot(pnps_file_1000GafricanpopsASD2, aes(ASD_status,DNDS))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)
p2 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 

#save as dsdn_adding1_1000genomes

#density graph
blou<-ggplot(pnps_file_1000GafricanpopsASD2, aes(DNDS, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)
blou + scale_fill_manual( values = c("green","yellow"))
#save as dnds_densitygraph_addingone_1000genomes
#not significatif
#Wilcoxon rank sum test with continuity correction
#data:  DNDS by ASD_status
#W = 85905, p-value = 0.2513
#alternative hypothesis: true location shift is not equal to 0

------------------------------------------------------------------------------------------------------------------------------------------
  #PSPN
------------------------------------------------------------------------------------------------------------------------------------------

#all data Nan -> zero
  wilcox.test(PNPS ~ ASD_status, data=pnps_file_1000GafricanpopsASD1)
  #not significatif
#Wilcoxon rank sum test with continuity correction
#data:  PNPS by ASD_status
#W = 77548, p-value = 0.982
#alternative hypothesis: true location shift is not equal to 0

  
#remove InF
inf_pb_pnps<-pnps_file_1000GafricanpopsASD_save$PNPS == "Inf"
pnps_file_1000GafricanpopsASD_save_noinfPNPS<-pnps_file_1000GafricanpopsASD_save_noinfPNPS[!inf_pb,]
#test
wilcox.test( PNPS ~ ASD_status, data= pnps_file_1000GafricanpopsASD_save_noinfPNPS)
#significant; W = 22838, p-value = 0.05449

#remove Inf and Nan
#test
zero<-is.na(pnps_file_1000GafricanpopsASD_save_noinfPNPS$PNPS)
no_zero_pnps_file_1000GafricanpopsASD_save_noinfPNPS<-pnps_file_1000GafricanpopsASD_save_noinfPNPS[!zero,]
wilcox.test( PNPS ~ ASD_status, data=no_zero_pnps_file_1000GafricanpopsASD_save_noinfPNPS)

#W = 22838, p-value = 0.05449 YES

#adding 1
#test
wilcox.test( PNPS ~ ASD_status, data= pnps_file_1000GafricanpopsASD2)
#yes significatif
#W = 88930.5, p-value = 0.037
#alternative hypothesis: true location shift is not equal to 

#plot
update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
p2<-ggplot(pnps_file_1000GafricanpopsASD2, aes(ASD_status,PNPS))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)
p2 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 

#density graph
blou<-ggplot(pnps_file_1000GafricanpopsASD2, aes(PNPS, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)
blou + scale_fill_manual( values = c("green","yellow"))
#saveas pnps_density_1000genomes
#save as pnps_adding1_1000genomes

#------------------------------------------------------------------------------------------------------------------------------------
 # McDonald–Kreitman test intrespecies with 1000 genomes
#------------------------------------------------------------------------------------------------------------------------------------
#can't do the equation 1- dnpn/dspn so remove all the NaN and InF
no_zero_pnps_file_1000GafricanpopsASD2 #no NaN anf zero for dnds, remove the pnps as well
#pnps_file_1000GafricanpopsASD[is.na(pnps_file_1000GafricanpopsASD)] <- 0
nonNanaPNPS<-no_zero_pnps_file_1000GafricanpopsASD2$PNPS == "Nan"
no_zero_pnps_file_1000GafricanpopsASD3<-no_zero_pnps_file_1000GafricanpopsASD2[!nonNanaPNPS,]
noInfPNPS<-no_zero_pnps_file_1000GafricanpopsASD3$PNPS == "Inf"
no_zero_pnps_file_1000GafricanpopsASD3<-no_zero_pnps_file_1000GafricanpopsASD3[!noInfPNPS,]

#lefy 427 genes
no_zero_pnps_file_1000GafricanpopsASD3$McDonald<-(1-(no_zero_pnps_file_1000GafricanpopsASD3$DNDS*no_zero_pnps_file_1000GafricanpopsASD3$PNPS))


#need to remove the PNPSand above too because I'm STUPID!!!!
wilcox.test( McDonald ~ ASD_status, data= no_zero_pnps_file_1000GafricanpopsASD3)

#yeah significant:W = 15920.5, p-value = 0.1061

p2<-ggplot(no_zero_pnps_file_1000GafricanpopsASD3, aes(ASD_status,McDonald))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)
p2 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 

#density graph
blou<-ggplot(no_zero_pnps_file_1000GafricanpopsASD3, aes(McDonald, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)
blou + scale_fill_manual( values = c("green","yellow"))



-------------------------------with everything
#------------------------------------------------------------------------------------------------------------------------------------
# McDonald–Kreitman test intrespecies with primate for dnds and 1000 genomes for pnps
#------------------------------------------------------------------------------------------------------------------------------------

#I need ot merge both df
head(primates)
head(pnps_file_1000GafricanpopsASD_save)
#merge 
primate2<-primates
colnames(primate2)<-c("ASD_status_onlyornot","Ensembl_Gene_ID","dnds_primate")
dim(primate2)
primates_1000gen <- merge(primate2,pnps_file_1000GafricanpopsASD_save,by="Ensembl_Gene_ID")
primates_1000save<-primates_1000gen

#cleaning up
primates_1000gen <-primates_1000gen [,-10]
primate2<-primates_save<-primate2<-primates
#more lcean up
noInfPNPS<-primates_1000gen$PNPS == "Inf"
primates_1000gen<-primates_1000gen[!noInfPNPS,]
noNaN<-primates_1000gen$PNPS == "NaN"
primates_1000gen<-primates_1000gen[!noNaN,]
dim(primates_1000gen) #439 genes 

#add McDonald–Kreitman for 
primates_1000gen$McDonald_primate<-(1-(primates_1000gen$dnds_primate*primates_1000gen$PNPS))


#now the test 
wilcox.test( McDonald_primate ~ ASD_status, data= primates_1000gen)

# ok significatif W = 17724.5, p-value = 0.01258

#quick fix
colnames(primates_1000gen)<-c("Ensembl_Gene_ID","ASD_status","dnds_primate","PN","PS","DN","DS","PNPS","DNDS","McDonald_primate")

#super significative #W = 24014, p-value = 0.01059
p2<-ggplot(primates_1000gen, aes(ASD_status,McDonald_primate))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)
p2 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 

#density graph 

#density graph
blou<-ggplot(primates_1000gen, aes(McDonald_primate, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)
blou + scale_fill_manual( values = c("green","yellow"))



#---------------------------------------let do the same but on everything
primates_1000save
primates_1000save$McDonald_primate<-(1-(primates_1000save$dnds_primate*primates_1000save$PNPS))


wilcox.test( McDonald_primate ~ ASD_status, data= primates_1000save)
#still pretty good significantE W = 17724.5, p-value = 0.01258
p2<-ggplot(primates_1000save, aes(ASD_status,McDonald_primate))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)
p2 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 

#density graph
blou<-ggplot(primates_1000save, aes(McDonald_primate, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)
blou + scale_fill_manual( values = c("green","yellow"))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#any genes above zero?
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(sort(primates_1000save$McDonald_primate))



