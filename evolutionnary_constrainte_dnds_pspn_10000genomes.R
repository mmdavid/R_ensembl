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

#------------------------------------------------------------------------------------------------------------------------------------
# McDonald–Kreitman test intrespecies with primate for dnds and 1000 genomes for pnps
#------------------------------------------------------------------------------------------------------------------------------------

#I need ot merge both df
head(primates)
head(pnps_file_1000GafricanpopsASD_save)
#merge 
colnames(primates)<-c("ASD_status_onlyornot","genes_names","dnds_primate")
primate2<-primates
dim(primate2)
primates_1000gen <- merge(primates,pnps_file_1000GafricanpopsASD_save,by="Ensembl_Gene_ID")
#cleaning up
primates_1000gen <-primates_1000gen [,-10]
#more lcean up
noInfPNPS<-primates_1000gen$PNPS == "Inf"
primates_1000gen<-primates_1000gen[!noInfPNPS,]
primates_1000gen[is.na(primates_1000gen)] <- 0

#add McDonald–Kreitman for 



#ok how about the polymorphism ------------------------------------------------------------------

#ok let's do the dnds boxplot


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


#--------- adding 1



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

--------------------------------------------------------------------------------------------------------------------------------------
  #dnds all I can do now
  wilcox.test( DNDS ~ ASD_status, data= pnps_file_1000GafricanpopsASD3)

#not significative 
#Wilcoxon rank sum test with continuity correction
#data:  DNDS by ASD_status
#W = 23366, p-value = 0.3547
#alternative hypothesis: true location shift is not equal to 0

p1<-ggplot(pnps_file_1000GafricanpopsASD3, aes(ASD_status,DNDS))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 5, outlier.shape = 1, fill = c("green", "yellow"),)

p1 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 

blou<-ggplot(pnps_file_1000GafricanpopsASD3, aes(PNPS, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow"))

-------------------------------------------------------------------------------------------------------- 
  #NI +1
  pnps_file_1000GafricanpopsASD2$alpha<-1-(pnps_file_1000GafricanpopsASD2$DNDS*pnps_file_1000GafricanpopsASD2$PNPS)
#test andplots
update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
p2<-ggplot(pnps_file_1000GafricanpopsASD2, aes(ASD_status,alpha))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)

p2 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 


blou<-ggplot(noinf_dnda_noinf_pnps_file_1000GafricanpopsASD, aes(NI, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow"))

#cool looking good
#test
wilcox.test( alpha ~ ASD_status, data= pnps_file_1000GafricanpopsASD2)

#Wilcoxon rank sum test with continuity correction
#data:  alpha by ASD_status
#W = 80070.5, p-value = 0.1324
#alternative hypothesis: true location shift is not equal to 0



p2<-ggplot(supzero, aes(ASD_status,alpha))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)

p2 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 

blou<-ggplot(supzero, aes(alpha, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow")) 
+ 
  scale_y_continuous(trans = 'log10')

#same but without adding one and removing the one missing infomration
noinf_dnda_noinf_pnps_file_1000GafricanpopsASD$alpha<-1-(noinf_dnda_noinf_pnps_file_1000GafricanpopsASD$NI)

wilcox.test( alpha ~ ASD_status, data= pnps_file_1000GafricanpopsASD2)
#so if I put eveything its significant, and if I remove the ones than are under zero: 



supzeronoinf<-noinf_dnda_noinf_pnps_file_1000GafricanpopsASD[noinf_dnda_noinf_pnps_file_1000GafricanpopsASD$alpha>0,]
wilcox.test( NI ~ ASD_status, data= supzeronoinf)

p2<-ggplot(supzeronoinf, aes(ASD_status,alpha))+
  geom_boxplot(aes(colour = factor(ASD_status)), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1, fill = c("green", "yellow"),)

p2 + geom_point(aes(colour = factor(ASD_status)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) +
  scale_color_manual(values=c("darkgreen","orange")) 


#it's significatif if you don't add +1 !! which is correct: because under it means that there is violation of the model:
#such as the segregation of slightly deleterious amino acid mutations
#Wilcoxon rank sum test with continuity correction
#data:  NI by ASD_status
#W = 63317.5, p-value = 0.07625
#alternative hypothesis: true location shift is not equal to 0

blou<-ggplot(supzeronoinf, aes(NI, fill = ASD_status)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)

blou + scale_fill_manual( values = c("green","yellow")) + 
  scale_y_continuous(trans = 'log10')
