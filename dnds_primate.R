--------------------------------------------------------------------------------------------------------------------------------------------
#6 PRIMATES ANALYSIS: files from David Enar, Petrov's lab
--------------------------------------------------------------------------------------------------------------------------------------------
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

colnames(primates)<-c("ASD_status_onlyornot","genes_names","dnds_primate")

#Other boxplot
update_geom_defaults("point", list(colour = NULL)) #to allo the dots to be the same color than the rest, need ot reset at the end update_geom_defaults("point", list(colour = "black"))
p1<-ggplot(primates, aes(ASD_status_onlyornot, dnds_primate))+
  geom_boxplot(aes(colour = factor(ASD_status_onlyornot)), fill = c("green", "yellow"), color = c("darkgreen","orange"), outlier.colour = NULL, outlier.size = 4, outlier.shape = 1)

p1 + geom_point(aes(colour = factor(ASD_status_onlyornot)), size = I(5), alpha = I(0.1), position = position_jitter(width = 0.4) ) + 
  scale_color_manual(values=c("darkgreen","orange")) 

#save as boxplot_dnds_density_graph

#test
wilcox.test( dnds_primate ~ ASD_status_onlyornot, data= primates)

#Wilcoxon rank sum test with continuity correction
#data:  dnds_primate by ASD_status_onlyornot
#W = 75151, p-value = 0.06528  #YEEAH! 
#alternative hypothesis: true location shift is not equal to 0

blou<-ggplot(primates, aes(dnds_primate, fill = ASD_status_onlyornot)) +
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.5)
blou + scale_fill_manual( values = c("green","yellow"))

#save as dnds_primares_density_graph
