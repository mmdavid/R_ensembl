------------------------------------------------------------------------------------------------------------------------------------------------
  #API ENSEMBL CONVERT SYMBOLS TO ENSEMBL IDs
  --------------------------------------------------------------------------------------------------------------------------------------------
  #OK so David sent me the value for several primate dn/ds and african american let's do this again
  #I need ot grab all the names of Ensembl into my 
  to_search<-attributes(mart.hs)
grep("ensembl", to_search)

to_search[1]
#let's try this attribute

getBM(attributes = "ensembl_gene_id", filters = "hgnc_symbol", values = "GABRA6", mart = mart.hs)

#that works! BIM! 
uniqueASD_emsblID=c()
uniqueASD_used=c()
for (i in 1:length(uniqueASD)){
  a<-getBM(attributes = "ensembl_gene_id", filters = "hgnc_symbol", values = uniqueASD[i], mart = mart.hs)
  uniqueASD_emsblID<-c(uniqueASD_emsblID,a)
  cat(i,"\t")
  b<-uniqueASD[i]
  uniqueASD_used<-c(uniqueASD,b)
}

names(uniqueASD_emsblID)<-uniqueASD
save(uniqueASD_emsblID,file="uniqueASD_emsblID.Rda")

length(uniqueASD)
#331

length(uniqueASD_emsblID) #but problem: several per list

tocorrectensembl=c()
for (i in 1:length(uniqueASD_emsblID)){
  a=c()
  if (is.character(uniqueASD_emsblID[[i]][1]) == 0) {a<-uniqueASD_emsblID[i]} 
  if (length(uniqueASD_emsblID[[i]]) > 1)  {a<-uniqueASD_emsblID[i]}
  tocorrectensembl<-c(tocorrectensembl,a)
}

#ok let's check the problems
tocorrectensembl

#$NOSTRIN
#[1] "ENSG00000275326" "ENSG00000163072"
# take ENSG00000163072 because th eother is less complete and don't have uniprot ID
uniqueASD_emsblID[names(uniqueASD_emsblID) == "NOSTRIN" ] [[1]] <- "ENSG00000163072"


#$CYP21A2
#[1] "ENSG00000235134" => ok "ENSG00000198457"=> ok "ENSG00000233151"=> ok  "ENSG00000206338"=> ok "ENSG00000232414" => ok "ENSG00000231852" => ok idea why there is so many IDs, just keep the first one
uniqueASD_emsblID[names(uniqueASD_emsblID) == "CYP21A2" ] [[1]] <- "ENSG00000235134"

#$PAX6
#[1] "ENSG00000007372" "LRG_720"  LRG even if it's the fixed transcript according to assembl can't be find in Deavid's file so I'll leep the ENS id      
uniqueASD_emsblID[names(uniqueASD_emsblID) == "PAX6" ] [[1]] <- "ENSG00000007372"

#$MAP2K2
#[1] "ENSG00000126934" "LRG_750"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "MAP2K2" ] [[1]] <- "ENSG00000126934"


#$AGRN
#[1] "ENSG00000188157" "LRG_198"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "AGRN" ] [[1]] <- "ENSG00000188157"

#$HSD17B10
#[1] "ENSG00000072506" "LRG_450"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "HSD17B10" ] [[1]] <- "ENSG00000072506"

#$KCNJ2
#[1] "ENSG00000123700" "LRG_328"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "KCNJ2" ] [[1]] <- "ENSG00000123700"

#$RS1
#[1] "ENSG00000102104" "LRG_702"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "RS1" ] [[1]] <- "ENSG00000102104"

#$CHD7
#[1] "ENSG00000171316" "LRG_176"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "CHD7" ] [[1]] <- "ENSG00000171316"

#$DOLK
#[1] "ENSG00000175283" "LRG_744"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "DOLK" ] [[1]] <- "ENSG00000175283"

#$DPYD
#[1] "ENSG00000188641" "LRG_722"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "DPYD" ] [[1]] <- "ENSG00000188641"

#$EPHB6
#[1] "ENSG00000275482" "ENSG00000106123"
uniqueASD_emsblID[names(uniqueASD_emsblID) == "EPHB6" ] [[1]] <- "ENSG00000275482"

#$EXT1
#[1] "ENSG00000182197" "LRG_493"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "EXT1" ] [[1]] <- "ENSG00000182197"

#$GAN
#[1] "ENSG00000261609" "LRG_242"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "GAN" ] [[1]] <- "ENSG00000261609"

#$PLN
#[1] "ENSG00000198523" "LRG_390"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "PLN" ] [[1]] <- "ENSG00000198523"

#$UBR7
#[1] "ENSG00000278787" "ENSG00000012963"
uniqueASD_emsblID[names(uniqueASD_emsblID) == "UBR7" ] [[1]] <- "ENSG00000278787"

#$XPC
#[1] "ENSG00000154767" "LRG_472"        
uniqueASD_emsblID[names(uniqueASD_emsblID) == "XPC" ] [[1]] <- "ENSG00000154767"

#$MIB1
#[1] "ENSG00000101752" "LRG_759"   
uniqueASD_emsblID[names(uniqueASD_emsblID) == "MIB1" ] [[1]] <- "ENSG00000101752"

tocorrectensembl2=c()
for (i in 1:length(uniqueASD_emsblID)){
  a=c()
  if (is.na(uniqueASD_emsblID[[i]][1])) {a<-uniqueASD_emsblID[i]} 
  tocorrectensembl2<-c(tocorrectensembl2,a)
}

#who are they? 
tocorrectensembl2

#$CAC1I ither name  CACNA1I, then ensembl gives me ENSG00000100346
#character(0)
uniqueASD_emsblID[names(uniqueASD_emsblID) == "CAC1I" ] [[1]] <- "ENSG00000100346"

#$DJC19 a couple of papers indicates it's DNADJ1 in yest not human, I'll removed it
#character(0)
uniqueASD_emsblID[names(uniqueASD_emsblID) == "DJC19"] <- NULL

#$JMJD7.PLA2G4B just changed for JMJD7.PLA2G4B and now ENSG00000168970
#character(0)
uniqueASD_emsblID[names(uniqueASD_emsblID) == "JMJD7.PLA2G4B" ] [[1]] <- "ENSG00000168970"

#$MSNP1AS psuedogene not in ensembl probably, will remove it
#character(0)
uniqueASD_emsblID[names(uniqueASD_emsblID) == "MSNP1AS"] <- NULL

#$MLL3 KEGG tells me it's KMT2C, ENSG00000055609
#character(0)
uniqueASD_emsblID[names(uniqueASD_emsblID) == "MLL3" ] [[1]] <- "ENSG00000055609"

#$CD42BPB can't find it anywhere need to remove it weird int he nature paper through
#character(0)

uniqueASD_emsblID[names(uniqueASD_emsblID) == "CD42BPB"] <- NULL

#$DSCAMF3 impossible to find either
#character(0)
uniqueASD_emsblID[names(uniqueASD_emsblID) == "DSCAMF3"] <- NULL

#$KAT another ghosts 
#character(0)
uniqueASD_emsblID[names(uniqueASD_emsblID) == "KAT"] <- NULL 

#$AL2 my last gene ghost ghost
#character(0)
uniqueASD_emsblID[names(uniqueASD_emsblID) == "AL2"] <- NULL

#looking into the list
uniqueASD_emsblID<-unlist(uniqueASD_emsblID)
names(uniqueASD_emsblID)<-c()


#ok involved in other co-morbid now

nonuniqueASD_emsblID=c()
for (i in 1:length(nonuniqueASD)){
  a<-getBM(attributes = "ensembl_gene_id", filters = "hgnc_symbol", values = nonuniqueASD[i], mart = mart.hs)
  nonuniqueASD_emsblID<-c(nonuniqueASD_emsblID,a)
  cat(i,"\t")
}

names(nonuniqueASD_emsblID)<-nonuniqueASD

tocorrectensembl2=c()
for (i in 1:length(nonuniqueASD_emsblID)){
  a=c()
  if (is.na(nonuniqueASD_emsblID[[i]][1])) {a<-nonuniqueASD_emsblID[i]} 
  tocorrectensembl2<-c(tocorrectensembl2,a)
}

#$AUTS5 : this paper http://www.ncbi.nlm.nih.gov/pubmed/19401682 just loci
#character(0)
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "AUTS5"] <- NULL

#$NOS2A A is a subunit: the NOS2 is ENSG00000007171
#character(0)
nonuniqueASD_emsblID[names(uniqueASD_emsblID) == "NOS2A" ] [[1]] <- "ENSG00000007171"

#$PFTK1 synonym to CDK14 ENSG00000058091
#character(0)
nonuniqueASD_emsblID[names(uniqueASD_emsblID) == "PFTK1" ] [[1]] <- "ENSG00000058091"

#$PRKCB1 same as PRKCB ENSG00000166501
#character(0)
nonuniqueASD_emsblID[names(uniqueASD_emsblID) == "PRKCB1" ] [[1]] <- "ENSG00000166501"

#$`PCDHA@` that means the whole gene cluster, too vague, I remove it 
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "PCDHA@"] <- NULL
#character(0)

#$`HOXB@` same gene cluster nothing precise
#character(0)
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "HOXB@"] <- NULL

#$`HOXD@`  same gene cluster nothing precise
#character(0)
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "HOXD@@"] <- NULL

#$CAC1C gene cards say same as CACNA1C  ENSG00000151067
#character(0)
nonuniqueASD_emsblID[names(uniqueASD_emsblID) == "CAC1C" ] [[1]] <- "ENSG00000151067"

#$CHR7 human alignment ? weird 
#character(0)
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "CHR7"] <- NULL

#$CNTP2 does not exist
#character(0)
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "CNTP2"] <- NULL

#$HLA.DRB1 ENSG00000196126
#character(0)
nonuniqueASD_emsblID[names(uniqueASD_emsblID) == "HLA.DRB1" ] [[1]] <- "ENSG00000196126"

#$IL8 CXCL8 ENSG00000169429
#character(0)
nonuniqueASD_emsblID[names(uniqueASD_emsblID) == "IL8" ] [[1]] <- "ENSG00000169429"

#$IV can;t find it
#character(0)
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "IV "] <- NULL

#$SP25 anohter ghost
#character(0)
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "SP25 "] <- NULL

tocorrectensembl=c()
for (i in 1:length(nonuniqueASD_emsblID)){
  a=c()
  if (length(nonuniqueASD_emsblID[[i]]) > 1)  {a<-nonuniqueASD_emsblID[i]}
  tocorrectensembl<-c(tocorrectensembl,a)
}


#pfiou 75 problems

#$MTHFR
#[1] "ENSG00000177000" "LRG_726"        
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "MTHFR" ] [[1]] <- "ENSG00000177000"

#$PTEN
#[1] "ENSG00000171862" "LRG_311"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "PTEN" ] [[1]] <- "ENSG00000171862"

#$UBE3A
#[1] "ENSG00000114062" "LRG_15"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "UBE3A" ] [[1]] <- "ENSG00000114062"

#$MET
#[1] "ENSG00000105976" "LRG_662"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "MET" ] [[1]] <- "ENSG00000105976"

#$`HLA-DRB1` the correct one is ENSG00000196126 the rest are subunit
#[1] "ENSG00000206306" "ENSG00000206240" "ENSG00000196126" "ENSG00000229074" "ENSG00000236884" "ENSG00000228080"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "HLA-DRB1" ] [[1]] <- "ENSG00000196126"

#$SLC6A3 ENSG00000142319 correct rest is subfamily of carrier
#[1] "ENSG00000142319" "ENSG00000276996"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "SLC6A3" ] [[1]] <- "ENSG00000142319"

#$CHRNA7 correct is ENSG00000175344 rest subunit
#[1] "ENSG00000175344" "ENSG00000274542"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "CHRNA7" ] [[1]] <- "ENSG00000175344"

#$ADA
#[1] "ENSG00000196839" "LRG_16"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "ADA" ] [[1]] <- "ENSG00000196839"

#$TNF
#[1] "ENSG00000223952":alternative "ENSG00000204490":alternative "ENSG00000232810": not alternative good one "ENSG00000228978":alternative "ENSG00000206439": can t find it "ENSG00000228849":alternative "ENSG00000230108":ENSG00000230108
#[8] "ENSG00000228321": alternative
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TNF" ] [[1]] <- "ENSG00000232810"

#$TCN2
#[1] "ENSG00000185339" "LRG_116"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TCN2" ] [[1]] <- "ENSG00000185339"

#$TSC1
#[1] "ENSG00000165699" "LRG_486"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TSC1" ] [[1]] <- "ENSG00000165699"

#$TSC2
#[1] "ENSG00000103197" "LRG_487"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TSC2" ] [[1]] <- "ENSG00000103197"

#$RAF1
#[1] "ENSG00000132155" "LRG_413"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "RAF1" ] [[1]] <- "ENSG00000132155"

#$RXRB
#[1] "ENSG00000235712" alternative "ENSG00000228333"alternative "ENSG00000227322"alternative "ENSG00000204231" good one "ENSG00000206289"alternative "ENSG00000231321"alternative
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "RXRB" ] [[1]] <-"ENSG00000204231"

#$NSD1
#[1] "ENSG00000165671" "LRG_512"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "NSD1" ] [[1]] <- "ENSG00000165671"

#$SCG5
#[1] "ENSG00000277614"alternative  "ENSG00000166922" good
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "SCG5" ] [[1]] <- "ENSG00000277614"

#$TBX1
#[1] "ENSG00000184058" "LRG_226"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TBX1" ] [[1]] <- "ENSG00000184058"

#$SLC12A6
[1] "ENSG00000140199" "LRG_270"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "SLC12A6" ] [[1]] <- "ENSG00000140199"

#$C4B
[1] "ENSG00000228454"alternative "ENSG00000228267" alternative"ENSG00000224389" good "ENSG00000224639"alternative"ENSG00000236625"alternative "LRG_138"        
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "C4B" ] [[1]] <- "ENSG00000224389"

#$CBS
[1] "ENSG00000274276" "ENSG00000160200":seems to be an archive version "LRG_777"        
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "CBS" ] [[1]] <- "ENSG00000274276"

#$SYNGAP1
[1] "ENSG00000227460" alternative "ENSG00000197283"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "SYNGAP1" ] [[1]] <- "ENSG00000197283"

#$NF1
#[1] "ENSG00000196712" "LRG_214"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "NF1" ] [[1]] <- "ENSG00000196712"

#$NGF
#[1] "ENSG00000134259" "LRG_260"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "NGF" ] [[1]] <- "ENSG00000134259"

#$OPRL1
#[1] "ENSG00000277044" alternative "ENSG00000125510"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "OPRL1" ] [[1]] <- "ENSG00000125510"

#$NTRK1
[1] "ENSG00000198400" "LRG_261"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "NTRK1" ] [[1]] <- "ENSG00000198400"

#$PIK3CA
#[1] "ENSG00000121879" "LRG_310"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "PIK3CA" ] [[1]] <- "ENSG00000121879"

#$MAP2K1
#[1] "ENSG00000169032" "LRG_725"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "MAP2K1" ] [[1]] <- "ENSG00000169032"

#$GSTP1
#[1] "ENSG00000084207" "LRG_723"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "GSTP1" ] [[1]] <- "ENSG00000084207"

#$HFE
#[1] "ENSG00000010704" "LRG_748"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "HFE" ] [[1]] <- "ENSG00000010704"

#$`HLA-A`: lots of subunit:ENSG00000206503 
#[1] "ENSG00000206505" "ENSG00000229215" "ENSG00000224320" "ENSG00000206503" "ENSG00000235657" "ENSG00000227715" "ENSG00000223980"
#[8] "ENSG00000231834"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "HLA-A" ] [[1]] <- "ENSG00000206503"

#$`HLA-B` ENSG00000234745
#[1] "ENSG00000232126" "ENSG00000228964" "ENSG00000234745" "ENSG00000224608" "ENSG00000223532" "ENSG00000206450"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "HLA-B" ] [[1]] <- "ENSG00000234745"

#$APBA2
#[1] "ENSG00000276495" "ENSG00000034053"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "APBA2" ] [[1]] <-"ENSG00000034053"

#$APC
#[1] "ENSG00000134982" "LRG_130"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "APC" ] [[1]] <- "ENSG00000134982"

#$HRAS
#[1] "ENSG00000276536" alternate "ENSG00000174775"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "HRAS" ] [[1]]<- "ENSG00000174775"

#$IL1RN
#[1] "ENSG00000136689" "LRG_188"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "IL1RN" ] [[1]] <- "ENSG00000136689"

#$MAPT
#[1] "ENSG00000186868" good "ENSG00000276155"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "MAPT" ] [[1]] <- "ENSG00000186868"

#$MOG
#[1] "ENSG00000232641" "ENSG00000236561" "ENSG00000204655" "ENSG00000234623" "ENSG00000237834" "ENSG00000137345" "ENSG00000234096"
#[8] "ENSG00000230885"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "MOG" ] [[1]] <-"ENSG00000204655"

#$KIR2DS1
#[1] "ENSG00000276327" "ENSG00000273603" "ENSG00000278304" "ENSG00000275421" "ENSG00000273517" "ENSG00000275921" "ENSG00000276387"
#[8] "ENSG00000278120" "ENSG00000275306" ENSG00000275306: only alternate ones, but the last one is the one given by default by ensembl
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "KIR2DS1" ] [[1]] <-"ENSG00000275306"

#$AVP
#[1] "ENSG00000101200" "LRG_715"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "AVP" ] [[1]] <- "ENSG00000101200"

#$BRD2
#[1] "ENSG00000230678" "ENSG00000234704" "ENSG00000204256" "ENSG00000234507" "ENSG00000215077" "ENSG00000236227" "ENSG00000235307"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "BRD2" ] [[1]] <- "ENSG00000204256"

#$CCL5
#[1] "ENSG00000274233": alternate "ENSG00000271503": right
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "CCL5" ] [[1]] <- "ENSG0000027150"

#$CD40LG
[1] "ENSG00000102245" "LRG_141"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "CD40LG" ] [[1]] <- "ENSG00000102245"

#$CYFIP1
#[1] "ENSG00000280618" "ENSG00000273749": good
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "CYFIP1" ] [[1]] <- "ENSG00000273749"

#$DMD
#[1] "ENSG00000198947" "LRG_199"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "DMD" ] [[1]] <- "ENSG00000198947"

$DRD4
[1] "ENSG00000276825" "ENSG00000069696"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "DRD4" ] [[1]]<- "ENSG00000069696"

#$GATA2
[1] "ENSG00000179348" "LRG_295"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "GATA2" ] [[1]] <- "ENSG00000179348"

#$NRAS
#[1] "ENSG00000213281" "LRG_92"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "NRAS" ] [[1]] <- "ENSG00000213281"

#$NSF
[1] "ENSG00000278174" "ENSG00000073969"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "NSF" ] [[1]] <- "ENSG00000073969"

#$SCN1A
#[1] "ENSG00000144285" "LRG_8"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "SCN1A" ] [[1]] <- "ENSG00000144285"

#$SOD1
#[1] "ENSG00000142168" "LRG_652"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "SOD1" ] [[1]] <- "ENSG00000142168"

#$TPO
[1] "ENSG00000277603" "ENSG00000115705"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TPO" ] [[1]] <- "ENSG00000115705"

#$TRIM32
[1] "ENSG00000119401" "LRG_211"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TRIM32" ] [[1]] <- "ENSG00000119401"

#$ADORA3
[1] "ENSG00000121933" "LRG_424"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "ADORA3" ] [[1]] <- "ENSG00000121933"

#$AFF2
[1] "ENSG00000281817" "ENSG00000155966"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "AFF2" ] [[1]] <- "ENSG00000155966"

#$AMT
#[1] "ENSG00000145020" "LRG_537"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "AMT" ] [[1]] <- "ENSG00000145020"

#$ANK2
#[1] "ENSG00000145362" "LRG_327"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "ANK2" ] [[1]] <- "ENSG00000145362"        

#$BRAF
#[1] "ENSG00000157764" "LRG_299"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "BRAF" ] [[1]] <- "ENSG00000157764"

#$BRCA2
#[1] "ENSG00000139618" "LRG_293"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "BRCA2" ] [[1]] <- "ENSG00000139618"

#$CACNA1C
[1] "ENSG00000151067" "LRG_334"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "CACNA1C" ] [[1]] <- "ENSG00000151067"

#$DHCR7
#[1] "ENSG00000172893" "LRG_340"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "DHCR7" ] [[1]] <- "ENSG00000172893"

#$EGR2
[1] "ENSG00000122877" "LRG_239"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "EGR2" ] [[1]] <- "ENSG00000122877"       

#$FGA
#[1] "ENSG00000171560" "LRG_557"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "FGA" ] [[1]] <- "ENSG00000171560"        

#$HERC2
[1] "ENSG00000128731" "ENSG00000277278"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "HERC2" ] [[1]] <- "ENSG00000128731" 

#$KCNQ2
#[1] "ENSG00000281151" "ENSG00000075043"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "KCNQ2" ] [[1]] <- "ENSG00000075043" 

#$KIF5C
#[1] "ENSG00000168280" "ENSG00000276734"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "KIF5C" ] [[1]] <- "ENSG00000168280"


#$KIT
#[1] "ENSG00000157404" "LRG_307"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "KIT" ] [[1]] <- "ENSG00000157404"        

#$PTPN11
#[1] "ENSG00000179295" "LRG_614"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "PTPN11" ] [[1]] <- "ENSG00000179295"        

#$PTPRC
#[1] "ENSG00000262418" "ENSG00000081237"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "PTPRC" ] [[1]] <- "ENSG00000081237" 

#$SNTG2
#[1] "ENSG00000281486" "ENSG00000172554" "ENSG00000281020"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "SNTG2" ] [[1]] <- "ENSG00000172554" 

#$TDO2
[1] "ENSG00000262635" "ENSG00000151790"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TDO2" ] [[1]] <- "ENSG00000151790" 

#$TRPM1
[1] "ENSG00000274965" "ENSG00000134160"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TRPM1" ] [[1]] <- "ENSG00000134160" 

#$TTN
[1] "ENSG00000155657" "LRG_391"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TTN" ] [[1]] <- "ENSG00000155657"        

#$TUBGCP5
#[1] "ENSG00000280807" "ENSG00000276856" "ENSG00000275835"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "TUBGCP5" ] [[1]] <- "ENSG00000275835"

#$VPS13B
#[1] "ENSG00000132549" "LRG_351"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "VPS13B" ] [[1]] <- "ENSG00000132549"        

#$YWHAE
[1] "ENSG00000108953" "ENSG00000274474"
nonuniqueASD_emsblID[names(nonuniqueASD_emsblID) == "YWHAE" ] [[1]] <- "ENSG00000108953"       


save(nonuniqueASD_emsblID,file="nonuniqueASD_emsblID.Rda")

nonuniqueASD_emsblID<- nonuniqueASD_emsblID[1:606]
nonuniqueASD_emsblID<-unlist(nonuniqueASD_emsblID)

#last check 
tocorrectensembl2=c()
for (i in 1:length(nonuniqueASD_emsblID)){
  a=c()
  if (is.na(nonuniqueASD_emsblID[[i]][1])) {a<-nonuniqueASD_emsblID[i]} 
  tocorrectensembl2<-c(tocorrectensembl2,a)
}

#ok empty
#tocorrectensembl=c()
tocorrectensembl=c()
for (i in 1:length(nonuniqueASD_emsblID)){
  a=c()
  if (length(nonuniqueASD_emsblID[[i]]) > 1)  {a<-nonuniqueASD_emsblID[i]}
  tocorrectensembl<-c(tocorrectensembl,a)
}
#empty it's clean