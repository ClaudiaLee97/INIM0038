##Packages required: bindrcpp, boot, datasets, dplyr, forcats, ggplot2, ggpubr, gplots, graphics, grDevices, magrittr, methods, purrr, readr, rlang, stas, stats4, stringr, survival, survminer, tibble, tidyr, tidyverse, utils, 


##Comparing 95th percentile growth rate in KTC and KPTC tumours
#Read file containing 95th percentile growth rates of KTC and KPTC tumours for each TSG
kptc.ktc.corr <- read.csv(file="Desktop/attemptingCAMP/kptc.ktc.corr.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

#Perform Pearson’s product-moment correlation
cor.test(kptc.ktc.corr$ktc.mean,kptc.ktc.corr$kptc.mean)

#Plot datapoints and separate context-dependent interactions with lines y=1, x=1
b <- ggplot(data=kptc.ktc.corr, aes(x=ktc.mean, y=kptc.mean, colour=sgID))+
  geom_point(data=kptc.ktc.corr, aes(x=ktc.mean, y=kptc.mean, colour=sgID))+
  #geom_smooth(data=kptc.ktc.corr,method="lm",colour="black")+
  geom_vline(xintercept = 1, aes(linetype="dotted"))+
  geom_hline(yintercept = 1, aes(linetype="dotted"))+
  labs(title="Comparing effect of Tp53 on 95th percentile growth rates\nR=0.8254419,p=5.29e-13", y="KPTC tumour growth rate", x="KTC tumour growth rate")+
  theme_classic()


##Organisng TRACERx patient data in preparation for correlation analyses
#Read table containing patient mutation data
mutTableAll <- read.table("Desktop/attemptingCAMP/mutTableAll.CSV", header = TRUE, sep=',')

#Sift out patients with driver mutations only, and patients with Kras mutations
mutTableAll.dr <- data.frame(mutTableAll[mutTableAll$DriverMut %in% 'TRUE',])
mutTableAll.dr.kras <- data.frame(mutTableAll.dr[mutTableAll.dr$Hugo_Symbol %in% 'KRAS',])

#Sift out patients with RAGE mutations 
rage.sgid <- c('ATM', 'NF1', 'PTEN', 'SETD2', 'STK11', 'RNF43') 
pat.rage <- data.frame(mutTableAll.dr[mutTableAll.dr$Hugo_Symbol %in% rage.sgid,]) pat.rage.1 <- data.frame(pat.rage[pat.rage$SampleID %in% sample.dr.kras,])

#Repeat above for patients with CALM, Context-dependent, Growth-attenuating and No significant effect mutations 

#Read table containing patient clinical data and sift out lung cancer patients with lung adenocarcinoma (LUAD) only
luad.clin.list <- read.table("Desktop/attemptingCAMP/luad.clin.list.CSV", header = TRUE, sep=',') 
pat.luad <- luad.clin.list$SampleID

#Create dataframe of LUAD patients with RAGE mutation and Kras background 
#rage.luad.kras.dr <- data.frame(pat.rage.1[pat.rage.1$SampleID %in% pat.luad,])

#Repeat above for patients with CALM, Context-dependent, Growth-attenuating and No significant effect mutations

#Create dataframe of mutation data of LUAD patients with Kras background and 48 candidate TSG mutations
clean.all <- rbind(rage.luad.kras.dr,calm.luad.kras.dr,cd.luad.kras.dr,tsg.luad.kras.dr,nse.luad.kras.dr)


##Organising region-specific patient data
#LUAD patients with driver candidate 48 TSG mutations in Kras-background data separated on a per region basis
mutTableAll_Region = clean.all %>%   separate_rows(RegionSum,sep=";") %>%   separate(RegionSum,into=c("region_id","region_var_ref_count"),sep=":") %>%   separate(region_var_ref_count,into=c("var_count","ref_count"),sep="/") %>%   mutate(SampleID=gsub("\\.","-",paste(site,SampleID,region_id,sep="_"))) mutTableAll_Region.driver<-mutTableAll_Region[mutTableAll_Region$DriverMut %in% TRUE,]


#LUAD patients with driver non-candidate 48 TSG mutations in Kras-background (driver wild-type patients) data separated on a per region basis
mutTableAll.cont.kras.sample <- data.frame(mutTableAll[mutTableAll$Hugo_Symbol %in% 'KRAS',]) 
mutTableAll.cont.kras <- data.frame(mutTableAll[mutTableAll$SampleID %in% mutTableAll.cont.kras.sample$SampleID,]) 
mutTableAll.cont.kras.luad <- data.frame(mutTableAll.cont.kras[mutTableAll.cont.kras$SampleID %in% pat.luad,]) 
'%ni%' <- Negate('%in%') mutTableAll.cont.kras.luad.nonsgid <- data.frame(mutTableAll.cont.kras.luad[mutTableAll.cont.kras.luad$Hugo_Symbol %ni% sgid,]) 
mutTableAll.cont.dr <- 
  data.frame(mutTableAll.cont.kras.luad.nonsgid[mutTableAll.cont.kras.luad.nonsgid$DriverMut %in% 'TRUE',])
mutTableAll_Region.cont.dr = mutTableAll.cont.dr %>%   separate_rows(RegionSum,sep=";") %>%   separate(RegionSum,into=c("region_id","region_var_ref_count"),sep=":") %>%   separate(region_var_ref_count,into=c("var_count","ref_count"),sep="/") %>%   mutate(SampleID=gsub("\\.","-",paste(site,SampleID,region_id,sep="_")))  


## Survival correlation from clinical data
#Read table containing patient overall survival in days data
surv <- as.data.frame(tracerx.surv500)
SampleID <- data.frame(surv$SampleID) 
os.OS <- data.frame(surv$OS_time_days) 
os <- data.frame(cbind(SampleID,os.OS)) 
colnames(os) <- c("SampleID","OS.time.days")

#Sift out survival data of LUAD patients with Kras background and 48 candidate TSG mutations
merged.os <- merge(merged.num,os,by="SampleID") merged.os<-merged.os[complete.cases(merged.os),]

#Classify survival data of patients based on their assigned growth phenotype
rage.os <- data.frame(merged.os.1$RAGE, merged.os.1$OS.time.days) 
names(rage.os) <- c("num.mut", "os.days") 
name <- rep("RAGE",95)
rage.os <- cbind(name,rage.os)

#Repeat above for each growth phenotype, to create combined dataframe of survival data for LUAD patients with Kras background sorted by the growth phenotype associated with their TSG mutation 

#Plot survival probability curve and derive risk table 
table(os.df.na.new$Growth.phenotype) 
surv.obj <- Surv(time=os.df.na.new$OS.days) 
fit.1 <- survfit(surv.obj~Growth.phenotype, data=os.df.na.new) summary(fit.1) 
a<-ggsurvplot(fit.1, data=os.df.na.new, risk.table = TRUE, pval=TRUE)

#Plot cox proportional hazards table 
os.df.na.new$Growth.phenotype <- factor(os.df.na.new$Growth.phenotype, levels = c("Control","RAGE","CALM","Context-dependent", "Growth attenuating", "No significant effect")) os.df.na.new$Growth.phenotype<-as.factor(os.df.na.new$Growth.phenotype) fit.coxph <-coxph(formula=Surv(time=os.df.na.new$OS.days)~ Growth.phenotype, os.df.na.new) ggforest(fit.coxph, data=os.df.na.new)


## Lesion size pathology from clinical data 
#Read table containing patient lesion size pathology data 
lsp.LSP <- data.frame(surv$Lesion1_size_pathology) 
lsp <- data.frame(cbind(SampleID,lsp.LSP)) 
colnames(lsp) <- c("SampleID","LesionSize")

#Sift out lesion size pathology data of LUAD patients with Kras background and 48 candidate TSG mutations
merged.lsp <- merge(merged.num,lsp,by="SampleID") 
merged.lsp<-merged.lsp[complete.cases(merged.lsp),]

#Classify lesion size pathology data of patients based on their assigned growth phenotype
rage.lsp <- data.frame(merged.lsp.1$RAGE, merged.lsp.1$LesionSize) 
names(rage.lsp) <- c("num.mut", "LesionSize") 
name <- rep("RAGE",99) 
rage.lsp <- cbind(name,rage.lsp)

#Repeat above for each growth phenotype, to create combined dataframe of lesion size pathology data for LUAD patients with Kras background sorted by the growth phenotype associated with their TSG mutation 
lsp.df <- data.frame(rbind(nse.lsp,rage.lsp,calm.lsp,cd.lsp,tsg.lsp)) 
lsp.df[lsp.df == 0] <- NA 
lsp.df.na<-lsp.df[complete.cases(lsp.df),] 

#Calculate tumour volume from lesion size pathology data
lsp.df.na$TumourVol <- (lsp.df.na$LesionSize)^3/2

#Plot violin plot of tumour volume
a <- ggplot()+   geom_violin(data=lsp.vol, aes(x=Growth.phenotype, y=TumourVol, color=Growth.phenotype))+   scale_y_log10()+   labs(title="Correlating tumour growth phenotype to tumour volume", y="Tumour volume", x="Tumour growth phenotype")+   theme_classic()


##MKI67, Mitosis-HPF and Percentage Necrosis from clinical data
#Plot boxplot for molecular growth indicators 
mki67 <- read.csv(file="Desktop/attemptingCAMP/mki67.new.csv", header=TRUE, sep=",", stringsAsFactors = FALSE) a <- ggplot(data=mki67, aes(x=Phenotype, y=MKI67, colour=Phenotype))+   geom_boxplot()+   labs(title="Mki67 RNA expression of LUAD tumours with Kras background and TSG mutations", y="Mki67 RNA expression", x="Growth phenotype")+   stat_summary(fun.y=mean)+   theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')

#Repeat above command for Mitosis-HPF and Percentage necrosis 

#Plot heatmap for molecular growth indiactors
rm(mean) mki67.means <- data.frame(aggregate(MKI67~Phenotype, mki67, mean)) mit.means <- data.frame(aggregate(Mitosis~Phenotype, nec.mit, mean)) nec.means <- data.frame(aggregate(Necrosis~Phenotype, nec.mit, mean)) gen.growth <- cbind(mki67.means,mit.means$Mitosis,nec.means$Necrosis) colnames(gen.growth) <- c("Phenotype", "Mki67", "Mitosis", "Necrosis") rownames(gen.growth) <- gen.growth$Phenotype gen.growth <- gen.growth[,2:4] a <- heatmap(as.matrix(gen.growth), scale="column",               Colv=NA, Rowv=NA,              col=cm.colors(256), xlab="Genomic growth indicators", ylab="TSG mutation", main="Heat map of genomic growth indicators",              margins=c(5,10), cexRow = 1, cexCol = 1)


##wGII and wFLOH from genomic instability data 
#Construct dataframes of non-clonal/clonal driver TSG mutations and wildtype mutations
my.mut.clonal <- data.frame(my.mut[my.mut$PyCloneClonal %in% c("C"),])
my.mut.subclonal <- data.frame(my.mut[my.mut$PyCloneClonal %in% c("S"),])
my.mut.short.clonal <- data.frame(cbind(my.mut.clonal$SampleID, my.mut.clonal$Hugo_Symbol))
my.mut.short.subclonal <- data.frame(cbind(my.mut.subclonal$SampleID, my.mut.subclonal$Hugo_Symbol))
colnames(my.mut.short.clonal) <- c("Sample.ID", "Mutation")
colnames(my.mut.short.subclonal) <- c("Sample.ID", "Mutation")
mutTableAll_Region.cont.dr <- read.csv(file="Desktop/attemptingCAMP/mutTableAll_Region.cont.dr.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
my.mut.long.0.clonal <- data.frame(mutTableAll_Region.cont.dr[mutTableAll_Region.cont.dr$PyCloneClonal %in% c("C"),])
my.mut.long.0.subclonal <- data.frame(mutTableAll_Region.cont.dr[mutTableAll_Region.cont.dr$PyCloneClonal %in% c("S"),])
my.mut.long.00.clonal <- data.frame(cbind(my.mut.long.0.clonal$SampleID, my.mut.long.0.clonal$Hugo_Symbol))
my.mut.long.00.subclonal <- data.frame(cbind(my.mut.long.0.subclonal$SampleID, my.mut.long.0.subclonal$Hugo_Symbol))
colnames(my.mut.long.00.clonal) <- c("Sample.ID", "Mutation")
colnames(my.mut.long.00.subclonal) <- c("Sample.ID", "Mutation")
my.mut.long.clonal <- data.frame(rbind(my.mut.short.clonal,my.mut.long.00.clonal))
my.mut.long.subclonal <- data.frame(rbind(my.mut.short.subclonal,my.mut.long.00.subclonal))
merged.gii.floh.clonal <- merge(gii.floh,my.mut.long.clonal,by="Sample.ID")
merged.gii.floh.subclonal <- merge(gii.floh,my.mut.long.subclonal,by="Sample.ID")


#Sift out wGII and wFLOH data for sample regions, and sort them into categories based on their clonal/subclonal growth phenotypes 
merged.gii.floh.clonal$Growth.phenotype<-NA
merged.gii.floh.clonal$Growth.phenotype[merged.gii.floh.clonal$Mutation %ni% c(nse.sgid,rage.sgid,calm.sgid,cd.sgid,tsg.sgid)] <- "Clonal Control"
merged.gii.floh.clonal$Growth.phenotype[merged.gii.floh.clonal$Mutation %in% rage.sgid] <- "Clonal RAGE"
merged.gii.floh.clonal$Growth.phenotype[merged.gii.floh.clonal$Mutation %in% calm.sgid] <- "Clonal CALM"
merged.gii.floh.clonal$Growth.phenotype[merged.gii.floh.clonal$Mutation %in% cd.sgid] <- "Clonal Context-dependent"
merged.gii.floh.clonal$Growth.phenotype[merged.gii.floh.clonal$Mutation %in% tsg.sgid] <- "Clonal Growth attenuating"
merged.gii.floh.clonal$Growth.phenotype[merged.gii.floh.clonal$Mutation %in% nse.sgid] <- "Clonal No significant effect"
merged.gii.floh.clonal<-merged.gii.floh.clonal[complete.cases(merged.gii.floh.clonal),]
#Boxplot for wGII scores 
comb.gii.floh %>%   mutate(Growth.phenotype = fct_relevel(Growth.phenotype, "Clonal Control", "Non-clonal Control", "Clonal RAGE", "Non-clonal RAGE", "Clonal CALM", "Non-clonal CALM", "Clonal Context-dependent", "Non-clonal Context-dependent", "Clonal Growth attenuator", "Non-clonal Growth attenuator","Clonal No significant effect","Non-clonal No significant effect")) %>%   ggplot(aes(x=Growth.phenotype, y=wGII, colour=Growth.phenotype))+   geom_boxplot()+   ggtitle("wGII for LUAD tumours with Kras background and clonal vs non-clonal TSG mutation")+   xlab("Growth phenotype")+   theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')

merged.gii.floh.subclonal$Growth.phenotype<-NA
merged.gii.floh.subclonal$Growth.phenotype[merged.gii.floh.subclonal$Mutation %ni% c(nse.sgid,rage.sgid,calm.sgid,cd.sgid,tsg.sgid)] <- "Sub-clonal Control"
merged.gii.floh.subclonal$Growth.phenotype[merged.gii.floh.subclonal$Mutation %in% rage.sgid] <- "Sub-clonal RAGE"
merged.gii.floh.subclonal$Growth.phenotype[merged.gii.floh.subclonal$Mutation %in% calm.sgid] <- "Sub-clonal CALM"
merged.gii.floh.subclonal$Growth.phenotype[merged.gii.floh.subclonal$Mutation %in% cd.sgid] <- "Sub-clonal Context-dependent"
merged.gii.floh.subclonal$Growth.phenotype[merged.gii.floh.subclonal$Mutation %in% tsg.sgid] <- "Sub-clonal Growth attenuating"
merged.gii.floh.subclonal$Growth.phenotype[merged.gii.floh.subclonal$Mutation %in% nse.sgid] <- "Sub-clonal No significant effect"
merged.gii.floh.subclonal<-merged.gii.floh.subclonal[complete.cases(merged.gii.floh.subclonal),]
comb.gii.floh.correct<-rbind(merged.gii.floh.clonal,merged.gii.floh.subclonal)

#plot boxplot for wGII scores
comb.gii.floh.correct %>%
  mutate(Growth.phenotype = fct_relevel(Growth.phenotype, "Clonal Control", "Sub-clonal Control", "Clonal RAGE", "Sub-clonal RAGE", "Clonal CALM", "Sub-clonal CALM", "Clonal Context-dependent", "Sub-clonal Context-dependent", "Clonal Growth attenuating", "Sub-clonal Growth attenuating","Clonal No significant effect","Sub-clonal No significant effect")) %>%
  ggplot(aes(x=Growth.phenotype, y=wGII, colour=Growth.phenotype))+
  geom_boxplot()+
  ggtitle("wGII for LUAD tumours with Kras background and clonal vs non-clonal TSG mutation")+
  xlab("Growth phenotype")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')

#Repeat above commands to plot boxplot for wFLOH scores

#Construct heat map for wGII and wFLOH scores
wgii.means <- data.frame(aggregate(wGII~Growth.phenotype, comb.gii.floh, mean)) wfloh.means <- data.frame(aggregate(wFLOH~Growth.phenotype, comb.gii.floh, mean)) instability <- merge(wgii.means,wfloh.means, by="Growth.phenotype") colnames(instability) <- c("Phenotype", "wGII", "wFLOH") rownames(instability) <- instability$Phenotype instability <- instability[,2:3] a <- heatmap.2(as.matrix(instability), scale="column", density.info = 'none',              Colv=NA, Rowv=NA,              col=cm.colors(256), xlab="Genomic instability indicators", ylab="TSG mutation", main="Heat map of genomic instability indicators",              margins=c(5,14), cexRow = 1, cexCol = 1)


## HLA LOH scores from genomic instability data 
#Calling HLA LOH was performed on Microsoft Excel. HLA LOH was called in a sample region if HLA allele copy number rounded to 0.
#Read files of list of regions with 48 TSG mutations and wildtype mutations, and list of regions with HLALOH and without HLALOH 
hlaloh <- read.csv(file="Desktop/new/hlaloh.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
non.hlaloh <- read.csv(file="Desktop/new/non-hlaloh.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
hlaloh.48 <- read.csv(file="Desktop/new/hlaloh.48.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
hlaloh.wt <- read.csv(file="Desktop/new/hlaloh.wt.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

#Create dataframe consisting of regions with 48 TSG mutations and HLALOH, regions with 48 TSG mutations without HLALOH, regions wildtype mutations and HLALOH, regions with wildtype mutations and nonHLALOH
hlaloh.and.48 <- merge(hlaloh.48,hlaloh,by="Sample.ID") 
hlaloh.and.48$Growth.phenotype<-NA
hlaloh.and.48$Growth.phenotype[hlaloh.and.48$Mutation %in% rage.sgid] <- "RAGE"
hlaloh.and.48$Growth.phenotype[hlaloh.and.48$Mutation %in% calm.sgid] <- "CALM"
hlaloh.and.48$Growth.phenotype[hlaloh.and.48$Mutation %in% cd.sgid] <- "Context-dependent"
hlaloh.and.48$Growth.phenotype[hlaloh.and.48$Mutation %in% tsg.sgid] <- "Growth attenuator"
hlaloh.and.48$Growth.phenotype[hlaloh.and.48$Mutation %in% nse.sgid] <- "No significant effect"
hlaloh.and.48.count <- table(hlaloh.and.48$Growth.phenotype) 
nonhlaloh.and.48 <- merge(hlaloh.48,non.hlaloh,by="Sample.ID") 
nonhlaloh.and.48$Growth.phenotype<-NA
nonhlaloh.and.48$Growth.phenotype[nonhlaloh.and.48$Mutation %in% rage.sgid] <- "RAGE"
nonhlaloh.and.48$Growth.phenotype[nonhlaloh.and.48$Mutation %in% calm.sgid] <- "CALM"
nonhlaloh.and.48$Growth.phenotype[nonhlaloh.and.48$Mutation %in% cd.sgid] <- "Context-dependent"
nonhlaloh.and.48$Growth.phenotype[nonhlaloh.and.48$Mutation %in% tsg.sgid] <- "Growth attenuator"
nonhlaloh.and.48$Growth.phenotype[nonhlaloh.and.48$Mutation %in% nse.sgid] <- "No significant effect"
nonhlaloh.and.48.count <- table(nonhlaloh.and.48$Growth.phenotype) 
hlaloh.and.wt <- merge(hlaloh.wt,hlaloh,by="Sample.ID") #wt with LOHHLA #453
nonhlaloh.and.wt <- merge(hlaloh.wt,non.hlaloh,by="Sample.ID") #wt without LOHHLA #2431

#Perform Fisher’s exact test to test for independence between wildtype/HLALOH, wildtype/non-HLALOH, TSG/HLALOH, TSG/non-HLALOH
fet.wt <-
  matrix(c(453, 453, 2431, 2431),
         nrow = 2,
         dimnames = list(Phenotype = c("WT", "Wild Type"),
                         HLALOH = c("HLALOH", "Non-HLALOH")))            
fisher.test(fet.wt, alternative = "two.sided")


## Immunology scores from immunology data
#Load immunology dataset, list of patient with Kras+48TSGs, and list of wild-type patients 
imm <- read.csv(file="Desktop/attemptingCAMP/immune.subsets.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
my.mut.short.edited <- read.csv(file="Desktop/attemptingCAMP/my.mut.short.edited.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
my.mut.short.cont <- read.csv(file="Desktop/attemptingCAMP/mutTableAll_Region.cont.dr.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

#Create dataframe of immunology scores sorted into growth phenotypes and control groups
merged.imm <- merge(imm,my.mut.short.edited,by="Sample.ID")
merged.imm$Growth.phenotype<-NA
merged.imm$Growth.phenotype[merged.imm$Mutation %in% rage.sgid] <- "RAGE"
merged.imm$Growth.phenotype[merged.imm$Mutation %in% calm.sgid] <- "CALM"
merged.imm$Growth.phenotype[merged.imm$Mutation %in% cd.sgid] <- "Context-dependent"
merged.imm$Growth.phenotype[merged.imm$Mutation %in% tsg.sgid] <- "Growth attenuator"
merged.imm$Growth.phenotype[merged.imm$Mutation %in% nse.sgid] <- "No significant effect"
merged.imm<-merged.imm[complete.cases(merged.imm),] 
merged.imm.cont <- merge(imm,my.mut.short.cont,by="Sample.ID")
merged.imm.cont$Growth.phenotype <- "Control"
comb.imm<-rbind(merged.imm,merged.imm.cont)

#Plot boxplot for immunology scores
imm.plot <- ggplot(merged.imm, aes(x=Growth.phenotype, y=treg.score.danaher, colour=Growth.phenotype))+ 
  geom_boxplot()+
  ggtitle("Treg score of patient samples with particular driver mutations")+
  ylab("Treg score")+
  xlab("Growth phenotype")

#Repeat above command for CD8, CD4, CD45, B-cell, T-cell, Th1, and TIL scores

#Find mean of immunology scores and plot normalised heat map
cd8.means <- data.frame(aggregate(cd8.score.danaher~Growth.phenotype, comb.imm, mean))
cd4.means <- data.frame(aggregate(cd4.score.danaher~Growth.phenotype, comb.imm, mean))
cd45.means <- data.frame(aggregate(cd45.score.danaher~Growth.phenotype, comb.imm, mean))
bcell.means <- data.frame(aggregate(bcell.score.danaher~Growth.phenotype, comb.imm, mean))
tcell.means <- data.frame(aggregate(tcells.score.danaher~Growth.phenotype, comb.imm, mean))
th1.means <- data.frame(aggregate(th1.score.danaher~Growth.phenotype, comb.imm, mean))
til.means <- data.frame(aggregate(total.til.score.danahaer~Growth.phenotype, comb.imm, mean))
treg.means <- data.frame(aggregate(treg.score.danaher~Growth.phenotype, comb.imm, mean))

imm.hm <- cbind(cd8.means,cd4.means$cd4.score.danaher,cd8.means$cd8.score.danaher,cd45.means$cd45.score.danaher,bcell.means$bcell.score.danaher,tcell.means$tcells.score.danaher,th1.means$th1.score.danaher,til.means$total.til.score.danahaer,treg.means$treg.score.danaher)
colnames(imm.hm) <- c("Phenotype", "CD8 score", "CD4 score", "CD45 score", "B-cell score", "T-cell score", "Th1 score", "Total TIL score", "T-reg Score")
rownames(imm.hm) <- imm.hm$Phenotype
imm.hm <- imm.hm[,2:9]
a <- heatmap(as.matrix(imm.hm), scale="column", 
             Colv=NA, Rowv=NA,
             col=cm.colors(256), xlab="Immunology scores", ylab="Growth phenotype",
             margins=c(6,10), cexRow = 1, cexCol = 0.7)


## Wnt signalling scores from RNASeq data 
#Read list of RNASeq patient IDs and their growth phenotype
rnaseq.list <- read.csv(file="Desktop/new/rnaseq.list.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

#Construct curated dataframe of selected Wnt genes 
rnaseq.all <- read.csv(file="Desktop/attemptingCAMP/rna.seq.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
rnaseq.curated <- data.frame(cbind(rnaseq.all$SampleID, rnaseq.all$FZR1, rnaseq.all$LRP5, rnaseq.all$LRP6, rnaseq.all$AXIN1, rnaseq.all$AXIN2, rnaseq.all$APC, rnaseq.all$GSK3A, rnaseq.all$GSK3B, rnaseq.all$BCAT1, rnaseq.all$BCAT2, rnaseq.all$LEF1, rnaseq.all$BCL9))
colnames(rnaseq.curated) <- c("SampleID", "FZR1","LRP5","LRP6","AXIN1","AXIN2","APC","GSK3A","GSK3B","BCAT1","BCAT2","LEF1","BCL9")

#Create dataframe of RNASeq patients and their Wnt gene expressions 
rnaseq.comb <- merge(rnaseq.list,rnaseq.curated,by='SampleID')

#Find mean of Wnt RNA expression for each growth phenotype group and plot normalised heat map
FZR1.means <- data.frame(aggregate(FZR1~Phenotype, rnaseq.comb.1, mean))
LRP5.means <- data.frame(aggregate(LRP5~Phenotype, rnaseq.comb.1, mean))
LRP6.means <- data.frame(aggregate(LRP6~Phenotype, rnaseq.comb.1, mean))
AXIN1.means <- data.frame(aggregate(AXIN1~Phenotype, rnaseq.comb.1, mean))
AXIN2.means <- data.frame(aggregate(AXIN2~Phenotype, rnaseq.comb.1, mean))
APC.means <- data.frame(aggregate(APC~Phenotype, rnaseq.comb.1, mean))
GSK3A.means <- data.frame(aggregate(GSK3A~Phenotype, rnaseq.comb.1, mean))
GSK3B.means <- data.frame(aggregate(GSK3B~Phenotype, rnaseq.comb.1, mean))
BCAT1.means <- data.frame(aggregate(BCAT1~Phenotype, rnaseq.comb.1, mean))
BCAT2.means <- data.frame(aggregate(BCAT2~Phenotype, rnaseq.comb.1, mean))
LEF1.means <- data.frame(aggregate(LEF1~Phenotype, rnaseq.comb.1, mean))
BCL9.means <- data.frame(aggregate(BCL9~Phenotype, rnaseq.comb.1, mean))
wnt.hm <- cbind(FZR1.means,LRP5.means$LRP5,LRP6.means$LRP6,AXIN1.means$AXIN1,AXIN2.means$AXIN2,APC.means$APC,GSK3A.means$GSK3A,GSK3B.means$GSK3B,BCAT1.means$BCAT1,BCAT2.means$BCAT2,LEF1.means$LEF1,BCL9.means$BCL9)
colnames(wnt.hm) <- c("Phenotype", "FZR1","LRP5","LRP6","AXIN1","AXIN2","APC","GSK3A","GSK3B","BCAT1","BCAT2","LEF1","BCL9")
rownames(wnt.hm) <- wnt.hm$Phenotype
wnt.hm <- wnt.hm[,2:13]
a <- heatmap(as.matrix(wnt.hm), scale="column", 
             Colv=NA, Rowv=NA,
             col=cm.colors(256), xlab="Wnt pathway molecular signatures", ylab="Growth phenotype", 
             margins=c(6,12), cexRow = 1, cexCol = 1)