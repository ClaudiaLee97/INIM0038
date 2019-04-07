#organising data
mutTableAll <- read.table("Desktop/attemptingCAMP/mutTableAll.CSV", header = TRUE, sep=',')
#driver mut only
mutTableAll.dr <- data.frame(mutTableAll[mutTableAll$DriverMut %in% 'TRUE',]) 
write.csv(mutTableAll.dr, "Desktop/attemptingCAMP/mutTableAll.dr.CSV")
#kras background only
mutTableAll.dr.kras <- data.frame(mutTableAll.dr[mutTableAll.dr$Hugo_Symbol %in% 'KRAS',])
sample.dr.kras <- mutTableAll.dr.kras$SampleID
#RAGE patients
rage.sgid <- c('ATM', 'NF1', 'PTEN', 'SETD2', 'STK11', 'RNF43')
pat.rage <- data.frame(mutTableAll.dr[mutTableAll.dr$Hugo_Symbol %in% rage.sgid,])
pat.rage.1 <- data.frame(pat.rage[pat.rage$SampleID %in% sample.dr.kras,])
write.csv(pat.rage.1, "Desktop/attemptingCAMP/patient.rage.CSV")
#CALM patients
calm.sgid <- c("DNMT3A", "CMTR2", "RBM10", "TSC1", "CDKN2C", "RB1", "STAG2")
pat.calm <- data.frame(mutTableAll.dr[mutTableAll.dr$Hugo_Symbol %in% calm.sgid,])
pat.calm.1 <- data.frame(pat.calm[pat.calm$SampleID %in% sample.dr.kras,])
write.csv(pat.calm.1, "Desktop/attemptingCAMP/patient.calm.CSV")
#context-dependent patients
cd.sgid <- c("DOT1L")
pat.cd <- data.frame(mutTableAll.dr[mutTableAll.dr$Hugo_Symbol %in% cd.sgid,])
pat.cd.1 <- data.frame(pat.cd[pat.cd$SampleID %in% sample.dr.kras,])
write.csv(pat.cd.1, "Desktop/attemptingCAMP/patient.cd.CSV")
#growth attenuator/tumour supressor patients 
tsg.sgid <- c("ATRX", "BAP1",'BRCA2','EP300',"ERCC4","FAT1","KDM6A","KEAP1","KMT2C","LRP1B","NF2","SMAD4","SMARCA4","UBR5","WRN")
pat.tsg <- data.frame(mutTableAll.dr[mutTableAll.dr$Hugo_Symbol %in% tsg.sgid,])
pat.tsg.1 <- data.frame(pat.tsg[pat.tsg$SampleID %in% sample.dr.kras,])
write.csv(pat.tsg.1, "Desktop/attemptingCAMP/patient.tsg.CSV")
#no significant effect patients
nse.sgid <- c('ARHGAP35','ARID1A','ARID1B','ARID2','ATF7IP','CDKN2A','CIC','FBXW7','GATA3','KMT2D','LATS1','MGA','NCOA6','NCOR1','POLE','PTPRD','RASA1','TET2', 'TP53')
pat.nse <- data.frame(mutTableAll.dr[mutTableAll.dr$Hugo_Symbol %in% nse.sgid,])
pat.nse.1 <- data.frame(pat.nse[pat.nse$SampleID %in% sample.dr.kras,])
write.csv(pat.nse.1, "Desktop/attemptingCAMP/patient.nse.CSV")
#missing ptprd

#luad only
luad.clin.list <- read.table("Desktop/attemptingCAMP/luad.clin.list.CSV", header = TRUE, sep=',')
pat.luad <- luad.clin.list$SampleID

rage.luad.kras.dr <- data.frame(pat.rage.1[pat.rage.1$SampleID %in% pat.luad,])
write.csv(rage.luad.kras.dr, "Desktop/attemptingCAMP/clean.rage.CSV")
calm.luad.kras.dr <- data.frame(pat.calm.1[pat.calm.1$SampleID %in% pat.luad,])
write.csv(calm.luad.kras.dr, "Desktop/attemptingCAMP/clean.calm.CSV")
cd.luad.kras.dr <- data.frame(pat.cd.1[pat.cd.1$SampleID %in% pat.luad,])
write.csv(cd.luad.kras.dr, "Desktop/attemptingCAMP/clean.cd.SV")
tsg.luad.kras.dr <- data.frame(pat.tsg.1[pat.tsg.1$SampleID %in% pat.luad,])
write.csv(tsg.luad.kras.dr, "Desktop/attemptingCAMP/clean.tsg.CSV")
nse.luad.kras.dr <- data.frame(pat.nse.1[pat.nse.1$SampleID %in% pat.luad,])
write.csv(nse.luad.kras.dr, "Desktop/attemptingCAMP/clean.nse.CSV")
clean.all <- rbind(rage.luad.kras.dr,calm.luad.kras.dr,cd.luad.kras.dr,tsg.luad.kras.dr,nse.luad.kras.dr)

write.csv(clean.all, "Desktop/attemptingCAMP/clean.all.CSV")
#genomic data of samples which contain driver mutations within the 48 genes AND Kras background

#survival correlation
surv <- as.data.frame(tracerx.surv500)
write.csv(surv, "Desktop/attemptingCAMP/surv.CSV")
#patients who have adenocarcinoma, who contain driver mutations within the 48 genes and have Kras background
surv <- read.csv(file="Desktop/new/csvs/surv.kev.clean.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
SampleID <- data.frame(surv$SampleID)
os.OS <- data.frame(surv$OS_time_days)
os <- data.frame(cbind(SampleID,os.OS))
colnames(os) <- c("SampleID","OS.time.days")
merged.os <- merge(merged.num,os,by="SampleID")
merged.os<-merged.os[complete.cases(merged.os),] #removed NAs
write.csv(merged.os, "Desktop/attemptingCAMP/merged.os.new.CSV") #CORRECT

merged.os.1 <- read.csv(file="Desktop/attemptingCAMP/merged.os.new.CSV", header=TRUE, sep=",", stringsAsFactors = FALSE)
rage.os <- data.frame(merged.os.1$RAGE, merged.os.1$OS.time.days)
names(rage.os) <- c("num.mut", "os.days")
name <- rep("RAGE",95)
rage.os <- cbind(name,rage.os)
calm.os <- data.frame(merged.os.1$CALM, merged.os.1$OS.time.days)
names(calm.os) <- c("num.mut", "os.days")
name <- rep("CALM",95)
calm.os <- cbind(name,calm.os)
cd.os <- data.frame(merged.os.1$CD, merged.os.1$OS.time.days)
names(cd.os) <- c("num.mut", "os.days")
name <- rep("Context-dependent",95)
cd.os <- cbind(name,cd.os)
tsg.os <- data.frame(merged.os.1$TSG, merged.os.1$OS.time.days)
names(tsg.os) <- c("num.mut", "os.days")
name <- rep("Growth attenuating",95)
tsg.os <- cbind(name,tsg.os)
nse.os <- data.frame(merged.os.1$NSE, merged.os.1$OS.time.days)
names(nse.os) <- c("num.mut", "os.days")
name <- rep("No significant effect",95)
nse.os <- cbind(name,nse.os)
os.df <- data.frame(rbind(nse.os,rage.os,calm.os,cd.os,tsg.os))
write.csv(os.df, "Desktop/new/csvs/os.df.new.csv")
os.df[os.df == 0] <- NA #change 0 to NA and remove rows
os.df.na<-os.df[complete.cases(os.df),]
colnames(os.df.na) <- c("Growth.phenotype", "Num.mutations", "OS.days")
write.csv(os.df.na, "Desktop/new/csvs/os.df.na.new.csv")
#manually added in control group!!!

#plotting survival curves
#overall survival
os.df.na.new <- read.csv(file="Desktop/new/csvs/os.df.na.new.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
table(os.df.na.new$Growth.phenotype)
surv.obj <- Surv(time=os.df.na.new$OS.days)
fit.1 <- survfit(surv.obj~Growth.phenotype, data=os.df.na.new)
summary(fit.1)
a<-ggsurvplot(fit.1, data=os.df.na.new, risk.table = TRUE, pval=TRUE)
a
ggsave("Desktop/new/comb.surv.pdf", dpi=600)
#individual comparison
os.df.0 <- os.df
colnames(os.df.0) <- c("Growth.phenotype", "Num.mutations", "OS.days")
#RAGE vs control
os.df.rage<-os.df.na.new[os.df.na.new$Growth.phenotype %in% c("RAGE","Control"),]
os.df.rage$OS.days<-as.numeric(as.character(os.df.rage$OS.days))
surv.obj <- Surv(time=os.df.rage$OS.days)
fit.1 <- survfit(surv.obj~Growth.phenotype, data=os.df.rage)
summary(fit.1)
a<-ggsurvplot(fit.1, data=os.df.rage, risk.table = TRUE, pval=TRUE)
a
ggsave("Desktop/attemptingCAMP/surv.rage.pdf", dpi=600)
#CALM vs control
os.df.calm<-os.df.na.new[os.df.na.new$Growth.phenotype %in% c("CALM","Control"),]
os.df.calm$OS.days<-as.numeric(as.character(os.df.calm$OS.days))
surv.obj <- Surv(time=os.df.calm$OS.days)
fit.1 <- survfit(surv.obj~Growth.phenotype, data=os.df.calm)
summary(fit.1)
a<-ggsurvplot(fit.1, data=os.df.calm, risk.table = TRUE, pval=TRUE)
a
ggsave("Desktop/attemptingCAMP/surv.calm.pdf", dpi=600)
#NSE vs control
os.df.nse<-os.df.na.new[os.df.na.new$Growth.phenotype %in% c("No significant effect","Control"),]
os.df.nse$OS.days<-as.numeric(as.character(os.df.nse$OS.days))
surv.obj <- Surv(time=os.df.nse$OS.days)
fit.1 <- survfit(surv.obj~Growth.phenotype, data=os.df.nse)
summary(fit.1)
a<-ggsurvplot(fit.1, data=os.df.nse, risk.table = TRUE, pval=TRUE)
a
ggsave("Desktop/attemptingCAMP/surv.prob.nse.new.pdf", dpi=600)
#CD vs control
os.df.cd<-os.df.na.new[os.df.na.new$Growth.phenotype %in% c("Context-dependent","Control"),]
os.df.cd$OS.days<-as.numeric(as.character(os.df.cd$OS.days))
surv.obj <- Surv(time=os.df.cd$OS.days)
fit.1 <- survfit(surv.obj~Growth.phenotype, data=os.df.cd)
summary(fit.1)
a<-ggsurvplot(fit.1, data=os.df.cd, risk.table = TRUE, pval=TRUE)
a
ggsave("Desktop/attemptingCAMP/surv.cd.pdf", dpi=600)
#TSG vs control
os.df.tsg<-os.df.na.new[os.df.na.new$Growth.phenotype %in% c("Growth attenuating","Control"),]
os.df.tsg$OS.days<-as.numeric(as.character(os.df.tsg$OS.days))
surv.obj <- Surv(time=os.df.tsg$OS.days)
fit.1 <- survfit(surv.obj~Growth.phenotype, data=os.df.tsg)
summary(fit.1)
a<-ggsurvplot(fit.1, data=os.df.tsg, risk.table = TRUE, pval=TRUE)
a
ggsave("Desktop/attemptingCAMP/surv.tsg.pdf", dpi=600)

#plotting coxph 
os.df.na.new$Growth.phenotype <- factor(os.df.na.new$Growth.phenotype, levels = c("Control","RAGE","CALM","Context-dependent", "Growth attenuating", "No significant effect"))
os.df.na.new$Growth.phenotype<-as.factor(os.df.na.new$Growth.phenotype)
fit.coxph <-coxph(formula=Surv(time=os.df.na.new$OS.days)~ Growth.phenotype, os.df.na.new)
ggforest(fit.coxph, data=os.df.na.new)
ggsave("Desktop/new/os.coxph.pdf", dpi=600)

#lesion pathology
surv <- read.csv(file="Desktop/new/csvs/surv.kev.clean.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
SampleID <- data.frame(surv$SampleID)
lsp.LSP <- data.frame(surv$Lesion1_size_pathology)
lsp <- data.frame(cbind(SampleID,lsp.LSP))
colnames(lsp) <- c("SampleID","LesionSize")
merged.lsp <- merge(merged.num,lsp,by="SampleID")
merged.lsp<-merged.lsp[complete.cases(merged.lsp),] #removed NAs
write.csv(merged.lsp, "Desktop/new/csvs/merged.lsp.csv") #CORRECT

merged.lsp.1 <- read.csv(file="Desktop/new/csvs/merged.lsp.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
rage.lsp <- data.frame(merged.lsp.1$RAGE, merged.lsp.1$LesionSize)
names(rage.lsp) <- c("num.mut", "LesionSize")
name <- rep("RAGE",99)
rage.lsp <- cbind(name,rage.lsp)
calm.lsp <- data.frame(merged.lsp.1$CALM, merged.lsp.1$LesionSize)
names(calm.lsp) <- c("num.mut", "LesionSize")
name <- rep("CALM",99)
calm.lsp <- cbind(name,calm.lsp)
cd.lsp <- data.frame(merged.lsp.1$CD, merged.lsp.1$LesionSize)
names(cd.lsp) <- c("num.mut", "LesionSize")
name <- rep("Context-dependent",99)
cd.lsp <- cbind(name,cd.lsp)
tsg.lsp <- data.frame(merged.lsp.1$TSG, merged.lsp.1$LesionSize)
names(tsg.lsp) <- c("num.mut", "LesionSize")
name <- rep("Growth attenuator",99)
tsg.lsp <- cbind(name,tsg.lsp)
nse.lsp <- data.frame(merged.lsp.1$NSE, merged.lsp.1$LesionSize)
names(nse.lsp) <- c("num.mut", "LesionSize")
name <- rep("No significant effect",99)
nse.lsp <- cbind(name,nse.lsp)

lsp.df <- data.frame(rbind(nse.lsp,rage.lsp,calm.lsp,cd.lsp,tsg.lsp))
lsp.df[lsp.df == 0] <- NA #change 0 to NA and remove rows
lsp.df.na<-lsp.df[complete.cases(lsp.df),]
colnames(lsp.df.na) <- c("Growth.phenotype", "Num.mutations", "LesionSize")
#calculate tumour volume
lsp.df.na$TumourVol <- (lsp.df.na$LesionSize)^3/2
write.csv(lsp.df.na, "Desktop/new/csvs/lsp.csv")
#manually added in control group!
lsp.vol <- read.csv(file="Desktop/new/csvs/lsp.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
a <- ggplot()+
  geom_violin(data=lsp.vol, aes(x=Growth.phenotype, y=TumourVol, color=Growth.phenotype))+
  scale_y_log10()+
  labs(title="Correlating tumour growth phenotype to tumour volume", y="Tumour volume", x="Tumour growth phenotype")+
  theme_classic()
a
ggsave("Desktop/new/lsp.violin.pdf", dpi=600, width=10, height=5)

#kptc.ktc.corr
kptc.ktc.corr <- read.csv(file="Desktop/attemptingCAMP/kptc.ktc.corr.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
head(kptc.ktc.corr)
cor.test(kptc.ktc.corr$ktc.mean,kptc.ktc.corr$kptc.mean)
b <- ggplot(data=kptc.ktc.corr, aes(x=ktc.mean, y=kptc.mean, colour=sgID))+
  geom_point(data=kptc.ktc.corr, aes(x=ktc.mean, y=kptc.mean, colour=sgID))+
  geom_smooth(data=kptc.ktc.corr,method="lm",colour="black")+
  labs(title="Comparing effect of Tp53 on 95th percentile growth rates\nR=0.8254419,p=5.29e-13", y="KPTC tumour growth rate", x="KTC tumour growth rate")+
  theme_classic()
b
ggsave("Desktop/new/kptc.ktc.corr.pdf", height=5, width=10)

#wGII and wFLOH
mutTableAll_Region = clean.all %>%
  separate_rows(RegionSum,sep=";") %>%
  separate(RegionSum,into=c("region_id","region_var_ref_count"),sep=":") %>%
  separate(region_var_ref_count,into=c("var_count","ref_count"),sep="/") %>%
  mutate(SampleID=gsub("\\.","-",paste(site,SampleID,region_id,sep="_")))
mutTableAll_Region.driver<-mutTableAll_Region[mutTableAll_Region$DriverMut %in% TRUE,] #only driver mutations allowed
write.csv(mutTableAll_Region.driver, "Desktop/attemptingCAMP/mutTableAll_Region.clean.csv")

my.mut <- read.csv(file="Desktop/attemptingCAMP/mutTableAll_Region.clean.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
my.mut.short <- data.frame(cbind(my.mut$SampleID, my.mut$Hugo_Symbol))
colnames(my.mut.short) <- c("Sample.ID", "Mutation")
gii.floh <- read.csv(file="Desktop/attemptingCAMP/gii.CSV", header=TRUE, sep=",", stringsAsFactors = FALSE)
merged.gii.floh <- merge(gii.floh,my.mut.short,by="Sample.ID")

merged.gii.floh$Growth.phenotype<-NA
merged.gii.floh$Growth.phenotype[merged.gii.floh$Mutation %in% rage.sgid] <- "Clonal RAGE"
merged.gii.floh$Growth.phenotype[merged.gii.floh$Mutation %in% calm.sgid] <- "Clonal CALM"
merged.gii.floh$Growth.phenotype[merged.gii.floh$Mutation %in% cd.sgid] <- "Clonal Context-dependent"
merged.gii.floh$Growth.phenotype[merged.gii.floh$Mutation %in% tsg.sgid] <- "Clonal Growth attenuator"
merged.gii.floh$Growth.phenotype[merged.gii.floh$Mutation %in% nse.sgid] <- "Clonal No significant effect"
merged.gii.floh<-merged.gii.floh[complete.cases(merged.gii.floh),] #removed NAs
write.csv(merged.gii.floh,"Desktop/attemptingCAMP/merged.gii.floh.new.csv")

my.mut.nondr <- read.csv(file="Desktop/attemptingCAMP/mutTableAll_Region.nondriver.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
my.mut.short.non.dr <- data.frame(cbind(my.mut.nondr$SampleID, my.mut.nondr$Hugo_Symbol))
colnames(my.mut.short.non.dr) <- c("Sample.ID", "Mutation")
gii.floh <- read.csv(file="Desktop/attemptingCAMP/gii.CSV", header=TRUE, sep=",", stringsAsFactors = FALSE)
merged.gii.floh.nondr <- merge(gii.floh,my.mut.short.non.dr,by="Sample.ID")

merged.gii.floh.nondr$Growth.phenotype<-NA
merged.gii.floh.nondr$Growth.phenotype[merged.gii.floh.nondr$Mutation %in% rage.sgid] <- "Non-clonal RAGE"
merged.gii.floh.nondr$Growth.phenotype[merged.gii.floh.nondr$Mutation %in% calm.sgid] <- "Non-clonal CALM"
merged.gii.floh.nondr$Growth.phenotype[merged.gii.floh.nondr$Mutation %in% cd.sgid] <- "Non-clonal Context-dependent"
merged.gii.floh.nondr$Growth.phenotype[merged.gii.floh.nondr$Mutation %in% tsg.sgid] <- "Non-clonal Growth attenuator"
merged.gii.floh.nondr$Growth.phenotype[merged.gii.floh.nondr$Mutation %in% nse.sgid] <- "Non-clonal No significant effect"
merged.gii.floh.nondr<-merged.gii.floh.nondr[complete.cases(merged.gii.floh),] #removed NAs

comb.gii.floh<-rbind(merged.gii.floh.nondr,merged.gii.floh)

comb.gii.floh %>%
  mutate(Growth.phenotype = fct_relevel(Growth.phenotype, "Clonal RAGE", "Non-clonal RAGE", "Clonal CALM", "Non-clonal CALM", "Clonal Context-dependent", "Non-clonal Context-dependent", "Clonal Growth attenuator", "Non-clonal Growth attenuator","Clonal No significant effect","Non-clonal No significant effect")) %>%
  ggplot(aes(x=Growth.phenotype, y=wGII, colour=Growth.phenotype))+
  geom_boxplot()+
  ggtitle("wGII for LUAD tumours with Kras background and clonal vs non-clonal TSG mutation")+
  xlab("Growth phenotype")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')
ggsave("Desktop/new/wgii.pdf", dpi=600, width = 8, height=5)

comb.gii.floh %>%
  mutate(Growth.phenotype = fct_relevel(Growth.phenotype, "Clonal RAGE", "Non-clonal RAGE", "Clonal CALM", "Non-clonal CALM", "Clonal Context-dependent", "Non-clonal Context-dependent", "Clonal Growth attenuator", "Non-clonal Growth attenuator","Clonal No significant effect","Non-clonal No significant effect")) %>%
  ggplot(aes(x=Growth.phenotype, y=wFLOH, colour=Growth.phenotype))+
  geom_boxplot()+
  ggtitle("wFLOH for LUAD tumours with Kras background and clonal vs non-clonal TSG mutation")+
  xlab("Growth phenotype")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')
ggsave("Desktop/new/wfloh.pdf", dpi=600, width = 8, height=5)

#mki67, necrosis, mit-hpf
mki67 <- read.csv(file="Desktop/attemptingCAMP/mki67.new.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
a <- ggplot(data=mki67, aes(x=Phenotype, y=MKI67, colour=Phenotype))+
  geom_boxplot()+
  labs(title="Mki67 RNA expression of LUAD tumours with Kras background and TSG mutations", y="Mki67 RNA expression", x="Growth phenotype")+
  stat_summary(fun.y=mean)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')
ggsave("Desktop/new/mki67.pdf", height=5, width=10)

nec.mit<- read.csv(file="Desktop/attemptingCAMP/nec_mit_data_manual.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
a <- ggplot(data=nec.mit, aes(x=Phenotype, y=Mitosis, colour=Phenotype))+
  geom_boxplot()+
  labs(title="Mitosis-HPF of Kras-background LUAD  tumours with Kras background and TSG mutations", y="Mitoses-HPF", x="Growth phenotype")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')
ggsave("Desktop/new/mit-hpf.pdf", height=5, width=10)

a <- ggplot(data=nec.mit, aes(x=Phenotype, y=Necrosis, colour=Phenotype))+
  geom_boxplot()+
  labs(title="Percentage necrosis of Kras-background LUAD  tumours with Kras background and TSG mutations", y="Percentage necrosis", x="Growth phenotype")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')
ggsave("Desktop/new/nec.pdf", height=5, width=10)

#try heat map for genomic growth indicators
rm(mean)
mki67.means <- data.frame(aggregate(MKI67~Phenotype, mki67, mean))
mit.means <- data.frame(aggregate(Mitosis~Phenotype, nec.mit, mean))
nec.means <- data.frame(aggregate(Necrosis~Phenotype, nec.mit, mean))
gen.growth <- cbind(mki67.means,mit.means$Mitosis,nec.means$Necrosis)
colnames(gen.growth) <- c("Phenotype", "Mki67", "Mitosis", "Necrosis")
rownames(gen.growth) <- gen.growth$Phenotype
gen.growth <- gen.growth[,2:4]
a <- heatmap(as.matrix(gen.growth), scale="column", 
             Colv=NA, Rowv=NA,
             col=cm.colors(256), xlab="Genomic growth indiactors", ylab="TSG mutation", main="Heat map of genomic growth indicators",
             margins=c(5,10), cexRow = 1, cexCol = 1)

#immunology
comb.imm <- read.csv(file="Desktop/attemptingCAMP/comb.imm.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
#plot cd8,cd4,b-cell,cd45,t-cell,th1,til,treg scores
imm.plot <- ggplot(comb.imm, aes(x=Growth.phenotype, y=total.til.score.danahaer, colour=Growth.phenotype))+ 
  geom_boxplot()+
  ggtitle("Total TIL score of LUAD tumours with Kras background and tumour suppressor mutations")+
  ylab("Total TIL score")+xlab("Growth phenotype")+theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')
imm.plot
ggsave("Desktop/new/imm.til.pdf", dpi=600,width=10,height=5)

is.numeric(comb.imm)
cd8.means <- data.frame(aggregate(cd8.score.danaher~Growth.phenotype, comb.imm, mean))
cd4.means <- data.frame(aggregate(cd4.score.danaher~Growth.phenotype, comb.imm, mean))
cd45.means <- data.frame(aggregate(cd45.score.danaher~Growth.phenotype, comb.imm, mean))
bcell.means <- data.frame(aggregate(bcell.score.danaher~Growth.phenotype, comb.imm, mean))
tcell.means <- data.frame(aggregate(tcells.score.danaher~Growth.phenotype, comb.imm, mean))
th1.means <- data.frame(aggregate(th1.score.danaher~Growth.phenotype, comb.imm, mean))
til.means <- data.frame(aggregate(total.til.score.danahaer~Growth.phenotype, comb.imm, mean))
treg.means <- data.frame(aggregate(treg.score.danaher~Growth.phenotype, comb.imm, mean))

imm.hm <- cbind(cd8.means,cd4.means$cd4.score.danaher,cd8.means$cd8.score.danaher,cd45.means$cd45.score.danaher,bcell.means$bcell.score.danaher,tcell.means$tcells.score.danaher,th1.means$th1.score.danaher,til.means$total.til.score.danahaer,treg.means$treg.score.danaher)
colnames(imm.hm) <- c("Phenotype", "CD8 score", "CD4 score", "CD45 score", "B-cell score", "T-cell score", "Th1 score", "Total TIL score", "Treg Score")
rownames(imm.hm) <- imm.hm$Phenotype
imm.hm <- imm.hm[,2:8]
a <- heatmap(as.matrix(imm.hm), scale="column", 
             Colv=NA, Rowv=NA,
             col=cm.colors(256), xlab="Immunology scores", ylab="TSG mutation", main="Heat map of immunology scores",
             margins=c(5,10), cexRow = 1, cexCol = 0.7)

#wnt signalling
#motifs: FZR1,LRP5,LRP6,AXIN1,AXIN2,APC,GSK3A,GSK3B,BCAT1,BCAT2,LEF1,BCL9
rnaseq.list <- read.csv(file="Desktop/new/rnaseq.list.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
rnaseq.all <- read.csv(file="Desktop/attemptingCAMP/rna.seq.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
rnaseq.curated <- data.frame(cbind(rnaseq.all$SampleID, rnaseq.all$FZR1, rnaseq.all$LRP5, rnaseq.all$LRP6, rnaseq.all$AXIN1, rnaseq.all$AXIN2, rnaseq.all$APC, rnaseq.all$GSK3A, rnaseq.all$GSK3B, rnaseq.all$BCAT1, rnaseq.all$BCAT2, rnaseq.all$LEF1, rnaseq.all$BCL9))
colnames(rnaseq.curated) <- c("SampleID", "FZR1","LRP5","LRP6","AXIN1","AXIN2","APC","GSK3A","GSK3B","BCAT1","BCAT2","LEF1","BCL9")
rnaseq.comb <- merge(rnaseq.list,rnaseq.curated,by='SampleID')
write.csv(rnaseq.comb,"Desktop/new/csvs/rnaseq.comb.csv")

rnaseq.comb.1 <- read.csv(file="Desktop/new/csvs/rnaseq.comb.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
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
             col=cm.colors(256), xlab="Wnt pathway molecular signatures", ylab="Growth phenotype", main="Heat map of Wnt pathway molecular signatures (RNA Seq)",
             margins=c(6,12), cexRow = 1, cexCol = 1)

#hlaloh
hlaloh <- read.csv(file="Desktop/new/hlaloh.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
non.hlaloh <- read.csv(file="Desktop/new/non-hlaloh.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
hlaloh.48 <- read.csv(file="Desktop/new/hlaloh.48.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
hlaloh.wt <- read.csv(file="Desktop/new/hlaloh.wt.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

hlaloh.and.48 <- merge(hlaloh.48,hlaloh,by="Sample.ID") #ts with LOHHLA
hlaloh.and.48$Growth.phenotype<-NA
hlaloh.and.48$Growth.phenotype[hlaloh.and.48$Mutation %in% rage.sgid] <- "RAGE"
hlaloh.and.48$Growth.phenotype[hlaloh.and.48$Mutation %in% calm.sgid] <- "CALM"
hlaloh.and.48$Growth.phenotype[hlaloh.and.48$Mutation %in% cd.sgid] <- "Context-dependent"
hlaloh.and.48$Growth.phenotype[hlaloh.and.48$Mutation %in% tsg.sgid] <- "Growth attenuator"
hlaloh.and.48$Growth.phenotype[hlaloh.and.48$Mutation %in% nse.sgid] <- "No significant effect"
hlaloh.and.48.count <- table(hlaloh.and.48$Growth.phenotype) #ts with LOHHLA

nonhlaloh.and.48 <- merge(hlaloh.48,non.hlaloh,by="Sample.ID") #ts without LOHHLA
nonhlaloh.and.48$Growth.phenotype<-NA
nonhlaloh.and.48$Growth.phenotype[nonhlaloh.and.48$Mutation %in% rage.sgid] <- "RAGE"
nonhlaloh.and.48$Growth.phenotype[nonhlaloh.and.48$Mutation %in% calm.sgid] <- "CALM"
nonhlaloh.and.48$Growth.phenotype[nonhlaloh.and.48$Mutation %in% cd.sgid] <- "Context-dependent"
nonhlaloh.and.48$Growth.phenotype[nonhlaloh.and.48$Mutation %in% tsg.sgid] <- "Growth attenuator"
nonhlaloh.and.48$Growth.phenotype[nonhlaloh.and.48$Mutation %in% nse.sgid] <- "No significant effect"
nonhlaloh.and.48.count <- table(nonhlaloh.and.48$Growth.phenotype) #ts without LOHHLA

hlaloh.and.wt <- merge(hlaloh.wt,hlaloh,by="Sample.ID") #wt with LOHHLA #453
nonhlaloh.and.wt <- merge(hlaloh.wt,non.hlaloh,by="Sample.ID") #wt without LOHHLA #2431

hlaloh.comb <- read.csv(file="Desktop/new/hlaloh.comb.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

fet.wt <-
  matrix(c(453, 453, 2431, 2431),
         nrow = 2,
         dimnames = list(Phenotype = c("WT", "Wild Type"),
                         HLALOH = c("HLALOH", "Non-HLALOH")))            
fisher.test(fet.wt, alternative = "two.sided")


#pycloneCCF
mutTableAll_Region.clean <- read.csv(file="Desktop/attemptingCAMP/mutTableAll_Region.clean.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
pyclone.dr = mutTableAll_Region.clean %>%
  separate_rows(PyCloneCCF,sep=";") %>%
  separate(PyCloneCCF,into=c("region_id","pycloneccf"),sep=":") 
write.csv(pyclone.dr, "Desktop/new/pycloneccf.dr.csv")
pyclone.dr <- read.csv(file="Desktop/new/pycloneccf.dr.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

pyclone.dr.cur <- data.frame(cbind(pyclone.dr$SampleID,pyclone.dr$Hugo_Symbol,pyclone.dr$pycloneccf))
colnames(pyclone.dr.cur) <- c("Sample ID", "Hugo_Symbol", "pycloneccf")
write.csv(pyclone.dr.cur, "Desktop/new/pycloneccf.dr.cur.csv")
pyclone.dr.cur <- read.csv(file="Desktop/new/pycloneccf.dr.cur.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

pyclone.dr.cur$Growth.phenotype<-NA
pyclone.dr.cur$Growth.phenotype[pyclone.dr.cur$Hugo_Symbol %in% rage.sgid] <- "Clonal RAGE"
pyclone.dr.cur$Growth.phenotype[pyclone.dr.cur$Hugo_Symbol %in% calm.sgid] <- "Clonal CALM"
pyclone.dr.cur$Growth.phenotype[pyclone.dr.cur$Hugo_Symbol %in% cd.sgid] <- "Clonal Context-dependent"
pyclone.dr.cur$Growth.phenotype[pyclone.dr.cur$Hugo_Symbol %in% tsg.sgid] <- "Clonal Growth attenuator"
pyclone.dr.cur$Growth.phenotype[pyclone.dr.cur$Hugo_Symbol %in% nse.sgid] <- "Clonal No significant effect"
pyclone.dr.means <- data.frame(aggregate(pycloneccf~Growth.phenotype, pyclone.dr.cur, mean))
pyclone.dr.var <- data.frame(aggregate(pycloneccf~Growth.phenotype, pyclone.dr.cur, var))

mutTableAll_Region.nondriver <- read.csv(file="Desktop/attemptingCAMP/mutTableAll_Region.nondriver.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
pyclone.nondr = mutTableAll_Region.nondriver %>%
  separate_rows(PyCloneCCF,sep=";") %>%
  separate(PyCloneCCF,into=c("region_id","pycloneccf"),sep=":") 
write.csv(pyclone.nondr, "Desktop/new/pycloneccf.nondr.csv")

pyclone.nondr.cur <- data.frame(cbind(pyclone.nondr$SampleID,pyclone.nondr$Hugo_Symbol,pyclone.nondr$pycloneccf))
colnames(pyclone.nondr.cur) <- c("Sample ID", "Hugo_Symbol", "pycloneccf")
write.csv(pyclone.nondr.cur, "Desktop/new/pycloneccf.nondr.cur.csv")
pyclone.nondr.cur <- read.csv(file="Desktop/new/pycloneccf.nondr.cur.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

pyclone.nondr.cur$Growth.phenotype<-NA
pyclone.nondr.cur$Growth.phenotype[pyclone.nondr.cur$Hugo_Symbol %in% rage.sgid] <- "Non-clonal RAGE"
pyclone.nondr.cur$Growth.phenotype[pyclone.nondr.cur$Hugo_Symbol %in% calm.sgid] <- "Non-clonal CALM"
pyclone.nondr.cur$Growth.phenotype[pyclone.nondr.cur$Hugo_Symbol %in% cd.sgid] <- "Non-clonal Context-dependent"
pyclone.nondr.cur$Growth.phenotype[pyclone.nondr.cur$Hugo_Symbol %in% tsg.sgid] <- "Non-clonal Growth attenuator"
pyclone.nondr.cur$Growth.phenotype[pyclone.nondr.cur$Hugo_Symbol %in% nse.sgid] <- "Non-clonal No significant effect"
pyclone.nondr.means <- data.frame(aggregate(pycloneccf~Growth.phenotype, pyclone.nondr.cur, mean))
pyclone.nondr.var <- data.frame(aggregate(pycloneccf~Growth.phenotype, pyclone.nondr.cur, var))

#ugly boxplot
pycloneccf.comb.means <- rbind(pyclone.nondr.means,pyclone.dr.means)
pycloneccf.comb.means %>%
  mutate(Growth.phenotype = fct_relevel(Growth.phenotype, "Clonal RAGE", "Non-clonal RAGE", "Clonal CALM", "Non-clonal CALM", "Clonal Context-dependent", "Non-clonal Context-dependent", "Clonal Growth attenuator", "Non-clonal Growth attenuator","Clonal No significant effect","Non-clonal No significant effect")) %>%
  ggplot(aes(x=Growth.phenotype, y=pycloneccf, color=Growth.phenotype))+
  geom_boxplot()+
  labs(title="PyCloneCCF of LUAD  tumours with Kras background and clonal vs non-clonal TSG mutations", y="PyCloneCCF", x="Growth phenotype")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')
ggsave("Desktop/new/pycloneccf.pdf", dpi=600, width = 8, height=5)

#need at least 2 rows for heatmap but only have 1 pycloneccf row
pycloneccf.hm <- rbind(pyclone.nondr.means,pyclone.dr.means)
colnames(pycloneccf.hm) <- "PyCloneCCF"
rownames(pycloneccf.hm) <- pycloneccf.hm$Growth.phenotype
pycloneccf.hm <- pycloneccf.hm[,2]
a <- heatmap(as.matrix(pycloneccf.hm), scale="column", 
             Colv=NA, Rowv=NA,
             col=cm.colors(256), xlab="PyCloneCCF", ylab="Growth phenotype", main="Heat map of PyCloneCCF for LUAD tumours with Kras background and clonal vs non-clonal TSG mutations",
             margins=c(6,12), cexRow = 1, cexCol = 1)

pyclone.comb.mean.std <- read.csv(file="Desktop/new/pycloneccf.comb.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
pyclone.comb.mean.std %>%
  mutate(Phenotype = fct_relevel(Phenotype, "Clonal RAGE", "Non-clonal RAGE", "Clonal CALM", "Non-clonal CALM", "Clonal Context-dependent", "Non-clonal Context-dependent", "Clonal Growth attenuator", "Non-clonal Growth attenuator","Clonal No significant effect","Non-clonal No significant effect")) %>%
  ggplot(aes(x=Phenotype, y=pycloneccf, color=Phenotype))+
  geom_boxplot()+
  geom_errorbar(aes(ymin=pyclone.comb.mean.std$std.low, ymax=pyclone.comb.mean.std$std.up, width=0.2))+
  labs(title="PyCloneCCF of LUAD  tumours with Kras background and clonal vs non-clonal TSG mutations", y="PyCloneCCF", x="Growth phenotype")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')
ggsave("Desktop/new/pycloneccf.mean.var.pdf", dpi=600, width = 8, height=5)

