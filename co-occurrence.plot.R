#tracerx
dat<-read.table("Desktop/Txt_R/CL_TRACERx_filtered_plot.txt",header=T)
cor.test(dat$Tubaseq_growth,dat$TRACERx_filtered_co_occurrence)
#cor.test(dat$Tubaseq_growth,dat$TRACERx_co.occurrence, method="spearman")
a<-ggplot(dat, aes(x=Tubaseq_growth, y = TRACERx_filtered_co_occurrence,color=GeneID))+
  geom_point(size = 3)+
  geom_smooth(method="lm",color="black")+
  labs(title="Gene interactions with TP53 to drive tumour growth\n(cor=0.06294855, p-value=0.6777)", y="TRACERx data - enrichment for co-occurrence with TP53", x="Tuba-seq data - growth advantage of being on KRAS+TP53 background")+
  theme_classic()
a
ggsave("Desktop/new.coc.tracerx.pdf",dpi=600)

#genie
dat.1<-read.table("Desktop/Txt_R/CL_GENIE_filtered_plot.txt",header=T) 
cor.test(dat.1$Tubaseq_growth,dat.1$GENIE_filtered_co_occurrence)
#cor.test(dat.1$Tubaseq_growth,dat.1$GENIE_co.occurrence, method="spearman")
a<-ggplot(dat.1, aes(x=Tubaseq_growth, y = GENIE_filtered_co_occurrence,color=GeneID))+
  geom_point(size = 3)+
  geom_smooth(method="lm",color="black")+
  labs(title="Gene interactions with TP53 to drive tumour growth\n(cor=-0.02297813, p-value=0.8881)", y="GENIE data - enrichment for co-occurrence with TP53", x="Tuba-seq data - growth advantage of being on KRAS+TP53 background")+
  theme_classic()
a
ggsave("Desktop/new.coc.genie.pdf",dpi=600)

#tcga
dat.2<-read.table("Desktop/Txt_R/TCGA_plot.txt",header=T) 
cor.test(dat.2$Tuba_Seq_Growth_Ix,dat.2$tcga_cooccurrence)
#cor.test(dat.2$Tuba_Seq_Growth_Ix,dat.2$tcga_cooccurrence, method="spearman")
a<-ggplot(dat.2, aes(x=Tuba_Seq_Growth_Ix, y = tcga_cooccurrence,color=GeneID))+
  geom_point(size = 3)+
  geom_smooth(method="lm",color="black")+
  labs(title="Gene interactions with TP53 to drive tumour growth\n(cor=0.3260421, p-value=0.02702)", y="TCGA data - enrichment for co-occurrence with TP53", x="Tuba-seq data - growth advantage of being on KRAS+TP53 background")+
  theme_classic()
a
ggsave("Desktop/new.coc.tcga.pdf",dpi=600)