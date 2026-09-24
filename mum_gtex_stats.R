basdir="data/"
'%ni%' <- Negate('%in%')
options("scipen"=100)
library(dplyr)
library(ggplot2)
library(forcats)

# load and process gtex data
gtexmet<-read.table(paste(basdir,"gtex/GTEx_Analysis_2025-08-22_v11_Annotations_SampleAttributesDS.txt",sep=""),sep="\t",header=T,fill=T,comment.char="",quote="")
gtex<-read.table(paste(basdir,"gtex/GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.gct.gz",sep=""),sep="\t",header=F,fill=T)
colnames(gtex)<-gtex[2,]
gtex<-gtex[-c(1:2),]
gnam<-gtex[,c(1:2)]
gtexn<-as.data.frame(gtex[,-c(1:2)])
for (j in (1:ncol(gtexn))) { gtexn[,j]<-as.numeric(gtexn[,j])}

save(gtexn,gnam,file=paste(basdir,"gtex/gtexn.Robject",sep=""))

# subset to quality samples (total reads >50M, RIN>7, mapping rate >0.9, not flagged, not tumor) & save
gtexsub<-gtexmet[which(gtexmet$SMRIN>7 & gtexmet$SMRDTTL>50000000 & gtexmet$SMMAPRT >0.9 & gtexmet$SMTORMVE!="FLAGGED" & gtexmet$SMSTYP=="Normal"),]
realsa<-colnames(gtexn)
gtexsub<-gtexsub[which(gtexsub$SAMPID%in%realsa),]

gtexq<-gtexn[,which(realsa%in%gtexsub$SAMPID)]
gnam[,1]<-do.call(rbind,strsplit(gnam[,1],split="\\."))[,1]
gtexsub<-gtexsub[which(gtexsub$SAMPID%in%colnames(gtexq)),]
save(gtexq,gnam,gtexsub,file=paste(basdir,"gtex/gtexn_qual.Robject",sep=""))
load(file=paste(basdir,"gtex/gtexn_qual.Robject",sep=""))

# use DESeq for normalization on genes expressed in most tissues
keep_genes <- rowSums(gtexq > 10) > 0.1 * ncol(gtexq)

library("DESeq2")

## add mummies 
meta<-read.table(paste(basdir,"/mummy_meta.txt",sep=""),sep="\t",header=F)
meta[20:22,6]<-paste(meta[20:22,6],"_a",sep="")
meta[23:nrow(meta),6]<-meta[23:nrow(meta),1]
inds=as.character(meta[,4])
mtc<-cbind(c(5,15,9),c(20,21,22))

## mummy data
mumdat<-list()
for (ind in inds) { mumdat[[ind]]<-unlist(read.table(paste(basdir,"/",ind,".gene.counts.txt",sep=""),sep="\t",header=T)[,7]) }
mumgen<-unlist(read.table(paste(basdir,"/",ind,".gene.counts.txt",sep=""),sep="\t",header=T)[,1]) 
mumdat<-do.call(cbind,mumdat)
mumgen<-do.call(rbind,strsplit(mumgen,split="\\."))[,1]

gnam[,1]<-do.call(rbind,strsplit(gnam[,1],split="\\."))[,1]
mdat<-mumdat[match(gnam[,1],mumgen),]
mdat[is.na(mdat)]<-0
mdatraw<-mdat
for (j in (1:nrow(mtc))) { mdat[,mtc[j,1]]<-mdat[,mtc[j,1]]+mdat[,mtc[j,2]] }

## gene lists for GO testing
goli<-list()
for (ind in c("286400","286411","368901")) {
#  coff<-ifelse(ind=="368901",50,10)
  goli[[ind]]<-cbind(gnam[which(mdat[keep_genes,ind]>0),],mdat[which(mdat[keep_genes,ind]>0),])
}

save(goli,gnam,mdat, keep_genes,meta,file=paste(basdir,"/bestgen.Robject",sep=""))

mdat<-mdat[,-c(20:22)]
colnames(mdat)<-NULL
selc<-which(colSums(mdat)>=1000)

## normalize gtex data
dds <- DESeqDataSetFromMatrix(countData = gtexq[keep_genes, ],
                              colData = data.frame(samples = colnames(gtexq), tissue = gtexsub$SMTS),
                              design = ~ tissue)
dds_gtex <- estimateSizeFactors(dds)  
dds_gtex <- estimateDispersions(dds_gtex)
vsd <- vst(dds_gtex, blind = TRUE)  
normalized_data <- assay(vsd)

disp_fun_gtex <- dispersionFunction(dds_gtex)
geo_means_gtex <- exp(rowMeans(log(gtexq[keep_genes, ] + 1))) - 1

## normalize mummy data
dds_lla <- DESeqDataSetFromMatrix(countData = mdat[keep_genes,selc],
                                  colData   = data.frame(samples = inds[selc], tissue = rep("mummy",length(selc))),
                                  design    = ~ 1)
dds_lla <- estimateSizeFactors(dds_lla, geoMeans = geo_means_gtex)
dispersionFunction(dds_lla) <- disp_fun_gtex

vsd_lla <- varianceStabilizingTransformation(dds_lla, blind = FALSE)
vst_lla_mat <- assay(vsd_lla)
save(vst_lla_mat, geo_means_gtex, normalized_data, file=paste(basdir,"gtex/normdat_mummy.Robject",sep=""))


load(file=paste(basdir,"gtex/normdat_mummy.Robject",sep=""))

### correlation to mean expression
tisues<-unique(gtexsub$SMTSD)

mdats<-vst_lla_mat
crlate<-list()
for (tis in tisues) {
  rome<-rowMeans(normalized_data[,which(colnames(normalized_data)%in%gtexsub[which(gtexsub$SMTSD==tis),]$SAMPID)])
  for (mm in 1:ncol(mdats)) { exv<-which(mdats[,mm]>0); crlate[[tis]][mm]<-cor.test(mdats[,mm],rome,method="spearman",exact=F)$estimate  }
}


smpl<-meta[which(meta[,4]%in%inds[selc]),6]
sitab4<-do.call(rbind,crlate)
colnames(sitab4)<-smpl
sitab4<-cbind(rownames(sitab4),sitab4)
write.table(sitab4,file=paste("/files/S4.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

crlp<-data.frame(tissue=rep(tisues,each=ncol(mdats)),value=unlist(crlate),sample=rep(smpl,length(crlate)))
besttis<-c()
for (mm in 1:ncol(mdats)) { sel<-which(crlp$sample==smpl[mm]); besttis[mm]<-crlp[sel,][which(crlp[sel,]$value==max(crlp[sel,]$value)),]$tissue }

skin_tissues  <- c("Skin - Sun Exposed (Lower leg)",
                   "Skin - Not Sun Exposed (Suprapubic)")
colon_tissues <- c("Colon - Sigmoid", "Colon - Transverse",
                   "Small Intestine - Terminal Ileum")

tiscol<-tisues
tiscol[-grep("Skin|Colon|Intestine",tiscol)]<-"grey70"
tiscol[grep("Skin",tiscol)]<-"darkcyan"
tiscol[grep("Colon|Intestine",tiscol)]<-"darkorange"
names(tiscol)<-tisues

max_by_sample <- crlp %>%
  group_by(sample) %>%
  summarise(max_val = max(value, na.rm = TRUE),
            top_tissue = tissue[which.max(value)],
            .groups = "drop") %>%
  mutate(violin_fill = ifelse(top_tissue %in% names(tiscol),
                              "#dddddd", "#dddddd"))

crlp2 <- crlp %>%
  left_join(max_by_sample, by = "sample") %>%
  mutate(sample = fct_reorder(sample, max_val, .desc = TRUE))

pt_cols <- setNames(rep("grey70", length(unique(crlp$tissue))),
                    sort(unique(crlp$tissue)))
overlap <- intersect(names(tiscol), names(pt_cols))
pt_cols[overlap] <- tiscol[overlap]


skin_tissues  <- c("Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)")
colon_tissues <- c("Colon - Sigmoid", "Colon - Transverse", "Small Intestine - Terminal Ileum")

# Classify tissues
crlp2 <- crlp2 %>% 
  mutate(tissue_class = case_when(
    tissue %in% colon_tissues ~ "Colon",
    tissue %in% skin_tissues  ~ "Skin",
    TRUE                      ~ "Other"
  ))

pal_pts  <- c(Colon = "darkorange", Skin = "darkcyan", Other = "grey60")


gtplot<-ggplot(crlp2, aes(x = sample, y = value)) +
  geom_violin(aes(fill = violin_fill),
              adjust = 1.0, scale = "width",
              draw_quantiles = 0.5,
              color = "grey35", linewidth = 0.25, na.rm = TRUE) +
  scale_fill_identity(guide = "none") +
  geom_jitter(aes(color = tissue_class),
              width = 0.22, height = 0,
              size = 1.2, alpha = 0.85, na.rm = TRUE) +
  scale_color_manual(values = pal_pts, name = NULL,
                     breaks = c("Colon", "Skin", "Other"),
                     labels = c("Colon/Ileum", "Skin", "Other")) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  labs(x = NULL, y = "Correlation coefficient", title = "Correlation with GTEx per sample") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.8), color = NA),
    legend.key.size = grid::unit(3, "mm"),
    legend.text = element_text(size = 11)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

save(gtplot,crlp2, crlp, crlate,file=paste(basdir,"/mummy/gtex_plot.Robject",sep=""))


## build combined plot for GTex and GO enrichment
load(file=paste(basdir,"/go_result.Robject",sep=""))
load(file=paste(basdir,"/gtex_plot.Robject",sep=""))
library(clusterProfiler)
library(patchwork)

ap=gtplot
bp=allplots[[1]]
cp=allplots[[3]]

design <- "
AAB
AAC
"

combined_plot <-
  ap + bp + cp +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )


pdf("files/Fig2_gtex.pdf", width = 10, height = 8)
combined_plot
dev.off()




## rna stats
qry<-c("Number of input reads", "Uniquely mapped reads number", "Uniquely mapped reads %", "Average mapped length", "Number of splices: Annotated","Number of splices: Non-canonical", "% of reads mapped to multiple loci","% of reads unmapped: other", "% of reads unmapped: too short")
qrc<-c("Assigned","Unassigned_NoFeatures","Unassigned_MultiMapping")

inds=as.character(meta[,4])
mumstats<-list()
for (ind in inds) {
  map<-read.table(paste(basdir,"/",ind,"Log.final.out",sep=""),sep="|",header=F,fill=T)
  st<-c()
  for (qr in qry) { st[qr]<-gsub("\t","",map[grep(qr,map[,1]),2]) }
  cnt<-read.table(paste(basdir,"/",ind,".gene.counts.txt.summary",sep=""),sep="\t",header=T,fill=T)
  for (qr in qrc) { st[qr]<-cnt[grep(qr,cnt[,1]),2] }
  if (ind %in% inds[c(20,22:35)]) { nm<-meta[which(meta[,4]==ind),6];mumstats[[nm]]<-c(st,0,0);next }
  cnm<-sum(read.table(paste(basdir,"/",ind,".mirna.counts",sep=""),sep=" ",header=F,fill=T)[,2]);names(cnm)<-"miRNAs"
  tnm<-sum(read.table(paste(basdir,"/",ind,".trna.counts",sep=""),sep=" ",header=F,fill=T)[,2]);names(tnm)<-"tRNAs"
  nm<-meta[which(meta[,4]==ind),6]
  mumstats[[nm]]<-c(st,cnm,tnm)
}

mumstats<-as.data.frame(do.call(rbind,mumstats))

## count rRNA
load(file=paste(basdir,"/bestgen.Robject",sep=""))
mdat<-mdat[,-c(20:22)]
colnames(mdat)<-NULL
selc<-which(colSums(mdat)>=1000)

colnames(mdat)<-inds[-c(20:22)]
pcod<-read.table(paste(basdir,"/ref_gen/STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh150/geneInfo.tab",sep""),sep="\t",header=F,fill=T)
pcodi<-do.call(rbind,strsplit(pcod[which(pcod[,3]=="protein_coding"),1],split="\\."))[,1]
pcd<-which(gnam[,1]%in%pcodi)
mdatc<-mdat[pcd,]
gnap<-gnam[pcd,]
rrna<-do.call(rbind,strsplit(pcod[which(pcod[,3]%in%c("rRNA","Mt_rRNA")),1],split="\\."))[,1]
rcd<-which(gnam[,1]%in%rrna);mdatc<-mdatraw[rcd,]
rrnc<-colSums(mdatc)
mumstats$rRNAs<-rrnc


##assigned over uniquely mapped
asi<-round(as.numeric(mumstats$Assigned)/as.numeric(mumstats$"Uniquely mapped reads number")*100,2)
## spliced of assigned
spl<-round(as.numeric(mumstats$"Number of splices: Annotated")/as.numeric(mumstats$Assigned)*100,2)

mumstats$Percentage_Assigned<-asi
mumstats$Percentage_Spliced<-spl
mumstats$Percentage_rRNA<-mumstats$rRNAs/as.numeric(mumstats$Assigned)*100

write.table(mumstats,file=paste(basdir,"/mumstats.tsv",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

## for best sample
map<-read.table(paste(basdir,"/368901Log.final.out",sep=""),sep="|",header=F,fill=T)
cnt<-read.table(paste(basdir,"/368901.gene.counts.txt.summary",sep=""),sep="\t",header=T,fill=T)
stats<-as.numeric(c(gsub("\t","",map[c(5,8,11,16,23,25,30,32),2]),cnt[c(1,9,12),2],9468))
names(stats)<-c("Total reads","Uniquely mapped","Splices","Non-canonical","Mapped to multiple","Mapped to many", "Unmapped too short","Unmapped other","Assigned","Unassigned MultiMapping","Unassigned NoFeatures","rRNA")
cats<-list(stats[c(2,5,6,7,8)],stats[c(9:11)],stats[c(9,3:4)],stats[c(9,12)])
cats[[4]][1]<-cats[[4]][1]-cats[[4]][2]
cats[[3]][1]<-cats[[3]][1]-sum(cats[[3]][2:3])
names(cats)<-c("Mapping","Assignment","Splicing","rRNA")
save(stats,cats,file=paste(basdir,"/mumstat_best.Robject",sep=""))


### highest expressed
load(file=paste(basdir,"/bestgen.Robject",sep=""))
mdatc<-mdat[pcd,]
highestgen<-list()
highestc<-list()
for (ind in inds) {
  exda<-mdatc[,which(colnames(mdatc)==ind)]
  highestgen[[ind]]<-paste(gnap[order(exda,decreasing=T)[1:10],2],collapse=",")
  highestc[[ind]]<-paste(exda[order(exda,decreasing=T)[1:10]],collapse=",")
}

mumstats$Highest_genes<-highestgen
mumstats$Highest_gene_count<-highestc

gcoun<-cbind(meta[,6],unlist(highestgen),unlist(highestc))
colnames(gcoun)<-c("Sample","Highest_genes","Highest_gene_count")

write.table(gcoun,file=paste(basdir,"/highcount.tsv",sep=""),sep="\t",row.names=F,col.names=T,quote=F)



### stats on mapping DNA with same pipeline
qry<-c("Number of input reads", "Uniquely mapped reads number", "Uniquely mapped reads %", "Average mapped length", "Number of splices: Annotated","Number of splices: Non-canonical", "% of reads mapped to multiple loci","% of reads unmapped: other", "% of reads unmapped: too short")
qrc<-c("Assigned","Unassigned_NoFeatures","Unassigned_MultiMapping")

inds=c(as.character(meta[1:18,6]),"Nino_Anus.deep")
inds<-gsub("La |El ","",inds)

mumstatsd<-list()
for (ind in inds) {
  map<-read.table(paste(basdir,"/",ind,"Log.final.out",sep=""),sep="|",header=F,fill=T)
  st<-c()
  for (qr in qry) { st[qr]<-gsub("\t","",map[grep(qr,map[,1]),2]) }
  cnt<-read.table(paste(basdir,"/",ind,".gene.counts.txt.summary",sep=""),sep="\t",header=T,fill=T)
  for (qr in qrc) { st[qr]<-cnt[grep(qr,cnt[,1]),2] }
  mumstatsd[[ind]]<-st
}

mumstatsd<-as.data.frame(do.call(rbind,mumstatsd))

## count rRNA
mumdatd<-list()
for (ind in inds) { mumdatd[[ind]]<-unlist(read.table(paste(basdir,"/",ind,".gene.counts.txt",sep=""),sep="\t",header=T)[,7]) }
mumgend<-unlist(read.table(paste(basdir,"/",ind,".gene.counts.txt",sep=""),sep="\t",header=T)[,1]) 
mumdatd<-do.call(cbind,mumdatd)
mumgend<-do.call(rbind,strsplit(mumgend,split="\\."))[,1]

mdatd<-mumdatd[match(gnam[,1],mumgend),]
mdatd[is.na(mdatd)]<-0
mdatrawd<-mdatd
rrna<-do.call(rbind,strsplit(pcod[which(pcod[,3]%in%c("rRNA","Mt_rRNA")),1],split="\\."))[,1]
rcd<-which(gnam[,1]%in%rrna);mdatc<-mdatrawd[rcd,]
rrnc<-colSums(mdatc)
mumstatsd$rRNAs<-rrnc

asi<-round(as.numeric(mumstatsd$Assigned)/as.numeric(mumstatsd$"Uniquely mapped reads number")*100,2)
spl<-round(as.numeric(mumstatsd$"Number of splices: Annotated")/as.numeric(mumstatsd$Assigned)*100,2)

mumstatsd$Percentage_Assigned<-asi
mumstatsd$Percentage_Spliced<-spl
mumstatsd$Percentage_rRNA<-mumstatsd$rRNAs/as.numeric(mumstatsd$Assigned)*100

write.table(mumstatsd,file=paste(basdir,"/mumstats_dna.tsv",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

tab1<-mumstats[c(1:22),]
tab1$extract_id<-meta[1:22,2]
tab1$sequence_id<-meta[1:22,4]
tab1$label<-meta[1:22,7]
tab1$label2<-meta[1:22,6]
tab2<-mumstatsd
rownames(tab2)<-c(meta[1:18,6],"El Nino_Anus_1_a")
colnames(tab2)<-paste(colnames(tab2),"_DNA",sep="")
tab2$extract_id<-meta[c(1:18,15),2]
tab3<-merge(tab1,tab2,by.x=0,by.y=0,all=T)
tab3<-as.matrix(tab3)
write.table(tab3,file=paste("/files/S1.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

## then subset by deep and shallow manually


## for best sample as DNA
map<-read.table(paste(basdir,"/Nino_Anus.deepLog.final.out",sep=""),sep="|",header=F,fill=T)
cnt<-read.table(paste(basdir,"/Nino_Anus.deep.gene.counts.txt.summary",sep=""),sep="\t",header=T,fill=T)
statsd<-as.numeric(c(gsub("\t","",map[c(5,8,11,16,23,25,30,32),2]),cnt[c(1,9,12),2],568))
names(statsd)<-c("Total reads","Uniquely mapped","Splices","Non-canonical","Mapped to multiple","Mapped to many", "Unmapped too short","Unmapped other","Assigned","Unassigned MultiMapping","Unassigned NoFeatures","rRNA")
catsd<-list(statsd[c(2,5,6,7,8)],statsd[c(9:11)],statsd[c(9,3:4)],statsd[c(9,12)])
catsd[[4]][1]<-catsd[[4]][1]-catsd[[4]][2]
catsd[[3]][1]<-catsd[[3]][1]-sum(catsd[[3]][2:3])
names(catsd)<-c("Mapping","Assignment","Splicing","rRNA")
save(statsd,catsd,file=paste(basdir,"/mumstat_bestDNA.Robject",sep=""))


save(mumstats,mumstatsd,tab3,file=paste(basdir,"/mumstats.Robject",sep=""))




################################################################################
## Assesssing microorganism  data
orgs<-t(matrix(c("GCF_003290485.1_ASM329048v1", "Malassezia",
"GCF_003812505.1_ASM381250v1", "Staphylococcus",
"GCF_006739385.1_ASM673938v1", "Cutibacterium",
"GCF_000002435.2_UU_WB_2.1", "Giardia",
"GCF_000154385.1_ASM15438v1", "Faecali",
"GCF_000210435.1_ASM21043v1", "Clostridioles"),ncol=6))

## get genes with highest count per organism
allmi<-list()
bestmi<-list()
for (j in (1:nrow(orgs))) {
  girn<-read.table(paste(basdir,"/",orgs[j,1],".gene.counts.txt",sep=""),sep="\t", header=T)
  girn$Start<-do.call(rbind,strsplit(as.character(girn$Start),split=";"))[,1]
  girn$End<-do.call(rbind,strsplit(as.character(girn$End),split=";"))[,1]
  colnames(girn)[7]<-orgs[j,2]
  allmi[[j]]<-girn
  bestmi[[j]]<-head(girn[order(girn[,7],decreasing=T),])
  }


