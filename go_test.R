basdir="data/"
'%ni%' <- Negate('%in%')
options("scipen"=100)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(biomaRt)

load(file=paste(basdir,"bestgen.Robject",sep=""))
counts<-mdat
rownames(counts)<-gnam[,1]
colnames(counts)<-meta[,6]

# Download & use Ensembl annotation
ensembl <- useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl",version = 113)

anno <- getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","gene_biotype"),
  mart = ensembl) |>
  rename(gene_id = ensembl_gene_id,symbol = external_gene_name,entrez = entrezgene_id,biotype = gene_biotype) |>
  filter(biotype == "protein_coding",!is.na(entrez),gene_id %in% rownames(counts)) |>
  distinct()

write.csv(
  anno,
  file = paste(basdir,"/human_gene_annotation_ensembl110_GRCh38.csv",sep=""),
  row.names = FALSE
)

## wrapper for GO test for each sample
get_go_cat<-function(counts, sample,anno) {
  print(sample)
  make_gsea_rank <- function(counts, sample, anno, min_count = 2) {
    dat <- data.frame(
      gene_id = sub("\\..*$", "", rownames(counts)),
      count = as.numeric(counts[, sample]),
      stringsAsFactors = FALSE
    ) |>
      left_join(
        anno |>
          mutate(gene_id = sub("\\..*$", "", gene_id)) |>
          dplyr::select(gene_id, entrez, biotype),
        by = "gene_id"
      ) |>
      filter(
        biotype == "protein_coding",
        !is.na(entrez),
        count >= min_count
      ) |>
      group_by(entrez) |>
      summarise(count = sum(count), .groups = "drop")
    
    # Within-sample abundance measure
    dat$logCPM <- log2(1e6 * dat$count / sum(dat$count) + 1)
    
    # Centered rank: positive = relatively abundant, negative = relatively low abundance
    dat$rank_stat <- rank(dat$logCPM, ties.method = "average") -
      (nrow(dat) + 1) / 2
    
    gene_rank <- dat$rank_stat
    names(gene_rank) <- dat$entrez
    
    sort(gene_rank, decreasing = TRUE)
  }
  
  my_rank <- make_gsea_rank(
    counts = counts,
    sample = sample,
    anno = anno,
    min_count = 5
  )

  set.seed( sample.int(1000,1) )
  
  my_go <- gseGO(
    geneList      = my_rank,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "ALL",
    minGSSize     = 10,
    maxGSSize     = 300,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    verbose       = FALSE,
    seed          = TRUE
  )

  # Convert Entrez IDs in leading-edge results to readable gene symbols
  my_go <- setReadable(
    my_go,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )
  
  return(my_go)
  }

## apply to samples with data (after first run removed those without data & shallower anal swab)
## second run: changed to at least 5 reads to reduce noise 
samples<-colnames(counts)[c(5,11,13,16,21)]
allgos<-lapply(samples,get_go_cat,counts = counts,anno = anno)
names(allgos)<-samples


# simplify for plotting & report lists
albp<-list()
alcc<-list()
for (sample in samples) {
  print(sample)
  simpset1<-allgos[[sample]] |> filter(NES > 0, ONTOLOGY == "BP", p.adjust < 0.05) |>  arrange(p.adjust)
  if(length(simpset1$Description)>0) {
    go_simple <- simpset1  |> 
    simplify(cutoff = 0.7,by = "p.adjust",select_fun = min  )
    albp[[sample]]<-go_simple
    }
  simpset2<-allgos[[sample]] |> filter(NES > 0, ONTOLOGY == "CC", p.adjust < 0.05) |>  arrange(p.adjust)
  if(length(simpset2$Description)>0) {
    go_simple <- simpset2  |> 
      simplify(cutoff = 0.7,by = "p.adjust",select_fun = min  )
    alcc[[sample]]<-go_simple
    }
  }

  
allplots<-list()
titl<-c("BP for La Doncella skin","CC for La Doncella skin","BP for El Nino anal swab","CC for El Nino anal swab")
j=1
for (sample in names(alcc)) {
  allplots[[j]]<-dotplot(albp[[sample]], showCategory = 10,font.size=10,title=titl[j]) + theme(plot.title = element_text(face = "bold",hjust = 0.5))
  j=j+1
  allplots[[j]]<-dotplot(alcc[[sample]], showCategory = 10,font.size=10,title=titl[j]) + theme(plot.title = element_text(face = "bold",hjust = 0.5))
  j=j+1
}

save(alcc, albp, allplots, file=paste(basdir,"/go_result.Robject",sep=""))

sitab5<-cbind(c(rep("BP for La Doncella skin",nrow(albp[[1]])),rep("CC for La Doncella skin",nrow(alcc[[1]])), rep("BP for El Nino anal swab",nrow(albp[[2]])),rep("CC for El Nino anal swab",nrow(alcc[[2]]))), rbind(as.data.frame(albp[[1]]),as.data.frame(alcc[[1]]),as.data.frame(albp[[2]]),as.data.frame(alcc[[2]])))
colnames(sitab5)[1]<-"Label"
write.table(sitab5,file=paste("/files/S5.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

  
