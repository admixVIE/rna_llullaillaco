# Packages
basdir="data/"
'%ni%' <- Negate('%in%')
options("scipen"=100)

library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)
library(GenomeInfoDb)
library(grid)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(forcats)
library(tidytext)

# a) MUC12 in Nino
# b) FLG in La Donc.

# Inputs
gff <- file.path(basdir, "/ref_gen/gencode.v47.annotation.gff3.gz")
ann <- import(gff)
# these bam files are an artifact, indeed they contain all reads mapped from the respective libraries
plist<-list(c( "chr7", 101008000, 101014159,file.path(basdir, "/MUCun.bam"),"ENST00000536621.6","MUC12 in El Nino"),
            c("chr1",152296396,152331008,file.path(basdir, "/FLGun.bam"),"ENST00000368799.2","FLG in La Doncella"))


plot_one <- function(entry) {
  chr=entry[1]
  from=as.numeric(entry[2])
  to=as.numeric(entry[3])
  bam <- entry[4]
  trsp <- entry[5]
  gen<-entry[6]
  
  region <- GRanges(chr, IRanges(from, to))
  
  # Gene/exon track
  exonM <- ann[ann$type=="exon" & seqnames(ann)==chr & start(ann)<=to & end(ann)>=from]
  exonM$gene <- exonM$gene_id; exonM$transcript <- exonM$transcript_id
  exonM<-exonM[which(exonM$transcript==trsp),]
  exonM$symbol <- if ("gene_name" %in% names(mcols(exonM))) exonM$gene_name else exonM$gene_id
  gtrack <- GeneRegionTrack(exonM, genome="hg38", chromosome=chr, name="GENCODE",
                            transcriptAnnotation="symbol")
  
  # Alignments track from BAM (requires index)
  #if (!file.exists(paste0(bam, ".bai"))) indexBam(bam)
  atrack <- AlignmentsTrack(bam, genome="hg38", isPaired=FALSE, name="RNA-seq")
  displayPars(atrack) <- list(lwd.sashimiMax = 5,sashimiHeight=0.15, sashimiNumbersCex = 0.5,coverageHeight=0.05)
  
  # 3) Axis track
  axis <- GenomeAxisTrack()
  
  # 4) Plot
  plotTracks(list(axis, gtrack, atrack),chromosome = chr, from = from, to = to,
             type = c("coverage", "pileup", "sashimi"),sashimiNumbers = TRUE,
             background.title = "white",main=gen,add=T, cex.main=0.85)
  
}


## plot most highly expressed genes
load(file=paste(basdir,"/mumstats.Robject",sep=""))
ta1<-cbind(tab1$label2,tab1$Highest_genes,tab1$Highest_gene_count)[c(5,16,21),]
tap<-list()
for (j in (1:nrow(ta1))) {
  sus<-as.numeric(unlist(strsplit(as.character(ta1[j,3]),split=",")))
  names(sus)<-unlist(strsplit(as.character(ta1[j,2]),split=","))
  tap[[unlist(ta1[j,1])]]<-sus
}

names(tap)<-c("La Doncella Skin","El Nino Ear", "El Nino anal swab")

tap_df <- imap_dfr(
  tap,
  ~ enframe(.x, name = "gene", value = "count") |>
    mutate(sample = .y)
)


tap_df <- tap_df |>
  mutate(
    gene_plot = reorder_within(gene, count, sample)
  )
tap_df

top_genes_plot <- ggplot(
  tap_df,
  aes(x = count, y = gene_plot)
) +
  geom_col(
    width = 0.75,
    fill = "#3C78A8"
  ) +
  geom_text(
    aes(label = count),
    hjust = -0.15,
    size = 3.2
  ) +
  scale_y_reordered() +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  facet_wrap(
    ~ sample,
    nrow = 1,
    scales = "free"
  ) +
  labs(
    x = "Uniquely assigned reads",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    strip.background = element_blank(),
    axis.text.y = element_text(
      face = "italic",
      size = 9
    ),
    panel.spacing = unit(1.5, "lines")
  )


# Open one PDF page and lay out two panels
pdf("/files/Fig3_RNAalign.pdf", width = 10, height = 11)
grid::grid.newpage()
grid::pushViewport(
  grid::viewport(
    layout = grid::grid.layout(
      nrow = 4,
      ncol = 1,
      heights = grid::unit.c(
        grid::unit(3, "null"),  # A: top genes
        grid::unit(4, "null"),  # B
        grid::unit(3, "mm"),    # spacer
        grid::unit(6, "null")   # C
      )
    )
  )
)

## Panel A: top genes
grid::pushViewport(
  grid::viewport(
    layout.pos.row = 1,
    layout.pos.col = 1
  )
)

print(top_genes_plot, newpage = FALSE)

grid::grid.text(
  "A",
  x = grid::unit(1, "mm"),
  y = grid::unit(1, "npc") - grid::unit(1, "mm"),
  just = c("left", "top"),
  gp = grid::gpar(fontsize = 13, fontface = 2)
)

grid::popViewport()

## Panel B
grid::pushViewport(
  grid::viewport(
    layout.pos.row = 2,
    layout.pos.col = 1
  )
)

plot_one(plist[[1]])

grid::grid.text(
  "B",
  x = grid::unit(1, "mm"),
  y = grid::unit(1, "npc") - grid::unit(1, "mm"),
  just = c("left", "top"),
  gp = grid::gpar(fontsize = 13, fontface = 2)
)

grid::popViewport()

## Panel C
grid::pushViewport(
  grid::viewport(
    layout.pos.row = 4,
    layout.pos.col = 1
  )
)

plot_one(plist[[2]])

grid::grid.text(
  "C",
  x = grid::unit(1, "mm"),
  y = grid::unit(1, "npc") - grid::unit(1, "mm"),
  just = c("left", "top"),
  gp = grid::gpar(fontsize = 13, fontface = 2)
)

grid::popViewport()

## Exit the outer layout viewport
grid::popViewport()

dev.off()




