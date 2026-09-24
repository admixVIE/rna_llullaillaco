# Analysis of transcriptomic data from Llullaillaco mummies

This repository serves to document the steps taken for analysis.



## Raw data processing

Raw data was trimmed with cutadapt and mapped to the human genome using STAR: `starmap.sh`

At the end of the same script, it is also shown how data was processed the same way for Giardia intestinalis and other microorganisms.

DNA sequencing data from the same extracts was processed the same way, documented in `starmap_dna.sh`

Mapping to tRNA and miRNA indexes is documented here: `bowtie_s.sh`


## Downstream analysis & plotting

Processing of GTEx data and obtaining summary statistics is found in the R script `mum_gtex_stats.R`

The visualizations for Fig. 1 and Fig. 2 are in the R script `mum_rna_sumplot.R`

Visualization of mapped reads is documented in `mum_rna_viz.R`

Calling of genotypes from RNA and overlap with genotype tables from genomic data follows procedures described in `genotype_comp.sh`

Gene Ontology testing is shown in `go_test.R`

Other things are in `other.txt`


## Files

The `files` directory contains the output tables from these scripts, which were then further processed in LibreOffice for the final Supplementary Tables, as well the pdf figures 1-3.
