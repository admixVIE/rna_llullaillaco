#!/bin/bash
#SBATCH --job-name=starHg
#SBATCH --cpus-per-task=32
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=./logs/%x_%a.out
#SBATCH --error=./logs/%x_%a.err
#SBATCH --array=15

module load cutadapt
module load STAR
module load Subread


basdir=data/

ID=${SLURM_ARRAY_TASK_ID}
name=$(sed -n "${ID}"p $basdir/mummy_meta.txt | sed -e "s/\r//g" | cut -f 6  )
name=$(echo $name |  sed 's/La //g')
name=$(echo $name |  sed 's/El //g')

##for ID=14, no shallow data available, instead downsampling the deep data to 1.5M reads
#module load SeqKit
#seqkit sample -p 0.02 $basdir/Nino_Hair_3_DEEP_SEQ_1.fastq.gz -o $basdir/Nino_Hair_subsamp.fastq.gz
#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 30 $basdir/Nino_Hair_subsamp.fastq.gz -o $basdir/${name}_1_cutadapt.fastq.gz

## main procedure for most single-end data
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 30 $basdir/${name}_SHALLOW_SEQ.fastq.gz -o $basdir/${name}_1_cutadapt.fastq.gz
STAR --runThreadN 4 --genomeDir $basdir/ref_gen/STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh150 --limitBAMsortRAM 38482792886 --outFilterMatchNmin 16 --outFilterMultimapNmax 10 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.9 --readFilesIn <( zcat $basdir/${name}_1_cutadapt.fastq.gz ) --outFileNamePrefix $basdir/${name} --outSAMtype BAM SortedByCoordinate --scoreGapNoncan -2 --scoreGapGCAG -1 --scoreGapATAC -2 --scoreDelOpen -1 --alignIntronMin 15 --alignIntronMax 1 --alignMatesGapMax 100000000 --alignSJoverhangMin 2 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignSJDBoverhangMin 2 --peOverlapMMp 0.2 --seedSearchStartLmax 12 --alignEndsType EndToEnd
echo "map done"
featureCounts -T 1 -O -t exon -g gene_id -a $basdir/ref_gen/gencode.v47.annotation.gff3.gz -o $basdir/${name}.gene.counts.txt  $basdir/${name}Aligned.sortedByCoord.out.bam
echo "count done"



## only for ID=15 (anal swab), also do paired-end (subsampled to ~25M reads)
module load SeqKit
seqkit sample2 -s 12 -p 0.1 $basdir/Nino_Anus_1_DEEP_SEQ_1.fastq.gz -o $basdir/anus_deep_sub_1.fastq.gz 
seqkit sample2 -s 12 -p 0.1 $basdir/Nino_Anus_1_DEEP_SEQ_2.fastq.gz -o $basdir/anus_deep_sub_2.fastq.gz 
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 30 $basdir/anus_deep_sub_1.fastq.gz $basdir/anus_deep_sub_2.fastq.gz -o $basdir/anus_deep_cutadapt_1.fastq.gz -p $basdir/anus_deep_cutadapt_2.fastq.gz
STAR --runThreadN 32 --genomeDir $basdir/ref_gen/STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh150 --limitBAMsortRAM 38482792886 --outFilterMatchNmin 16 --outFilterMultimapNmax 10 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.9 --readFilesIn <( zcat $basdir/anus_deep_cutadapt_1.fastq.gz ) <(zcat $basdir/anus_deep_cutadapt_2.fastq.gz) --outFileNamePrefix $basdir/$name.deep --outSAMtype BAM SortedByCoordinate --scoreGapNoncan -2 --scoreGapGCAG -1 --scoreGapATAC -2 --scoreDelOpen -1 --alignIntronMin 15 --alignIntronMax 1 --alignMatesGapMax 100000000 --alignSJoverhangMin 2 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignSJDBoverhangMin 2 --peOverlapMMp 0.2 --seedSearchStartLmax 12 --alignEndsType EndToEnd
echo "map done"
featureCounts -p -T 1 -O -t exon -g gene_id -a $basdir/ref_gen/gencode.v47.annotation.gff3.gz -o $basdir/$name.deep.gene.counts.txt  $basdir/${name}.deepAligned.sortedByCoord.out.bam
echo "count done"

exit

