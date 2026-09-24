#!/bin/bash
#SBATCH --job-name=starHg
#SBATCH --cpus-per-task=32
#SBATCH --mem=40G
#SBATCH --time=8:00:00
#SBATCH --output=./logs/%x_%a.out
#SBATCH --error=./logs/%x_%a.err
#SBATCH --array=21

module load cutadapt
basdir=data/

ID=${SLURM_ARRAY_TASK_ID}
name=$(sed -n "${ID}"p $basdir/mummy_meta.txt | sed -e "s/\r//g" | cut -f 4  )

#echo $name
if [[ "$ID" -lt 20 ]]; then
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 30 $basdir/${name}/${name}.fastq.gz -o $basdir/${name}_1_cutadapt.fastq.gz
fi


module load STAR
module load Subread

if [[ "$ID" -lt 20 ]]; then
    STAR --runThreadN 4 --genomeDir $basdir/ref_gen/STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh150 --limitBAMsortRAM 38482792886 --outFilterMatchNmin 16 --outFilterMultimapNmax 10 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.9 --readFilesIn <( zcat $basdir/${name}_1_cutadapt.fastq.gz ) --outFileNamePrefix $basdir/${name} --outSAMtype BAM SortedByCoordinate --scoreGapNoncan -2 --scoreGapGCAG -1 --scoreGapATAC -2 --scoreDelOpen -1 --alignIntronMin 15 --alignIntronMax 1 --alignMatesGapMax 100000000 --alignSJoverhangMin 2 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignSJDBoverhangMin 2 --peOverlapMMp 0.2 --seedSearchStartLmax 12 --alignEndsType EndToEnd
    echo "work done"
    exit
fi

if [[ "$ID" -gt 19 ]]; then
    dir=$basdir/${name}/
    STAR --runThreadN 32 --genomeDir $basdir/ref_gen/STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh150 --limitBAMsortRAM 38482792886 --outFilterMatchNmin 16 --outFilterMultimapNmax 10 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.9 --readFilesIn <( zcat ${dir}/${name}_1_cutadapt.fastq.gz ) <( zcat ${dir}/${name}_2_cutadapt.fastq.gz ) --outFileNamePrefix $basdir/${name} --outSAMtype BAM SortedByCoordinate --scoreGapNoncan -2 --scoreGapGCAG -1 --scoreGapATAC -2 --scoreDelOpen -1 --alignIntronMin 15 --alignIntronMax 1 --alignMatesGapMax 100000000 --alignSJoverhangMin 2 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignSJDBoverhangMin 2 --peOverlapMMp 0.2 --seedSearchStartLmax 12 --alignEndsType EndToEnd
fi

echo "work done"
exit


for ID in {1..19};do
    name=$(sed -n "${ID}"p $basdir/mummy_meta.txt | sed -e "s/\r//g" | cut -f 4  )
    featureCounts -T 1 -O -t exon -g gene_id -a $basdir/ref_gen/gencode.v47.annotation.gff3.gz -o $basdir/${name}.gene.counts.txt  $basdir/${name}Aligned.sortedByCoord.out.bam
done


## paired end
for ID in {20..22};do
    name=$(sed -n "${ID}"p $basdir/mummy_meta.txt | sed -e "s/\r//g" | cut -f 4  )
    featureCounts -T 2 -p -B -C -O -t exon -g gene_id -a $basdir/ref_gen/gencode.v47.annotation.gff3.gz -o $basdir/${name}.gene.counts.txt  $basdir/${name}Aligned.sortedByCoord.out.bam
done

echo "count done"
exit






### Once for microorganisms in La D.
# GCF_003290485.1_ASM329048v1 = Malazzesia
# GCF_003812505.1_ASM381250v1 = Staphylococcus 
# GCF_006739385.1_ASM673938v1 = Cutibacterium
for org in GCF_003290485.1_ASM329048v1 GCF_003812505.1_ASM381250v1 GCF_006739385.1_ASM673938v1; do
echo $org
STAR --runThreadN 1 --runMode genomeGenerate --genomeDir $basdir/$org --genomeFastaFiles $basdir/$org/"$org"_genomic.fna --sjdbGTFfile $basdir/$org/"$org"_genomic.gtf
STAR --runThreadN 32 --genomeDir $basdir/$org --limitBAMsortRAM 38482792886 --outFilterMatchNmin 16 --outFilterMultimapNmax 10 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.9 --readFilesIn <( zcat $basdir/286400_cutadapt_dedupe.fastq.gz ) --outFileNamePrefix $basdir/"$org" --outSAMtype BAM SortedByCoordinate --scoreGapNoncan -2 --scoreGapGCAG -1 --scoreGapATAC -2 --scoreDelOpen -1 --alignIntronMin 15 --alignIntronMax 1 --alignMatesGapMax 100000000 --alignSJoverhangMin 2 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignSJDBoverhangMin 2 --peOverlapMMp 0.2 --seedSearchStartLmax 12 --alignEndsType EndToEnd
featureCounts -T 2 -B -C -O -t CDS -g gene_id -a $basdir/$org/"$org"_genomic.gtf -o $basdir/$org.gene.counts.txt  $basdir/"$org"Aligned.sortedByCoord.out.bam
done


### Once for microorganisms in El N.
# GCF_000002435.2_UU_WB_2.1 = Giardia
# GCF_000154385.1_ASM15438v1 = Faecali
# GCF_000953275.1_CD630DERM = Clostridioles
for org in GCF_000002435.2_UU_WB_2.1 GCF_000154385.1_ASM15438v1 GCF_000953275.1_CD630DERM; do
echo $org
STAR --runThreadN 1 --runMode genomeGenerate --genomeDir $basdir/$org --genomeFastaFiles $basdir/$org/"$org"_genomic.fna --sjdbGTFfile $basdir/$org/"$org"_genomic.gtf
STAR --runThreadN 32 --genomeDir $basdir/$org --limitBAMsortRAM 38482792886 --outFilterMatchNmin 16 --outFilterMultimapNmax 10 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.9 --readFilesIn <( zcat $basdir/368901_1_cutadapt.fastq.gz ) <( zcat $basdir/368901_2_cutadapt.fastq.gz ) --outFileNamePrefix $basdir/"$org" --outSAMtype BAM SortedByCoordinate --scoreGapNoncan -2 --scoreGapGCAG -1 --scoreGapATAC -2 --scoreDelOpen -1 --alignIntronMin 15 --alignIntronMax 1 --alignMatesGapMax 100000000 --alignSJoverhangMin 2 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignSJDBoverhangMin 2 --peOverlapMMp 0.2 --seedSearchStartLmax 12 --alignEndsType EndToEnd
featureCounts -T 2 -p -B -C -O -t CDS -g gene_id -a $basdir/$org/"$org"_genomic.gtf -o $basdir/$org.gene.counts.txt  $basdir/"$org"Aligned.sortedByCoord.out.bam
done

exit


