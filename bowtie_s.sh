#!/bin/bash
#SBATCH --job-name=bowtie
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --output=./logs/%x_%a.out
#SBATCH --error=./logs/%x_%a.err
#SBATCH --array=1-21

## run once
#module load Bowtie SAMtools
#bowtie-build hsa-mature.fas mirgenedb_mature
#awk '/^>/ {print; next} {print $0 "CCA"}' hg38-mature-tRNAs.fa > hg38-mature-tRNAs-CCA.fa
#bowtie-build hg38-mature-tRNAs-CCA.fa tRNA_mature

basdir=data/
mirdex=$basdir/ref_gen/smallrna/mirgenedb_mature
tdex=$basdir/ref_gen/smallrna/tRNA_mature

ID=${SLURM_ARRAY_TASK_ID}
name=$(sed -n "${ID}"p $basdir/mummy_meta.txt | sed -e "s/\r//g" | cut -f 4  )

module load cutadapt fastp

echo $name
if [[ "$ID" -lt 20 ]]; then
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 15 -q 20,20 --max-n 0 --trim-n $basdir/${name}/${name}.fastq.gz -o $basdir/${name}_mi_cutadapt.fastq.gz
fi

if [[ "$ID" -gt 19 &&  "$ID" -lt 23 ]]; then
    dir=$basdir/${name}/
fi


# double index
if [[ "$ID" -gt 19 ]]; then
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 15 -q 20,20 --max-n 0 --trim-n $basdir/${name}/${name}*R1*.fastq.gz $basdir/${name}/${name}*R2*.fastq.gz  -o $basdir/${name}_mi_cutadapt_1.fastq.gz -p $basdir/${name}_mi_cutadapt_2.fastq.gz
    fastp -i $basdir/${name}_mi_cutadapt_1.fastq.gz -I $basdir/${name}_mi_cutadapt_2.fastq.gz -m --merged_out $basdir/${name}_mi_cutadapt.fastq.gz -w 8
fi


module load Bowtie SAMtools

if [[ "$ID" -lt 23 ]]; then
    ## mirna
    bowtie -v 1 -a --best --strata -S $mirdex $basdir/${name}_mi_cutadapt_1.fastq.gz | samtools view -b -q 1 | samtools sort -o $basdir/${name}.mirna.bam
    samtools index $basdir/${name}.mirna.bam
    # Count per miRNA (fractional assignment for multi-mappers)
    samtools idxstats $basdir/${name}.mirna.bam | awk '{print $1,$3}' > $basdir/${name}.mirna.counts
    # miRNA: RPM (reads per million miRNA-mapped reads)
    total_miRNA=$(awk '{sum+=$2} END {print sum}' $basdir/${name}.mirna.counts)
    awk -v total=$total_miRNA '{print $1,$2/total*1e6}' $basdir/${name}.mirna.counts > $basdir/${name}.mirna.rpm
    
    ## trna
    bowtie -v 1 -a --best --strata -k 50 -S $tdex $basdir/${name}_mi_cutadapt_1.fastq.gz | samtools view -b -q 1 | samtools sort -o $basdir/${name}.trna.bam
    samtools index $basdir/${name}.trna.bam
    # Count at anticodon/isotype level (use tRNAscan-SE anticodon annotation)
    samtools idxstats $basdir/${name}.trna.bam | awk '{print $1,$3}' > $basdir/${name}.trna.counts
    # tRNA: RPM (reads per million tRNA-mapped reads)
    total_tRNA=$(awk '{sum+=$2} END {print sum}' $basdir/${name}.trna.counts)
    awk -v total=$total_tRNA '{print $1,$2/total*1e6}' $basdir/${name}.trna.counts > $basdir/${name}.trna.rpm
    echo "work done"
    exit
fi



exit
