#!/bin/bash
#SBATCH --job-name=gecomp2
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=0-10:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=./log/%x.out
#SBATCH --error=./log/%x.err
#SBATCH --array=1-19

## array only for the individual libraries

basdir=data/

# 1) extract homozygous alternative sites per indvidual
module load bcftools
bcftools view -m2 $basdir/Ind3_real.vcf.gz | bcftools query -f '%CHROM\t%POS\tALT[\t%GT]'| awk '$3=="A" || $3=="T" || $3=="G" || $3=="C"'| grep '1/1' > $basdir/ind3
bcftools view -m2 $basdir/Ind1_wo_file6.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ALT[\t%GT]'| awk '$3=="A" || $3=="T" || $3=="G" || $3=="C"'| grep '1/1' > $basdir/ind1
bcftools view -m2 $basdir/Ind2_real.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ALT[\t%GT]'| awk '$3=="A" || $3=="T" || $3=="G" || $3=="C"'| grep '1/1' > $basdir/ind2


# extract homozygous reference sites for those positions
bcftools view -R ind1_exon $basdir/Ind3_real.vcf.gz | bcftools query -f '%CHROM\t%POS[\t%GT]' | grep '0/0' > ind1_0_in_in3
bcftools view -R ind1_exon $basdir/Ind2_real.vcf.gz | bcftools query -f '%CHROM\t%POS[\t%GT]' | grep '0/0' > ind1_0_in_in2

bcftools view -R ind2_exon $basdir/Ind3_real.vcf.gz | bcftools query -f '%CHROM\t%POS[\t%GT]' | grep '0/0' > ind2_0_in_in3
bcftools view -R ind2_exon $basdir/Ind1_wo_file6.vcf.gz | bcftools query -f '%CHROM\t%POS[\t%GT]' | grep '0/0' > ind2_0_in_in1

bcftools view -R ind3_exon $basdir/Ind1_wo_file6.vcf.gz | bcftools query -f '%CHROM\t%POS[\t%GT]' | grep '0/0' > ind3_0_in_in1
bcftools view -R ind3_exon $basdir/Ind2_real.vcf.gz | bcftools query -f '%CHROM\t%POS[\t%GT]' | grep '0/0' > ind3_0_in_in2


# 2) extract exonic regions 
grep 'exon' $basdir/Homo_sapiens.GRCh37.87.gtf | cut -f 1,4,5 | sort -k 1,1 -k2,2n | mergeBed -i stdin > $basdir/Homo_sapiens.GRCh37.87.exon.bed
bedtools intersect -a <(awk '{print ,,,,}' OFS='\t' ind1) -b $basdir/Homo_sapiens.GRCh37.87.exon.bed | cut -f 1,2,4,5 > $basdir/ind1_exon
awk 'NR==FNR{a[,];next} (,) in a' $basdir/ind1_0_in_in2 $basdir/ind1_0_in_in3 > $basdir/ind1_0_in_in2_in3


# 3) call genotypes from RNA seq data (for each individual library)
name=$(sed -n "${SLURM_ARRAY_TASK_ID}"p $basdir/list | cut -f 2 | sed 's/\//\t/g' | cut -f 8)
out="$basdir/mummies"
in="$basdir/mapped"

module load samtools bcftools htslib
refgenome="$basdir/human_g1k_v37.fasta"

samtools index ${in}/${name}_Hg19Aligned.sortedByCoord.out.bam
for chr in {1..22}; do
bcftools mpileup -f ${refgenome} ${in}/${name}_Hg19Aligned.sortedByCoord.out.bam -r $chr --threads 2 -Oz -o ${out}/${name}_chr${chr}.mpileup.vcf.gz
tabix -f ${out}/${name}_chr${chr}.mpileup.vcf.gz
bcftools call -a GQ --ploidy 2 -mv --threads 2 -Oz ${out}/${name}_chr${chr}.mpileup.vcf.gz -r $chr -o ${out}/${name}_chr${chr}.calls.vcf.gz
tabix -f ${out}/${name}_chr${chr}.calls.vcf.gz
done

# 4) extract overlap (for each individual library)
for i in {1..22}; do bcftools view -R $basdir/ind1_exon -m2 $basdir/286407_chr${i}.calls.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ALT[\t%GT]'| awk '$3=="A" || $3=="T" || $3=="G" || $3=="C"' >> $basdir/286407_exon; done


