#!/bin/bash
#SBATCH --mail-user="yeast_genome"
#SBATCH --job-name=genome_assemble
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mem=100G
#SBATCH -p compute
#SBATCH -c 10
#SBATCH -n 1
#SBATCH -t 2-20

bwa index ../refer/yeast_ref_genome.fna # Index the reference genome (.fna) for BWA alignment

id="P1-A01 P1-A02 P1-A03 P1-A04 P1-A05 P1-A06 P1-A07 P1-A08 P1-A09 P1-A10 P1-A11 P1-A12 P1-B01 P1-B02 P1-B03 P1-B04 P1-B05 P1-B06 P1-B07 P1-B08 P1-B09 P1-B10 P1-B11 P1-B12 P1-C01 P1-C02 P1-C03 P1-C05 P1-C06 P1-C07 P1-C08 P1-C09 P1-C10 P1-C11 P1-C12 P1-D01 P1-D02 P1-D03 P1-D04 P1-D05 P1-D06 P1-D07 P1-D08 P1-D09 P1-D10 P1-D11 P1-D12 P1-E01 P1-E02 P1-E03 P1-E04 P1-E05 P1-E06 P1-E07 P1-E08 P1-E09 P1-E10 P1-E11 P1-E12 P1-F01 P1-F02 P1-F03 P1-F04 P1-F05 P1-F06 P1-F07 P1-F08 P1-F09 P1-F10 P1-F11 P1-F12 P1-G01 P1-G02 P1-G03 P2-A01 P2-A02 P2-A03 P2-A04 P2-A05 P2-A06 P2-A07 P2-A08 P2-A09 P2-A10 P2-A11 P2-A12 P2-B01 P2-B02 P2-B03 P2-B04 P2-B05 P2-B06 P2-B07 P2-B08 P2-B09 P2-B10 P2-B11 P2-B12 P2-C01 P2-C02 P2-C03 P2-C04 P2-C05 P2-C06 P2-C07 P2-C08 P2-C09 P2-C10 P2-C11 P2-C12 P2-D01 P2-D02 P2-D03 P2-D04 P2-D05 P2-D07 P2-D09 P2-D10 P2-D11 P2-D12 P2-E01 P2-E03 P2-E04 P2-E05 P2-E06 P2-E07 P2-E08 P2-E09 P2-E10 P2-E11 P2-E12 P2-F01 P2-F02 P2-F03 P2-F04 P2-F05 P2-F06 P2-F07 P2-F08 P2-F09 P2-F10 P2-F11 P2-F12 P2-G01 P2-G02 P2-G03 P2-G04 P2-G05 P2-G06 P2-G07 P2-G08 P2-G09 P2-G10 P2-G11 P2-G12 P2-H01 P2-H03 P2-H04 P2-H05 P2-H06 P2-H07 P2-H09 P2-H10 P2-H11 P2-H12"

for i in $id
do
        echo "starting genome assemlby of" $i
        echo "aligment"
        bwa mem ../refer/yeast_ref_genome.fna "../trimed_data/"$i* > $i"_aligned.sam"
        samtools view -S -b $i"_aligned.sam" > $i"_aligned.bam"
        samtools sort $i"_aligned.bam" > $i"_aligned_sort.bam"
        samtools index $i"_aligned_sort.bam"
        echo "variant calling"
        bcftools mpileup -O b -o $i"_raw_variants.bcf" -f ../refer/yeast_ref_genome.fna $i"_aligned_sort.bam"
        bcftools call -m -v -p 2 -o $i"_final_variants.vcf" $i"_raw_variants.bcf"
        bcftools view $i"_final_variants.vcf" -Oz -o $i"_final_variants.vcf.gz"
        bcftools index  $i"_final_variants.vcf.gz"
        "assemlby and consensus"
        bcftools consensus -f ../refer/yeast_ref_genome.fna  $i"_final_variants.vcf.gz" > $i"_genome.fa"        
        echo "fone for" $i
done
