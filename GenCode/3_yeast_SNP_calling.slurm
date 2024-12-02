#!/bin/bash
#SBATCH --mail-user="yeast_genome"
#SBATCH --job-name=snp_calling
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mem=100G
#SBATCH -p compute
#SBATCH -c 10
#SBATCH -n 1
#SBATCH -t 2-20

strain="1039 1051 1063 1069"
# Define reference genomes for each strain: those are the ancestral strains.
declare -A ref_genomes
ref_genomes["1039"]="P2-H05_genome.fa"
ref_genomes["1051"]="P2-H04_genome.fa"
ref_genomes["1063"]="P2-H03_genome.fa"
ref_genomes["1069"]="P2-H09_genome.fa"

for s in $strain
do
        cd $s
        # Set the reference genome for the current strain
        ref_genome=${ref_genomes[$s]}
        # Index the reference genome (if not already indexed)
        if [ ! -f "${ref_genome}.bwt" ]; then
        echo "Indexing the reference genome for strain $s: $ref_genome"
        bwa index "$ref_genome"
        fi
        echo "Align Assembled Genomes to the Reference"
        for genome in "./"*.fa
        do
                genome_name=$(basename "$genome" .fa)
                bwa mem $ref_genome "$genome" > "$genome_name".sam
        done
        echo "Convert SAM to BAM, Sort, and Index"
        for samfile in "./"*.sam
        do
        bamfile=$(basename "$samfile" .sam).bam
        sorted_bamfile=$(basename "$samfile" .sam)_sorted.bam
        samtools view -S -b "$samfile" > "$bamfile"
        samtools sort "$bamfile" -o "$sorted_bamfile"
        samtools index "$sorted_bamfile"
        done
        echo "Variant Calling"
        for bamfile in "./"*_sorted.bam
        do
        vcf_file=$(basename "$bamfile" _sorted.bam).vcf
        bcftools mpileup -Ou -f $ref_genome "$bamfile" | bcftools call --ploidy 1 -mv -Ov -o "$vcf_file"
        bcftools view $vcf_file -Oz -o $vcf_file".gz"
        bcftools index  $vcf_file".gz"
        done
        echo "Merge VCF Files"
        bcftools merge "./"*.vcf.gz -Ov -o $s"_merged_variants.vcf"
        cd ..
done
echo "all vcf calling done"
