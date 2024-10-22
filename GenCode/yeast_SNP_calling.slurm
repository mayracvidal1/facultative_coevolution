#!/bin/bash
strain="1063 1069 1051"
ref_genome="../../refer/yeast_ref_genome.fna"
for s in $strain
do
        cd $s
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
        bcftools mpileup -Ou -f $ref_genome "$bamfile" | bcftools call -mv -Ov -o "$vcf_file"
        bcftools view $vcf_file -Oz -o $vcf_file".gz"
        bcftools index  $vcf_file".gz"
        done
        echo "Merge VCF Files"
        bcftools merge "./"*.vcf.gz -Ov -o $s"_merged_variants.vcf"
        cd ..
done
echo "all vcf calling done"
