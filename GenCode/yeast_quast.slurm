#!/bin/bash

# Directory of assembled genomes
genome_directory="/home/l/liu-cong/work/y_genome/assembles"
reference_genome="/home/l/liu-cong/work/y_genome/refer/yeast_ref_genome.fna"


# Loop through each genome file
for genome in "$genome_directory"/*.fa; do
    #Extract filename without path and extension.
    genome_name=$(basename "$genome" .fa)

    #Run Quast
    quast -o "$genome_name" -r "$reference_genome" "$genome"

    echo "Quast completed for $genome_name"
done

echo "All Quast are completed."
