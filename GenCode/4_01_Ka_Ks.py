'''
Calculating Ka/Ks (Non-synonymous/synonymous substitution rate) ratio
'''
import allel
import numpy as np
import pandas as pd
# Load the annotated VCF
callset = allel.read_vcf("variants_annotated.vcf", fields=['samples', 'calldata/GT','ANN'])

# Extract necessary information
samples = callset['samples']
genotypes = allel.GenotypeArray(callset['calldata/GT'])
annotations = callset['variants/ANN']

# Initialize a dictionary to store counts
# Structure: mutations[sample][gene]['nonsynonymous' or 'synonymous']
mutations = {sample: {} for sample in samples}

for i, ann in enumerate(annotations):
    for effect in ann.split(','):
        # Parse the effect details
        details = effect.split('|')
        effect_type = details[1] 
        gene_name = details[3]  

        # Determine mutation type
        mutation_type = None
        if 'missense_variant' in effect_type or 'stop_gained' in effect_type:
            mutation_type = 'nonsynonymous'
        elif 'synonymous_variant' in effect_type:
            mutation_type = 'synonymous'
        
        # Skip if not a target mutation type
        if not mutation_type:
            continue

        # Check each sample for the variant
        for sample_idx, sample in enumerate(samples):
            if np.any(genotypes[i, sample_idx] > 0):  # Variant present in this sample
                if gene_name not in mutations[sample]:
                    mutations[sample][gene_name] = {'nonsynonymous': 0, 'synonymous': 0}
                mutations[sample][gene_name][mutation_type] += 1

# Calculate Ka/Ks ratios and prepare the table
data = []
for sample, genes in mutations.items():
    for gene, counts in genes.items():
        ka = counts['nonsynonymous']
        ks = counts['synonymous']
        ratio = ka / ks if ks > 0 else 'inf'  # Handle division by zero
        data.append([sample, gene, ka, ks, ratio])

# Create a DataFrame
df = pd.DataFrame(data, columns=['Sample', 'Gene', 'Ka', 'Ks', 'Ka/Ks Ratio'])
print(df)

# Save the table to a CSV file
df.to_csv('ka_ks_ratios_output.csv', index=False)
