'''output the unique genes under positive selection from each group'''
import pandas as pd

# Load the dataset
df = pd.read_csv("ka_ks_ratios_output.csv")
df_filtered = df[df['Ka/Ks Ratio'] > 1]
# Get unique groups
groups = df_filtered['Group'].unique()

# Prepare a dictionary to hold unique genes for each group
unique_genes_info = {group: None for group in groups}

for group in groups:
    # Genes in the current group
    current_genes = set(df_filtered[df_filtered['Group'] == group]['Gene'])
    
    # Genes in other groups
    other_genes = set(df_filtered[df_filtered['Group'] != group]['Gene'])
    
    # Find genes unique to the current group
    unique_genes = current_genes - other_genes
    
    # Extract rows for unique genes
    unique_genes_df = df_filtered[(df_filtered['Gene'].isin(unique_genes)) & (df_filtered['Group'] == group)]
    
    # Store in dict
    unique_genes_info[group] = unique_genes_df[['Gene', 'Ka/Ks Ratio']]
    
    # Print or save to CSV
    print(f"Unique genes for {group}:")
    print(unique_genes_info[group], "\n")

for group, unique_genes_df in unique_genes_info.items():
    output_path = f'unique_genes_{group}.csv'
    unique_genes_df.to_csv(output_path, index=False)
    print(f"Saved unique genes for {group} to {output_path}")
