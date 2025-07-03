import pandas as pd
import os

# File path - update this to your actual file path
file_path = "Combined_Basic_Results.tsv"

# Read the TSV file
# Using low_memory=False to avoid the DtypeWarning
df = pd.read_csv(file_path, sep='\t', header=None, low_memory=False)

# If the file actually has headers, you should modify to:
# df = pd.read_csv(file_path, sep='\t', low_memory=False)

# If there's no header, assign column names
if len(df.columns) == 7:  # Assuming 7 columns based on your example
    df.columns = ['name', 'sample', 'Length', 'Chromosome', 'Start', 'End', 'Source']
else:
    # Adjust this part if the actual file has a different number of columns
    print(f"Warning: Expected 7 columns but found {len(df.columns)}. Please check the file structure.")

# 1. Count occurrences of each sample
sample_counts = df['sample'].value_counts().reset_index()
sample_counts.columns = ['Sample', 'Count']
print("\nSample counts:")
print(sample_counts)

# 2. Count length distributions for each sample
# Create length categories
def length_category(length):
    try:
        # Convert length to integer if it's a string
        length_val = int(length) if isinstance(length, str) else length
        
        if length_val <= 1000:
            return '0-1000'
        elif length_val <= 10000:
            return '1000-10000'
        else:
            return '>10000'
    except (ValueError, TypeError):
        # In case of conversion error, return 'Unknown'
        return 'Unknown'

# Apply the length category function
df['length_category'] = df['Length'].apply(length_category)

# Group by sample and length category
length_distribution = df.groupby(['sample', 'length_category']).size().reset_index()
length_distribution.columns = ['Sample', 'Length Range', 'Count']

# Calculate total counts per sample for percentage calculation
sample_totals = df.groupby('sample').size().reset_index()
sample_totals.columns = ['Sample', 'Total']

# Merge the length distribution with sample totals
length_distribution = pd.merge(length_distribution, sample_totals, on='Sample')

# Calculate percentages
length_distribution['Percentage'] = (length_distribution['Count'] / length_distribution['Total'] * 100).round(2)

# Sort by sample and length range for better readability
length_distribution = length_distribution.sort_values(['Sample', 'Length Range'])

# Export results to CSV
length_distribution.to_csv('sample_length_distribution.csv', index=False)
sample_counts.to_csv('sample_counts.csv', index=False)

# Print the results
print("\nLength distribution by sample:")
for sample in df['sample'].unique():
    print(f"\nSample: {sample}")
    sample_data = length_distribution[length_distribution['Sample'] == sample]
    for _, row in sample_data.iterrows():
        print(f"  {row['Length Range']}: {row['Count']} ({row['Percentage']}%)")

print("\nResults have been saved to 'sample_length_distribution.csv' and 'sample_counts.csv'")
