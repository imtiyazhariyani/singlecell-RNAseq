import pandas as pd
import anndata
import loompy
with loompy.connect("s_fca_biohub_antenna_10x.loom", validate=False) as ds:
	print(ds.shape)
	print(ds.ra.keys())
	print(ds.ca.keys())
	genes = ds.ra["Gene"]
	cells = ds.ca["CellID"]
	df = pd.DataFrame(ds[:, :], index=genes, columns=cells).T

df['Total_Mapped_Reads'] = df.sum(axis=1)

# Assuming 'umi_cell_info.csv' contains the cell type annotations with the 'UMI' column
# Replace the file path with the actual path to your CSV file
cell_annotations = pd.read_csv('umi_cell_info.csv')

# Extract the UMI from the index and create a new column 'UMI' in the DataFrame
df['UMI'] = df.index.to_series().apply(lambda x: x.split('-')[0])

# Merge cell annotations with the DataFrame based on the 'UMI' column
merged_df = df.merge(cell_annotations, on='UMI', how='left')

# Convert 'Total_Mapped_Reads' to integers
merged_df['Total_Mapped_Reads'] = merged_df['Total_Mapped_Reads'].astype(int)

# Group the data by cell type and calculate the counts for each cell type
cell_type_counts = merged_df.groupby('Cell_Type')['Total_Mapped_Reads'].sum()

# Print the counts for each cell type
print(cell_type_counts)

