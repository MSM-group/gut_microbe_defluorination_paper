# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 14:49:55 2023

@author: felderfl
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

csv_filepath = 'filepath'

overalltitle = 'C-Alpha_variance'
measure = 'Variance'

def create_heatmap(csv_filepath, save_filepath, dpi=600, custom_order=None):
    # Read the CSV file into a Pandas DataFrame
    data = pd.read_csv(csv_filepath)  # Do not use the first column as the index

    # Extract data for the heatmap
    heatmap_data = data.values[:, 1:]  # Exclude the first column ('File')
    row_labels = data['File'].values  # Use 'File' column as row labels (strings)
    column_labels = data.columns[1:]  # Use the header row as column names

    # Convert the data to a NumPy array with float data type
    heatmap_data = heatmap_data.astype(float)

    # Normalize the colors for each column
    norm = Normalize(vmin=heatmap_data.min(), vmax=heatmap_data.max())

    # Reorder rows based on the custom order
    if custom_order:
        custom_order_indices = [row_labels.tolist().index(row) for row in custom_order]
        heatmap_data = heatmap_data[custom_order_indices]
        row_labels = row_labels[custom_order_indices]

    # Create a heatmap using Matplotlib
    plt.figure(figsize=(12, 10), dpi =600)  # Increase the figure size for higher resolution
    plt.imshow(heatmap_data, cmap='viridis', interpolation='nearest', aspect='auto', norm=norm)

    # Display individual values as text annotations
    for i in range(heatmap_data.shape[0]):
        for j in range(heatmap_data.shape[1]):
            plt.text(j, i, f'{heatmap_data[i, j]:.2f}', ha='center', va='center', color='white')

    plt.colorbar()
    plt.title(f'{overalltitle}' ' Distances [nm]')
    plt.xlabel('Conserved Residue Pair')
    plt.ylabel('Non-Defluorinating Proteins / Defluorinating Proteins')

    # Set row and column labels
    plt.xticks(range(len(column_labels)), column_labels, rotation=90)
    plt.yticks(range(len(row_labels)), row_labels)

    # Draw horizontal lines to separate defluorinating and non-defluorinating rows
    defluorinating_row_count = 4
    plt.axhline(defluorinating_row_count - 0.5, color='white', linewidth=3)

    # Save the heatmap with higher resolution
    plt.savefig(save_filepath, dpi=dpi)
    plt.show()

# Custom order of rows
custom_order = ['WP178618037', 'GKG91607', 'WP087189991', 'Rha0230', 'WP118709078', 'PA0810', 'RSc1362', 'DehlVa', 'DhlB']


save_filepath = f'{overalltitle}''_Heatmap.png'  # You can specify the output image file path here
create_heatmap(csv_filepath, save_filepath, dpi=600, custom_order=custom_order)  # Adjust the DPI as needed



from scipy import stats

# Read the CSV file into a Pandas DataFrame
data = pd.read_csv(csv_filepath)

# Extract data for the boxplots, excluding the first column and the first row
boxplot_data = data.iloc[:, 1:].values

# Create two groups: rows 1-4 and rows 5-9
group1 = boxplot_data[(0,1,3,5,7), :]
group2 = boxplot_data[(2,4,6,8), :]

# Create independent boxplots for each column and perform a t-test
num_columns = boxplot_data.shape[1]

for i in range(num_columns):
    column_title = data.columns[i + 1]  # Get the value from the first row
    plt.figure(figsize=(6, 4), dpi=300)
    plt.boxplot([group1[:, i], group2[:, i]], labels=["Non-Defluorinating n=5", "Defluorinating n=4"])
    
    # Perform the Student's t-test
    t_stat, p_value = stats.ttest_ind(group1[:, i], group2[:, i])
    
    # Add the p-value to the title
    plt.title(f'{overalltitle}' ' 'f'{column_title} (p-value: {p_value:.4f})')
    plt.ylabel(f'{measure}'' (nm)')
    plt.savefig(f'{overalltitle}' f'{column_title}' '.png', dpi=300)
    plt.show()





