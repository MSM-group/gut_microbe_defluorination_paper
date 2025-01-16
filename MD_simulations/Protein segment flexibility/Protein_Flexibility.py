# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 15:25:59 2024

@author: felderfl
"""
import os
import matplotlib.pyplot as plt
import pandas as pd

# Path to the folder containing the CSV files
folder_path = 'filepath'

# Get a list of all CSV files in the folder
csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]

# Colors for the boxplots
colors = ["#c5b100", "#1E90FF", "#FF00FF", "#32CD32", "#999999", "#B03060", "#FFFFFF"]

# Create a list to hold RMSD data for each CSV file
rmsd_data = []

# Read RMSD data from each CSV file
for file in csv_files:
    file_path = os.path.join(folder_path, file)
    data = pd.read_csv(file_path)
    
    # Assuming the second column is RMSD
    rmsd = data.iloc[:, 1]
    
    # Exclude the first 10,000 and last 1,000 measurement points
    rmsd = rmsd[10000:-1000]
    
    rmsd_data.append(rmsd)

# Change default font to Arial
plt.rcParams['font.family'] = 'Arial'

# Create a figure and axes
fig, ax = plt.subplots(figsize=(8, 12))

# Create a boxplot for the RMSD data
boxprops = dict(linestyle='-', linewidth=1.5)
medianprops = dict(linestyle='-', linewidth=3, color='black')
flierprops = dict(marker='o', markersize=5, linestyle='none', markerfacecolor='gray', markeredgecolor='gray')

# Plot each boxplot with a specific color
for i in range(len(rmsd_data)):
    box = ax.boxplot(rmsd_data[i], positions=[i], widths=0.6, vert=True, patch_artist=True, boxprops=boxprops, medianprops=medianprops, flierprops=flierprops, showfliers=False)
    for patch in box['boxes']:
        patch.set_facecolor(colors[i % len(colors)])

# Adding labels and title with increased font size
ax.set_xticklabels(['1', '2', '3', '4', '5', '6', '1-6'], fontsize=18)
ax.set_xlabel('Protein segment', fontsize=20)
ax.set_ylabel('RMSD (Flexibility)', fontsize=20)

# Set y-axis tick labels fontsize
ax.tick_params(axis='y', labelsize=18)

# Rotate x-axis labels to 0 degrees
plt.xticks(rotation=0)

# Set higher DPI for better resolution
fig.set_dpi(1000)

# Show plot
plt.show()






