import os
import matplotlib.pyplot as plt
import pandas as pd

# Path to the folder containing the CSV files
folder_path = 'filepath'

# Get a list of all CSV files in the folder
csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]

# Check if there is at least one CSV file
if len(csv_files) == 0:
    raise ValueError("No CSV files found in the specified folder.")

# Use the first CSV file found
file_path = os.path.join(folder_path, csv_files[0])
data = pd.read_csv(file_path, sep=';')

# Assuming the first column is time and the second column is RMSD
rmsf = data.iloc[:, 1]

# Generating x values as count of data points
x_values = range(1, len(rmsf) + 1)

# Create figure and axis with specified size
fig, ax = plt.subplots(figsize=(8, 3))  # width, height

# Plotting RMSF data
ax.plot(x_values, rmsf, marker='o', markersize=0, linestyle='-', color='tab:red')

# Labeling axes and title
ax.set_xlabel('C-Alpha Index')
ax.set_ylabel('RMSF')
plt.title('Fluctuation of the C-Alpha Atoms in WP178')

# Set higher DPI for better resolution
plt.gcf().set_dpi(300)

# Show plot
plt.grid(True)
plt.show()
