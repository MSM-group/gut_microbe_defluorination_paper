import csv

def normalize_csv_values(csv_file):
    """
    Reads and normalizes the values from the CSV file using the updated normalization formula.
    Returns a dictionary mapping residue indices (1-based) to normalized values.
    """
    all_values = []
    normalized_values = {}

    # Read CSV values
    with open(csv_file, 'r') as file:
        reader = csv.reader(file, delimiter=';')
        for row in reader:
            row_values = []
            for value in row:
                try:
                    row_values.append(float(value))
                except ValueError:
                    row_values.append(0.0)  # Handle non-numeric values
            all_values.append(row_values)

    if not all_values:
        return {}

    # Find the minimum and maximum values
    min_value = min(min(row) for row in all_values if row)
    max_value = max(max(row) for row in all_values if row)

    # Normalize each value and assign to sequential residue indices starting from residue 2
    current_index = 2  # Start from residue 2
    normalized_values[1] = "00.00"  # Set residue 1's B-factor to 00.00

    for row in all_values:
        for value in row:
            normalized_value = ((value - min_value) / (max_value - min_value)) * 99.99 if (max_value - min_value) != 0 else 0
            normalized_values[current_index] = f"{normalized_value:05.2f}"
            current_index += 1

    return normalized_values

def modify_pdb_file_by_residue(input_file, output_file, normalized_csv):
    """
    Modifies the PDB file by injecting the normalized values from the CSV
    into the B-factor field (columns 61–66), matching by residue index.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):  # Modify only ATOM/HETATM lines
                try:
                    # Extract the residue index from columns 23–26
                    residue_index = int(line[22:26].strip())
                    # Find the corresponding normalized value for the residue
                    if residue_index in normalized_csv:
                        value_for_residue = normalized_csv[residue_index]
                        # Inject the normalized value into the B-factor field
                        line = line[:60] + f"{value_for_residue:>6}" + line[66:]
                except (ValueError, IndexError):
                    pass  # Skip if residue index is invalid or out of range
            outfile.write(line)

if __name__ == "__main__":
    # File paths
    csv_file = r"\\eawag\userdata\felderfl\My Documents\GitHub\gut_microbe_defluorination_paper\Dehalogenation_colorimetric_assays\Visualization on Protein Backbone\Dechlorination-Defluorination.csv"
    pdb_file = r"\\eawag\userdata\felderfl\My Documents\GitHub\gut_microbe_defluorination_paper\Dehalogenation_colorimetric_assays\Visualization on Protein Backbone\WP_1786180371.pdb"
    output_pdb_file = r"\\eawag\userdata\felderfl\My Documents\GitHub\gut_microbe_defluorination_paper\Dehalogenation_colorimetric_assays\Visualization on Protein Backbone\WP_1786180371_Dechlorination-Defluorination.pdb"

    # Normalize values from CSV
    normalized_csv = normalize_csv_values(csv_file)

    # Modify the PDB file using normalized CSV values, matched by residue index
    modify_pdb_file_by_residue(pdb_file, output_pdb_file, normalized_csv)

    print(f"Modified PDB file saved to {output_pdb_file}")
