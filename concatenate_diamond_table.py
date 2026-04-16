import csv
from collections import defaultdict

# **Block 1: Processing the Input File**

# Define file names for the first block
input_file = 'diamond.txt'
processed_output_file = 'processed_table.txt'

# Define the index for the staxids column
STAXIDS_COLUMN_INDEX = 3  # Replace with the actual column index

# Read the input file and process data
data = defaultdict(lambda: defaultdict(list))

with open(input_file, 'r') as infile:
    reader = csv.reader(infile, delimiter='\t')  # Adjust delimiter if necessary
    for row in reader:
        sample_id = row[0]
        for i, value in enumerate(row):
            if i == STAXIDS_COLUMN_INDEX:
                try:
                    # Convert each taxid to an integer if possible
                    numeric_values = [str(int(float(v.strip()))) for v in value.split(';')]
                    value = ';'.join(numeric_values)
                except ValueError:
                    pass  # Handle the case where conversion to int is not possible
            data[sample_id][i].append(value)

# Write the processed data to the processed output file
with open(processed_output_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile, delimiter='\t')  # Adjust delimiter if necessary
    for sample_id, columns in data.items():
        concatenated_row = [sample_id] + [';'.join(columns[i]) for i in range(len(columns))]
        writer.writerow(concatenated_row)


# **Block 2: Removing the Second Column**

# Define file names for the second block
output_file = 'contigs_LCA_OTU.txt'

# Open the processed output file for reading
with open(processed_output_file, 'r') as infile:
    # Open the final output file for writing
    with open(output_file, 'w') as outfile:
        # Iterate through each line in the processed output file
        for line in infile:
            # Split the line into columns using the tab delimiter
            columns = line.strip().split('\t')
            
            # Remove the second column (index 1)
            # Rejoin the remaining columns into a new line with tab separator
            new_line = '\t'.join([columns[i] for i in range(len(columns)) if i != 1])
            
            # Write the new line to the final output file
            outfile.write(new_line + '\n')

print(f"Final output saved to {output_file}")
