#!/bin/bash

# Paths to input files
diamond_output="contigs_LCA_OTU.txt"
viral_list="viral_list.txt"
filtered_output="contigs_LCA_OTU_corrected.txt"

# Read the viral list into an array
mapfile -t viral_taxids < "$viral_list"

# Create an associative array for quick lookup
declare -A viral_taxid_map
for taxid in "${viral_taxids[@]}"; do 
    viral_taxid_map["$taxid"]=1
done

# Process the diamond output
while IFS=$'\t' read -r qseqid sseqid pident length mismatch gapopen evalue bitscore staxid qcovhsp; do
    # Split the staxid field into an array using ';' as delimiter
    IFS=';' read -ra taxid_array <<< "$staxid"
    
    # Initialize an array to hold selected viral taxids
    selected_taxids=()

    # Check if any of the taxids are in the viral list and collect them
    for taxid in "${taxid_array[@]}"; do 
        if [[ -n "${viral_taxid_map[$taxid]}" ]]; then 
            selected_taxids+=("$taxid")
        fi 
    done

    # If there are no viral taxids found, default to the first taxid in the original list
    if [ ${#selected_taxids[@]} -eq 0 ]; then 
        selected_taxids=("${taxid_array[0]}")
    fi

    # Join the selected viral taxids back into a string separated by ';'
    selected_taxids_str=$(IFS=';'; echo "${selected_taxids[*]}")

    # Print the filtered line with the selected taxids
    echo -e "$qseqid\t$sseqid\t$pident\t$length\t$mismatch\t$gapopen\t$evalue\t$bitscore\t$selected_taxids_str\t$qcovhsp" 
done < "$diamond_output" > "$filtered_output"

echo "Filtered output written to $filtered_output"
