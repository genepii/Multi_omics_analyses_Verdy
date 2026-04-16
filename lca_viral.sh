#!/bin/bash

# Input files
diamond_output="contigs_LCA_OTU_corrected.txt"
viral_list="viral_list.txt"
viral_taxonomy="viral_taxonomy.txt"  # This table should map taxid to full lineage

# Output file
output_with_lca="contigs_LCA_OTU_corrected_final.txt"

# Step 1: Prepare viral taxid list and taxonomy information
mapfile -t viral_taxids < "$viral_list"

# Create an associative array for quick viral taxid lookup
declare -A viral_taxid_map
for taxid in "${viral_taxids[@]}"; do 
    viral_taxid_map["$taxid"]=1
done

# Step 2: Load taxonomy information into an associative array
declare -A taxid_to_lineage
while IFS=$'\t' read -r taxid lineage; do 
    taxid_to_lineage["$taxid"]="$lineage"
done < "$viral_taxonomy"

# Step 3: Process the diamond output to compute LCA
while IFS=$'\t' read -r qseqid sseqid pident length mismatch gapopen evalue bitscore staxids qcovhsp; do 

    # Split the staxids field into an array using ';' as a delimiter
    IFS=';' read -ra taxid_array <<< "$staxids"

    selected_taxids=()
    
    # Collect viral taxids that are present in the list
    for taxid in "${taxid_array[@]}"; do 
        if [[ -n "${viral_taxid_map[$taxid]}" ]]; then 
            selected_taxids+=("$taxid")
        fi 
    done

    # If multiple viral taxids are found, compute LCA
    if [[ ${#selected_taxids[@]} -gt 1 ]]; then
        # Get the lineages for the selected taxids
        lineages=() 
        for taxid in "${selected_taxids[@]}"; do 
            lineages+=("${taxid_to_lineage[$taxid]}")
        done

        # Split the lineages into arrays by their taxonomic levels
        IFS=';' read -ra first_lineage <<< "${lineages[0]}" 
        lca=("${first_lineage[@]}")
        
        for lineage in "${lineages[@]:1}"; do 
            IFS=';' read -ra current_lineage <<< "$lineage"
            # Compare each level of the taxonomy
            for i in "${!lca[@]}"; do 
                if [[ "${lca[$i]}" != "${current_lineage[$i]}" ]]; then 
                    lca=("${lca[@]:0:$i}")
                    break
                fi 
            done 
        done
    else
        IFS=';' read -ra lca <<< "${taxid_to_lineage[${selected_taxids[0]}]}"
    fi

    # Ensure that the LCA has exactly eight ranks
    while [ ${#lca[@]} -lt 8 ]; do 
        lca+=("__")
    done

    # If LCA has more than eight ranks, truncate it to the first eight
    lca=("${lca[@]:0:8}")

    # Join the LCA taxonomy back into a string with "; " as a separator
    lca=$(IFS='; '; echo "${lca[*]}")

    # Replace occurrences of ";__" with "; __"
    lca=$(echo "$lca" | sed 's/;__/; __/g')

    # Prepare the output line
    lca_taxid=$(IFS=';' echo "${selected_taxids[*]}")
    echo -e "$qseqid\t$sseqid\t$pident\t$length\t$mismatch\t$gapopen\t$evalue\t$bitscore\t$staxids\t$qcovhsp\t$lca\t$lca_taxid" 

done < "$diamond_output" > "$output_with_lca"

echo "LCA analysis complete. Results are in $output_with_lca"
