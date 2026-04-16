#### Metatranscriptomic reads processing and filtering ###

#1 Interleave

    // Convert raw paired-end fastq to an interleaved fastq

    reformat.sh -Xmx${memory}G threads=$task.cpus in1="${R1}" in2="${R2}" out="interleaved/${sampleId}.fastq.gz" qin=auto qout=33
    
#2 Dehosting

	// Remove huamn reads from the interleaved fastq
	
	db = "/opt/scrubber/data/human_filter.db"
	 
    unpigz -p $task.cpus -c "${sampleId}.fastq.gz" | /opt/scrubber/scripts/scrub.sh -p $task.cpus -x -d ["db"] -o - -i - | repair.sh -Xmx${memory}G threads=$task.cpus in=stdin out=stdout | pigz -6 -p $task.cpus -c > "dehosted/${sampleId}.fastq.gz"
	
#3 QC
	
	// QC for dehosted reads
	
    fastqc "work/${sampleId}_R1.fastq.gz" "work/${sampleId}_R2.fastq.gz" --outdir fastqc -t $task.cpus
    unzip -d work "fastqc/${sampleId}_R1_fastqc.zip"
    unzip -d work "fastqc/${sampleId}_R2_fastqc.zip"
    cat \$(find work/ -type f -name fastqc_data.txt) >> "work/${sampleId}.txt"
    echo -e "RAW\\tfastq_readcount\\t\$(grep "Total Sequences" "work/${sampleId}.txt" | cut -f 2 | paste -sd+ | bc)\\t${sampleId}" > "fastqc/${sampleId}.tsv"

#4 Trimming

    // Trim adapters from each read, drop any read pair if at least one read is shorter than a specified length, multiple adapters can be provided

	--> params.trim_cutadapt : 
	"trim_cutadapt" : {
    "todo" : 1,
    "minimum_length" : 30,
    "error_rate" : 0,
    "quality" : 20,
    "overlap" : 1,
    "adapter_f" : [],
    "adapter_r" : [],
    "front_f" : ["ACTTTGTGTTTGA","TGCTCTTCCGATC"],
    "front_r" : ["ACTTTGTGTTTGA","TGCTCTTCCGATC"],
    "anywhere_f" : [],
    "anywhere_r" : []
}

    if [[ ! -s "${sampleId}.fastq.gz" ]]; then
        touch "trimmed/${sampleId}.fastq.gz" "trimmed/${sampleId}.txt"
        
        exit 0
    fi
    
    if [[ "${params.trim_cutadapt["adapter_f"]}" != [] || "${params.trim_cutadapt["adapter_r"]}" != [] ]]; then
        adapter_f=\$(for item in \$(echo ${params.trim_cutadapt["adapter_f"]} | tr -d '[\\[\\],\\"]'); do echo -n "-a \${item} "; done)
        adapter_r=\$(for item in \$(echo ${params.trim_cutadapt["adapter_r"]} | tr -d '[\\[\\],\\"]'); do echo -n "-A \${item} "; done)
    elif [[ "${params.readtype}" == "paired-end" ]]; then
        BBMERGE_ADAPTER_F="\$(bbmerge.sh in="${sampleId}.fastq.gz" outa=stdout | grep -A 1 'Read1' | grep -v '>' | grep -v '^--\$' | perl -pe 'chomp if eof' | tr '\\n' ',')"
        BBMERGE_ADAPTER_R="\$(bbmerge.sh in="${sampleId}.fastq.gz" outa=stdout | grep -A 1 'Read2' | grep -v '>' | grep -v '^--\$' | perl -pe 'chomp if eof' | tr '\\n' ',')"
        if [[ \${BBMERGE_ADAPTER_F} == "N" ]]; then BBMERGE_ADAPTER_F=""; fi
        if [[ \${BBMERGE_ADAPTER_R} == "N" ]]; then BBMERGE_ADAPTER_R=""; fi
        adapter_f=\$(for item in \$(echo \${BBMERGE_ADAPTER_F} | tr -d '[\\[\\],\\"]'); do echo -n "-a \${item} "; done)
        adapter_r=\$(for item in \$(echo \${BBMERGE_ADAPTER_R} | tr -d '[\\[\\],\\"]'); do echo -n "-A \${item} "; done)
    else
        adapter_f=""
        adapter_r=""
    fi

    front_f=\$(for item in \$(echo ${params.trim_cutadapt["front_f"]} | tr -d '[\\[\\],\\"]'); do echo -n "-g \${item} "; done)
    front_r=\$(for item in \$(echo ${params.trim_cutadapt["front_r"]} | tr -d '[\\[\\],\\"]'); do echo -n "-G \${item} "; done)
    anywhere_f=\$(for item in \$(echo ${params.trim_cutadapt["anywhere_f"]} | tr -d '[\\[\\],\\"]'); do echo -n "-b \${item} "; done)
    anywhere_r=\$(for item in \$(echo ${params.trim_cutadapt["anywhere_r"]} | tr -d '[\\[\\],\\"]'); do echo -n "-B \${item} "; done)

    if [[ ${params.readtype} == "paired-end" ]]; then
        ARGS="--interleaved --pair-filter=any \${adapter_f}\${adapter_r}\${front_f}\${front_r}\${anywhere_f}\${anywhere_r}"
    else
        ARGS="\${adapter_f}\${front_f}\${anywhere_f}"
    fi

    echo "ARGS : \$ARGS"
    echo "BBMERGE_ADAPTER_F : \$BBMERGE_ADAPTER_F"
    echo "BBMERGE_ADAPTER_R : \$BBMERGE_ADAPTER_R"

    cutadapt \${ARGS} --cores $task.cpus --minimum-length ${params.trim_cutadapt["minimum_length"]} --error-rate ${params.trim_cutadapt["error_rate"]} --nextseq-trim ${params.trim_cutadapt["quality"]} --trim-n --overlap ${params.trim_cutadapt["overlap"]} -o "trimmed/${sampleId}.fastq" "${sampleId}.fastq.gz" > "trimmed/${sampleId}.txt"
    pigz -p $task.cpus -6 "trimmed/${sampleId}.fastq"

    rm -rf work
	
#5 Kraken2 for rapid assignation and Krakentools to extract remained human and bacterial reads

##  a) print samples into a txt file
mkdir ${workfile}
cd ${workfile}
samples=$(
    for i in *; do echo $i;done)
printf '%s\n' "$samples"|awk -F'.' '{print $1}' > "sampleID.txt"

#   b) reformat trimmed reads into R1 and R2 reads
cd ${workfile}
mkdir ${workfile}/reformat 
while IFS= read -r sampleId; do
singularity exec seqmet/singularity/denovo.sif reformat.sh threads=60 in="${sampleId}.fastq" out1="${sampleId}_R1.fastq" out2="${sampleId}_R2.fastq" int=t qin=auto qout=33
done < "sampleID.txt"

#   c) kraken2 for trimmed reads
mkdir krak2
mkdir krak2/un

input_directory="trimmed_reads_R1/R2.fastq"
kraken2_db="/kraken2/"
while IFS= read -r sample_id; do
    # Run the Kraken2 command for each sample ID
    singularity exec seqmet/singularity/denovo.sif \
    kraken2 --db "$kraken2_db" \
            --confidence 0.1 \
            --minimum-hit-groups 2 \
            --classified-out "krak2/${sample_id}#.fastq.gz" \
            --output "krak2/${sample_id}.out" \
            --report-minimizer-data \
            --report "krak2/${sample_id}.tsv" \
            --unclassified-out "krak2/un/${sample_id}#_unclassified.fastq.gz" \
            --paired "${input_directory}/${sample_id}_R1.fastq.gz" "${input_directory}/${sample_id}_R2.fastq.gz" \
            --threads 64
done < "sampleID.txt"

#   d) krakentools to filter out remaining humand, bacterial, fungal or archaeal reads

mkdir krak2/krak_filter_human_bacteria

while IFS= read -r sample_id; do
mv "krak2/${sample_id}_1.fastq" "krak2/${sample_id}_R1.fastq"
mv "krak2/${sample_id}_2.fastq" "krak2/${sample_id}_R2.fastq"
done < "sampleID.txt"

while IFS= read -r sample_id; do
mv "krak2/un/${sample_id}_1.fastq" "krak2/un/${sample_id}_R1.fastq"
mv "krak2/un/${sample_id}_2.fastq" "krak2/un/${sample_id}_R2.fastq"
done < "sampleID.txt"

# Loop over each sample ID
while IFS= read -r sample_id; do
    # Set the input and output file paths
    input_out="krak2/${sample_id}.out"
    input_fastq="krak2/${sample_id}_R1.fastq"
    output_fastq="krak2/krak_filter_human_bacteria/${sample_id}_R1.fastq"
    input_tsv="krak2/${sample_id}.tsv"

    # Run the extraction script with the appropriate arguments
    singularity exec /ananihu/krakentools.sif \
    extract_kraken_reads.py --max 1000000000 -t 2 2759 4751 2157 -k "$input_out" -s "$input_fastq" -o "$output_fastq" -r "$input_tsv" --include-children --exclude --fastq-output
    echo "Extraction complete for sample: ${sample_id}"
done < "sampleID.txt"

while IFS= read -r sample_id; do
    # Set the input and output file paths
    input_out="krak2/${sample_id}.out"
    input_fastq="krak2/${sample_id}_R2.fastq"
    output_fastq="krak2/krak_filter_human_bacteria/${sample_id}_R2.fastq"
    input_tsv="krak2/${sample_id}.tsv"

    # Run the extraction script with the appropriate arguments
    singularity exec ananihu/krakentools.sif \
    extract_kraken_reads.py --max 1000000000 -t 2 2759 4751 2157 -k "$input_out" -s "$input_fastq" -o "$output_fastq" -r "$input_tsv" --include-children --exclude --fastq-output
    echo "Extraction complete for sample: ${sample_id}"
done < "sampleID.txt"


#6 Concatenate classified and unclassified reads for assembly

mkdir concatenated_reads

while IFS= read -r sampleId; do
cat "krak2/krak_filter_human_bacteria/${sampleId}_R1.fastq" "krak2/un/${sampleId}_R1.fastq" > "concatenated_reads/${sampleId}_R1.fastq"
done < "sampleID.txt"

while IFS= read -r sampleId; do
cat "krak2/krak_filter_human_bacteria/${sampleId}_R2.fastq" "krak2/un/${sampleId}_R2.fastq" > "concatenated_reads/${sampleId}_R2.fastq"
done < "sampleID.txt"
