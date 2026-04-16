#### Remapping concatenated reads on vOTUs ### 

(we followed the mapping protocol from: Kaelin et al., Nat Microbiol. 2022 May;7(5):653-662)

# 1) remap concatenated reads on vOTUs

mkdir mapping
mkdir count

input_fastq="filtered_reads"

while IFS= read -r sampleId; do
singularity exec /seqmet/singularity/denovo.sif \
bwa index vOTUs.fasta
singularity exec /seqmet/singularity/denovo.sif \
bwa mem -t 100 -M -L 97,97 vOTUs.fasta ${input_fastq}/${sampleId}_R1.fastq ${input_fastq}/${sampleId}_R2.fastq > mapping/${sampleId}-readsMapped.sam;
singularity exec /seqmet/singularity/denovo.sif \
samtools view -@ 30 -h -F 0x900 mapping/${sampleId}-readsMapped.sam > mapping/${sampleId}-secondaryRemoved.sam;
singularity exec /seqmet/singularity/denovo.sif \
samtools view -@ 30 -h -F 0x4 mapping/${sampleId}-secondaryRemoved.sam > mapping/${sampleId}-secondaryUnMappedRemoved.sam;
singularity exec /seqmet/singularity/denovo.sif \
samtools view -@ 30 -S -b mapping/${sampleId}-secondaryUnMappedRemoved.sam > mapping/${sampleId}-secondaryUnMappedRemoved.bam;
singularity exec /seqmet/singularity/denovo.sif \
samtools sort -@ 30 mapping/${sampleId}-secondaryUnMappedRemoved.bam > mapping/${sampleId}-secondaryUnMappedRemoved_sorted.bam;
singularity exec /seqmet/singularity/denovo.sif \
samtools index -@ 30 mapping/${sampleId}-secondaryUnMappedRemoved_sorted.bam;
singularity exec /seqmet/singularity/denovo.sif \
samtools idxstats -@ 30 mapping/${sampleId}-secondaryUnMappedRemoved_sorted.bam > count/${sampleId}-counts.txt;
done < "sampleID.txt"


# 2) remap concatenated reads on the internal control

###MS2 : NC_001417.2 phage MS2 genome

mkdir mapping_MS2
mkdir count_ms2

input_fastq="filtered_reads"

while IFS= read -r sampleId; do
singularity exec /seqmet/singularity/denovo.sif \
bwa index MS2.fasta
singularity exec /seqmet/singularity/denovo.sif \
bwa mem -t 40 -M -L 97,97 MS2.fasta ${input_fastq}/${sampleId}_R1.fastq ${input_fastq}/${sampleId}_R2.fastq > mapping_MS2/${sampleId}-readsMapped.sam;
singularity exec /seqmet/singularity/denovo.sif \
samtools view -@ 30 -h -F 0x900 mapping_MS2/${sampleId}-readsMapped.sam > mapping_MS2/${sampleId}-secondaryRemoved.sam;
singularity exec /seqmet/singularity/denovo.sif \
samtools view -@ 30 -h -F 0x4 mapping_MS2/${sampleId}-secondaryRemoved.sam > mapping_MS2/${sampleId}-secondaryUnMappedRemoved.sam;
singularity exec /seqmet/singularity/denovo.sif \
samtools view -@ 30 -S -b mapping_MS2/${sampleId}-secondaryUnMappedRemoved.sam > mapping_MS2/${sampleId}-secondaryUnMappedRemoved.bam;
singularity exec /seqmet/singularity/denovo.sif \
samtools sort -@ 30 mapping_MS2/${sampleId}-secondaryUnMappedRemoved.bam > mapping_MS2/${sampleId}-secondaryUnMappedRemoved_sorted.bam;
singularity exec /seqmet/singularity/denovo.sif \
samtools index -@ 30 mapping_MS2/${sampleId}-secondaryUnMappedRemoved_sorted.bam;
singularity exec /seqmet/singularity/denovo.sif \
samtools idxstats -@ 30 mapping_MS2/${sampleId}-secondaryUnMappedRemoved_sorted.bam > count_ms2/${sampleId}-counts.txt;
done < "sampleID.txt"