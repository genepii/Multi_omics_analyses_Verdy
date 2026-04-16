#### Metagenomics contigs building, dereplication and taxonomic assignment ###

#   1) denovo contigs with spades

mkdir denovo
mkdir denovo/spades

while IFS= read -r sampleId; do
    singularity exec seqmet/singularity/denovo.sif \
    spades.py -t 100 -1 concatenated_reads/${sampleId}_R1.fastq -2 concatenated_reads/${sampleId}_R2.fastq -o "denovo/spades/${sampleId}/" --meta
done < "/sampleID.txt"

#   2) copy all contigs file into one file
mkdir denovo/contigs
input_fasta="/denovo/spades"
while IFS= read -r sampleId; do
cp "${input_fasta}/${sampleId}/contigs.fasta" "denovo/contigs/${sampleId}_contigs.fasta"
done < "sampleID.txt"

#   3) dereplicate all contigs
mkdir denovo/cdhit
cat denovo/contigs/*.fasta > denovo/cdhit/allcontigs.fasta
singularity exec /seqmet/singularity/denovo.sif \
reformat.sh minlength=$length_threshold in=denovo/cdhit/allcontigs.fasta out=denovo/cdhit/allcontigs_filtered.fasta
cd-hit-est -i denovo/cdhit/allcontigs_filtered.fasta -o denovo/cdhit/OTUs.fasta -c 0.95 -G 0 -aS 0.95 -g 1 -r 1 -M 0 -d 0
awk '/^>/{$0=">OTU_"++i}1' denovo/cdhit/OTUs.fasta > denovo/cdhit/OTUs.fasta

#   4) assign dereplicated contigs to diamond (NR database from NCBI)

mkdir diamond

db="db/refseq_protein_nonredund_diamond"

./diamond blastx --query /denovo/cdhit/OTUs.fasta -p 114 -e 0.001 -k 1 --max-hsps 1 --db $db --outfmt 102 --out denovo/cdhit/contigs_LCA_OTU.txt

tax2lin="seqmet/db/ncbitax2lin/ncbitax2lin-221121_hosttax_sorted.tsv"
cat /denovo/cdhit/contigs_LCA_OTU.txt | awk -F'\t' '{print $2 "\t" $1}' > /denovo/cdhit/contigs_LCA_OTU_sorted.txt
sort /denovo/cdhit/contigs_LCA_OTU_sorted.txt >> /denovo/cdhit/contigs_LCA_OTU_sorted2.txt
join -a1 -o 2.2,1.2 -t$'\t' /denovo/cdhit/contigs_LCA_OTU_sorted2.txt $tax2lin | grep "Viruses" > /denovo/cdhit/viral_contigs.tsv
cat denovo/cdhit/OTUs.fasta |awk '/^>/ {if(N>0) printf("\n"); printf("%s ",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' > denovo/cdhit/vOTUs.fasta
cat diamond/viral_contigs.tsv | awk -F'\t' '{print $2}'|uniq > viral_contigs_list.txt
grep -f viral_contigs_list.txt vOTUs.fasta | sed 's/ /\n/'g > vOTUs.fasta
