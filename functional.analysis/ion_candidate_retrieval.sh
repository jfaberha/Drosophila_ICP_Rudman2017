### This was the first pass, but it was discovered the allele frequency files were generated using dm5.39
### instead of dm6, so use 'candidate_retrieval_dm5_fixed.sh' after this script to get the properly
### converted dm5 coordinates.

#Pull candidate genes from gff file to develop a bed file
sed 's/,/\t/g' DrosICP_CandidateGenes_190619.csv | sed 1d > DrosICP_CandidateGenes_190619.tsv

awk -F'\t' '{print "gene=" $1 ";"}' DrosICP_CandidateGenes_190619.tsv | uniq > geneNames.txt
filename='geneNames.txt'

while read i; do
	grep -i $i ./GCF_000001215.4/genomic.gff | awk '$3=="gene"' > temp
	echo "$i" | paste - temp >> temp2
done < $filename

mv temp2 ion_candidates.gff
rm temp

#convert gff to bed file and swap refseq names for chrom number
awk '{print $2 "\t" $5 "\t" $6 "\t" $1}' ion_candidates.gff > ion_candidates.bed

sed -i'' -e 's/NC_004354.4/X/g' ion_candidates.bed
sed -i'' -e 's/NT_033779.5/2L/g' ion_candidates.bed
sed -i'' -e 's/NT_033778.4/2R/g' ion_candidates.bed
sed -i'' -e 's/NT_037436.4/3L/g' ion_candidates.bed
sed -i'' -e 's/NT_033777.3/3R/g' ion_candidates.bed
sed -i'' -e 's/NC_004353.4/4/g' ion_candidates.bed
sed -i'' -e 's/NC_024512.1/Y/g' ion_candidates.bed
sed -i'' -e 's/NC_024511.2/MT/g' ion_candidates.bed

#add 500bp boundary regions to capture snps in immediate upstream/downstream of genes
awk '$2<$3' ion_candidates.bed  | awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4 "\t+"}'> temp
awk '$3<$2' ion_candidates.bed | awk '{print $1 "\t" $3-500 "\t" $2+500 "\t" $4 "\t-"}'>> temp
#awk '$2<$3' ion_candidates.bed  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t+"}'> temp
#awk '$3<$2' ion_candidates.bed | awk '{print $1 "\t" $3 "\t" $2 "\t" $4 "\t-"}'>> temp
sort -k1,1 -k2,2 -k3,3 temp > temp2
echo "chrom	start	stop	name	strand" | cat - temp2 > ion_candidates_sorted.bed
rm temp*

#sites.17.bed exported from r and was the sites file converted to a bed file
#by setting position as the start and position +1 as the stop
bedtools intersect -a sites.17.bed -b ion_candidates_sorted.bed -wao > ion_candidate_snp_overlap.txt

#remove any row duplicates and add headers
awk '{print $0 "\t" $1 "\t" $2}' ion_candidate_snp_overlap.txt | uniq -f 11 | cut -f1-10 > temp
awk '{print $0 "\t" $1 "\t" $2}' ion_candidate_snp_overlap.txt | uniq -c -f 11 | awk '{print $1}' > temp2
paste temp temp2 > temp3
paste sites.17.bed ion_candidates_sorted.bed | head -1 | awk '{print $0 "\toverlapcount"}' | cat - temp3 > ion_candidate_snp_overlap_nr.txt
rm temp*