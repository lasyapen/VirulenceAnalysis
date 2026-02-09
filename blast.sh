ml BLAST+

makeblastdb -in VFDB_setA_nt.fas -dbtype nucl -out vfdb_nucleotide_db

QUERY="/scratch/lp27068/Olivia1_Data_Reads/virulence_analysis/concatenated_query.fasta"
DB="vfdb_nucleotide_db"
OUT_TAB="vfdb_blast_results_tabular.txt"
OUT_STD="vfdb_blast_results_standard.txt"

blastn -query $QUERY -db $DB -evalue 1e-5 -num_threads 4 -out $OUT_TAB -outfmt 6
blastn -query $QUERY -db $DB -evalue 1e-5 -num_threads 4 -out $OUT_STD -outfmt 0

#convert vfdb_blast_results_tabular.tx into a file called "annotated_blast_result.txt" that contains the full MAG name as column 1 (including the k141 part of the header)


