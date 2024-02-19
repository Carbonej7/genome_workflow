# Set of Codes to preprocess data for c/u genomic analysis


Python enviorment requirements in requirements.txt
Copy Initial Data:  <virus_name>_sequences.fasta,
                    <virus_name_reference.fasta,
                    <virus_name>_dates.csv
      To <Run_genome__processing/_internal/data/
      Change the names:
                    <virus_name>_sequences.fasta > sequences.fasta,
                    <virus_name_reference.fasta >  ref_sequence.fasta
                    <virus_name>_dates.csv > metadata.csv
                  
To Run:
  "python Run_genome_processing.py <target_rule>"
