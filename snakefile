sequences_input_file = "sequences.fasta"
sequence_ref_file = "ref_sequence.fasta"
sequences_aln_file = "sequences_aln.fasta"
sequences_aln_cleaned_file = "sequences_aln_clean.fasta"
metadata_file = "metadata.csv"
cleaned_metadata_file = "metadata_clean.csv"
removed_duplicates_fasta = "sequences_aln_clean_removed_duplicates.fasta"
removed_duplicates_metadata = "metadata_clean_removed_duplicates.csv"
preprocessing_done = "preprocessing_done.flag"
treetime_output_dir = "treetime"
ancestral_sequences = f"{treetime_output_dir}/ancestral_sequences.fasta"
dates_input = f"{treetime_output_dir}/dates.tsv"
treetime_done_flag = f"{treetime_output_dir}/treetime_done.flag"
treetime_done = "treetime_done.flag"
ancestral_sequences_cleaned = f"{treetime_output_dir}/ancestral_cleaned.fasta"
root_sequence = f"{treetime_output_dir}/root.fasta"
cleaned_dates = f"{treetime_output_dir}/dates_cleaned.tsv"
root_dates_output = f"{treetime_output_dir}/root.tsv"
root_filled_ref = f"{treetime_output_dir}/root_filled.fasta"
ancestral_sequences_filled_root = f"{treetime_output_dir}/ancestral_sequences_filled.fasta"
dates_sorted = f"{treetime_output_dir}/dates_sorted.tsv"
sequences_sorted = f"{treetime_output_dir}/ancestral_sequences_sorted.fasta"

rule all:
    input: f"{treetime_output_dir}/ancestral_sequences.fasta"

rule mafft_alignment:
    input:
        sequences=sequences_input_file,
        ref_sequence=sequence_ref_file
    output:
        sequences_aln=sequences_aln_file
    shell:
        "python3 Run_mafft.py {input.sequences} {input.ref_sequence} {output.sequences_aln}"

rule clean_fasta_version_numbers:
    input:
        sequences_aln=sequences_aln_file
    output:
        sequences_aln_clean=sequences_aln_cleaned_file
    shell:
        "python Run_clean_fasta.py {input.sequences_aln} {output.sequences_aln_clean}"

rule clean_metadata:
    input: metadata=metadata_file
    output: cleaned_metadata=cleaned_metadata_file
    shell: "python Run_clean_metadata.py {input.metadata} {output.cleaned_metadata}"

rule remove_duplicates:
    input:
        input_fasta=sequences_aln_cleaned_file,
        input_metadata=cleaned_metadata_file
    output:
        fasta_output_removed_duplicates= removed_duplicates_fasta,
        metadata_output_removed_duplicates= removed_duplicates_metadata
    shell:
        "python Remove_duplicates.py {input.input_fasta} {input.input_metadata} {output.fasta_output_removed_duplicates} {output.metadata_output_removed_duplicates}"

rule preprocessing_done_flag:
    output: touch("preprocessing.flag")

rule run_treetime:
    input:
        processed_fasta= removed_duplicates_fasta,
        processed_metadata= removed_duplicates_metadata,
        preprocessing_done_flag= "preprocessing.flag"
    output:
        ancestral_sequences=ancestral_sequences,
        dates=dates_input
    shell:
        "treetime --aln {input.processed_fasta} --dates {input.processed_metadata} --stochastic-resolve --outdir {treetime_output_dir}"

rule treetime_done_flag:
    output: touch(f"{treetime_output_dir}/treetime_done.flag")

rule export_root_remove_nodes:
    input:
        ancestral_sequences_input=ancestral_sequences
    output:
        root_seq=root_sequence,
        ancestral_cleaned=ancestral_sequences_cleaned
    shell:
        "python Extract_root_remove_nodes.py {input.ancestral_sequences_input} {output.ancestral_cleaned} {output.root_seq}"

rule dates_remove_nodes:
    input:
        dates=dates_input
    output:
        root_dates_out= root_dates_output,
        cleaned_dates_out= cleaned_dates
    shell:
        "python Extract_root_remove_nodes_tsv.py {input.dates} {output.cleaned_dates_out} {output.root_dates_out}"

rule replace_insertions_root:
    input:
        root_file=root_sequence,
        reference_sequence=sequence_ref_file
    output:
        root_file_out=root_filled_ref
    shell:
        "python Replace_insertions_root.py {input.root_file} {input.reference_sequence} {output.root_file_out}"

rule root_insertions_replaced_flag:
    output:
        touch(f"{treetime_output_dir}/root_insertions_replaced.flag")

rule replace_insertions_sequences:
    input:
        ancestral_sequences_in=ancestral_sequences_cleaned ,
        ref_root_sequence= root_filled_ref ,
        root_insertions_check=f"{treetime_output_dir}/root_insertions_replaced.flag"
    output:
        ancestral_sequences_filled=ancestral_sequences_filled_root
    shell:
        "python Replace_insertions.py {input.ancestral_sequences_in} {input.ref_root_sequence} {output.ancestral_sequences_filled}"
rule ancestral_sequences_replaced_flag:
    output:
        touch(f"{treetime_output_dir}/ancestral_insertions_replaced.flag")

rule check_sort_dates:
    input:
        ancestral_inertions_check=f"{treetime_output_dir}/ancestral_insertions_replaced.flag",
        ancestral_sequences_new=ancestral_sequences_filled_root,
        dates_cleaned=cleaned_dates
    output:
        dates_sorted=dates_sorted,
        ancestral_seq_sorted=sequences_sorted
    shell:
        "python Run_sort_dates_2.py  {input.dates_cleaned} {input.ancestral_sequences_new} {output.dates_sorted} {output.ancestral_seq_sorted}"



ruleorder: mafft_alignment > clean_fasta_version_numbers > clean_metadata > remove_duplicates > preprocessing_done_flag > run_treetime > treetime_done_flag > export_root_remove_nodes > dates_remove_nodes > replace_insertions_root > root_insertions_replaced_flag > replace_insertions_sequences > ancestral_sequences_replaced_flag > check_sort_dates

