data_dir = "data"
sequences_input_file = f"{data_dir}/sequences.fasta"
sequence_ref_file = f"{data_dir}/ref_sequence.fasta"
sequences_aln_file = f"{data_dir}/sequences_aln.fasta"
sequences_aln_cleaned_file = f"{data_dir}/sequences_aln_clean.fasta"
metadata_file = f"{data_dir}/metadata.csv"
cleaned_metadata_file = f"{data_dir}/metadata_clean.csv"
removed_duplicates_fasta = f"{data_dir}/sequences_aln_clean_removed_duplicates.fasta"
removed_duplicates_metadata = f"{data_dir}/metadata_clean_removed_duplicates.csv"
preprocessing_done = "preprocessing_done.flag"
treetime_output_dir = "treetime"
ancestral_sequences = f"{treetime_output_dir}/ancestral_sequences.fasta"
dates_input = f"{treetime_output_dir}/dates.tsv"
treetime_done_flag = f"{treetime_output_dir}/treetime_done.flag"
treetime_done = "treetime_done.flag"
ancestral_sequences_cleaned = f"{treetime_output_dir}/ancestral_cleaned.fasta"
outliers = f"{treetime_output_dir}/outliers.tsv"
root_sequence = f"{treetime_output_dir}/root.fasta"
cleaned_dates = f"{treetime_output_dir}/dates_cleaned.tsv"
root_filled_ref = f"{treetime_output_dir}/root_filled.fasta"
ancestral_sequences_filled_root = f"{treetime_output_dir}/ancestral_sequences_filled.fasta"
dates_sorted = f"{treetime_output_dir}/dates_sorted.tsv"
sequences_sorted = f"{treetime_output_dir}/ancestral_sequences_sorted.fasta"
nodes_add_dates_sorted = f"{treetime_output_dir}/nodes_add_dates_sorted.tsv"
dates_csv = f"{treetime_output_dir}/nodes_add_dates_sorted.csv"
nodes_add_fasta_sorted = f"{treetime_output_dir}/nodes_add_dates_sorted.fasta"
processing_done_output_dir = "processed"
processed_fasta = f"{processing_done_output_dir}/processed.fasta"
processed_csv = f"{processing_done_output_dir}/processed.csv"
average_mutations = f"{processing_done_output_dir}/average_mutations.csv"
genome_c_u = f"{processing_done_output_dir}/genome_c_u_all.csv"
genome_c_u_avg = f"{processing_done_output_dir}/genome_c_u_all_avg.csv"
genome_mut_L_all = f"{processing_done_output_dir}/genome_mut_L_all.csv"



ruleorder: mafft_alignment > clean_fasta_version_numbers > clean_metadata > remove_duplicates > preprocessing_done_flag > run_treetime > treetime_done_flag > export_root_remove_nodes > dates_remove_nodes > replace_insertions_root > root_insertions_replaced_flag > replace_insertions_sequences > ancestral_sequences_replaced_flag > check_sort_dates > sort_dates_flag > add_nodes_back > convert_to_csv > nodes_added_csv_sorted_flag > add_nodes_fasta > check_lengths > count_mutations > L_dist


rule mafft_alignment:
    input:
        sequences=sequences_input_file,
        ref_sequence=sequence_ref_file
    output:
        sequences_aln=sequences_aln_file
    group: 1
    shell:
        "python3 scripts/Run_mafft.py {input.sequences} {input.ref_sequence} {output.sequences_aln}"

rule clean_fasta_version_numbers:
    input:
        sequences_aln=sequences_aln_file
    output:
        sequences_aln_clean=sequences_aln_cleaned_file
    group: 1
    shell:
        "python scripts/Run_clean_fasta.py {input.sequences_aln} {output.sequences_aln_clean}"

rule clean_metadata:
    input: metadata=metadata_file
    output: cleaned_metadata=cleaned_metadata_file
    group: 1
    shell: "python scripts/Run_clean_metadata.py {input.metadata} {output.cleaned_metadata}"

rule remove_duplicates:
    input:
        input_fasta=sequences_aln_cleaned_file,
        input_metadata=cleaned_metadata_file
    output:
        fasta_output_removed_duplicates= removed_duplicates_fasta,
        metadata_output_removed_duplicates= removed_duplicates_metadata
    group: 1
    shell:
        "python scripts/Remove_duplicates.py {input.input_fasta} {input.input_metadata} {output.fasta_output_removed_duplicates} {output.metadata_output_removed_duplicates}"

rule preprocessing_done_flag:
    output: touch("preprocessing.flag")
    group: 1

rule run_treetime:
    input:
        processed_fasta= removed_duplicates_fasta,
        processed_metadata= removed_duplicates_metadata,
        preprocessing_done_flag= "preprocessing.flag"
    output:
        ancestral_sequences=ancestral_sequences,
        dates=dates_input,
        outliers_dates=outliers
    group: 2
    shell:
        "treetime --aln {input.processed_fasta} --dates {input.processed_metadata} --stochastic-resolve --outdir {treetime_output_dir}"

rule treetime_done_flag:
    output: touch(f"{treetime_output_dir}/treetime_done.flag")
    group: 2


rule export_root_remove_nodes:
    input:
        ancestral_sequences_input=ancestral_sequences,
        treetime_flag = f"{treetime_output_dir}/treetime_done.flag"
    output:
        root_seq=root_sequence,
        ancestral_cleaned=ancestral_sequences_cleaned
    group: 3
    shell:
        "python scripts/Extract_root_remove_nodes.py {input.ancestral_sequences_input} {output.ancestral_cleaned} {output.root_seq}"

rule dates_remove_nodes:
    input:
        dates=dates_input
    output:
        cleaned_dates_output= cleaned_dates
    group: 3
    shell:
        "python scripts/Extract_root_remove_nodes_tsv.py {input.dates} {output.cleaned_dates_output}"

rule replace_insertions_root:
    input:
        root_file=root_sequence,
        reference_sequence=sequence_ref_file
    output:
        root_file_out=root_filled_ref
    group: 3
    shell:
        "python scripts/Replace_insertions_root.py {input.root_file} {input.reference_sequence} {output.root_file_out}"

rule root_insertions_replaced_flag:
    output:
        touch(f"{treetime_output_dir}/root_insertions_replaced.flag")
    group: 3

rule replace_insertions_sequences:
    input:
        ancestral_sequences_in=ancestral_sequences_cleaned ,
        ref_root_sequence= root_filled_ref ,
        root_insertions_check=f"{treetime_output_dir}/root_insertions_replaced.flag"
    output:
        ancestral_sequences_filled=ancestral_sequences_filled_root
    group: 4
    shell:
        "python scripts/Replace_insertions.py {input.ancestral_sequences_in} {input.ref_root_sequence} {output.ancestral_sequences_filled}"
rule ancestral_sequences_replaced_flag:
    output:
        touch(f"{treetime_output_dir}/ancestral_insertions_replaced.flag")
    group: 4

rule check_sort_dates:
    input:
        ancestral_inertions_check=f"{treetime_output_dir}/ancestral_insertions_replaced.flag",
        ancestral_sequences_new=ancestral_sequences_filled_root,
        dates_cleaned=cleaned_dates
    output:
        dates_sorted=dates_sorted,
        ancestral_seq_sorted=sequences_sorted
    group: 5
    shell:
        "python scripts/Run_sort_dates_2.py  {input.dates_cleaned} {input.ancestral_sequences_new} {output.dates_sorted} {output.ancestral_seq_sorted}"

rule sort_dates_flag:
    output:
        touch(f"{treetime_output_dir}/sort_dates.flag")
    group: 5

rule add_nodes_back:
    input:
        sort_dates_flag = f"{treetime_output_dir}/sort_dates.flag",
        reference_dates = dates_input,
        dates_to_add_sorted = dates_sorted,
        outlier_dates = outliers
    output:
        filtered_nodes_dates = nodes_add_dates_sorted
    group: 6
    shell:
        "python scripts/Add_nodes_tsv.py {input.reference_dates} {input.dates_to_add_sorted} {input.outlier_dates} {output.filtered_nodes_dates}"

rule convert_to_csv:
    input:
        input_tsv = nodes_add_dates_sorted
    output:
        output_csv = dates_csv
    group: 6
    shell:
        "python scripts/tsv_to_csv.py {input.input_tsv} {output.output_csv}"

rule nodes_added_csv_sorted_flag:
    group: 6
    output:
        touch(f"{treetime_output_dir}/nodes_added_csv_sorted.flag")

rule add_nodes_fasta:
    input:
        csv_done_flag= f"{treetime_output_dir}/nodes_added_csv_sorted.flag",
        ref_sequence= ancestral_sequences,
        input_sequences= sequences_sorted,
        dates= dates_csv
    output:
        sorted_fasta=nodes_add_fasta_sorted
    group: 7
    shell:
        "python scripts/Sort_dates_3.py {input.ref_sequence} {input.input_sequences} {input.dates} {output.sorted_fasta}"
rule check_lengths:
    input:
        sorted_fasta= nodes_add_fasta_sorted,
        dates= dates_csv
    output:
        fasta_output=processed_fasta,
        dates_output=processed_csv
    group: 8
    shell:
        "python scripts/Check_sequence_length.py {input.sorted_fasta} {input.dates} {output.fasta_output} {output.dates_output}"

rule count_mutations:
    input:
        genome_sequences= processed_fasta,
        genome_dates= processed_csv
    output:
        average_mutations= average_mutations
    group: 9
    shell:
        "python scripts/Count_genome_time_4.py {input.genome_sequences} {input.genome_dates} {output.average_mutations}"

rule L_dist:
    input:
        genome=processed_fasta
    output:
        c_u=genome_c_u,
        c_u_avg=genome_c_u_avg,
        L_dist=genome_mut_L_all
    group: 10
    shell:
        "python scripts/L_distribution.py {input.genome} {output.c_u} {output.c_u_avg} {output.L_dist}"

rule all:
    input: genome_mut_L_all
