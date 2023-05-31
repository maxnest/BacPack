# Snakemake fules for in silico taxonomic classification of Bacillus cereus group isolates using assembled genome

rule btyper3_classification:
    input:
        selected_seqs= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    output:
        btyper3_report= OUTPUT_DIR + f"/{SPECIES_TAG}_btyper3/btyper3_final_results/{SPECIES_TAG}.selected_seqs_final_results.txt"
    conda: "../envs/BTyper3_env.yaml"
    shell:
        "{BTYPER3} -i {input.selected_seqs} -o {OUTPUT_DIR}/{SPECIES_TAG}_btyper3"
