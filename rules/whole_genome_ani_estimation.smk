# Snakemake rules for whole-genome similarity (ANI) estimation between strain analyzed and references assemblies

rule ref_fasta_prep:
    input:
        ref_nucl_fasta_path=REF_NUCL_FASTA_PATH
    output:
        ref_nucl_fasta_list= OUTPUT_DIR + f"/References/Ref_nucl_fasta_list.txt"
    shell:
        "find '{input.ref_nucl_fasta_path}' -name '*.fna' > {output.ref_nucl_fasta_list}"

rule whole_genome_ani_estimation:
    input:
        selected_seqs= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta",
        ref_nucl_fasta_list= OUTPUT_DIR + f"/References/Ref_nucl_fasta_list.txt"
    output:
        fastani_out= OUTPUT_DIR + f"/{SPECIES_TAG}_fastANI/{SPECIES_TAG}_selected_seqs_vs_ref.fastANI.out"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_selected_seqs_fastANI.log"
    threads: config['threads']
    message: "Executing fastANI with {threads} threads on the following file {input}."
    shell:
        "{FASTANI} -q {input.selected_seqs} --rl {input.ref_nucl_fasta_list} --kmer 16 --threads {threads} --output {output.fastani_out} 2> {log}"

rule top_references_selecting:
    input:
        fastani_out= OUTPUT_DIR + f"/{SPECIES_TAG}_fastANI/{SPECIES_TAG}_selected_seqs_vs_ref.fastANI.out"
    output:
        fastani_top= OUTPUT_DIR + f"/{SPECIES_TAG}_fastANI/{SPECIES_TAG}_selected_seqs_vs_ref.fastANI.top_10_best_hits.tsv"
    script: "../scripts/fastANI_top_references.py"
