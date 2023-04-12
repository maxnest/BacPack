# Snakemake rules for de novo assembling of short paired-end reads using the SPAdes program

rule paired_reads_assembling:
    input:
        tmm_r1= OUTPUT_DIR + f"/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp_tmm.R1.fastq",
        tmm_r2= OUTPUT_DIR + f"/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp_tmm.R2.fastq"
    output:
        contigs= OUTPUT_DIR + f"/{SPECIES_TAG}_spades/contigs.fasta",
        scaffolds= OUTPUT_DIR + f"/{SPECIES_TAG}_spades/scaffolds.fasta"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_spades.log"
    threads: config['threads']
    message: "Executing SPAdes with {threads} threads on the following files {input}."
    shell:
        "{PYTHON} {SPADES} -1 {input.tmm_r1} -2 {input.tmm_r2} --careful --threads {threads} -o {OUTPUT_DIR}/{SPECIES_TAG}_spades 2> {log}"

rule short_seqs_filtering_after_assembly:
    input:
        input_fasta= OUTPUT_DIR + f"/{SPECIES_TAG}_spades/scaffolds.fasta"
    output:
        output_fasta= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    script: "../scripts/filtering_assembed_seqs_by_length.py"
