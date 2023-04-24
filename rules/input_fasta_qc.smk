# Snakemake rules for quality control of input assembly using the QUAST and BUSCO programs

rule short_seqs_filtering_after_assembly:
    input:
        input_fasta=f"{FASTA}"
    output:
        output_fasta= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    script: "../scripts/filtering_assembed_seqs_by_length.py"

rule selected_seqs_quast_qc:
    input:
        selected_seqs= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    output:
        quast_report= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_selected_seqs_quast/report.tsv"
    params:
        "--min-contig 500 --gene-finding"
    threads: config['threads']
    shell:
        "{PYTHON} {QUAST} --output-dir {OUTPUT_DIR}/QC/{SPECIES_TAG}_selected_seqs_quast {params} --threads {threads} {input.selected_seqs}"

rule selected_seqs_busco_qc_bacillales_odb10:
    input:
        selected_seqs= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    output:
        short_summary= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_selected_seqs_vs_Bacillales_odb10/run_bacillales_odb10/short_summary.txt"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_selected_seqs_busco_qc_bacillales_obd10.log"
    threads: config['threads']
    conda: "../envs/BUSCO_env.yaml"
    params:
        "--mode genome --offline"
    shell:
        "{BUSCO} -i {input.selected_seqs} -o {SPECIES_TAG}_selected_seqs_vs_Bacillales_odb10 -l {BACILLALES_ODB} {params} --cpu {threads} --out_path {OUTPUT_DIR}/QC -f 2> {log}"

rule selected_seqs_busco_qc_bacilli_odb10:
    input:
        selected_seqs= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    output:
        short_summary= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_selected_seqs_vs_Bacilli_odb10/run_bacilli_odb10/short_summary.txt"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_selected_seqs_busco_qc_bacilli_obd10.log"
    threads: config['threads']
    conda: "../envs/BUSCO_env.yaml"
    params:
        "--mode genome --offline"
    shell:
        "{BUSCO} -i {input.selected_seqs} -o {SPECIES_TAG}_selected_seqs_vs_Bacilli_odb10 -l {BACILLI_ODB} {params} --cpu {threads} --out_path {OUTPUT_DIR}/QC -f 2> {log}"
