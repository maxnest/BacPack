# Snakemake rules for quality control of long reads contigs using the BUSCO and consensus identification with the Medaka program

rule long_read_contigs_busco_qc_bacillales_odb10:
    input:
        flye_contigs= OUTPUT_DIR + f"/{SPECIES_TAG}_flye/assembly.fasta"
    output:
        short_summary= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_flye_contigs_vs_Bacillales_odb10/run_bacillales_odb10/short_summary.txt"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_flye_contigs_busco_qc_bacillales_obd10.log"
    threads: config['threads']
    message: "Executing BUSCO with {threads} threads on the following files {input}."
    conda: "../envs/BUSCO_env.yaml"
    params:
        "--mode genome --offline"
    shell:
        "{BUSCO} -i {input.flye_contigs} -o {SPECIES_TAG}_flye_contigs_vs_Bacillales_odb10 -l {BACILLALES_ODB} {params} --cpu {threads} --out_path {OUTPUT_DIR}/QC -f 2> {log}"

rule flye_medaka_consensus:
    input:
        long_fq={LONG_FQ},
        flye_contigs= OUTPUT_DIR + f"/{SPECIES_TAG}_flye/assembly.fasta"
    output:
        medaka_consensus= OUTPUT_DIR + f"/{SPECIES_TAG}_medaka/consensus.fasta"
    container: f"{MEDAKA_DOCKER}"
    log: OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_medaka.log"
    threads: config['threads']
    shell:
        """
        {DOCKER} run --gpus all -it --rm \
             -v={input.long_fq}:{input.long_fq} \
             -v={OUTPUT_DIR}/{SPECIES_TAG}_flye:{OUTPUT_DIR}/{SPECIES_TAG}_flye \
             -v={OUTPUT_DIR}/{SPECIES_TAG}_medaka:{OUTPUT_DIR}/{SPECIES_TAG}_medaka \
             -v={OUTPUT_DIR}/Logs:{OUTPUT_DIR}/Logs \
             {MEDAKA_DOCKER} medaka_consensus -i {input.long_fq} -d {input.flye_contigs} -o {OUTPUT_DIR}/{SPECIES_TAG}_medaka -m {MEDAKA_MODEL} -t {threads}
        """

rule flye_medaka_consensus_busco_qc_bacillales_odb10:
    input:
        flye_medaka_consensus= OUTPUT_DIR + f"/{SPECIES_TAG}_medaka/consensus.fasta"
    output:
        short_summary= OUTPUT_DIR + f"/QC/Polishing/{SPECIES_TAG}_flye_medaka_consensus_vs_Bacillales_odb10/run_bacillales_odb10/short_summary.txt"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_flye_medaka_consensus_busco_qc_bacillales_obd10.log"
    threads: config['threads']
    message: "Executing BUSCO with {threads} threads on the following files {input}."
    conda: "../envs/BUSCO_env.yaml"
    params:
        "--mode genome --offline"
    shell:
        "{BUSCO} -i {input.flye_medaka_consensus} -o {SPECIES_TAG}_flye_medaka_consensus_vs_Bacillales_odb10 -l {BACILLALES_ODB} {params} --cpu {threads} --out_path {OUTPUT_DIR}/QC/Polishing -f 2> {log}"
