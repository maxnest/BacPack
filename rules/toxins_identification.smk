# Snakemake rules for toxins identification for strain using the CryProcessor and BtToxin_Digger programs

rule run_cryprocessor:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    output:
        cryprocessor_report= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka_cryprocessor/logs/diamond_matches_{SPECIES_TAG}_prokka.txt"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_cry_processor.log"
    shell:
        "{PYTHON} {CRYPROCESSOR} --force -fi {input.prokka_proteins} -od {OUTPUT_DIR}/{SPECIES_TAG}_prokka_cryprocessor -r do 2> {log}"

rule run_bttoxin_digger:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    output:
        bttoxin_digger_report= OUTPUT_DIR + f"/{SPECIES_TAG}_bttoxin_digger/Results/Toxins/All_Toxins.txt"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_BtToxin_Digger.log"
    threads: config['threads']
    conda: "../envs/bttoxin_digger_env.yaml"
    shell:
        "mkdir -p {OUTPUT_DIR}/{SPECIES_TAG}_bttoxin_digger && cd {OUTPUT_DIR}/{SPECIES_TAG}_bttoxin_digger && cp {input.prokka_proteins} . && {BTTOXIN_DIGGER} --SeqPath {OUTPUT_DIR}/{SPECIES_TAG}_bttoxin_digger --SequenceType prot --prot_suffix .faa --threads {threads} 2> {log}"

