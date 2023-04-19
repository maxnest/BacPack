# Snakemake rules for toxins identification for strain using the CryProcessor and IDOPs programs

rule run_cryprocessor:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_cry_processor.log"
    shell:
        "{PYTHON} {CRYPROCESSOR} -fi {input.prokka_proteins} -od {OUTPUT_DIR}/{SPECIES_TAG}_prokka_cryprocessor -r do 2> {log}"

rule run_idops:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_idops.log"
    conda: "../envs/IDOPS_env.yaml"
    shell:
        "{IDOPS} --verbosity --output {OUTPUT_DIR}/{SPECIES_TAG}_prokka_idops -t {input.prokka_proteins} 2> {log}"
