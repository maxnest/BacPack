# Snakemake rule for assembling long reads using the Flye program

rule long_reads_assembling:
    input: 
        long_fq={LONG_FQ}
    output:
        long_fasta = OUTPUT_DIR + f"/{SPECIES_TAG}_flye/assembly.fasta"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_flye.log"
    threads: config['threads']
    shell:
        "{PYTHON} {FLYE} --nano-hq {input.long_fq} --threads {threads} --out-dir {OUTPUT_DIR}/{SPECIES_TAG}_flye --iterations 3 2> {log}"
