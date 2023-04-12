# Snakemake rule for quality control of long reads using the RabbitQCPlus program

rule rabbit_qc:
    input:
        long_fq={LONG_FQ}
    output:
        long_fq_html = OUTPUT_DIR + f"/{SPECIES_TAG}_rabbitqc/" + f"{long_reads_tag}_RabbitQCPlus.html"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_long_reads_rabbitqc.log"
    threads: config['threads']
    shell:
        "{RABBITQC} -i {input.long_fq} --overWrite --threadNum {threads} --TGS --doOverrepresentation --printORPSeqs && mv {long_reads_tag}_RabbitQCPlus.html {OUTPUT_DIR}/{SPECIES_TAG}_rabbitqc/ 2> {log}"
