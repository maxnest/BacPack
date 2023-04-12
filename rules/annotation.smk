# Snakemake rules for sequence annotation using the DeepBGC and antiSMASH programs

rule run_deepbgc:
    input:
        selected_seqs= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_DeepBGC.log"
    conda: "../envs/DeepBGC_env.yaml"
    params:
        "--score 0.7 --prodigal-meta-mode"
    shell:
        """
        mkdir {OUTPUT_DIR}/{SPECIES_TAG}_DeepBGC
        cd {OUTPUT_DIR}/{SPECIES_TAG}_DeepBGC

        mkdir {OUTPUT_DIR}/{SPECIES_TAG}_DeepBGC/deepdata        
        export DEEPBGC_DOWNLOADS_DIR={OUTPUT_DIR}/{SPECIES_TAG}_DeepBGC/deepdata ;

        {DEEPBGC} download
        {DEEPBGC} pipeline -o {OUTPUT_DIR}/{SPECIES_TAG}_DeepBGC {params} {input.selected_seqs} 2> {log}
        """

rule run_antismash:
    input:
        selected_seqs= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta",
        prokka_gff3= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.gff"
    output:
        antismash_gbk= OUTPUT_DIR + f"/{SPECIES_TAG}_antiSMASH/{SPECIES_TAG}_antismash.gbk"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_antismash.log"
    conda: "../envs/antiSMASH_env.yaml"
    threads: config['threads']
    shell:
        "{ANTISMASH} --cb-knownclusters --fullhmmer --cpus {threads} --output-basename {SPECIES_TAG}_antismash --output-dir {OUTPUT_DIR}/{SPECIES_TAG}_antiSMASH --genefinding-gff3 {input.prokka_gff3} {input.selected_seqs}"
