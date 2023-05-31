# Snakemake rules for sequence annotation using the DeepBGC, antiSMASH, and MMseq2 (versus VFDB) programs

rule run_deepbgc:
    input:
        selected_seqs= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    output:
        deepbgc_report= OUTPUT_DIR + f"/{SPECIES_TAG}_DeepBGC/{SPECIES_TAG}_DeepBGC.bgc.tsv"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_DeepBGC.log"
    conda: "../envs/DeepBGC_env.yaml"
    params:
        "--score 0.7 --prodigal-meta-mode"
    shell:
        "{DEEPBGC} pipeline -o {OUTPUT_DIR}/{SPECIES_TAG}_DeepBGC {params} {input.selected_seqs} 2> {log}"

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

rule run_mmseqs_vs_vfdb:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    output:
        mmseqs_vfdb= OUTPUT_DIR + f"/{SPECIES_TAG}_mmseqs_vs_VFDB/{SPECIES_TAG}_mmseqs_vs_VFDB.tsv"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_mmseqs_vs_VFDB.log"
    threads: config['threads']
    shell:
        "{MMSEQS} easy-search {input.prokka_proteins} {VFDB} {output.mmseqs_vfdb} tmp_mmseqs --threads {threads} --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov'"

rule mmseqs_hits_filtering:
    input:
        mmseqs_vfdb= OUTPUT_DIR + f"/{SPECIES_TAG}_mmseqs_vs_VFDB/{SPECIES_TAG}_mmseqs_vs_VFDB.tsv"
    output:
        mmseqs_vfdb_filtered= OUTPUT_DIR + f"/{SPECIES_TAG}_mmseqs_vs_VFDB/{SPECIES_TAG}_mmseqs_vs_VFDB.filtered.tsv"
    script: "../scripts/mmseqs_hits_vs_VFDB_filtering.py"
