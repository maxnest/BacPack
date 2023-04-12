# Snakemake rules for protein-coding genes identification and both completeness and quality checking

rule prodigal_prep:
    input:
        fastani_top= OUTPUT_DIR + f"/{SPECIES_TAG}_fastANI/{SPECIES_TAG}_selected_seqs_vs_ref.fastANI.top_10_best_hits.tsv"
    output:
        ref_fna_merged= OUTPUT_DIR + f"/{SPECIES_TAG}_fastANI/fastANI_top_ref.merged.fna"
    shell:
        "awk -F' ' 'NR>1 {{print $2}}' {input.fastani_top} | tr '\n' ' ' | sed 's/^/cat /' | sed -e 's#$#> {output.ref_fna_merged}#' > {OUTPUT_DIR}/{SPECIES_TAG}_fastANI/prodigal_prep.sh && bash {OUTPUT_DIR}/{SPECIES_TAG}_fastANI/prodigal_prep.sh"

rule prodigal_model_building:
    input:
        ref_fna_merged= OUTPUT_DIR + f"/{SPECIES_TAG}_fastANI/fastANI_top_ref.merged.fna"
    output:
        ref_fna_merged_model= OUTPUT_DIR + f"/{SPECIES_TAG}_prodigal/{SPECIES_TAG}.fastANI_top_ref.trn"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_prodigal_model_building.log"
    shell:
        "{PRODIGAL} -i {input.ref_fna_merged} -t {output.ref_fna_merged_model} 2> {log}"

rule selected_seqs_annotation_with_prokka:
    input:
        selected_seqs= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta",
        ref_fna_merged_model= OUTPUT_DIR + f"/{SPECIES_TAG}_prodigal/{SPECIES_TAG}.fastANI_top_ref.trn"
    output:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa",
        prokka_gff3= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.gff"
    params:
        "--force --addgenes --gffver '3' --genus 'Bacillus' --usegenus"
    threads: config['threads']
    shell:
        "{PROKKA} --outdir {OUTPUT_DIR}/{SPECIES_TAG}_prokka --prefix {SPECIES_TAG}_prokka {params} --prodigaltf {input.ref_fna_merged_model} --cpus {threads} {input.selected_seqs} --proteins {PROKKA_DB}"

rule prokka_proteins_busco_qc_bacillales_odb10:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    output:
        short_summary= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_proteins_vs_Bacillales_odb10/run_bacillales_odb10/short_summary.txt"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_prokka_proteins_busco_qc_bacillales_obd10.log"
    threads: config['threads']
    message: "Executing BUSCO with {threads} threads on the following file {input}."
    conda: "../envs/BUSCO_env.yaml"
    params:
        "--mode proteins --offline"
    shell:
        "{BUSCO} -i {input.prokka_proteins} -o {SPECIES_TAG}_prokka_proteins_vs_Bacillales_odb10 -l {BACILLALES_ODB} {params} --cpu {threads} --out_path {OUTPUT_DIR}/QC -f 2> {log}"

rule prokka_proteins_busco_qc_bacilli_odb10:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    output:
        short_summary= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_proteins_vs_Bacilli_odb10/run_bacilli_odb10/short_summary.txt"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_prokka_proteins_busco_qc_bacilli_obd10.log"
    threads: config['threads']
    message: "Executing BUSCO with {threads} threads on the following files {input}."
    conda: "../envs/BUSCO_env.yaml"
    params:
        "--mode proteins --offline"
    shell:
        "{BUSCO} -i {input.prokka_proteins} -o {SPECIES_TAG}_prokka_proteins_vs_Bacilli_odb10 -l {BACILLI_ODB} {params} --cpu {threads} --out_path {OUTPUT_DIR}/QC -f 2> {log}"

rule prokka_checkm_qc:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    output:
        prokka_proteins_checkm= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_checkm/{SPECIES_TAG}_prokka_checkm"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_prokka_checkm.log"
    threads: config['threads']
    conda: "../envs/CheckM_env.yaml"
    shell:
        "mkdir {OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm/Prokka_faa && cp {input.prokka_proteins} {OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm/Prokka_faa/{SPECIES_TAG}_prokka.faa && cd {OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm/ && {CHECKM} lineage_wf --genes -t {threads} --pplacer_threads {threads} --tab_table -x faa -f {SPECIES_TAG}_prokka_checkm {OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm/Prokka_faa {OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm 2> {log}"
