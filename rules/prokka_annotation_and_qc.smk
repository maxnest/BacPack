# Snakemake rules for protein-coding genes identification and both completeness and quality checking

rule selected_seqs_annotation_with_prokka:
    input:
        fastani_out= OUTPUT_DIR + f"/{SPECIES_TAG}_fastANI/{SPECIES_TAG}_selected_seqs_vs_ref.fastANI.out",
        selected_seqs= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    output:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa",
        prokka_gff3= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.gff",
        prokka_report= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.tsv"
    params:
        "--force --addgenes --gffver '3' --genus 'Bacillus' --usegenus"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_prodigal_model_building.log"
    threads: config['threads']
    run:
        import os
        import pandas as pd

        if os.stat(input.fastani_out).st_size != 0:
            # Top 10 reference genomes
            fastANI_output_df = pd.read_csv(input.fastani_out, sep='\t', header=None).rename(
                columns={0: 'Query', 1: 'Ref', 2: 'ANI', 3: 'Orthologous_matches', 4: 'Query_fragments'})
            fastani_top_refs = fastANI_output_df.sort_values('ANI', ascending=False).head(10)
            fastani_top_refs.to_csv(f"{OUTPUT_DIR}/{SPECIES_TAG}_fastANI/{SPECIES_TAG}_selected_seqs_vs_ref.fastANI.top_10_best_hits.tsv", sep='\t', index=False)
            # Prodigal model training
            shell("awk -F' ' 'NR>1 {{print $2}}' {OUTPUT_DIR}/{SPECIES_TAG}_fastANI/{SPECIES_TAG}_selected_seqs_vs_ref.fastANI.top_10_best_hits.tsv | tr '\n' ' ' | sed 's/^/cat /' | sed -e 's#$#> {OUTPUT_DIR}/{SPECIES_TAG}_fastANI/fastANI_top_ref.merged.fna#' > {OUTPUT_DIR}/{SPECIES_TAG}_fastANI/prodigal_prep.sh && bash {OUTPUT_DIR}/{SPECIES_TAG}_fastANI/prodigal_prep.sh && {PRODIGAL} -i {OUTPUT_DIR}/{SPECIES_TAG}_fastANI/fastANI_top_ref.merged.fna -t {OUTPUT_DIR}/{SPECIES_TAG}_fastANI/{SPECIES_TAG}.fastANI_top_ref.trn 2> {log}")
            # Prokka annotation with model with references
            shell("{PROKKA} --outdir {OUTPUT_DIR}/{SPECIES_TAG}_prokka --prefix {SPECIES_TAG}_prokka {params} --prodigaltf {OUTPUT_DIR}/{SPECIES_TAG}_fastANI/{SPECIES_TAG}.fastANI_top_ref.trn --cpus {threads} {input.selected_seqs} --proteins {PROKKA_DB}")
        else:
            # Prokka annotation with model without references
            shell("{PROKKA} --outdir {OUTPUT_DIR}/{SPECIES_TAG}_prokka --prefix {SPECIES_TAG}_prokka {params} --cpus {threads} {input.selected_seqs} --proteins {PROKKA_DB}")


rule prokka_proteins_busco_qc_bacillales_odb10:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    output:
        short_summary= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_proteins_vs_Bacillales_odb10/run_bacillales_odb10/short_summary.txt",
        busco_bacillales_report= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_proteins_vs_Bacillales_odb10/short_summary.specific.bacillales_odb10.{SPECIES_TAG}_prokka_proteins_vs_Bacillales_odb10.json"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_prokka_proteins_busco_qc_bacillales_obd10.log"
    threads: config['threads']
    conda: "../envs/BUSCO_env.yaml"
    params:
        "--mode proteins --offline"
    shell:
        "{BUSCO} -i {input.prokka_proteins} -o {SPECIES_TAG}_prokka_proteins_vs_Bacillales_odb10 -l {BACILLALES_ODB} {params} --cpu {threads} --out_path {OUTPUT_DIR}/QC -f 2> {log}"

rule prokka_proteins_busco_qc_bacilli_odb10:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    output:
        short_summary= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_proteins_vs_Bacilli_odb10/run_bacilli_odb10/short_summary.txt",
        busco_bacilli_report= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_proteins_vs_Bacilli_odb10/short_summary.specific.bacilli_odb10.{SPECIES_TAG}_prokka_proteins_vs_Bacilli_odb10.json"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_prokka_proteins_busco_qc_bacilli_obd10.log"
    threads: config['threads']
    conda: "../envs/BUSCO_env.yaml"
    params:
        "--mode proteins --offline"
    shell:
        "{BUSCO} -i {input.prokka_proteins} -o {SPECIES_TAG}_prokka_proteins_vs_Bacilli_odb10 -l {BACILLI_ODB} {params} --cpu {threads} --out_path {OUTPUT_DIR}/QC -f 2> {log}"

rule prokka_checkm_qc:
    input:
        prokka_proteins= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa"
    output:
        checkm_report= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_checkm/{SPECIES_TAG}_prokka_checkm"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_prokka_checkm.log"
    threads: config['threads']
    conda: "../envs/CheckM_env.yaml"
    shell:
        "mkdir {OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm/Prokka_faa && cp {input.prokka_proteins} {OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm/Prokka_faa/{SPECIES_TAG}_prokka.faa && cd {OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm/ && {CHECKM} lineage_wf --genes -t {threads} --pplacer_threads {threads} --tab_table -x faa -f {SPECIES_TAG}_prokka_checkm {OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm/Prokka_faa {OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm 2> {log}"
