# Snakemake rules for pipeline results summarizing

rule general_description:
    input:
        btyper3_report= OUTPUT_DIR + f"/{SPECIES_TAG}_btyper3/btyper3_final_results/{SPECIES_TAG}.selected_seqs_final_results.txt",
        quast_report= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_selected_seqs_quast/report.tsv",
        prokka_report= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.tsv",
        checkm_report= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_checkm/{SPECIES_TAG}_prokka_checkm",
        busco_bacillales_report= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_proteins_vs_Bacillales_odb10/short_summary.specific.bacillales_odb10.{SPECIES_TAG}_prokka_proteins_vs_Bacillales_odb10.json",
        busco_bacilli_report= OUTPUT_DIR + f"/QC/{SPECIES_TAG}_prokka_proteins_vs_Bacilli_odb10/short_summary.specific.bacilli_odb10.{SPECIES_TAG}_prokka_proteins_vs_Bacilli_odb10.json"
    output:
        pipeline_general_description= OUTPUT_DIR + f"/{SPECIES_TAG}_summary/{SPECIES_TAG}_general_description.tsv"
    script: "../scripts/pipeline_results_general_summary.py"

rule toxins_description:
    input:
        bttoxin_digger_report= OUTPUT_DIR + f"/{SPECIES_TAG}_bttoxin_digger/Results/Toxins/All_Toxins.txt", 
        cryprocessor_report= OUTPUT_DIR + f"/{SPECIES_TAG}_prokka_cryprocessor/logs/diamond_matches_{SPECIES_TAG}_prokka.txt",
        mmseqs_vfdb_filtered= OUTPUT_DIR + f"/{SPECIES_TAG}_mmseqs_vs_VFDB/{SPECIES_TAG}_mmseqs_vs_VFDB.filtered.tsv",
    output:
        pipeline_toxins_description= OUTPUT_DIR + f"/{SPECIES_TAG}_summary/{SPECIES_TAG}_toxins_description.tsv"
    script: "../scripts/pipeline_toxins_identification_results_summary.py"

rule bgc_description:
    input:
        deepbgc_report= OUTPUT_DIR + f"/{SPECIES_TAG}_DeepBGC/{SPECIES_TAG}_DeepBGC.bgc.tsv"
    output:
        pipeline_bgc_description= OUTPUT_DIR + f"/{SPECIES_TAG}_summary/{SPECIES_TAG}_bgc_description.tsv"
    script: "../scripts/pipeline_bgc_identification_results_summary.py"
