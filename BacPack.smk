configfile: 'config.yaml'

# Input: #
SPECIES_TAG = config['species_tag']
FASTA = config['fasta']
SAMPLE_R1_FQ = config['r1_fq']
SAMPLE_R2_FQ = config['r2_fq']
LONG_FQ = config['long_fq']
OUTPUT_DIR = config['output_dir']
REF_NUCL_FASTA_PATH = config['ref_nucl_fasta_path']
MEDAKA_MODEL = config['medaka_model']

# TEMP:
ROUNDS = 2

# Resources: #
PROKKA_DB = config['prokka_db']
BACILLALES_ODB = config['bacillales_odb']
BACILLI_ODB = config['bacilli_odb']
VFDB = config['vfdb']

# Soft: #
PYTHON = config['python']
FASTQC = config['fastqc']
FASTP = config['fastp']
SPADES = config['spades']
QUAST = config['quast']
BUSCO = config['busco']
FASTANI = config['fastani']
PRODIGAL = config['prodigal']
PROKKA = config['prokka']
CHECKM = config['checkm']
CRYPROCESSOR = config['cryprocessor']
BTTOXIN_DIGGER = config['bttoxin_digger']
RABBITQC = config['rabbitqc']
FLYE = config['flye']
MEDAKA_CONSENSUS = config['medaka_consensus']
BWA = config['bwa']
SAMTOOLS = config['samtools']
PILON = config['pilon']
DEEPBGC = config['deepbgc']
ANTISMASH = config['antismash']
MMSEQS = config['mmseqs']

# Extracting suffixes: #
if SAMPLE_R1_FQ and SAMPLE_R2_FQ:
    r1_tag=SAMPLE_R1_FQ.split("/")[-1].split(".fastq")[0]
    r2_tag=SAMPLE_R2_FQ.split("/")[-1].split(".fastq")[0]

if LONG_FQ: long_reads_tag=LONG_FQ.split('/')[-1].split('.fastq')[0]

# Main results: #
quast_report = f"{OUTPUT_DIR}/QC/{SPECIES_TAG}_selected_seqs_quast/report.tsv",
selected_seqs_vs_bacillales_odb10 = f"{OUTPUT_DIR}/QC/{SPECIES_TAG}_selected_seqs_vs_Bacillales_odb10/run_bacillales_odb10/short_summary.txt",
selected_seqs_vs_bacilli_odb10 = f"{OUTPUT_DIR}/QC/{SPECIES_TAG}_selected_seqs_vs_Bacilli_odb10/run_bacilli_odb10/short_summary.txt",
selected_seqs_fastani = f"{OUTPUT_DIR}/{SPECIES_TAG}_fastANI/{SPECIES_TAG}_selected_seqs_vs_ref.fastANI.out",
selected_seqs_fastani_top_hits = f"{OUTPUT_DIR}/{SPECIES_TAG}_fastANI/{SPECIES_TAG}_selected_seqs_vs_ref.fastANI.top_10_best_hits.tsv",
selected_seqs_fastani_top_hits_merged = f"{OUTPUT_DIR}/{SPECIES_TAG}_fastANI/fastANI_top_ref.merged.fna",
prodigal_trn = f"{OUTPUT_DIR}/{SPECIES_TAG}_prodigal/{SPECIES_TAG}.fastANI_top_ref.trn",
prokka_proteins = f"{OUTPUT_DIR}/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.faa",
prokka_gff3 = f"{OUTPUT_DIR}/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.gff",
prokka_report = f"{OUTPUT_DIR}/{SPECIES_TAG}_prokka/{SPECIES_TAG}_prokka.tsv",
checkm_report = f"{OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_checkm/{SPECIES_TAG}_prokka_checkm",
busco_bacillales_report = f"{OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_proteins_vs_Bacillales_odb10/short_summary.specific.bacillales_odb10.{SPECIES_TAG}_prokka_proteins_vs_Bacillales_odb10.json",
busco_bacilli_report = f"{OUTPUT_DIR}/QC/{SPECIES_TAG}_prokka_proteins_vs_Bacilli_odb10/short_summary.specific.bacilli_odb10.{SPECIES_TAG}_prokka_proteins_vs_Bacilli_odb10.json",
prokka_cryprocessor_log = f"{OUTPUT_DIR}/Logs/{SPECIES_TAG}_cry_processor.log",
cryprocessor_report = f"{OUTPUT_DIR}/{SPECIES_TAG}_prokka_cryprocessor/logs/diamond_matches_{SPECIES_TAG}_prokka.txt",
bttoxin_digger_report= f"{OUTPUT_DIR}/{SPECIES_TAG}_bttoxin_digger/Results/Toxins/All_Toxins.txt",
deepbgc_report = f"{OUTPUT_DIR}/{SPECIES_TAG}_DeepBGC/{SPECIES_TAG}_DeepBGC.bgc.tsv",
antismash_gbk = f"{OUTPUT_DIR}/{SPECIES_TAG}_antiSMASH/{SPECIES_TAG}_antismash.gbk", 
mmseqs_vs_vfdb = f"{OUTPUT_DIR}/{SPECIES_TAG}_mmseqs_vs_VFDB/{SPECIES_TAG}_mmseqs_vs_VFDB.tsv",
mmseqs_vfdb_filtered = f"{OUTPUT_DIR}/{SPECIES_TAG}_mmseqs_vs_VFDB/{SPECIES_TAG}_mmseqs_vs_VFDB.filtered.tsv",
pipeline_general_description = f"{OUTPUT_DIR}/{SPECIES_TAG}_summary/{SPECIES_TAG}_general_description.tsv",
pipeline_toxins_description = f"{OUTPUT_DIR}/{SPECIES_TAG}_summary/{SPECIES_TAG}_toxins_description.tsv",
pipeline_bgc_description = f"{OUTPUT_DIR}/{SPECIES_TAG}_summary/{SPECIES_TAG}_bgc_description.tsv"

rule_all_results = [quast_report, selected_seqs_vs_bacillales_odb10, selected_seqs_vs_bacilli_odb10, selected_seqs_fastani_top_hits, prokka_proteins, prokka_gff3, prokka_report, busco_bacillales_report, busco_bacilli_report, checkm_report, prokka_cryprocessor_log, cryprocessor_report, bttoxin_digger_report, deepbgc_report, antismash_gbk, mmseqs_vs_vfdb, mmseqs_vfdb_filtered, pipeline_general_description, pipeline_toxins_description, pipeline_bgc_description]

# Main: #
if FASTA:
    # Results:
    output_nucl_fasta = f"{OUTPUT_DIR}/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    # Add results wanted to common list
    rule_all_results.append(output_nucl_fasta)
    # Snakemake files:
    include: "rules/input_fasta_qc.smk"
else:
    if LONG_FQ:
        # Results:
        short_reads_qc = [f"{OUTPUT_DIR}/{SPECIES_TAG}_fastqc/{lib}_fastqc.html" for lib in [r1_tag, r2_tag]]
        tmm_reads = [f"{OUTPUT_DIR}/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp_tmm.{pair}.fastq" for pair in ["R1", "R2"]]
        long_reads_qc = [f"{OUTPUT_DIR}/{SPECIES_TAG}_rabbitqc/{lib}_RabbitQCPlus.html" for lib in [long_reads_tag]]
        flye_contigs = f"{OUTPUT_DIR}/{SPECIES_TAG}_flye/assembly.fasta"
        flye_contigs_vs_bacillales_odb10 = f"{OUTPUT_DIR}/QC/{SPECIES_TAG}_flye_contigs_vs_Bacillales_odb10/run_bacillales_odb10/short_summary.txt"
        flye_medaka_consensus = f"{OUTPUT_DIR}/{SPECIES_TAG}_medaka/consensus.fasta"
        flye_medaka_consensus_vs_bacillales_odb10 = f"{OUTPUT_DIR}/QC/Polishing/{SPECIES_TAG}_flye_medaka_consensus_vs_Bacillales_odb10/run_bacillales_odb10/short_summary.txt"
        aln_bam = f"{OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_{ROUNDS}/bwa.round_{ROUNDS}.bam"
        pilon_fasta = f"{OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_{ROUNDS}/pilon.round_{ROUNDS}.fasta"
        output_nucl_fasta = f"{OUTPUT_DIR}/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
        # Add results wanted to common list
        rule_all_results.extend([short_reads_qc, tmm_reads, long_reads_qc, flye_contigs, flye_contigs_vs_bacillales_odb10, flye_medaka_consensus, flye_medaka_consensus_vs_bacillales_odb10, aln_bam, pilon_fasta, output_nucl_fasta])
        # Snakemake files: 
        include: "rules/short_reads_qc_and_trimming.smk"
        include: "rules/long_reads_qc.smk"
        include: "rules/long_reads_assembly.smk"
        include: "rules/long_reads_medaka_consensus.smk"
        include: "rules/medaka_consensus_polishing.smk"
        include: "rules/selected_sequences_qc.smk"
    elif SAMPLE_R1_FQ and SAMPLE_R2_FQ:
        # Results:
        short_reads_qc = [f"{OUTPUT_DIR}/{SPECIES_TAG}_fastqc/{lib}_fastqc.html" for lib in [r1_tag, r2_tag]]
        tmm_reads = [f"{OUTPUT_DIR}/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp_tmm.{pair}.fastq" for pair in ["R1", "R2"]]
        spades_contigs = f"{OUTPUT_DIR}/{SPECIES_TAG}_spades/contigs.fasta"
        spades_scaffolds = f"{OUTPUT_DIR}/{SPECIES_TAG}_spades/scaffolds.fasta"
        output_nucl_fasta = f"{OUTPUT_DIR}/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
        # Add results wanted to common list
        rule_all_results.extend([short_reads_qc, tmm_reads, spades_contigs, spades_scaffolds, output_nucl_fasta])
        # Snakemake files:
        include: "rules/short_reads_qc_and_trimming.smk"
        include: "rules/short_reads_assembly.smk"
        include: "rules/selected_sequences_qc.smk"

rule all:
    input:
        all_results = rule_all_results

##### Rules: #####
include: "rules/whole_genome_ani_estimation.smk"
include: "rules/prokka_annotation_and_qc.smk"
include: "rules/toxins_identification.smk"
include: "rules/annotation.smk"
include: "rules/results_summary.smk"
