# Snakemake rules for quality control and trimming of short paired-end reads using FastQC and FastP programs

rule fastq_qc:
    input:
        raw_r1={SAMPLE_R1_FQ},
        raw_r2={SAMPLE_R2_FQ}
    output:
        raw_r1_html= OUTPUT_DIR + f"/{SPECIES_TAG}_fastqc/{r1_tag}_fastqc.html",
        raw_r2_html= OUTPUT_DIR + f"/{SPECIES_TAG}_fastqc/{r2_tag}_fastqc.html"    
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_fastqc.log"
    threads: config['threads']
    message: "Executing FastQC with {threads} threads on the following files {input}."    
    shell:
        "{FASTQC} --outdir {OUTPUT_DIR}/{SPECIES_TAG}_fastqc --noextract --format fastq --threads {threads} {input.raw_r1} {input.raw_r2} 2> {log}" 

rule fastq_trimming:
    input:
        raw_r1={SAMPLE_R1_FQ},
        raw_r2={SAMPLE_R2_FQ}
    output:
        tmm_r1= OUTPUT_DIR + f"/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp_tmm.R1.fastq",
        tmm_r2= OUTPUT_DIR + f"/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp_tmm.R2.fastq",
        failed= OUTPUT_DIR + f"/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp_failed.fastq",
        html= OUTPUT_DIR + f"/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp.html",
        json= OUTPUT_DIR + f"/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp.json"
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_fastp.log"
    params:
        "--cut_right --cut_window_size 4 --cut_mean_quality 20 --qualified_quality_phred 20 --length_required 25"
    threads: config['threads']
    message: "Executing FastP with {threads} threads on the following files {input}."    
    shell:
        "{FASTP} --in1 {input.raw_r1} --out1 {output.tmm_r1} --in2 {input.raw_r2} --out2 {output.tmm_r2} --failed_out {output.failed} {params} --html {output.html} --json {output.json} --thread {threads} 2> {log}"
