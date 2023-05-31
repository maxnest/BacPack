# Snakemake rules for polishing of long reads contigs using short reads libraries and the PILON program

rule polishing_using_pilon:
    input:
        flye_medaka_consensus= OUTPUT_DIR + f"/{SPECIES_TAG}_medaka/consensus.fasta",
        tmm_r1= OUTPUT_DIR + f"/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp_tmm.R1.fastq",
        tmm_r2= OUTPUT_DIR + f"/{SPECIES_TAG}_fastp/{SPECIES_TAG}.fastp_tmm.R2.fastq"
    output:
        aln_bam= OUTPUT_DIR + f"/{SPECIES_TAG}_pilon/Round_{ROUNDS}/bwa.round_{ROUNDS}.bam",
        pilon_fasta= OUTPUT_DIR + f"/{SPECIES_TAG}_pilon/Round_{ROUNDS}/pilon.round_{ROUNDS}.fasta"
    threads: config['threads']
    log:
        OUTPUT_DIR + f"/Logs/{SPECIES_TAG}_polishing_using_pilon.log"
    conda: "../envs/BUSCO_env.yaml" # BUSCO
    params:
        "--mode genome --offline" # BUSCO
    shell:
        """
        mkdir {OUTPUT_DIR}/{SPECIES_TAG}_pilon/prev_fasta

        cp {input.flye_medaka_consensus} {OUTPUT_DIR}/{SPECIES_TAG}_pilon/prev_fasta/prev.fasta

        for i in `seq 1 {ROUNDS}`; do
            if [ $i -eq 1 ]
            then
                printf '***** Polishing round %i / %i ***** \n' $i {ROUNDS} >> {log}
                printf 'Input fasta: %s\n' '{input.flye_medaka_consensus}' >> {log}
                {BUSCO} -i {OUTPUT_DIR}/{SPECIES_TAG}_pilon/prev_fasta/prev.fasta -o {SPECIES_TAG}_before_polishing.Bacillales_odb10 -l {BACILLALES_ODB} {params} --cpu {threads} --out_path {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Polishing_QC -f
                printf 'BUSCO before polishing:' >> {log}
                fgrep 'C:' {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Polishing_QC/{SPECIES_TAG}_before_polishing.Bacillales_odb10/short_summary.specific.bacillales_odb10.{SPECIES_TAG}_before_polishing.Bacillales_odb10.txt >> {log}
            else
                printf '***** Polishing round %i / %i ***** \n' $i {ROUNDS} >> {log}
                prev_round=$(expr $i - 1)
                printf 'Input fasta: %s%s%s%s%s\n' '{OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_' $prev_round '/pilon.round_' $prev_round '.fasta' >> {log} 
            fi
            printf 'Output fasta: %s%i%s%i%s\n' '{OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_' $i '/pilon.round_' $i '.fasta' >> {log}
            
            mkdir -p {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i
            # Step 1: Short paired-end reads alignment
            {BWA} index -a bwtsw {OUTPUT_DIR}/{SPECIES_TAG}_pilon/prev_fasta/prev.fasta
            {BWA} mem {OUTPUT_DIR}/{SPECIES_TAG}_pilon/prev_fasta/prev.fasta {input.tmm_r1} {input.tmm_r2} -t {threads} > {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/bwa.round_$i.sam
            {SAMTOOLS} view -S -b {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/bwa.round_$i.sam > {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/bwa.round_$i.bam
            {SAMTOOLS} sort -@ {threads} -o {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/bwa.round_$i.bam {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/bwa.round_$i.bam
            {SAMTOOLS} index {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/bwa.round_$i.bam && rm {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/bwa.round_$i.sam
            # Step 2: Polishing
            java -jar {PILON} --genome {OUTPUT_DIR}/{SPECIES_TAG}_pilon/prev_fasta/prev.fasta --frags {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/bwa.round_$i.bam --outdir {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i --output pilon.round_$i
            sed -i 's/_pilon//g' {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/pilon.round_$i.fasta
            rm {OUTPUT_DIR}/{SPECIES_TAG}_pilon/prev_fasta/prev.fasta*
            cp {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/pilon.round_$i.fasta {OUTPUT_DIR}/{SPECIES_TAG}_pilon/prev_fasta/prev.fasta
            # Step 3: Completeness control using BUSCO 
            {BUSCO} -i {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Round_$i/pilon.round_$i.fasta -o {SPECIES_TAG}_pilon.round_$i.Bacillales_odb10 -l {BACILLALES_ODB} {params} --cpu {threads} --out_path {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Polishing_QC -f
            printf 'BUSCO after %i round of polishing:' $i >> {log}
            fgrep 'C:' {OUTPUT_DIR}/{SPECIES_TAG}_pilon/Polishing_QC/{SPECIES_TAG}_pilon.round_$i.Bacillales_odb10/short_summary.specific.bacillales_odb10.{SPECIES_TAG}_pilon.round_$i.Bacillales_odb10.txt >> {log}
        done

        cp {OUTPUT_DIR}/{SPECIES_TAG}_pilon/prev_fasta/prev.fasta {output.pilon_fasta}
        """

rule short_seqs_filtering_after_polishing:
    input:
        input_fasta= OUTPUT_DIR + f"/{SPECIES_TAG}_pilon/Round_{ROUNDS}/pilon.round_{ROUNDS}.fasta"
    output:
        output_fasta= OUTPUT_DIR + f"/{SPECIES_TAG}_selected_sequences/{SPECIES_TAG}.selected_seqs.fasta"
    script: "../scripts/filtering_assembed_seqs_by_length.py"
