from Bio import SeqIO

def seqs_filtering(assembly, min_len, output_fasta):
    fasta_seqs = SeqIO.parse(assembly, "fasta")
    with open(output_fasta, "a") as outfasta:
        for fasta in fasta_seqs:
            name, seq = fasta.id, fasta.seq
            if len(seq) >= int(min_len):
                outfasta.write(f">{name}\n{seq}\n")


if __name__ == "__main__":
    seqs_filtering(snakemake.input.input_fasta, snakemake.config['min_len'], snakemake.output.output_fasta)
