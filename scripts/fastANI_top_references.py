import pandas as pd

def fastANI_top_references(fastani_output, top, output_tag):
    fastANI_output_df = pd.read_csv(fastani_output, sep='\t', header=None).rename(
        columns={0: 'Query', 1: 'Ref', 2: 'ANI', 3: 'Orthologous_matches', 4: 'Query_fragments'})
    top_refs = fastANI_output_df.sort_values('ANI', ascending=False).head(top)
    top_refs.to_csv(output_tag, sep='\t', index=False)


if __name__ == "__main__":
    fastANI_top_references(snakemake.input.fastani_out, 10, snakemake.output.fastani_top)
