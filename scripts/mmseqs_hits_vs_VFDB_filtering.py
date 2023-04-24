import pandas as pd


def filter_mmseqs_results(vir_pd, id_threshold, cov_threshold):
    filtered_vir_pd = vir_pd[(vir_pd['qcov'] * 100 >= cov_threshold) &
                             (vir_pd['tcov'] * 100 >= cov_threshold) &
                             (vir_pd['pident'] >= id_threshold)]
    
    filtered_vir_pd_grouped = filtered_vir_pd.groupby('query', as_index=False).apply(
        lambda group: group[(group['pident'] == group['pident'].max()) &
                            (group['evalue'] == group['evalue'].min())])

    return filtered_vir_pd_grouped



if __name__ == "__main__":
    vir_pd = pd.read_csv(snakemake.input.mmseqs_vfdb, sep='\t',
                         names=['query', 'target', 'pident', 'alnlen',
                                                     'mismatch', 'gapopen', 'qstart', 'qend',
                                                     'tstart', 'tend', 'evalue', 'bits',
                                                     'qcov', 'tcov'],
                         dtype={'query': str, 'target': str, 'pident': float, 'alnlen': int,
                                                     'mismatch': int, 'gapopen': int, 'qstart': int, 'qend': int,
                                                     'tstart': int, 'tend': int, 'evalue': float, 'bits': float,
                                                     'qcov': float, 'tcov': float})
    filtered_vir_pd = filter_mmseqs_results(vir_pd, 70, 70)

    if len(filtered_vir_pd.index) != 0:
        filtered_vir_pd.to_csv(snakemake.output.mmseqs_vfdb_filtered, sep='\t', header=True, index=False)
