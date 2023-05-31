import pandas as pd


def deepbgc_tsv_parsing(deepbgc_tsv, summary_dict):
    deepbgc_pd = pd.read_csv(deepbgc_tsv,
                             dtype={'sequence_id': str, 'detector': str, 'detector_version': str, 'detector_label': str,
                                    'bgc_candidate_id': str, 'nucl_start': int,	'nucl_end': int, 'nucl_length': int,
                                    'num_proteins': int, 'num_domains': int, 'num_bio_domains': int,
                                    'deepbgc_score': float, 'product_activity': str, 'antibacterial': float,
                                    'cytotoxic': float, 'inhibitor': float, 'antifungal': float, 'product_class': str,
                                    'Alkaloid':	float, 'NRP': float, 'Other': float, 'Polyketide': float, 'RiPP': float,
                                    'Saccharide': float, 'Terpene': float,	'protein_ids': str, 'bio_pfam_ids': str,
                                    'pfam_ids': str}, sep='\t')

    deepbgc_pd_subset = deepbgc_pd[(deepbgc_pd['product_activity'].notnull()) |
                                   (deepbgc_pd['product_class'].notnull())]

    for seqid, method, method_version, score, activity, product_class, domains in zip(
            deepbgc_pd_subset['sequence_id'], deepbgc_pd_subset['detector'], deepbgc_pd_subset['detector_version'],
            deepbgc_pd_subset['deepbgc_score'], deepbgc_pd_subset['product_activity'],
            deepbgc_pd_subset['product_class'], deepbgc_pd_subset['pfam_ids']):
        summary_dict[len(summary_dict) + 1] = {'seqid': seqid, 'method': f'{method} (v{method_version})',
                                               'score': score, 'activity': activity, 'product_class': product_class,
                                               'domains': domains}


def output_writing(outfile, summary_pd):
    summary_pd.to_csv(outfile, sep='\t', index=False,
                      header=['Sequence_ID', 'Method', 'Score', 'Activity', 'Product_class', 'PfamA domains'])


if __name__ == '__main__':
    summary_dict = {}
    # Input files parsing #
    deepbgc_tsv_parsing(snakemake.input.deepbgc_report, summary_dict)
    # Summary: from dict to pandas df #
    summary_pd = pd.DataFrame.from_dict(summary_dict, orient='index')
    # Output writing #
    output_writing(snakemake.output.pipeline_bgc_description, summary_pd)