import pandas as pd


def quast_report_parsing(quast, summary_dict):
    with open(quast) as quast_report:
        species_tag = quast_report.readline().split("\t")[1].split(".selected_seqs")[0]
        summary_dict['Species_tag'] = species_tag
        for line in quast_report:
            description = line.strip().split("\t")
            summary_dict[f'Quast: {description[0]}'] = description[1]


def prokka_tsv_parsing(prokka, summary_dict):
    prokka_pd = pd.read_csv(prokka, sep='\t',
                            dtype={'locus_tag': str, 'ftype': str, 'length_bp': int, 'gene': str, 'EC_number': str,
                                   'COG': str, 'product': str})
    only_cds, only_genes = prokka_pd[(prokka_pd['ftype'] == 'CDS')], prokka_pd[(prokka_pd['ftype'] == 'gene')]
    summary_dict['Prokka: # genes'] = only_genes.shape[0]
    summary_dict['Prokka: # CDS'] = only_cds.shape[0]
    summary_dict['Prokka: # hypothetical protein'] = len([protein for protein in only_cds['product'].to_list()
                                                          if 'hypothetical protein' in protein.lower()])


def checkm_tab_parsing(checkm, summary_dict):
    checkm_pd = pd.read_csv(checkm, sep='\t',
                            dtype={'Bin Id': str, 'Marker lineage': str, '# genomes': int, '# markers': int,
                                  '# marker sets': int, '0': int, '1': int, '2': int, '3': int,	'4': int, '5+': int,
                                  'Completeness': str, 'Contamination': str, 'Strain heterogeneity': str})

    summary_dict['CheckM: Completeness'] = checkm_pd['Completeness'].to_list()[0]
    summary_dict['CheckM: Contamination'] = checkm_pd['Contamination'].to_list()[0]
    summary_dict['CheckM: Strain heterogeneity'] = checkm_pd['Strain heterogeneity'].to_list()[0]


def busco_json_parsing(busco, summary_dict):
    busco_pd = pd.read_json(busco)
    lineage_info = busco_pd.loc['lineage_dataset']['parameters'].split('/')[-1]
    one_line_summary = busco_pd.loc['one_line_summary']['results']
    summary_dict[f'BUSCO: {lineage_info}'] = one_line_summary


def summary_output(outfile, summary_dict):
    with open(outfile, 'a') as output_file:
        for key, value in summary_dict.items():
            output_file.write(f'{key}\t{value}\n')


if __name__ == '__main__':
    summary_dict = {}
    # Quast #
    quast_report_parsing(snakemake.input.quast_report, summary_dict)
    # Prokka #
    prokka_tsv_parsing(snakemake.input.prokka_report, summary_dict)
    # CheckM #
    checkm_tab_parsing(snakemake.input.checkm_report, summary_dict)
    # BUSCO #
    busco_json_parsing(snakemake.input.busco_bacillales_report, summary_dict)
    busco_json_parsing(snakemake.input.busco_bacilli_report, summary_dict)
    # Output writing #
    summary_output(snakemake.output.pipeline_general_description, summary_dict)