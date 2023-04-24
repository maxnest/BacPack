import pandas as pd


def cryprocessor_txt_parsing(cryprocessor, summary_dict):
    with open(cryprocessor) as cryprocessor_tab:
        for line in cryprocessor_tab:
            description = line.strip().split('\t')
            qseqid, sseqid, pident, evalue, bitscore = description[0], description[1], description[2], \
                                                       description[-2], description[-1]
            summary_dict[len(summary_dict) + 1] = {'qseqid': qseqid, 'sseqid': sseqid, 'pident': pident,
                                                   'evalue': evalue, 'bitscore': bitscore, 'method': 'CryProcessor'}


def mmseqs_tsv_parsing(mmseqs, summary_dict):
    with open(mmseqs) as mmseqs_tab:
        header = mmseqs_tab.readline()
        for line in mmseqs_tab:
            description = line.strip().split('\t')
            qseqid, sseqid, pident, evalue, bitscore = description[0], description[1], description[2], \
                                                       description[-4], description[-3]
            summary_dict[len(summary_dict) + 1] = {'qseqid': qseqid, 'sseqid': sseqid, 'pident': pident,
                                                   'evalue': evalue, 'bitscore': bitscore, 'method': 'VFDB: MMseq2'}


def bttoxin_digger_txt_parsing(bttoxin_digger, summary_dict):
    bttoxin_pd = pd.read_csv(
        bttoxin_digger, dtype={'Strain': str, 'Protein_id': str, 'Protein_len': int, 'Strand': str,
                               'Gene location on scaffold': str, 'SVM': str, 'BLAST': str, 'HMM': str, 'Hit_id': str,
                               'Hit_length': str, 'Aln_length': str, 'Query start-end': str,	'Hit stard-end': str,
                               'Identity': str,	'Evalue of blast': str,	'Hmm hit': str, 'Hmm hit length': str,
                               'Evalue of Hmm': str, 'Nomenclature': str, 'Endotoxin_N': str, 'Endotoxin_M': str,
                               'Endotoxin_C': str,	'Endotoxin_mid': str, 'Toxin_10': str,	'ETX_MTX2': str,
                               'Gene sequence': str,	'Protein sequence': str, 'Scaffold sequence': str}, sep='\t')

    bttoxin_pd_blast, bttoxin_pd_hmm = bttoxin_pd[(bttoxin_pd['BLAST'] == 'YES')], \
                                       bttoxin_pd[(bttoxin_pd['HMM'] == 'YES')]

    if len(bttoxin_pd_blast.index) != 0:
        for qseqid, sseqid, pident, evalue in zip(bttoxin_pd_blast['Protein_id'], bttoxin_pd_blast['Hit_id'],
                                                 bttoxin_pd_blast['Identity'], bttoxin_pd_blast['Evalue of blast']):
            summary_dict[len(summary_dict) + 1] = {'qseqid': qseqid, 'sseqid': sseqid, 'pident': pident,
                                                   'evalue': evalue, 'bitscore': '-', 'method': 'BtToxin_Digger: BLAST'}

    if len(bttoxin_pd_hmm.index) != 0:
        for qseqid, sseqid, evalue in zip(bttoxin_pd_hmm['Protein_id'], bttoxin_pd_hmm['Hmm hit'],
                                          bttoxin_pd_hmm['Evalue of Hmm']):
            summary_dict[len(summary_dict) + 1] = {'qseqid': qseqid, 'sseqid': sseqid, 'pident': '-',
                                                   'evalue': evalue, 'bitscore': '-', 'method': 'BtToxin_Digger: HMM'}


def vfdb_tab_parsing(vfdb):
    vfdb_dict = {}

    with open(vfdb) as vfdb_tab:
        header = vfdb_tab.readline()
        for line in vfdb_tab:
            description = line.strip().split("\t")
            vfdb_id, vfdb_description = description[0], description[1]
            vfdb_dict[vfdb_id] = vfdb_description

    return vfdb_dict


def bpprc_toxins_dp_parsing(bpprc_toxins_db):
    bpprc_toxins_db_pd = pd.read_csv(
        bpprc_toxins_db,
        dtype={'name': str, 'partnerprotein': str, 'partnerprotein_textbox': str, 'target_order': str,
               'target_species': str, 'activity': str,	'taxonid': str, 'lc50': str, 'units': str,
               'percentage_mortality': str, 'publication': str, 'other_citations': str,	'life_stage': str,
               'instar': str, 'assay_material': str, 'assay_method': str, 'comment': str},
        sep='\t'
    )
    return bpprc_toxins_db_pd


def crop_factor_name(sseqid):
    end_ind = len(sseqid)
    for symb in list(sseqid)[::-1]:
        if symb.isnumeric():
            end_ind -= 1
        else:
            break

    return sseqid[:end_ind]


def add_vfdb_info(summary_pd, vfdb_dict, output_dict):
    for qseqid, sseqid, method, evalue, pident, bitscore in \
            zip(summary_pd['qseqid'], summary_pd['sseqid'], summary_pd['method'], summary_pd['evalue'],
                summary_pd['pident'], summary_pd['bitscore']):

        if sseqid in vfdb_dict:
            gene_id = f"({vfdb_dict[sseqid].split(') (')[1].split('[')[0]}"
            output_dict[len(output_dict) + 1] = {'qseqid': qseqid, 'sseqid': gene_id, 'method': method, 'evalue': evalue,
                                                 'pident': pident, 'bitscore': bitscore, 'database': 'VFDB',
                                                 'name': '-', 'target_order': '-',
                                                 'target_species': '-', 'activity': '-',
                                                 'percentage_mortality': '-', 'other': vfdb_dict[sseqid]}


def add_bpprc_toxicity_info(summary_pd, bpprc_toxins_db_pd, output_dict):
    for qseqid, sseqid, method, evalue, pident, bitscore in \
            zip(summary_pd['qseqid'], summary_pd['sseqid'], summary_pd['method'], summary_pd['evalue'],
                summary_pd['pident'], summary_pd['bitscore']):

        # Split factor name:
        factor_name_cropped = crop_factor_name(sseqid)
        # Subset database:
        bpprc_toxins_db_pd_subset = bpprc_toxins_db_pd[(bpprc_toxins_db_pd['name'] == sseqid) |
                                                       (bpprc_toxins_db_pd['name'] == factor_name_cropped)]

        for name, target_order, target_species, activity, percentage_mortality in zip(
                bpprc_toxins_db_pd_subset['name'], bpprc_toxins_db_pd_subset['target_order'],
                bpprc_toxins_db_pd_subset['target_species'], bpprc_toxins_db_pd_subset['activity'],
                bpprc_toxins_db_pd_subset['percentage_mortality']):
            output_dict[len(output_dict) + 1] = {'qseqid': qseqid, 'sseqid': sseqid, 'method': method, 'evalue': evalue,
                                                 'pident': pident, 'bitscore': bitscore, 'database': 'BPPRC',
                                                 'name': name, 'target_order': target_order,
                                                 'target_species': target_species, 'activity': activity,
                                                 'percentage_mortality': percentage_mortality, 'other': '-'}


def output_writing(outfile, output_pd):
    output_pd.to_csv(outfile, sep='\t', index=False, header=['Protein_ID', 'Factor', 'Method', 'E-value',
                                                             'Pident', 'Bitscore', 'Database', 'Name in Database',
                                                             'Target_order', 'Target_species', 'Activity',
                                                             'Percentage_mortality', 'Other'])


if __name__ == '__main__':
    summary_dict, output_dict = {}, {}
    # Input files parsing #
    cryprocessor_txt_parsing(snakemake.input.cryprocessor_report, summary_dict)
    mmseqs_tsv_parsing(snakemake.input.mmseqs_vfdb_filtered, summary_dict)
    bttoxin_digger_txt_parsing(snakemake.input.bttoxin_digger_report, summary_dict)
    # Summary: from dict to pandas df #
    summary_pd = pd.DataFrame.from_dict(summary_dict, orient='index')
    # Database parsing #
    bpprc_toxins_db_pd = bpprc_toxins_dp_parsing(snakemake.config['bpprc_tab'])
    vfdb_dict = vfdb_tab_parsing(snakemake.config['vfdb_tab'])
    # Add description from databases to hits#
    add_bpprc_toxicity_info(summary_pd, bpprc_toxins_db_pd, output_dict)
    add_vfdb_info(summary_pd, vfdb_dict, output_dict)
    # Output: from dict to pandas df #
    output_pd = pd.DataFrame.from_dict(output_dict, orient='index')
    # Output writing #
    output_writing(snakemake.output.pipeline_toxins_description, output_pd)

