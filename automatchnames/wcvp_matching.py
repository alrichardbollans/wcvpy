import numpy as np
import pandas as pd
from tqdm import tqdm

from automatchnames import clean_urn_ids
from taxa_lists import get_all_taxa


def _get_dict_from_wcvp_record(record: pd.DataFrame, taxa_list: pd.DataFrame) -> dict:
    """
    Formats a record from wcvp into a dictionary to integrate into other data
    :param record:
    :return:
    """
    taxonomic_status = record['taxonomic_status'].values[0]
    Accepted_Name = record['accepted_name'].values[0]
    Accepted_ID = record['accepted_kew_id'].values[0]
    if taxonomic_status == 'Accepted':
        accepted_taxon = record
    else:
        accepted_taxon = taxa_list[(taxa_list['taxon_name'] == Accepted_Name) & (taxa_list['kew_id'] == Accepted_ID)]
        if len(accepted_taxon.index) == 0:
            print('Warning: id given for a synonym whose corresponding accepted taxon is not given in given taxa list')
            print('Consider using larger set of taxa')
            return {'Accepted_Name': np.nan, 'Accepted_ID': np.nan, 'Accepted_Rank': np.nan,
                    'Accepted_Species': np.nan, 'Accepted_Species_ID': np.nan,
                    'taxonomic_status_of_submitted_name': taxonomic_status}
    Accepted_Rank = accepted_taxon['rank'].values[0]
    if Accepted_Rank in ["Synonym", "Homotypic_Synonym", "Species"]:
        Accepted_Species = Accepted_Name
        Accepted_Species_ID = Accepted_ID
    elif Accepted_Rank == "Genus":
        Accepted_Species = np.nan
        Accepted_Species_ID = np.nan

    else:
        # When subspecies and varieties are not accepted we need to find their parent

        Accepted_Species = accepted_taxon['parent_name'].values[0]
        Accepted_Species_ID = accepted_taxon['parent_kew_id'].values[0]

    return {'Accepted_Name': Accepted_Name, 'Accepted_ID': Accepted_ID, 'Accepted_Rank': Accepted_Rank,
            'Accepted_Species': Accepted_Species, 'Accepted_Species_ID': Accepted_Species_ID,
            'taxonomic_status_of_submitted_name': taxonomic_status}


def id_lookup_wcvp(all_taxa: pd.DataFrame, given_id: str) -> dict:
    """
    Looks for id in list of taxa, returns a dictionary of accepted information
    :param all_taxa:
    :param given_id:
    :return:
    """
    given_id = str(given_id)
    if "urn:lsid:ipni.org:names:" in given_id:
        given_id = clean_urn_ids(given_id)
    record = all_taxa[all_taxa['kew_id'] == given_id]
    nan_dict = {'Accepted_Name': np.nan, 'Accepted_ID': np.nan, 'Accepted_Rank': np.nan,
                'Accepted_Species': np.nan, 'Accepted_Species_ID': np.nan, 'taxonomic_status_of_submitted_name': np.nan}
    if len(record.index) == 0:
        print(f"Can't find id: {given_id} in given wcvp taxa data")
        return nan_dict
    if len(record.index) > 1:
        print(f"Multiple id matches found in given wcvp data for id: {given_id}")
        return nan_dict
    return _get_dict_from_wcvp_record(record, all_taxa)


def get_wcvp_info_for_names_in_column(df: pd.DataFrame, name_col: str, all_taxa: pd.DataFrame = None):
    """
    Appends accepted info columns to df from list of taxa, based on names in name_col
    :param df:
    :param name_col:
    :param all_taxa:
    :return:
    """
    if all_taxa is None:
        all_taxa = get_all_taxa()

    taxa_in_df = all_taxa[all_taxa['taxon_name'].isin(df[name_col])]

    # First get df of wcvp matches for each sample
    dict_for_matches = {name_col: [], 'Accepted_Name': [], 'Accepted_ID': [], 'Accepted_Rank': [],
                        'Accepted_Species': [], 'Accepted_Species_ID': [], 'taxonomic_status_of_submitted_name': []}
    for sample in df[name_col].values:

        matching_records = taxa_in_df[taxa_in_df['taxon_name'] == sample]
        for record_id in matching_records['kew_id'].values:
            dict_for_matches[name_col].append(sample)
            dict_for_record = id_lookup_wcvp(all_taxa, record_id)

            for k in dict_for_record:
                dict_for_matches[k].append(dict_for_record[k])

    match_df = pd.DataFrame(dict_for_matches)

    # Remove duplicates in match_df based on priority
    status_priority = ["Accepted", "Synonym", "Homotypic_Synonym"]
    for r in match_df["taxonomic_status_of_submitted_name"].unique():
        if r not in status_priority:
            raise ValueError(f'Status priority list does not contain {r} and needs updating.')
    match_df['taxonomic_status_of_submitted_name'] = pd.Categorical(match_df['taxonomic_status_of_submitted_name'],
                                                                    status_priority)
    match_df.sort_values('taxonomic_status_of_submitted_name', inplace=True)

    match_df.drop_duplicates(subset=[name_col], inplace=True, keep='first')

    return match_df
