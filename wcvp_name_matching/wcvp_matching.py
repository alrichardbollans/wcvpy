import pandas as pd

from wcvp_download import get_all_taxa, wcvp_columns, wcvp_accepted_columns
from wcvp_name_matching import clean_urn_ids, acc_info_col_names

status_priority = ['Accepted', 'Artificial Hybrid', 'Synonym', 'Illegitimate', 'Invalid', 'Local Biotype',
                   'Misapplied', 'Orthographic', 'Unplaced']


def id_lookup_wcvp(all_taxa: pd.DataFrame, given_id: str) -> pd.DataFrame:
    """
    Looks for id in list of taxa, returns a dictionary of accepted information
    :param all_taxa:
    :param given_id:
    :return:
    """

    clean_id = str(given_id)
    if "urn:lsid:ipni.org:names:" in clean_id:
        clean_id = clean_urn_ids(clean_id)
    record = all_taxa[all_taxa[wcvp_columns['id']] == clean_id]
    if len(record.index) == 0:
        print(f"Can't find id: {clean_id} in given wcvp taxa data")

    if len(record.index) > 1:
        print(f"Multiple id matches found in given wcvp data for id: {clean_id}")
        print(record[wcvp_columns['name']].unique())

    return record[[wcvp_columns['name']] + acc_info_col_names]


def get_family_specific_resolutions(resolution_df: pd.DataFrame, family_column: str = None) -> pd.DataFrame:
    if family_column is not None:
        match_df = resolution_df.loc[
            (resolution_df[wcvp_columns['family']] == resolution_df[family_column]) | (
                    resolution_df[wcvp_accepted_columns['family']] == resolution_df[family_column]) |
            (resolution_df[family_column].isna())]
    else:
        match_df = resolution_df

    return match_df


def get_wcvp_info_for_names_in_column(df: pd.DataFrame, matching_name_col: str, unique_submission_id_col: str,
                                      all_taxa: pd.DataFrame = None, family_column: str = None):
    """
    Appends accepted info columns to df from list of taxa, based on names in matching_name_col
    :param df:
    :param matching_name_col:
    :param all_taxa:
    :return:
    """
    if all_taxa is None:
        all_taxa = get_all_taxa()

    merged_with_wcvp = pd.merge(df, all_taxa, how='left', left_on=matching_name_col,
                                right_on=wcvp_columns['name'])

    match_df = get_family_specific_resolutions(merged_with_wcvp, family_column=family_column)
    match_df = match_df[[unique_submission_id_col] + acc_info_col_names]

    # Some items in WCVP have no accepted info, remove these
    match_df = match_df.dropna(subset=['accepted_name'])
    # Remove duplicates in match_df based on priority

    for r in match_df["taxon_status"].unique():
        if r not in status_priority:
            raise ValueError(f'Status priority list does not contain {r} and needs updating.')
    match_df['taxon_status'] = pd.Categorical(match_df['taxon_status'],
                                              status_priority)
    match_df.sort_values('taxon_status', inplace=True)

    match_df.drop_duplicates(subset=[unique_submission_id_col], inplace=True, keep='first')
    match_df['taxon_status'] = match_df['taxon_status'].astype(object)
    match_df['matched_by'] = 'direct_wcvp'
    return match_df
