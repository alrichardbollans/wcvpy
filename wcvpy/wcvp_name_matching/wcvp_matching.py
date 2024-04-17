from typing import List

import numpy as np
import pandas as pd

from wcvpy.wcvp_download import get_all_taxa, wcvp_columns, wcvp_accepted_columns, clean_whitespaces_in_names
from wcvpy.wcvp_name_matching import clean_urn_ids, output_record_col_names, lowercase_name_col, \
    remove_fullstop, tidied_taxon_authors_col, tidy_authors, \
    status_priority


def lookup_ipni_id_in_wcvp(all_taxa: pd.DataFrame, given_id: str) -> pd.DataFrame:
    """
    Looks for id in list of taxa, returns a dictionary of accepted information
    :param all_taxa:
    :param given_id:
    :return:
    """

    clean_id = str(given_id)
    if "urn:lsid:ipni.org:names:" in clean_id:
        clean_id = clean_urn_ids(clean_id)
    record = all_taxa[all_taxa[wcvp_columns['ipni_id']] == clean_id]
    if len(record.index) == 0:
        print(f"Can't find id: {clean_id} in given wcvp taxa data")

    if len(record.index) > 1:
        print(f"Multiple id matches found in given wcvp data for id: {clean_id}")
        print(record[wcvp_columns['name']].unique())

    return record[[wcvp_columns['name']] + output_record_col_names]


def get_family_specific_resolutions(resolution_df: pd.DataFrame, family_column: str = None) -> pd.DataFrame:
    if family_column is not None:
        match_df = resolution_df.loc[
            (resolution_df[wcvp_columns['family']] == resolution_df[family_column]) | (
                    resolution_df[wcvp_accepted_columns['family']] == resolution_df[family_column]) |
            (resolution_df[family_column].isna())]
    else:
        match_df = resolution_df

    return match_df


def tidy_value_for_matching(given_value: str) -> str:
    lower = given_value.lower()
    without_fullstop = remove_fullstop(lower)
    rmved_whitespace = clean_whitespaces_in_names(without_fullstop)
    return rmved_whitespace


def match_name_to_concatenated_columns(df: pd.DataFrame, matching_name_col: str, all_taxa: pd.DataFrame,
                                       columns: List[str]):
    # Remove taxa without author information as these are matched later without authors
    taxa_to_use = all_taxa.dropna(subset=columns).copy()
    column_series = [all_taxa[c].fillna('') for c in columns]
    taxa_to_use['taxon_name_with_extra_columns'] = taxa_to_use[wcvp_columns['name']].str.cat(column_series,
                                                                                             sep=' ')
    taxa_to_use['taxon_name_with_extra_columns'] = taxa_to_use['taxon_name_with_extra_columns'].apply(
        tidy_value_for_matching)

    df[lowercase_name_col] = df[matching_name_col].apply(tidy_value_for_matching)
    # Match with taxon authors
    author_merged = pd.merge(df, taxa_to_use, how='left', left_on=lowercase_name_col,
                             right_on='taxon_name_with_extra_columns')
    author_merged = author_merged.dropna(subset=[wcvp_columns['wcvp_id']])

    unmatched_with_authors_df = df[~df[lowercase_name_col].isin(author_merged[lowercase_name_col].values)].copy()

    # Repeat but with 'tidied' authors
    unmatched_with_authors_df[tidied_taxon_authors_col] = unmatched_with_authors_df[matching_name_col].apply(
        tidy_authors)
    unmatched_with_authors_df[tidied_taxon_authors_col] = unmatched_with_authors_df[
        tidied_taxon_authors_col].apply(tidy_value_for_matching)

    tidy_author_merged = pd.merge(unmatched_with_authors_df, taxa_to_use, how='left',
                                  left_on=tidied_taxon_authors_col,
                                  right_on='taxon_name_with_extra_columns')
    tidy_author_merged = tidy_author_merged.dropna(subset=[wcvp_columns['wcvp_id']])

    matched = pd.concat([author_merged, tidy_author_merged], ignore_index=True)

    matched['matched_name'] = matched[wcvp_columns['name']].str.cat([matched[c].fillna('') for c in columns],
                                                                    sep=' ')

    unmatched = unmatched_with_authors_df[
        ~unmatched_with_authors_df[tidied_taxon_authors_col].isin(
            tidy_author_merged[tidied_taxon_authors_col].values)]

    return matched, unmatched


def get_wcvp_info_for_names_in_column(df: pd.DataFrame, matching_name_col: str, unique_submission_id_col: str,
                                      all_taxa: pd.DataFrame = None, family_column: str = None, wcvp_version: str = None):
    """
    Appends accepted info columns to df from list of taxa, based on names in matching_name_col
    :param df:
    :param matching_name_col:
    :param unique_submission_id_col:
    :param all_taxa:
    :param family_column:
    :param wcvp_version:
    :return:
    """
    if all_taxa is None:
        all_taxa = get_all_taxa(version=wcvp_version)

    # First try with author info i.e. taxon name + taxon_authors and then
    # taxon name + parenthetical_author + primary_author then taxon name + primary author
    author_merged, unmatched_with_authors_df = match_name_to_concatenated_columns(df, matching_name_col,
                                                                                  all_taxa,
                                                                                  [wcvp_columns['authors']])

    author_merged['matched_by'] = 'direct_wcvp_w_author'

    # Try with paranthetical and primary author columns - slight difference with taxon authors as no parantheses around paranthet author
    paranthet_author_merged, unmatched_with_paranthet_authors_df = match_name_to_concatenated_columns(
        unmatched_with_authors_df, matching_name_col,
        all_taxa,
        [wcvp_columns['paranthet_author'], wcvp_columns['primary_author']])

    paranthet_author_merged['matched_by'] = 'direct_wcvp_w_author'

    # Try with primary author columns
    primary_author_merged, unmatched_with_primary_author_df = match_name_to_concatenated_columns(
        unmatched_with_paranthet_authors_df, matching_name_col,
        all_taxa,
        [wcvp_columns['primary_author']])

    primary_author_merged['matched_by'] = 'direct_wcvp_w_author'

    # Match with just name
    all_taxa['tidied_taxon_name'] = all_taxa[wcvp_columns['name']].apply(tidy_value_for_matching)
    unmatched_with_primary_author_df[lowercase_name_col] = unmatched_with_primary_author_df[
        matching_name_col].apply(tidy_value_for_matching)
    just_name_merged = pd.merge(unmatched_with_primary_author_df, all_taxa, how='left',
                                left_on=lowercase_name_col,
                                right_on='tidied_taxon_name')
    just_name_merged = just_name_merged.dropna(subset=[wcvp_columns['wcvp_id']])
    just_name_merged['matched_by'] = 'direct_wcvp'
    just_name_merged['matched_name'] = just_name_merged[wcvp_columns['name']]

    merged_with_wcvp = pd.concat([author_merged, paranthet_author_merged, primary_author_merged, just_name_merged])
    match_df = get_family_specific_resolutions(merged_with_wcvp, family_column=family_column)
    match_df = match_df[[unique_submission_id_col] + output_record_col_names + ['matched_by', 'matched_name']].copy()

    # Remove duplicates in match_df based on priority

    for r in match_df["taxon_status"].unique():
        if r not in status_priority:
            raise ValueError(f'Status priority list does not contain {r} and needs updating.')
    match_df['taxon_status'] = pd.Categorical(match_df['taxon_status'],
                                              status_priority)
    match_df = match_df.sort_values('taxon_status')
    # Appropriately label unique matches
    match_df['matched_by'] = np.where(match_df.duplicated(keep=False, subset=[unique_submission_id_col]),
                                      match_df['matched_by'],
                                      match_df['matched_by'] + '_unique')
    match_df = match_df.drop_duplicates(subset=[unique_submission_id_col], keep='first')
    match_df['taxon_status'] = match_df['taxon_status'].astype(object)
    return match_df
