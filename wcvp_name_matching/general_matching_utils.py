from typing import List

import numpy as np
import pandas as pd

from wcvp_download import wcvp_columns, wcvp_accepted_columns
from wcvp_name_matching import output_record_col_names

status_priority = ['Accepted', 'Artificial Hybrid', 'Synonym', 'Illegitimate', 'Invalid', 'Local Biotype',
                   'Misapplied', 'Orthographic', 'Unplaced']
rank_priority = ["nothof.", "Form", "Subspecies", "Subvariety", "Variety", "Species", "Genus"]


def resolve_matches_by_priorities(match_df: pd.DataFrame, submission_id_col: str,
                                  priority_order: List[str]) -> pd.DataFrame:
    out_df = match_df.copy(deep=True)
    # priority order should be ['rank', 'status'] or ['status', 'rank']
    for p in priority_order:
        if p == 'status':
            specific_order = status_priority
        elif p == 'rank':
            specific_order = rank_priority
        else:
            raise ValueError("priority order should be ['rank', 'status'] or ['status', 'rank']")
        out_df[wcvp_columns[p]] = pd.Categorical(out_df[wcvp_columns[p]], specific_order)

    out_df.sort_values([wcvp_columns[x] for x in priority_order], inplace=True)
    out_df.drop_duplicates(subset=[submission_id_col], inplace=True,
                           keep='first')

    return out_df


# def resolve_matches_by_status_priority(match_df: pd.DataFrame, submission_id_col: str) -> pd.DataFrame:
#     # Remove duplicate matches with worse status
#     for r in match_df["taxon_status"].unique():
#         if r not in status_priority:
#             raise ValueError(f'Status priority list does not contain {r} and needs updating.')
#     match_df[wcvp_columns['status']] = pd.Categorical(match_df[wcvp_columns['status']],
#                                                       status_priority)
#     match_df.sort_values(wcvp_columns['status'], inplace=True)
#     # Drop duplicated submissions of the same rank with worse taxonomic status
#     match_df.drop_duplicates(subset=[submission_id_col, wcvp_accepted_columns['rank']], inplace=True,
#                              keep='first')
#     match_df[wcvp_columns['status']] = match_df[wcvp_columns['status']].astype(object)
#     return match_df
#
#
# def resolve_matches_by_rank_priority(match_df: pd.DataFrame, submission_id_col: str) -> pd.DataFrame:
#     # Remove duplicate matches with worse specificity
#     for r in match_df[wcvp_accepted_columns['rank']].unique():
#         if r not in rank_priority:
#             print(match_df[match_df[wcvp_accepted_columns['rank']] == r][wcvp_columns['name']].values)
#             raise ValueError(f'Rank priority list does not contain {r} and needs updating.')
#     match_df[wcvp_accepted_columns['rank']] = pd.Categorical(match_df[wcvp_accepted_columns['rank']],
#                                                              rank_priority)
#     match_df.sort_values(wcvp_accepted_columns['rank'], inplace=True)
#
#     # Drop duplicated submissions of the same status with higher rank
#     match_df.drop_duplicates(subset=[submission_id_col, wcvp_columns['status']], keep='first', inplace=True)
#     match_df[wcvp_accepted_columns['rank']] = match_df[wcvp_accepted_columns['rank']].astype(object)
#
#     return match_df


def get_accepted_wcvp_info_from_ipni_ids_in_column(df: pd.DataFrame, id_col_name: str,
                                                   all_taxa: pd.DataFrame) -> \
        pd.DataFrame:
    """
    Appends accepted info columns to df from list of taxa, based on ids in id_col_name
    :param all_taxa:
    :param df:
    :param id_col_name:
    :return:
    """

    in_df = df.copy(deep=True)
    in_df[id_col_name] = in_df[id_col_name].fillna('no_id_placeholder')
    match_df = pd.merge(in_df,
                        all_taxa[[wcvp_columns['ipni_id'], wcvp_columns['family'],
                                  wcvp_columns['rank']] + output_record_col_names],
                        how='left',
                        left_on=id_col_name,
                        right_on=wcvp_columns['ipni_id'])
    if wcvp_columns['ipni_id'] not in in_df.columns:
        match_df = match_df.drop(columns=[wcvp_columns['ipni_id']])

    if len(match_df.index) != len(in_df.index):
        raise ValueError('Generating accepted info is mismatched')

    match_df[id_col_name] = match_df[id_col_name].replace('no_id_placeholder', np.nan)
    match_df['matched_by'] = 'ipni_id'

    return match_df
