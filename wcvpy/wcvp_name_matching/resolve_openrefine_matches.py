from typing import List

import pandas as pd

from wcvpy.OpenRefineMatching import reco_submitted_name_col_id
from wcvpy.wcvp_download import wcvp_accepted_columns, wcvp_columns
from wcvpy.wcvp_name_matching import resolve_matches_by_priorities, get_accepted_wcvp_info_from_ipni_ids_in_column


def resolve_openrefine_to_best_matches(reco_df: pd.DataFrame, all_taxa: pd.DataFrame, families_of_interest: List[str] = None):
    # There shouldn't be any repeated reco_ids for the same submitted names so check this first
    problems = reco_df[reco_df.duplicated(subset=['reco_id', reco_submitted_name_col_id], keep=False)]

    if len(problems.index) > 0:
        raise ValueError(
            f'Repeated conflicting matches in reconciled data: {problems[reco_submitted_name_col_id]}')

    out_df = get_accepted_wcvp_info_from_ipni_ids_in_column(reco_df, 'reco_id', all_taxa)

    # with info from ipni ids
    out_df = out_df[~out_df[wcvp_columns['status']].isna()]
    #
    # out_df.to_csv('reco_find_with_info.csv')
    # Within family (and synonym families)
    if families_of_interest is not None:
        # Begin by removing incorrectly matched families
        out_df = out_df.loc[
            (out_df[wcvp_columns['family']].isin(families_of_interest)) | (
                out_df[wcvp_accepted_columns['family']].isin(families_of_interest))]

    # Unique matches
    unique_matches = out_df.drop_duplicates(subset=[reco_submitted_name_col_id], keep=False).copy()
    unique_matches['matched_by'] = 'openrefine_unique'

    non_unique_reco_matches = out_df[out_df[reco_submitted_name_col_id].duplicated(keep=False)]

    # Get matches where accepted name is unique
    unique_acc_names = non_unique_reco_matches.drop_duplicates(
        subset=[reco_submitted_name_col_id, wcvp_accepted_columns['ipni_id']], keep='first')
    submitted_names_with_single_accepted_match = unique_acc_names.drop_duplicates(
        subset=[reco_submitted_name_col_id],
        keep=False).copy()
    submitted_names_with_single_accepted_match['matched_by'] = 'openrefine_unique_accepted_name'

    non_unique_matches = non_unique_reco_matches[
        ~non_unique_reco_matches[reco_submitted_name_col_id].isin(
            submitted_names_with_single_accepted_match[reco_submitted_name_col_id].values)].copy()

    # Best scoring matches
    non_unique_matches['reco_score_max'] = non_unique_matches.groupby([reco_submitted_name_col_id])[
        'reco_score'].transform('max')

    top_scorers = non_unique_matches[non_unique_matches['reco_score'] == 'reco_score_max']
    top_unique_scorers = top_scorers[~top_scorers[reco_submitted_name_col_id].duplicated(keep=False)]
    top_unique_scorers['matched_by'] = 'openrefine_unique_top_score'

    unmatched = non_unique_matches[~non_unique_matches[reco_submitted_name_col_id].isin(
        top_unique_scorers[reco_submitted_name_col_id].values)]

    # Matches based on priority
    best_priority_matches = resolve_matches_by_priorities(unmatched, reco_submitted_name_col_id,
                                                          ['status', 'rank'])
    best_priority_matches['matched_by'] = 'openrefine_best_priority'

    resolved_matches = pd.concat(
        [unique_matches, submitted_names_with_single_accepted_match, top_unique_scorers,
         best_priority_matches])
    resolved_matches = resolved_matches.drop(columns=['reco_score_max'])

    return resolved_matches
