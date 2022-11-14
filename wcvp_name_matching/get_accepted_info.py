import hashlib
import os

import numpy as np
import pandas as pd
# Add progress bar to apply method
from tqdm import tqdm

tqdm.pandas()

from typing import List
from pkg_resources import resource_filename

from wcvp_name_matching import get_wcvp_info_for_names_in_column, \
    get_knms_name_matches, id_lookup_wcvp, clean_urn_ids, acc_info_col_names, temp_outputs_dir, \
    tidy_names_in_column, \
    hybrid_characters, recapitalised_name_col, submitted_name_col_id
from wcvp_download import get_all_taxa, wcvp_columns

matching_data_path = resource_filename(__name__, 'matching data')

rank_priority = ["Subspecies", "Subvariety", "Variety", "Species", "Genus"]
status_priority = ["Accepted", "Synonym"]


def _temp_output(df: pd.DataFrame, tag: str, warning: str = None):
    df_str = df.to_string()
    str_to_hash = str(df_str).encode()
    temp_basename_csv = str(hashlib.md5(str_to_hash).hexdigest()) + ".csv"
    try:
        os.mkdir(temp_outputs_dir)
    except FileExistsError as error:
        pass
    outfile = os.path.join(temp_outputs_dir, tag + temp_basename_csv)
    if warning is not None:
        print(f'{warning}. Check tempfile: {outfile}.')
    else:
        print(f'Temp file for this run: {outfile}')
    df.to_csv(outfile)


def _autoresolve_missing_matches(unmatched_submissions_df: pd.DataFrame, matching_name_col: str,
                                 name_id_col: str,
                                 all_taxa: pd.DataFrame,
                                 family_column: str = None) -> pd.DataFrame:
    """

    :param unmatched_submissions_df:
    :return:
    """

    if len(unmatched_submissions_df.index) > 0:
        _temp_output(unmatched_submissions_df, 'unmatched_to_autoresolve',
                     "Resolving submitted names which weren't initially matched using KNMS.")
        print(
            "This may take some time... This can be sped up by specifying families of interest (if you "
            "haven't already done so) or checking the temp file for misspelled submissions.")

        # For each submission, check if any name is contained in the name, then take the lowest rank of
        # matches
        # Create a dict of submissions with possible matches
        match_df = pd.DataFrame(
            {matching_name_col: [], wcvp_columns['name']: [], 'accepted_name': [], 'accepted_ipni_id': [],
             'accepted_rank': [],
             'accepted_parent': [], 'accepted_parent_ipni_id': [],
             'accepted_species': [], 'accepted_species_ipni_id': [],
             'accepted_family': [],
             'taxon_status': []})

        # Get more precise list of taxa which possibly matches submissions
        wcvp_name_containment_df = all_taxa[
            all_taxa.progress_apply(
                lambda x: any(
                    x[wcvp_columns['name']] in y for y in unmatched_submissions_df[matching_name_col].values),
                axis=1)]

        for i in tqdm(range(len(unmatched_submissions_df[matching_name_col].values)),
                      desc="Searching automated matches",
                      ascii=False, ncols=72):
            s = unmatched_submissions_df[matching_name_col].values[i]
            for index, row in wcvp_name_containment_df.iterrows():
                record = row.to_frame().transpose()[acc_info_col_names + [wcvp_columns['name']]]

                taxa = row[wcvp_columns['name']]
                if taxa in s:
                    record[matching_name_col] = [s]
                    match_df = pd.concat([match_df, record])

        match_df = match_df.dropna(subset=['accepted_name'])

        # Remove genera matches where submitted name has more than one word.
        # This to avoid matching mispelt species to genera
        if len(match_df.index) > 0:
            if family_column is None:
                unique_family_genera_pairs = all_taxa[
                    [wcvp_columns['family'], wcvp_columns['genus']]].drop_duplicates(keep='first')
                genera_unique_to_family = unique_family_genera_pairs[wcvp_columns['genus']].unique()

                match_df = match_df[
                    (~((match_df['accepted_rank'] == 'Genus') & match_df[matching_name_col].str.contains(
                        " "))) | (match_df[wcvp_columns['name']].isin(genera_unique_to_family))]

        # Remove duplicate matches with worse priority

        for r in match_df["taxon_status"].unique():
            if r not in status_priority:
                raise ValueError(f'Status priority list does not contain {r} and needs updating.')
        match_df['taxon_status'] = pd.Categorical(match_df['taxon_status'],
                                                  status_priority)
        match_df.sort_values('taxon_status', inplace=True)
        # Drop duplicated submissions of the same rank with worse taxonomic status
        match_df.drop_duplicates(subset=[matching_name_col, 'accepted_rank'], inplace=True, keep='first')
        match_df['taxon_status'] = match_df['taxon_status'].astype(object)
        # Remove duplicate matches with worse specificity
        for r in match_df['accepted_rank'].unique():
            if r not in rank_priority:
                raise ValueError(f'Rank priority list does not contain {r} and needs updating.')
        match_df['accepted_rank'] = pd.Categorical(match_df['accepted_rank'], rank_priority)
        match_df.sort_values('accepted_rank', inplace=True)

        # Get the most precise match by dropping duplicate submissions
        match_df.drop_duplicates(subset=[matching_name_col], keep='first', inplace=True)
        match_df['accepted_rank'] = match_df['accepted_rank'].astype(object)
        # Merge with original data
        matches = pd.merge(unmatched_submissions_df[[name_id_col, matching_name_col]], match_df,
                           on=matching_name_col,
                           sort=False)
        matches = matches.dropna(subset=['accepted_name'])
        return matches
    else:
        return unmatched_submissions_df


def _get_knms_matches_and_accepted_info_from_names_in_column(df: pd.DataFrame, matching_name_col: str,
                                                             all_taxa: pd.DataFrame) -> pd.DataFrame:
    """
    Matches names in df using knms and gets corresponding accepted info from wcvp.
    :param df:
    :param matching_name_col:
    :return:
    """
    if len(df.index) > 0:
        match_records = get_knms_name_matches(df[matching_name_col].values)
        match_records['ipni_id'] = match_records['ipni_id'].apply(clean_urn_ids)

        single_matches = match_records[match_records['match_state'] == 'true'].copy()
        single_matches_with_info = get_accepted_info_from_ids_in_column(single_matches, 'ipni_id', all_taxa)

        single_matches_with_info['matched_by'] = 'knms_single'

        multiple_matches = match_records[match_records['match_state'] == 'multiple_matches'].copy()
        if len(multiple_matches.index) > 0:
            best_matches = _find_best_matches_from_multiples(multiple_matches, all_taxa=all_taxa)

            resolved_df = pd.concat([single_matches_with_info, best_matches], axis=0)
        else:
            resolved_df = single_matches_with_info
        resolved_df = resolved_df.rename(columns={'submitted': matching_name_col})

        resolved_df = resolved_df.drop(columns=['match_state', 'ipni_id', 'matched_name'])

        # This trick allows merging on columns with duplicates and matches duplicates rather than repeating
        # them
        df['cc'] = df.groupby(matching_name_col).cumcount()
        resolved_df['cc'] = resolved_df.groupby(matching_name_col).cumcount()
        merged_df = df.merge(resolved_df, how='outer').drop(columns='cc', axis=1)
        merged_df = merged_df.dropna(subset=['accepted_name'])
        return merged_df
    else:
        return df


def _find_best_matches_from_multiples(multiple_match_records: pd.DataFrame,
                                      all_taxa: pd.DataFrame = None) -> pd.DataFrame:
    """
    Gets best matching record for each set of multiple matches, and returns a dataframe of single matches
    :param multiple_match_records:
    :return:
    """
    if all_taxa is None:
        all_taxa = get_all_taxa()
    # First find accepted info for the multiple matches
    multiple_match_records['ipni_id'] = multiple_match_records['ipni_id'].apply(clean_urn_ids)
    multiple_matches_with_accepted_ids = get_accepted_info_from_ids_in_column(multiple_match_records,
                                                                              'ipni_id', all_taxa)

    # First, use matches where accepted name is the same as the submitted name
    accepted_names_matching_submitted_names = multiple_matches_with_accepted_ids[
        multiple_matches_with_accepted_ids['submitted'] == multiple_matches_with_accepted_ids[
            'accepted_name']]

    unmatched_containment_df = multiple_matches_with_accepted_ids[
        ~multiple_matches_with_accepted_ids['submitted'].isin(
            accepted_names_matching_submitted_names['submitted'].values)]

    # reduce list to remove essentially repeated matches
    unique_accepted_matches = unmatched_containment_df.drop_duplicates(
        subset=['submitted', 'accepted_ipni_id'], keep='first')

    #  Next use matches where the submitted name has a unique match
    submitted_names_with_single_accepted_match = unique_accepted_matches.drop_duplicates(subset=['submitted'],
                                                                                         keep=False)

    # Where neither of the above apply, sort by accepted rank
    unresolved_submissions = unique_accepted_matches[
        ~unique_accepted_matches["submitted"].isin(submitted_names_with_single_accepted_match["submitted"])]

    # matches where the accepted name is contained in the submitted name
    # In case of duplicates, these are sorted by specifity of rank
    unresolved_submissions = unresolved_submissions.dropna(subset=['accepted_name'])
    accepted_names_in_submitted_names = unresolved_submissions[
        unresolved_submissions.apply(lambda x: x['accepted_name'] in x['submitted'],
                                     axis=1)]

    for r in accepted_names_in_submitted_names['accepted_rank'].unique():
        if r not in rank_priority:
            raise ValueError(f'Rank priority list does not contain {r} and needs updating.')
    accepted_names_in_submitted_names['accepted_rank'] = pd.Categorical(
        accepted_names_in_submitted_names['accepted_rank'], rank_priority)
    accepted_names_in_submitted_names.sort_values('accepted_rank', inplace=True)

    # Get the most precise match by dropping duplicate submissions
    accepted_names_in_submitted_names.drop_duplicates(subset=['submitted'], keep='first', inplace=True)
    accepted_names_in_submitted_names['accepted_rank'] = accepted_names_in_submitted_names[
        'accepted_rank'].astype(object)

    accepted_names_matching_submitted_names['matched_by'] = 'knms_multiple_1'
    submitted_names_with_single_accepted_match['matched_by'] = 'knms_multiple_2'
    accepted_names_in_submitted_names['matched_by'] = 'knms_multiple_3'
    matches_to_use = pd.concat([accepted_names_matching_submitted_names,
                                submitted_names_with_single_accepted_match,
                                accepted_names_in_submitted_names])
    matches_to_use.drop_duplicates(subset=['submitted'], keep='first', inplace=True)

    unmatched_df = multiple_matches_with_accepted_ids[
        ~multiple_matches_with_accepted_ids['submitted'].isin(matches_to_use['submitted'].values)]
    if len(unmatched_df.index) > 0:
        _temp_output(unmatched_df, 'unmatched_samples_with_multiple_knms_hits',
                     'Warning: Some samples have multiple matches in knms but havent been automatically '
                     'resolved. These may still be resolved in later steps.')

    return matches_to_use


def get_accepted_info_from_ids_in_column(df: pd.DataFrame, id_col_name: str, all_taxa: pd.DataFrame) -> \
        pd.DataFrame:
    """
    Appends accepted info columns to df from list of taxa, based on ids in id_col_name
    :param all_taxa:
    :param df:
    :param id_col_name:
    :return:
    """

    match_df = pd.merge(df, all_taxa[[wcvp_columns['id']] + acc_info_col_names], how='left',
                        left_on=id_col_name,
                        right_on=wcvp_columns['id'])
    if wcvp_columns['id'] not in df.columns:
        match_df = match_df.drop(columns=[wcvp_columns['id']])

    if len(match_df.index) != len(df.index):
        raise ValueError('Generating accepted info is mismatched')

    return match_df


def get_accepted_info_from_names_in_column(in_df: pd.DataFrame, name_col: str,
                                           families_of_interest: List[str] = None,
                                           family_column: str = None,
                                           manual_resolution_csv: str = None,
                                           match_level: str = 'full') -> pd.DataFrame:
    """
    First tries to match names in df to wcvp directly to obtain accepted info and then
    matches names in df using knms and gets corresponding accepted info from wcvp
    :param family_column:
    :param in_df:
    :param families_of_interest:
    :param manual_resolution_csv:
    :param match_level:
    :param name_col:
    :return:
    """
    match_levels = ['full', 'weak', 'knms']
    if match_level not in match_levels:
        raise ValueError(f'match_level should be one of {match_levels}')
    problem_columns = [x for x in in_df.columns if
                       x in [submitted_name_col_id, recapitalised_name_col] + list(wcvp_columns.keys())]
    if len(problem_columns) > 0:
        raise ValueError(
            f'Column names used in input data will be confused in matching process: {problem_columns}')

    if len(in_df.index) > 0:
        if family_column is not None and families_of_interest is not None:
            fams_in_family_column = in_df[family_column].unique().dropna()
            for f in fams_in_family_column:
                if f not in families_of_interest:
                    raise ValueError(
                        f'Family {f} given in family column but not in families of interest')
        if family_column is not None and families_of_interest is None:
            families_of_interest = in_df[family_column].unique().dropna()
        all_taxa = get_all_taxa(families_of_interest=families_of_interest)

        df = in_df.copy(deep=True)

        # Standardise input names
        # Non standard captialisation causes issues
        # TODO: fix capitalisation by parsing. Then add epithet columns
        # TODO: Adding suggested epithets and ranks will help with later parsing

        df = df.drop_duplicates(subset=[name_col])
        df = df.dropna(subset=[name_col])
        tidy_names_in_column(df, name_col)

        # First get manual matches
        if manual_resolution_csv is not None:
            manual_match_df = pd.read_csv(manual_resolution_csv)
            manual_match_df = manual_match_df[
                manual_match_df['submitted'].isin(df[submitted_name_col_id].values.tolist())]
            man_matches_with_accepted_info = get_accepted_info_from_ids_in_column(manual_match_df,
                                                                                  'resolution_id', all_taxa)
            man_matches_with_accepted_info = man_matches_with_accepted_info.dropna(subset=['accepted_name'])
            manual_matches = pd.merge(df, man_matches_with_accepted_info, left_on=submitted_name_col_id,
                                      right_on='submitted',
                                      sort=False)
            manual_matches['matched_by'] = 'manual'
            unmatched_manual_df = df[
                ~df[submitted_name_col_id].isin(manual_matches[submitted_name_col_id].values)]
        else:
            manual_matches = pd.DataFrame()
            unmatched_manual_df = df

        # Then match with exact matches in wcvp
        wcvp_exact_name_match_df = get_wcvp_info_for_names_in_column(unmatched_manual_df,
                                                                     recapitalised_name_col,
                                                                     submitted_name_col_id,
                                                                     all_taxa=all_taxa)
        wcvp_exact_name_match_df['matched_by'] = 'direct_wcvp'
        wcvp_resolved_df = pd.concat([wcvp_exact_name_match_df, manual_matches], axis=0)
        unmatched_name_df = df[
            ~df[submitted_name_col_id].isin(wcvp_resolved_df[submitted_name_col_id].values)]

        if match_level in ['full', 'knms']:
            # If exact matches aren't found in wcvp, use knms
            matches_with_knms = _get_knms_matches_and_accepted_info_from_names_in_column(unmatched_name_df,
                                                                                         recapitalised_name_col,
                                                                                         all_taxa)
            knms_resolved_df = pd.concat([wcvp_resolved_df, matches_with_knms], axis=0)
            unmatched_df = df[~df[submitted_name_col_id].isin(knms_resolved_df[submitted_name_col_id].values)]

            if match_level == 'full':
                # Get autoresolved matches
                unmatched_resolutions = _autoresolve_missing_matches(unmatched_df, recapitalised_name_col,
                                                                     submitted_name_col_id,
                                                                     all_taxa,
                                                                     family_column=family_column)

                unmatched_resolutions['matched_by'] = 'autoresolution'

                final_resolved_df = pd.concat([unmatched_resolutions, knms_resolved_df], axis=0)

            else:
                final_resolved_df = knms_resolved_df
        else:
            final_resolved_df = wcvp_resolved_df

        # Provide temp outputs
        unmatched_final_df = df[
            ~df[submitted_name_col_id].isin(final_resolved_df[submitted_name_col_id].values)]
        if len(unmatched_final_df.index) > 0:
            _temp_output(unmatched_final_df, 'unmatched_samples',
                         'WARNING: some submissions have not been resolved and must be manually resolved. '
                         'Consider fixing names in your original data.')
            final_resolved_df = pd.concat([final_resolved_df, unmatched_final_df])

        _temp_output(final_resolved_df, 'final_resolutions')

        final_resolved_df = final_resolved_df[[submitted_name_col_id] + acc_info_col_names + ['matched_by']]
        out_df = pd.merge(in_df, final_resolved_df, left_on=name_col, right_on=submitted_name_col_id,
                          how='left')
        # Rename columns such that name_col column is original submitted names (without cleaning)
        out_df = out_df.drop(columns=[name_col])
        out_df = out_df.rename(columns={submitted_name_col_id: name_col})
        # reorder
        out_df = out_df[in_df.columns.tolist() + acc_info_col_names + ['matched_by']]
        return out_df
    else:
        out_copy = in_df.copy()
        for a in acc_info_col_names + ['matched_by']:
            out_copy[a] = np.nan
        return out_copy


def post_checks():
    #
    pass


if __name__ == '__main__':
    taxa = get_all_taxa()
    sarc = taxa[taxa['genus'] == 'Sarcorhiza']
    sarc.to_csv('test.csv')