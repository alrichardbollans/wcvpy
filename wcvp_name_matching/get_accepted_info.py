import hashlib
import os
from itertools import combinations

import numpy as np
import pandas as pd
# Add progress bar to apply method
from tqdm import tqdm

tqdm.pandas()

from typing import List
from pkg_resources import resource_filename

from wcvp_name_matching import get_wcvp_info_for_names_in_column, \
    get_knms_name_matches, clean_urn_ids, output_record_col_names, temp_outputs_dir, \
    tidy_names_in_column, recapitalised_name_col, submitted_name_col_id, \
    tidy_families_in_column, status_priority, submitted_family_name_col_id, unique_submission_index_col, \
    lowercase_name_col, tidied_taxon_authors_col, get_word_combinations, \
    remove_whitespace_at_beginning_and_end
from wcvp_download import get_all_taxa, wcvp_columns, wcvp_accepted_columns

matching_data_path = resource_filename(__name__, 'matching data')

rank_priority = ["nothof.", "Form", "Subspecies", "Subvariety", "Variety", "Species", "Genus"]


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
        print(f'{warning}.')
        print(f'Check tempfile: {outfile}.')
    else:
        print(f'Temp file for this run:')
        print(f'{outfile}')
    df.to_csv(outfile)


def _autoresolve_missing_matches(unmatched_submissions_df: pd.DataFrame, matching_name_col: str,
                                 submission_id_col: str,
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
        # Create a df of submissions with possible matches
        match_df = pd.DataFrame()

        if family_column is not None:
            match_df[family_column] = []

        # Get more precise list of taxa which possibly matches submissions
        substrings_of_unmatched_names = []
        for sub_name in unmatched_submissions_df[matching_name_col].values:
            substrings_of_unmatched_names += [sub_name[x:y] for x, y in
                                              combinations(range(len(sub_name) + 1), r=2)]

        substrings_of_unmatched_names = list(set(substrings_of_unmatched_names))

        wcvp_name_containment_df = all_taxa[
            all_taxa[wcvp_columns['name']].isin(substrings_of_unmatched_names)]

        wcvp_names = wcvp_name_containment_df[wcvp_columns['name']].tolist()
        mask = unmatched_submissions_df[matching_name_col].str.contains(
            '|'.join(wcvp_names))
        unmatched_submissions_worth_checking = unmatched_submissions_df[mask]

        for i in tqdm(range(len(unmatched_submissions_worth_checking[matching_name_col].values)),
                      desc="Searching automated matches",
                      ascii=False, ncols=72):
            submitted_name = unmatched_submissions_worth_checking[matching_name_col].values[i]
            submitted_id = unmatched_submissions_worth_checking[submission_id_col].values[i]
            word_combinations_to_check = get_word_combinations(submitted_name)
            if family_column is not None:
                submitted_family = unmatched_submissions_worth_checking[family_column].values[i]

                for index, row in wcvp_name_containment_df.iterrows():
                    record = row.to_frame().transpose()[
                        output_record_col_names + [wcvp_columns['family'], wcvp_columns['name']]]

                    if submitted_family == row[wcvp_columns['family']] or submitted_family == row[
                        wcvp_accepted_columns['family']]:
                        taxa = row[wcvp_columns['name']]

                        if taxa in word_combinations_to_check:
                            record[matching_name_col] = [submitted_name]
                            record[family_column] = [submitted_family]
                            record['matched_name'] = taxa
                            record[submission_id_col] = submitted_id
                            match_df = pd.concat([match_df, record])
            else:
                for index, row in wcvp_name_containment_df.iterrows():
                    record = row.to_frame().transpose()[
                        output_record_col_names + [wcvp_columns['family'], wcvp_columns['name']]]

                    taxa = row[wcvp_columns['name']]

                    if taxa in word_combinations_to_check:
                        record[matching_name_col] = [submitted_name]
                        record[submission_id_col] = submitted_id
                        record['matched_name'] = taxa
                        match_df = pd.concat([match_df, record])


        if len(match_df.index) > 0:
            match_df = match_df.dropna(subset=[wcvp_accepted_columns['name']])
            if family_column is None:
                # Remove genera matches where genus appears in different families
                unique_family_genera_pairs = all_taxa[
                    [wcvp_columns['family'], wcvp_columns['genus']]].drop_duplicates(keep='first')
                genera_unique_to_family = unique_family_genera_pairs[wcvp_columns['genus']].unique()

                match_df = match_df[
                    (~(match_df[wcvp_accepted_columns['rank']] == 'Genus')) | (
                        match_df[wcvp_columns['name']].isin(genera_unique_to_family))]

            match_df['matched_by'] = 'autoresolution'
            # Appropriately label unique matches
            match_df['matched_by'] = np.where(
                match_df.duplicated(keep=False, subset=[submission_id_col]),
                match_df['matched_by'],
                match_df['matched_by'] + '_unique')

            # Remove duplicate matches with worse status
            for r in match_df["taxon_status"].unique():
                if r not in status_priority:
                    raise ValueError(f'Status priority list does not contain {r} and needs updating.')
            match_df[wcvp_columns['status']] = pd.Categorical(match_df[wcvp_columns['status']],
                                                              status_priority)
            match_df.sort_values(wcvp_columns['status'], inplace=True)
            # Drop duplicated submissions of the same rank with worse taxonomic status
            match_df.drop_duplicates(subset=[submission_id_col, wcvp_accepted_columns['rank']], inplace=True,
                                     keep='first')
            match_df[wcvp_columns['status']] = match_df[wcvp_columns['status']].astype(object)
            # Remove duplicate matches with worse specificity
            for r in match_df[wcvp_accepted_columns['rank']].unique():
                if r not in rank_priority:
                    print(match_df[match_df[wcvp_accepted_columns['rank']] == r][wcvp_columns['name']].values)
                    raise ValueError(f'Rank priority list does not contain {r} and needs updating.')
            match_df[wcvp_accepted_columns['rank']] = pd.Categorical(match_df[wcvp_accepted_columns['rank']],
                                                                     rank_priority)
            match_df.sort_values(wcvp_accepted_columns['rank'], inplace=True)

            # Get the most precise match by dropping duplicate submissions
            match_df.drop_duplicates(subset=[submission_id_col], keep='first', inplace=True)
            match_df[wcvp_accepted_columns['rank']] = match_df[wcvp_accepted_columns['rank']].astype(object)
            # Merge with original data
            matches = pd.merge(unmatched_submissions_df[[submission_id_col]], match_df,
                               on=submission_id_col,
                               sort=False)
            matches = matches.dropna(subset=[wcvp_accepted_columns['name']])

            return matches
        else:
            return match_df
    else:
        return unmatched_submissions_df


def _get_knms_matches_and_accepted_info_from_names_in_column(df: pd.DataFrame, matching_name_col: str,
                                                             unique_submission_id_col: str,
                                                             all_taxa: pd.DataFrame,
                                                             family_column: str = None
                                                             ) -> pd.DataFrame:
    """
    Matches names in df using knms and gets corresponding accepted info from wcvp.
    :param df:
    :param matching_name_col:
    :return:
    """
    if len(df.index) > 0:
        match_records = get_knms_name_matches(df[matching_name_col].unique())
        match_records = pd.merge(match_records, df, left_on='submitted', right_on=matching_name_col)
        match_records['ipni_id'] = match_records['ipni_id'].apply(clean_urn_ids)
        match_records = get_accepted_wcvp_info_from_ipni_ids_in_column(match_records, 'ipni_id', all_taxa)
        match_records = match_records[match_records['match_state'] != 'false']
        if family_column is not None:

            # Begin by removing incorrectly matched families
            match_records = match_records.loc[
                (match_records[wcvp_columns['family']] == match_records[family_column]) | (
                        match_records[wcvp_accepted_columns['family']] == match_records[family_column]) |
                (match_records[family_column].isna())]
            if len(match_records.index) > 0:
                new_uniques = match_records[
                    ~match_records.duplicated(subset=[unique_submission_id_col], keep=False)]

                match_records['match_state'] = match_records.apply(
                    lambda row: 'true' if row[unique_submission_id_col] in new_uniques[
                        unique_submission_id_col].values else row[
                        'match_state'], axis=1)

        single_matches = match_records[match_records['match_state'] == 'true'].copy()
        single_matches['matched_by'] = 'knms_single'

        multiple_matches = match_records[match_records['match_state'] == 'multiple_matches'].copy()
        if len(multiple_matches.index) > 0:
            best_matches = _find_best_matches_from_multiple_knms_matches(multiple_matches,
                                                                         unique_submission_id_col)

            resolved_df = pd.concat([single_matches, best_matches], axis=0)
        else:
            resolved_df = single_matches

        resolved_df = resolved_df[
            [unique_submission_id_col] + output_record_col_names + ['matched_by', 'matched_name']]
        resolved_df = resolved_df.dropna(subset=[wcvp_accepted_columns['name']])
        resolved_df = resolved_df.drop_duplicates(subset=[unique_submission_id_col], keep='first')

        return resolved_df
    else:
        return df


def _find_best_matches_from_multiple_knms_matches(multiple_match_records: pd.DataFrame,
                                                  unique_submission_id_col: str) -> pd.DataFrame:
    """
    Find the best match from 'multiple' matches in KNMS. unique_submission_id_col specifies unique
    submissions, this can just be the submitted name or in cases where a list of families is included this must
    distinguish cases with the same name but different families
    :param multiple_match_records:
    :param unique_submission_id_col:
    :return:
    """

    # First, use matches where accepted name is the same as the submitted name
    accepted_names_matching_submitted_names = multiple_match_records[
        multiple_match_records['submitted'] == multiple_match_records[
            wcvp_accepted_columns['name']]]

    unmatched_containment_df = multiple_match_records[
        ~multiple_match_records[unique_submission_id_col].isin(
            accepted_names_matching_submitted_names[unique_submission_id_col].values)]

    # reduce list to remove essentially repeated matches
    unique_accepted_matches = unmatched_containment_df.drop_duplicates(
        subset=[unique_submission_id_col, wcvp_accepted_columns['ipni_id']], keep='first')

    #  Next use matches where the submitted name has a unique match
    submitted_names_with_single_accepted_match = unique_accepted_matches.drop_duplicates(
        subset=[unique_submission_id_col],
        keep=False)

    # Where neither of the above apply, sort by accepted rank
    unresolved_submissions = unique_accepted_matches[
        ~unique_accepted_matches[unique_submission_id_col].isin(
            submitted_names_with_single_accepted_match[unique_submission_id_col])]

    # matches where the accepted name is contained in the submitted name
    # In case of duplicates, these are sorted by specifity of rank
    unresolved_submissions = unresolved_submissions.dropna(subset=[wcvp_accepted_columns['name']])
    accepted_names_in_submitted_names = unresolved_submissions[
        unresolved_submissions.apply(lambda x: x[wcvp_accepted_columns['name']] in x['submitted'],
                                     axis=1)]

    for r in accepted_names_in_submitted_names[wcvp_accepted_columns['rank']].unique():
        if r not in rank_priority:
            raise ValueError(f'Rank priority list does not contain {r} and needs updating.')
    accepted_names_in_submitted_names[wcvp_accepted_columns['rank']] = pd.Categorical(
        accepted_names_in_submitted_names[wcvp_accepted_columns['rank']], rank_priority)
    accepted_names_in_submitted_names.sort_values(wcvp_accepted_columns['rank'], inplace=True)

    # Get the most precise match by dropping duplicate submissions
    accepted_names_in_submitted_names.drop_duplicates(subset=[unique_submission_id_col], keep='first',
                                                      inplace=True)
    accepted_names_in_submitted_names[wcvp_accepted_columns['rank']] = accepted_names_in_submitted_names[
        wcvp_accepted_columns['rank']].astype(object)

    accepted_names_matching_submitted_names['matched_by'] = 'knms_multiple_1'
    submitted_names_with_single_accepted_match['matched_by'] = 'knms_multiple_2'
    accepted_names_in_submitted_names['matched_by'] = 'knms_multiple_3'
    matches_to_use = pd.concat([accepted_names_matching_submitted_names,
                                submitted_names_with_single_accepted_match,
                                accepted_names_in_submitted_names])
    matches_to_use.drop_duplicates(subset=[unique_submission_id_col], keep='first', inplace=True)

    unmatched_df = multiple_match_records[
        ~multiple_match_records[unique_submission_id_col].isin(
            matches_to_use[unique_submission_id_col].values)]
    if len(unmatched_df.index) > 0:
        _temp_output(unmatched_df, 'unmatched_samples_with_multiple_knms_hits',
                     'Warning: Some samples have multiple matches in knms but havent been automatically '
                     'resolved. These may still be resolved in later steps.')

    return matches_to_use


def get_accepted_wcvp_info_from_ipni_ids_in_column(df: pd.DataFrame, id_col_name: str, all_taxa: pd.DataFrame) -> \
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
                        all_taxa[[wcvp_columns['ipni_id'], wcvp_columns['family']] + output_record_col_names],
                        how='left',
                        left_on=id_col_name,
                        right_on=wcvp_columns['ipni_id'])
    if wcvp_columns['ipni_id'] not in in_df.columns:
        match_df = match_df.drop(columns=[wcvp_columns['ipni_id']])

    if len(match_df.index) != len(in_df.index):
        raise ValueError('Generating accepted info is mismatched')

    match_df[id_col_name] = match_df[id_col_name].replace('no_id_placeholder', np.nan)

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
    # Check for bad inputs
    match_levels = ['full', 'direct', 'knms']
    if match_level not in match_levels:
        raise ValueError(f'match_level should be one of {match_levels}')

    reserved_column_names = [submitted_name_col_id, submitted_family_name_col_id, recapitalised_name_col,
                             lowercase_name_col,
                             unique_submission_index_col, 'submitted', 'matched_by', 'matched_name',
                             'resolution_id',
                             'taxon_name_with_taxon_authors', tidied_taxon_authors_col
                             ] + list(wcvp_columns.values()) + output_record_col_names
    problem_columns = [x for x in in_df.columns if
                       x in reserved_column_names]
    if len(problem_columns) > 0:
        raise ValueError(
            f'Column names used in input data will be confused in matching process: {problem_columns}. '
            f'Following column names are reserved: {reserved_column_names}')

    if len(in_df.index) > 0:
        df = in_df.copy(deep=True)
        # Standardise inputs
        df = df.dropna(subset=[name_col])
        tidy_names_in_column(df, name_col)
        if family_column is not None:
            tidy_families_in_column(df, family_column)

        # Check families of interest aligns with names in family column
        if family_column is not None:
            if families_of_interest is not None:
                fams_in_family_column = in_df[family_column].unique().dropna()
                for f in fams_in_family_column:
                    if f not in families_of_interest:
                        raise ValueError(
                            f'Family {f} given in family column but not in families of interest')

            if df[family_column].isnull().values.any():
                na_families = df[df[family_column].isna()]
                print(f'Warning: Missing values in specified family column')
                print(na_families)
            elif families_of_interest is None:
                families_of_interest = df[family_column].unique()

        # Check families of interest and in family column are in wcvp, and remove if not
        if families_of_interest is not None or family_column is not None:
            wcvp_families_df = get_all_taxa()
            wcvp_families = list(wcvp_families_df[wcvp_columns['family']].unique())
            wcvp_acc_families = list(wcvp_families_df[wcvp_accepted_columns['family']].unique())
            wcvp_all_families = wcvp_families + wcvp_acc_families
            if families_of_interest is not None:
                problem_fams = [f for f in families_of_interest if f not in wcvp_all_families]
                if len(problem_fams) > 0:
                    print(f'WARNING: Given families not in wcvp: {problem_fams}')
                    families_of_interest = None
            if family_column is not None:
                problem_fams = []
                for f in df[family_column].unique():
                    if pd.notna(f):
                        if f not in wcvp_all_families:
                            problem_fams.append(f)
                if len(problem_fams) > 0:
                    print(f'WARNING: Given families in family column  not in wcvp: {problem_fams},')
                    print('Setting these families to Nan')
                    for f in problem_fams:
                        df[family_column].replace(f, np.nan, inplace=True)

        # Create a unique identifier for submissions
        if family_column is None:
            df[unique_submission_index_col] = df[submitted_name_col_id].values.tolist()
            in_df[unique_submission_index_col] = in_df[name_col].values.tolist()

        else:
            df[unique_submission_index_col] = df[
                [submitted_family_name_col_id, submitted_name_col_id]].values.tolist()
            in_df[unique_submission_index_col] = in_df[
                [family_column, name_col]].values.tolist()

        df[unique_submission_index_col] = df[unique_submission_index_col].astype(str)
        in_df[unique_submission_index_col] = in_df[unique_submission_index_col].astype(str)
        df = df.drop_duplicates(subset=[unique_submission_index_col])

        all_taxa = get_all_taxa(families_of_interest=families_of_interest)
        # First get manual matches using given ipni ids
        if manual_resolution_csv is not None:
            manual_match_df = pd.read_csv(manual_resolution_csv)
            manual_match_df = manual_match_df[
                manual_match_df['submitted'].isin(df[submitted_name_col_id].values.tolist())]
            man_matches_with_accepted_info = get_accepted_wcvp_info_from_ipni_ids_in_column(manual_match_df,
                                                                                       'resolution_id',
                                                                                            all_taxa)
            man_matches_with_accepted_info = man_matches_with_accepted_info.dropna(
                subset=[wcvp_accepted_columns['name']])
            manual_matches = pd.merge(df, man_matches_with_accepted_info, left_on=submitted_name_col_id,
                                      right_on='submitted',
                                      sort=False)
            manual_matches['matched_by'] = 'manual'
            manual_matches['matched_name'] = np.nan
            unmatched_manual_df = df[
                ~df[unique_submission_index_col].isin(manual_matches[unique_submission_index_col].values)]
        else:
            manual_matches = pd.DataFrame()
            unmatched_manual_df = df

        # Then match with exact matches in wcvp
        wcvp_exact_name_match_df = get_wcvp_info_for_names_in_column(unmatched_manual_df,
                                                                     recapitalised_name_col,
                                                                     unique_submission_index_col,
                                                                     family_column=family_column,
                                                                     all_taxa=all_taxa)

        wcvp_resolved_df = pd.concat([wcvp_exact_name_match_df, manual_matches], axis=0)
        unmatched_name_df = df[
            ~df[unique_submission_index_col].isin(wcvp_resolved_df[unique_submission_index_col].values)]

        if match_level in ['full', 'knms']:
            # If exact matches aren't found in wcvp, use knms

            matches_with_knms = _get_knms_matches_and_accepted_info_from_names_in_column(unmatched_name_df,
                                                                                         recapitalised_name_col,
                                                                                         unique_submission_index_col,
                                                                                         all_taxa,
                                                                                         family_column=family_column)
            knms_resolved_df = pd.concat([wcvp_resolved_df, matches_with_knms], axis=0)
            unmatched_df = df[
                ~df[unique_submission_index_col].isin(knms_resolved_df[unique_submission_index_col].values)]

            if match_level == 'full':
                # Get autoresolved matches
                unmatched_resolutions = _autoresolve_missing_matches(unmatched_df, recapitalised_name_col,
                                                                     unique_submission_index_col,
                                                                     all_taxa,
                                                                     family_column=family_column)

                final_resolved_df = pd.concat([unmatched_resolutions, knms_resolved_df], axis=0)

            else:
                final_resolved_df = knms_resolved_df
        else:
            final_resolved_df = wcvp_resolved_df

        # Provide temp outputs
        unmatched_final_df = df[
            ~df[unique_submission_index_col].isin(final_resolved_df[unique_submission_index_col].values)]
        if len(unmatched_final_df.index) > 0:
            _temp_output(unmatched_final_df, 'unmatched_samples',
                         'WARNING: some submissions have not been resolved and must be manually resolved. '
                         'Consider fixing names in your original data.')
            final_resolved_df = pd.concat([final_resolved_df, unmatched_final_df])

        _temp_output(final_resolved_df, 'final_resolutions')

        final_resolved_df = final_resolved_df[
            [unique_submission_index_col] + output_record_col_names + ['matched_by', 'matched_name']]
        out_df = pd.merge(in_df, final_resolved_df, on=unique_submission_index_col,
                          how='left')
        in_df.drop(columns=[unique_submission_index_col], inplace=True)
        out_df = out_df.drop(columns=[unique_submission_index_col])
        # reorder
        out_df = out_df[in_df.columns.tolist() + output_record_col_names + ['matched_by', 'matched_name']]
        out_df['matched_name'] = out_df['matched_name'].apply(remove_whitespace_at_beginning_and_end)
        return out_df
    else:
        out_copy = in_df.copy()
        for a in output_record_col_names + ['matched_by', 'matched_name']:
            out_copy[a] = np.nan
        return out_copy
