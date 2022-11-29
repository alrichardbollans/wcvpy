from typing import List

import numpy as np
import pandas as pd

from wcvp_download import get_up_to_date_wcvp_zip, get_all_taxa, wcvp_columns, wcvp_accepted_columns

native_code_column = 'native_tdwg3_codes'
introduced_code_column = 'intro_tdwg3_codes'


def get_distributions_for_ipni_ids(ipni_ids: List[str], include_doubtful: bool = False,
                                   include_extinct: bool = False):
    wcvp_data = get_all_taxa()
    relevant_data = wcvp_data[wcvp_data[wcvp_columns['id'].isin(ipni_ids)]]

    return add_distribution_list_to_wcvp(relevant_data,
                                         include_doubtful=include_doubtful,
                                         include_extinct=include_extinct)


def _sorted_tuple(iterable):
    return tuple(sorted(iterable))


def add_accepted_dist_info_to_rows(taxa_df: pd.DataFrame, all_accepted: pd.DataFrame) -> pd.DataFrame:
    all_accepted = all_accepted.assign(accepted_native_code_column=all_accepted[wcvp_columns['id']])

    all_accepted = all_accepted.rename(columns={wcvp_columns['plant_name_id']: 'plant_name_id_acc'})
    all_accepted = all_accepted[
        ['plant_name_id_acc', 'accepted_ipni_id', 'accepted_name', 'accepted_family', 'accepted_rank',
         'accepted_parent', 'accepted_parent_ipni_id', 'accepted_parent_rank']]
    taxa_df_with_accepted_id = pd.merge(all_accepted, taxa_df, left_on='plant_name_id_acc',
                                        right_on='accepted_plant_name_id', how='right')
    taxa_df_with_accepted_id = taxa_df_with_accepted_id.drop(columns=['plant_name_id_acc'])

    return taxa_df_with_accepted_id


def add_distribution_list_to_wcvp(include_doubtful: bool = False,
                                  include_extinct: bool = False):
    """
    Gets a copy of WCVP with distribution data for all taxa
    :param include_doubtful:
    :param include_extinct:
    :return:
    """
    all_wcvp_data = get_all_taxa()
    wcvp_zip = get_up_to_date_wcvp_zip()

    csv_file = wcvp_zip.open('wcvp_distribution.csv')
    all_dist_data = pd.read_csv(csv_file, encoding='utf-8', sep='|')
    all_dist_data = all_dist_data.dropna(subset=['area_code_l3'])
    csv_file.close()

    all_dist_data.head(100).to_csv('dist_example.csv')

    merged = pd.merge(all_wcvp_data, all_dist_data, on='plant_name_id', how='left')
    if include_doubtful and include_extinct:
        natives = merged[merged['introduced'] == 0]
        intros = merged[merged['introduced'] == 1]

    elif include_extinct and not include_doubtful:
        natives = merged[(merged['introduced'] == 0) & (merged['location_doubtful'] == 0)]
        intros = merged[(merged['introduced'] == 1) & (merged['location_doubtful'] == 0)]

    elif include_doubtful and not include_extinct:
        natives = merged[(merged['introduced'] == 0) & (merged['extinct'] == 0)]
        intros = merged[(merged['introduced'] == 1) & (merged['extinct'] == 0)]

    else:
        natives = merged[
            (merged['introduced'] == 0) & (merged['extinct'] == 0) & (merged['location_doubtful'] == 0)]
        intros = merged[
            (merged['introduced'] == 1) & (merged['extinct'] == 0) & (merged['location_doubtful'] == 0)]

    grouped_natives = natives.groupby('plant_name_id')['area_code_l3'].apply(_sorted_tuple).reset_index(
        name=native_code_column)

    grouped_intros = intros.groupby('plant_name_id')['area_code_l3'].apply(_sorted_tuple).reset_index(
        name=introduced_code_column)

    taxa_with_natives = pd.merge(all_wcvp_data, grouped_natives, how='left', on='plant_name_id')
    taxa_with_natives_intros = pd.merge(taxa_with_natives, grouped_intros, how='left', on='plant_name_id')

    # Update taxa distributions with distributions from accepted taxa
    ### Native
    accepted_native_dists = taxa_with_natives_intros.dropna(
        subset=[native_code_column, wcvp_columns['acc_plant_name_id']], how='any')
    accepted_native_dists = accepted_native_dists[[native_code_column, wcvp_columns['acc_plant_name_id']]]

    wcvp_data_with_distributions = pd.merge(taxa_with_natives_intros, accepted_native_dists, how='left',
                                            on=wcvp_columns['acc_plant_name_id'])

    if native_code_column + '_x' in wcvp_data_with_distributions.columns:
        wcvp_data_with_distributions[native_code_column] = np.where(
            wcvp_data_with_distributions[native_code_column + '_x'].isna(),
            wcvp_data_with_distributions[native_code_column + '_y'],
            wcvp_data_with_distributions[native_code_column + '_x'])

        check_df = wcvp_data_with_distributions[
            wcvp_data_with_distributions[wcvp_columns['status']] == 'Accepted']
        problem_df = check_df.loc[
            ~(check_df[native_code_column + '_x'] == check_df[native_code_column + '_y'])].dropna(
            subset=[native_code_column + '_x', native_code_column + '_y'], how='all')
        assert len(problem_df.index) == 0

        check_df = wcvp_data_with_distributions[
            (wcvp_data_with_distributions[wcvp_columns['status']] != 'Accepted') & (
                wcvp_data_with_distributions[wcvp_columns['acc_plant_name_id']].isna())].dropna(
            subset=[native_code_column + '_x', native_code_column + '_y'], how='all')
        pd.testing.assert_series_equal(check_df[native_code_column + '_x'], check_df[native_code_column],
                                       check_names=False)
        assert check_df[native_code_column + '_y'].isnull().all()

        check_df = wcvp_data_with_distributions[
            (wcvp_data_with_distributions[wcvp_columns['status']] != 'Accepted') & (
                ~wcvp_data_with_distributions[wcvp_columns['acc_plant_name_id']].isna())].dropna(
            subset=[native_code_column + '_x', native_code_column + '_y'], how='all')
        pd.testing.assert_series_equal(check_df[native_code_column + '_y'], check_df[native_code_column],
                                       check_names=False)
        assert check_df[native_code_column + '_x'].isnull().all()

    ### Introduced
    accepted_intro_dists = taxa_with_natives_intros.dropna(
        subset=[introduced_code_column, wcvp_columns['acc_plant_name_id']], how='any')
    accepted_intro_dists = accepted_intro_dists[[introduced_code_column, wcvp_columns['acc_plant_name_id']]]

    wcvp_data_with_distributions = pd.merge(wcvp_data_with_distributions, accepted_intro_dists, how='left',
                                            on=wcvp_columns['acc_plant_name_id'])

    if introduced_code_column + '_x' in wcvp_data_with_distributions.columns:
        wcvp_data_with_distributions[introduced_code_column] = np.where(
            wcvp_data_with_distributions[introduced_code_column + '_x'].isna(),
            wcvp_data_with_distributions[introduced_code_column + '_y'],
            wcvp_data_with_distributions[introduced_code_column + '_x'])

        check_df = wcvp_data_with_distributions[
            wcvp_data_with_distributions[wcvp_columns['status']] == 'Accepted']
        problem_df = check_df.loc[
            ~(check_df[introduced_code_column + '_x'] == check_df[introduced_code_column + '_y'])].dropna(
            subset=[introduced_code_column + '_x', introduced_code_column + '_y'], how='all')
        assert len(problem_df.index) == 0

        check_df = wcvp_data_with_distributions[
            (wcvp_data_with_distributions[wcvp_columns['status']] != 'Accepted') & (
                wcvp_data_with_distributions[wcvp_columns['acc_plant_name_id']].isna())].dropna(
            subset=[introduced_code_column + '_x', introduced_code_column + '_y'], how='all')
        pd.testing.assert_series_equal(check_df[introduced_code_column + '_x'],
                                       check_df[introduced_code_column],
                                       check_names=False)
        assert check_df[introduced_code_column + '_y'].isnull().all()

        check_df = wcvp_data_with_distributions[
            (wcvp_data_with_distributions[wcvp_columns['status']] != 'Accepted') & (
                ~wcvp_data_with_distributions[wcvp_columns['acc_plant_name_id']].isna())].dropna(
            subset=[introduced_code_column + '_x', introduced_code_column + '_y'], how='all')
        pd.testing.assert_series_equal(check_df[introduced_code_column + '_y'],
                                       check_df[introduced_code_column],
                                       check_names=False)
        assert check_df[introduced_code_column + '_x'].isnull().all()

    return wcvp_data_with_distributions
