import time

import pandas as pd

from wcvp_download import get_up_to_date_wcvp_zip, get_all_taxa, wcvp_columns

native_code_column = 'native_tdwg3_codes'
introduced_code_column = 'intro_tdwg3_codes'


def get_distributions_for_taxa(df: pd.DataFrame, wcvp_id_col: str, include_doubtful: bool = False,
                               include_extinct: bool = False):
    start = time.time()
    wcvp_with_dists = add_distribution_list_to_wcvp(include_doubtful, include_extinct)
    wcvp_with_dists = wcvp_with_dists.dropna(subset=wcvp_columns['plant_name_id'])
    wcvp_with_dists = wcvp_with_dists[[wcvp_columns['plant_name_id'], native_code_column, introduced_code_column]]
    # relevant_data = wcvp_with_dists[wcvp_with_dists[wcvp_columns['plant_name_id'].isin(df[wcvp_id_col].values)]]
    output = pd.merge(df, wcvp_with_dists, how='left', left_on=wcvp_id_col, right_on=wcvp_columns['plant_name_id'])
    if wcvp_columns['plant_name_id'] not in df.columns:
        output = output.drop(columns=[wcvp_columns['plant_name_id']])

    end = time.time()
    print(f'Time elapsed for getting taxa distributions: {end - start}s')
    return output


def _sorted_tuple(iterable):
    return tuple(sorted(iterable))


def add_distribution_list_to_wcvp(include_doubtful: bool = False,
                                  include_extinct: bool = False):
    """
    Gets a copy of WCVP with distribution data for all taxa
    :param include_doubtful:
    :param include_extinct:
    :return:
    """
    # Only use accepted taxa for distributions as everything else is unreliable
    accepted_wcvp_data = get_all_taxa(accepted=True)
    wcvp_zip = get_up_to_date_wcvp_zip()

    csv_file = wcvp_zip.open('wcvp_distribution.csv')
    all_dist_data = pd.read_csv(csv_file, encoding='utf-8', sep='|',
                                dtype={wcvp_columns['plant_name_id']: object,
                                       'plant_locality_id': object})
    all_dist_data = all_dist_data.dropna(subset=['area_code_l3'])
    csv_file.close()

    merged = pd.merge(accepted_wcvp_data, all_dist_data, on='plant_name_id', how='left')
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

    accepted_taxa_with_natives = pd.merge(accepted_wcvp_data, grouped_natives, how='left', on='plant_name_id')
    accepted_taxa_with_natives_intros = pd.merge(accepted_taxa_with_natives, grouped_intros, how='left',
                                                 on='plant_name_id')
    accepted_taxa_with_natives_intros = accepted_taxa_with_natives_intros[
        [introduced_code_column, native_code_column, wcvp_columns['acc_plant_name_id']]]
    accepted_taxa_with_natives_intros = accepted_taxa_with_natives_intros.dropna(
        subset=[wcvp_columns['acc_plant_name_id']])
    # Update taxa list with distributions from accepted taxa
    all_wcvp = get_all_taxa()
    wcvp_data_with_distributions = pd.merge(all_wcvp, accepted_taxa_with_natives_intros,
                                            on=wcvp_columns['acc_plant_name_id'], how='left')
    wcvp_data_with_distributions = wcvp_data_with_distributions.dropna(
        subset=[introduced_code_column, native_code_column], how='all')

    return wcvp_data_with_distributions
