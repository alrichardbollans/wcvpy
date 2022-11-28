import datetime
import os
import time
import zipfile
from typing import List

import numpy as np
import pandas as pd
import requests
from dateutil.parser import parse as parsedate
from pkg_resources import resource_filename

_inputs_path = resource_filename(__name__, 'inputs')

wcvp_columns = {'family': 'family',
                'rank': 'taxon_rank',
                'genus': 'genus',
                'name': 'taxon_name',
                'id': 'ipni_id',
                'status': 'taxon_status',
                'parent_name': 'parent_name',
                'parent_ipni_id': 'parent_ipni_id',
                'authors': 'taxon_authors',
                'paranthet_author': 'parenthetical_author',
                'primary_author': 'primary_author',
                'plant_name_id':'plant_name_id'
                }

wcvp_accepted_columns = {'family': 'accepted_family',
                         'id': 'accepted_ipni_id',
                         'name': 'accepted_name'

                         }


def add_accepted_info_to_rows(taxa_df: pd.DataFrame, all_accepted: pd.DataFrame) -> pd.DataFrame:
    all_accepted = all_accepted.assign(accepted_ipni_id=all_accepted[wcvp_columns['id']])
    all_accepted = all_accepted.assign(accepted_name=all_accepted[wcvp_columns['name']])
    all_accepted = all_accepted.assign(accepted_family=all_accepted[wcvp_columns['family']])
    all_accepted = all_accepted.assign(accepted_rank=all_accepted[wcvp_columns['rank']])
    all_accepted = all_accepted.assign(accepted_parent=all_accepted[wcvp_columns['parent_name']])
    all_accepted = all_accepted.assign(accepted_parent_ipni_id=all_accepted[wcvp_columns['parent_ipni_id']])
    all_accepted = all_accepted.assign(accepted_parent_rank=all_accepted['parent_rank'])

    all_accepted = all_accepted.rename(columns={wcvp_columns['plant_name_id']: 'plant_name_id_acc'})
    all_accepted = all_accepted[
        ['plant_name_id_acc', 'accepted_ipni_id', 'accepted_name', 'accepted_family', 'accepted_rank',
         'accepted_parent', 'accepted_parent_ipni_id', 'accepted_parent_rank']]
    taxa_df_with_accepted_id = pd.merge(all_accepted, taxa_df, left_on='plant_name_id_acc',
                                        right_on='accepted_plant_name_id', how='right')
    taxa_df_with_accepted_id = taxa_df_with_accepted_id.drop(columns=['plant_name_id_acc'])

    return taxa_df_with_accepted_id


def get_parent_names_and_ipni_ids(taxa_df: pd.DataFrame, all_data: pd.DataFrame) -> pd.DataFrame:
    parent_data = all_data.drop(columns=['parent_plant_name_id'])
    parent_data = parent_data.rename(columns={wcvp_columns['plant_name_id']: 'parent_plant_name_id',
                                              wcvp_columns['name']: wcvp_columns['parent_name'],
                                              wcvp_columns['id']: wcvp_columns['parent_ipni_id'],
                                              wcvp_columns['rank']: 'parent_rank'})
    parent_data = parent_data[
        ['parent_plant_name_id', wcvp_columns['parent_name'], wcvp_columns['parent_ipni_id'], 'parent_rank']]
    taxa_df_with_parent_info = pd.merge(parent_data, taxa_df, left_on='parent_plant_name_id',
                                        right_on='parent_plant_name_id', how='right')

    return taxa_df_with_parent_info


def get_species_names_and_ipni_ids(taxa_df: pd.DataFrame):
    taxa_df['accepted_species'] = np.where(taxa_df['accepted_rank'] == 'Species', taxa_df['accepted_name'],
                                           np.nan)
    taxa_df['accepted_species'] = np.where(taxa_df['accepted_parent_rank'] == 'Species',
                                           taxa_df['accepted_parent'],
                                           taxa_df['accepted_species'])

    taxa_df['accepted_species_ipni_id'] = np.where(taxa_df['accepted_rank'] == 'Species',
                                                   taxa_df[wcvp_accepted_columns['id']],
                                                   np.nan)
    taxa_df['accepted_species_ipni_id'] = np.where(taxa_df['accepted_parent_rank'] == 'Species',
                                                   taxa_df['accepted_parent_ipni_id'],
                                                   taxa_df['accepted_species_ipni_id'])

    return taxa_df


def get_up_to_date_wcvp_zip(force_use_existing: bool = False):
    wcvp_path = 'http://sftp.kew.org/pub/data-repositories/WCVP'
    wcvp_link = '/'.join([wcvp_path, 'wcvp.zip'])

    input_zip_file = os.path.join(_inputs_path, 'wcvp.zip')
    if not force_use_existing:
        # Download if doesn't exist
        if not os.path.exists(input_zip_file):
            print('Downloading latest WCVP version...')
            r = requests.get(wcvp_link, stream=True)
            with open(input_zip_file, 'wb') as fd:
                for chunk in r.iter_content(chunk_size=128):
                    fd.write(chunk)

        else:
            # Download if online version is newer
            r = requests.head(wcvp_link)
            url_time = r.headers['last-modified']
            url_date = parsedate(url_time).astimezone()
            file_time = datetime.datetime.fromtimestamp(os.path.getmtime(input_zip_file)).astimezone()
            if url_date > file_time:
                print('Downloading latest WCVP version...')
                r = requests.get(wcvp_link, stream=True)
                with open(input_zip_file, 'wb') as fd:
                    for chunk in r.iter_content(chunk_size=128):
                        fd.write(chunk)

    return zipfile.ZipFile(input_zip_file)


def get_all_taxa(families_of_interest: List[str] = None, ranks: List[str] = None, genera: List[str] = None,
                 species: List[str] = None,
                 specific_taxa: List[str] = None,
                 accepted: bool = False, statuses_to_drop=None, output_csv: str = None,
                 force_use_existing: bool = False) -> pd.DataFrame:
    start = time.time()

    if output_csv is not None:
        new_output_dir = os.path.dirname(output_csv)
        if not os.path.isdir(new_output_dir) and new_output_dir != '':
            os.mkdir(new_output_dir)

    zf = get_up_to_date_wcvp_zip(force_use_existing=force_use_existing)
    csv_file = zf.open('wcvp_names.csv')
    all_wcvp_data = pd.read_csv(csv_file, encoding='utf-8', sep='|', dtype={'homotypic_synonym': object})

    csv_file.close()
    all_wcvp_data['accepted_plant_name_id'] = all_wcvp_data['accepted_plant_name_id'].astype(float)
    all_wcvp_data[wcvp_columns['plant_name_id']] = all_wcvp_data[wcvp_columns['plant_name_id']].astype(float)
    all_accepted = all_wcvp_data[all_wcvp_data[wcvp_columns['status']] == 'Accepted']

    wcvp_data = all_wcvp_data.copy(deep=True)

    if statuses_to_drop is not None:
        # Remove taxa with poor status
        wcvp_data = all_wcvp_data[~all_wcvp_data[wcvp_columns['status']].isin(statuses_to_drop)]

    if genera is not None:
        wcvp_data = wcvp_data.loc[wcvp_data['genus'].isin(genera)]

    if species is not None:
        wcvp_data = wcvp_data.loc[wcvp_data['species'].isin(species)]

    if accepted:
        wcvp_data = wcvp_data[wcvp_data[wcvp_columns['status']] == 'Accepted']

    if ranks is not None:
        wcvp_data = wcvp_data[wcvp_data[wcvp_columns['rank']].isin(ranks)]

    if specific_taxa is not None:
        wcvp_data = wcvp_data[wcvp_data[wcvp_columns['name']].isin(specific_taxa)]
    wcvp_data = add_accepted_info_to_rows(wcvp_data,
                                          get_parent_names_and_ipni_ids(all_accepted, all_wcvp_data))
    wcvp_data = get_species_names_and_ipni_ids(wcvp_data)

    if families_of_interest is not None:
        wcvp_data = wcvp_data.loc[(wcvp_data[wcvp_columns['family']].isin(families_of_interest)) | (
            wcvp_data[wcvp_accepted_columns['family']].isin(families_of_interest))]

    if output_csv is not None:
        wcvp_data.to_csv(output_csv)

    end = time.time()
    print(f'Time elapsed for (down)loading WCVP: {end - start}s')
    return wcvp_data
