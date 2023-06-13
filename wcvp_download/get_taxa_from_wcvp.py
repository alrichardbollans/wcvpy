import datetime
import hashlib
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
                'ipni_id': 'ipni_id',
                'status': 'taxon_status',
                'parent_name': 'parent_name',
                'parent_ipni_id': 'parent_ipni_id',
                'authors': 'taxon_authors',
                'paranthet_author': 'parenthetical_author',
                'primary_author': 'primary_author',
                'wcvp_id': 'plant_name_id',
                'parent_plant_name_id': 'parent_plant_name_id',
                'acc_plant_name_id': 'accepted_plant_name_id',
                'lifeform': 'lifeform_description'
                }

wcvp_accepted_columns = {'family': 'accepted_family',
                         'ipni_id': 'accepted_ipni_id',
                         'wcvp_id': 'accepted_plant_name_id',
                         'name': 'accepted_name',
                         'name_w_author': 'accepted_name_w_author',
                         'species': 'accepted_species',
                         'species_w_author': 'accepted_species_w_author',
                         'species_ipni_id': 'accepted_species_ipni_id',
                         'species_wcvp_id': 'accepted_species_id',
                         'rank': 'accepted_rank',
                         'parent_name': 'accepted_parent',
                         'parent_rank': 'accepted_parent_rank'
                         }

wcvp_columns_used_in_direct_matching = [wcvp_columns['genus'],
                                        wcvp_columns['family'],
                                        wcvp_columns['name'], wcvp_columns['authors'],
                                        wcvp_columns['paranthet_author'], wcvp_columns['primary_author']]

hybrid_characters = ["×", "+"]

infraspecific_chars = ['agamosp.', 'convar.', 'ecas.', 'f.', 'grex', 'group', 'lusus', 'microf.', 'microgene',
                       'micromorphe', 'modif.', 'monstr.', 'mut.', 'nid', 'nothof.', 'nothosubsp.',
                       'nothovar.', 'positio', 'proles', 'provar.', 'psp.', 'stirps', 'subf.', 'sublusus',
                       'subproles', 'subsp.', 'subspecioid', 'subvar.', 'unterrasse', 'var.']


def _clean_whitespaces(given_str: str):
    if pd.isnull(given_str):
        return given_str
    else:
        stripped = given_str.strip()
        out = " ".join(stripped.split())
        return out


def add_accepted_info_to_rows(taxa_df: pd.DataFrame, all_accepted: pd.DataFrame) -> pd.DataFrame:
    all_accepted = all_accepted.assign(accepted_ipni_id=all_accepted[wcvp_columns['ipni_id']])
    all_accepted = all_accepted.assign(accepted_name=all_accepted[wcvp_columns['name']])
    all_accepted = all_accepted.assign(accepted_family=all_accepted[wcvp_columns['family']])
    all_accepted = all_accepted.assign(accepted_rank=all_accepted[wcvp_columns['rank']])
    all_accepted = all_accepted.assign(accepted_parent=all_accepted[wcvp_columns['parent_name']])
    all_accepted = all_accepted.assign(accepted_parent_w_author=all_accepted['parent_name_w_author'])
    all_accepted = all_accepted.assign(accepted_parent_ipni_id=all_accepted[wcvp_columns['parent_ipni_id']])
    all_accepted = all_accepted.assign(accepted_parent_id=all_accepted[wcvp_columns['parent_plant_name_id']])
    all_accepted = all_accepted.assign(accepted_parent_rank=all_accepted['parent_rank'])

    all_accepted['accepted_name_w_author'] = all_accepted[wcvp_columns['name']].str.cat(
        all_accepted[wcvp_columns['authors']].fillna(''),
        sep=' ').str.strip()

    all_accepted = all_accepted.rename(columns={wcvp_columns['wcvp_id']: 'plant_name_id_acc'})
    all_accepted = all_accepted[
        ['plant_name_id_acc', 'accepted_ipni_id', 'accepted_name', 'accepted_name_w_author',
         'accepted_family', 'accepted_rank',
         'accepted_parent', 'accepted_parent_w_author', 'accepted_parent_id', 'accepted_parent_ipni_id',
         'accepted_parent_rank']]
    taxa_df_with_accepted_id = pd.merge(all_accepted, taxa_df, left_on='plant_name_id_acc',
                                        right_on=wcvp_columns['acc_plant_name_id'], how='right')
    taxa_df_with_accepted_id = taxa_df_with_accepted_id.drop(columns=['plant_name_id_acc'])

    return taxa_df_with_accepted_id


def get_parent_names_and_ipni_ids(taxa_df: pd.DataFrame, all_data: pd.DataFrame) -> pd.DataFrame:
    parent_data = all_data.drop(columns=['parent_plant_name_id'])
    parent_data['parent_name_w_author'] = parent_data[wcvp_columns['name']].str.cat(
        parent_data[wcvp_columns['authors']].fillna(''),
        sep=' ').str.strip()

    parent_data = parent_data.rename(columns={wcvp_columns['wcvp_id']: 'parent_plant_name_id',
                                              wcvp_columns['name']: wcvp_columns['parent_name'],
                                              wcvp_columns['ipni_id']: wcvp_columns['parent_ipni_id'],
                                              wcvp_columns['rank']: 'parent_rank'})

    parent_data = parent_data[
        ['parent_plant_name_id', wcvp_columns['parent_name'], 'parent_name_w_author',
         wcvp_columns['parent_ipni_id'], 'parent_rank']]
    taxa_df_with_parent_info = pd.merge(parent_data, taxa_df, left_on='parent_plant_name_id',
                                        right_on='parent_plant_name_id', how='right')

    return taxa_df_with_parent_info


def get_species_names_and_ipni_ids(taxa_df: pd.DataFrame):
    taxa_df['accepted_species'] = np.where(taxa_df['accepted_rank'] == 'Species', taxa_df['accepted_name'],
                                           np.nan)

    taxa_df['accepted_species_w_author'] = np.where(taxa_df['accepted_rank'] == 'Species',
                                                    taxa_df['accepted_name_w_author'],
                                                    np.nan)

    taxa_df['accepted_species'] = np.where(taxa_df['accepted_parent_rank'] == 'Species',
                                           taxa_df['accepted_parent'],
                                           taxa_df['accepted_species'])

    taxa_df['accepted_species_w_author'] = np.where(taxa_df['accepted_parent_rank'] == 'Species',
                                                    taxa_df['accepted_parent_w_author'],
                                                    taxa_df['accepted_species_w_author'])

    # IPNI ids
    taxa_df['accepted_species_ipni_id'] = np.where(taxa_df['accepted_rank'] == 'Species',
                                                   taxa_df[wcvp_accepted_columns['ipni_id']],
                                                   np.nan)
    taxa_df['accepted_species_ipni_id'] = np.where(taxa_df['accepted_parent_rank'] == 'Species',
                                                   taxa_df['accepted_parent_ipni_id'],
                                                   taxa_df['accepted_species_ipni_id'])

    # WCVP ids
    taxa_df['accepted_species_id'] = np.where(taxa_df['accepted_rank'] == 'Species',
                                              taxa_df[wcvp_accepted_columns['wcvp_id']],
                                              np.nan)
    taxa_df['accepted_species_id'] = np.where(taxa_df['accepted_parent_rank'] == 'Species',
                                              taxa_df['accepted_parent_id'],
                                              taxa_df['accepted_species_id'])

    return taxa_df


def get_wcvp_zip(get_new_version: bool = False):
    wcvp_path = 'http://sftp.kew.org/pub/data-repositories/WCVP'
    wcvp_link = '/'.join([wcvp_path, 'wcvp.zip'])

    input_zip_file = os.path.join(_inputs_path, 'wcvp.zip')

    if not os.path.exists(_inputs_path):
        os.mkdir(_inputs_path)

    def download_newest():
        print('Downloading latest WCVP version...')
        r = requests.get(wcvp_link, stream=True)
        with open(input_zip_file, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

    print(f'Loading WCVP...')
    if get_new_version:

        print(f'The latest version will be downloaded if not already available at {input_zip_file}')
        # Download if doesn't exist
        if not os.path.exists(input_zip_file):
            download_newest()

        else:
            # Download if online version is newer
            r = requests.head(wcvp_link)
            url_time = r.headers['last-modified']
            url_date = parsedate(url_time).astimezone()
            file_time = datetime.datetime.fromtimestamp(os.path.getmtime(input_zip_file)).astimezone()
            if url_date > file_time:
                download_newest()

    elif not os.path.exists(input_zip_file):
        download_newest()
    else:
        r = requests.head(wcvp_link)
        url_time = r.headers['last-modified']
        url_date = parsedate(url_time).astimezone()
        file_time = datetime.datetime.fromtimestamp(os.path.getmtime(input_zip_file)).astimezone()
        if url_date > file_time:
            print(
                f'WARNING: Loading your existing version of WCVP which was is out of date. Downloaded at: {file_time}')
            print(f'A new checklist version was released at: {url_date}')
            print('To up date the WCVP version, run get_all_taxa(get_new_version=True)')

        else:
            print('Using up to date WCVP.')

    file_time = datetime.datetime.fromtimestamp(os.path.getmtime(input_zip_file)).astimezone()
    return file_time, zipfile.ZipFile(input_zip_file)


def get_all_taxa(families_of_interest: List[str] = None, ranks: List[str] = None, genera: List[str] = None,
                 species: List[str] = None,
                 specific_taxa: List[str] = None,
                 accepted: bool = False, statuses_to_drop=None, output_csv: str = None,
                 get_new_version: bool = False, clean_strings: bool = True) -> pd.DataFrame:
    start = time.time()

    if output_csv is not None:
        new_output_dir = os.path.dirname(output_csv)
        if not os.path.isdir(new_output_dir) and new_output_dir != '':
            os.mkdir(new_output_dir)

    filetime, zf = get_wcvp_zip(get_new_version=get_new_version)
    csv_file = zf.open('wcvp_names.csv')

    reading_dtypes = {'homotypic_synonym': object, wcvp_columns['wcvp_id']: object,
                      wcvp_columns['acc_plant_name_id']: object,
                      'parent_plant_name_id': object,
                      'basionym_plant_name_id': object}
    all_wcvp_data = pd.read_csv(csv_file, encoding='utf-8', sep='|', quotechar='"', quoting=3,
                                dtype=reading_dtypes)

    csv_file.close()

    str_to_hash = str([filetime]).encode()
    filepath_for_zipped_parsed_version = os.path.join(_inputs_path, "parsed_wcvp" + str(
        hashlib.md5(str_to_hash).hexdigest()) + ".csv.zst")

    if os.path.isfile(filepath_for_zipped_parsed_version):
        print(f'Using preparsed file:{filepath_for_zipped_parsed_version}')
        # need to correctly specify dtypes of extra columns
        reading_dtypes[wcvp_accepted_columns['species_wcvp_id']] = object
        reading_dtypes['accepted_parent_id'] = object
        reading_dtypes['genus_hybrid'] = object

        parsed_wcvp_data = pd.read_csv(filepath_for_zipped_parsed_version, encoding='utf-8', sep='|',
                                       quotechar='"', quoting=3,
                                       dtype=reading_dtypes)
    else:
        print(f'Parsing the checklist. Parsed file will be saved to:{filepath_for_zipped_parsed_version}')

        if clean_strings:
            # Clean strings
            for col in wcvp_columns_used_in_direct_matching:
                all_wcvp_data[col] = all_wcvp_data[col].apply(_clean_whitespaces)

        all_accepted = all_wcvp_data[
            all_wcvp_data[wcvp_columns['status']].isin(['Accepted', 'Artificial Hybrid'])]

        parsed_wcvp_data = add_accepted_info_to_rows(all_wcvp_data,
                                                     get_parent_names_and_ipni_ids(all_accepted,
                                                                                   all_wcvp_data))
        parsed_wcvp_data = get_species_names_and_ipni_ids(parsed_wcvp_data)
        if clean_strings:
            # only write if strings have been cleaned
            parsed_wcvp_data.to_csv(filepath_for_zipped_parsed_version, index=False, compression='zstd',
                                    encoding='utf-8', sep='|', quotechar='"', quoting=3)

    if statuses_to_drop is None:
        statuses_to_drop = ['Local Biotype']

    parsed_wcvp_data = parsed_wcvp_data[~all_wcvp_data[wcvp_columns['status']].isin(statuses_to_drop)]

    if genera is not None:
        for g in genera:
            if g not in all_wcvp_data[wcvp_columns['genus']].values:
                raise ValueError(f'Given genus: {g} not in WCVP. CASE SENSITIVE')
        parsed_wcvp_data = parsed_wcvp_data.loc[parsed_wcvp_data['genus'].isin(genera)]

    if species is not None:
        for s in species:
            if s not in all_wcvp_data['species'].values:
                raise ValueError(f'Given species: {s} not in WCVP. CASE SENSITIVE')
        parsed_wcvp_data = parsed_wcvp_data.loc[parsed_wcvp_data['species'].isin(species)]

    if accepted:
        parsed_wcvp_data = parsed_wcvp_data[parsed_wcvp_data[wcvp_columns['status']] == 'Accepted']

    known_ranks_in_wcvp = ['Species', 'nothosubsp.', 'Subspecies', 'Form', 'Variety', 'microgene', 'Genus',
                           'proles', 'nothof.', 'Subvariety', 'nothovar.', 'Subform', 'lusus', 'monstr.',
                           '[**]', '[*]', 'sublusus', 'Convariety', 'psp.', 'subspecioid', 'group', 'grex',
                           'stirps', 'mut.', 'subproles', 'nid', 'provar.', 'positio', 'micromorphe',
                           'modif.',
                           'ecas.', 'microf.', 'agamosp.']
    if ranks is not None:
        for r in ranks:
            if r not in known_ranks_in_wcvp:
                raise ValueError(f'Given rank: {r} not in wcvp ranks: {known_ranks_in_wcvp}. CASE SENSITIVE')
        parsed_wcvp_data = parsed_wcvp_data[parsed_wcvp_data[wcvp_columns['rank']].isin(ranks)]

    if specific_taxa is not None:
        for f in specific_taxa:
            if f not in all_wcvp_data[wcvp_columns['name']].values:
                raise ValueError(f'Given specific taxa: {f} not in WCVP. CASE SENSITIVE')
        parsed_wcvp_data = parsed_wcvp_data[parsed_wcvp_data[wcvp_columns['name']].isin(specific_taxa)]

    if families_of_interest is not None:
        for f in families_of_interest:
            if f not in all_wcvp_data[wcvp_columns['family']].values:
                raise ValueError(f'Given family: {f} not in WCVP. CASE SENSITIVE')
        parsed_wcvp_data = parsed_wcvp_data.loc[
            (parsed_wcvp_data[wcvp_columns['family']].isin(families_of_interest)) | (
                parsed_wcvp_data[wcvp_accepted_columns['family']].isin(families_of_interest))]

    if output_csv is not None:
        parsed_wcvp_data.to_csv(output_csv)

    end = time.time()
    print(f'Time elapsed for (down)loading WCVP: {end - start}s')
    return parsed_wcvp_data
