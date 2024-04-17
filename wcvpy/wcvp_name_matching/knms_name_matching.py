import hashlib
import json
import os
import unicodedata as ud

import numpy as np
import pandas as pd
import requests
from typing import List

from pkg_resources import resource_filename

inputs_path = resource_filename(__name__, 'inputs')
temp_outputs_dir = 'name matching temp outputs'
knms_outputs_dir = os.path.join(temp_outputs_dir, 'knms matches')

latin_letters = {}


def is_latin(uchr):
    # https://stackoverflow.com/a/3308844/8633026
    try:
        return latin_letters[uchr]
    except KeyError:
        return latin_letters.setdefault(uchr, 'LATIN' in ud.name(uchr))


def only_roman_chars(unistr):
    # https://stackoverflow.com/a/3308844/8633026
    return all(is_latin(uchr)
               for uchr in unistr
               if uchr.isalpha())


def get_knms_name_matches(names: List[str]) -> pd.DataFrame:
    """
    Searches knms for matching names
    :param names:
    :return:
    """
    try:
        os.mkdir(temp_outputs_dir)
    except FileExistsError as error:
        pass

    try:
        os.mkdir(knms_outputs_dir)
    except FileExistsError as error:
        pass

    unique_name_list = list(np.unique(names))

    temp_file_tag = 'knms_matches_cache_'
    existing_df = None
    existing_info = []
    for temp_cl_file in os.listdir(knms_outputs_dir):
        if temp_cl_file.startswith(temp_file_tag):
            # Pandas will read TRUE/true as bools and therefore as True rather than true
            existing_info.append(pd.read_csv(os.path.join(knms_outputs_dir, temp_cl_file), dtype={'match_state': str}, index_col=0))
    if len(existing_info) > 0:
        existing_df = pd.concat(existing_info)

        already_known_names = existing_df['submitted'].tolist()
        # remove the item for all its occurrences
        for alread_known in already_known_names:
            c = unique_name_list.count(alread_known)
            for i in range(c):
                unique_name_list.remove(alread_known)
    if len(unique_name_list) > 0:
        knms_url = "http://namematch.science.kew.org/api/v2/powo/match"
        res = requests.post(knms_url, json=unique_name_list)
        headings = ['submitted', 'match_state', 'ipni_id', 'matched_name']

        if res.status_code == 500:
            print('Possibly from non-latin scripts in names')
            print(unique_name_list)
            raise ValueError('Internal Server error from KNMS.')
        elif res.status_code == 429:
            raise ConnectionRefusedError('KNMS Rate limiting')

        elif res.status_code == 504:
            raise TimeoutError('Network Timedout (possibly due to lots of names)')

        else:
            content = json.loads(res.content.decode('utf-8'))
            try:
                if all(len(content["records"][x]) == 2 for x in range(len(content["records"]))):
                    shortened_headings = ['submitted', 'match_state']
                    records = pd.DataFrame(content["records"], columns=shortened_headings)
                    records['ipni_id'] = np.nan
                    records['matched_name'] = np.nan
                else:
                    records = pd.DataFrame(content["records"], columns=headings)
            except ValueError:
                raise requests.ConnectionError('records not retrieved due to server error')
            records.replace('', np.nan, inplace=True)
            records['submitted'].ffill(inplace=True)
            records['match_state'].ffill(inplace=True)

            if (records['match_state'] == 'false').all():
                print('All KNMS records return false. Not saving these records as this sometimes indicates server issues.')
            else:
                str_to_hash = str(unique_name_list).encode()
                temp_output_knms_csv = os.path.join(knms_outputs_dir, temp_file_tag + str(hashlib.md5(str_to_hash).hexdigest()) + ".csv")
                records.to_csv(temp_output_knms_csv)
            if existing_df is not None:
                records = pd.concat([records, existing_df])
    else:
        print(f'Already searched for these name in KNMS. Returning records from cache directory: {knms_outputs_dir}')
        records = existing_df
    return records


def clear_stored_knms_matches():
    for f in os.listdir(knms_outputs_dir):
        os.remove(os.path.join(knms_outputs_dir, f))
