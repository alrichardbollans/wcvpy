import os
import zipfile
import io

import numpy as np
import pandas as pd
import requests
from pkg_resources import resource_filename
from typing import List

_inputs_path = resource_filename(__name__, 'inputs')
_outputs_path = resource_filename(__name__, 'outputs')


# Standardise rank names
def capitalize_first_letter_of_rank(g: str):
    try:
        l = g.lower()

        return l.capitalize()
    except AttributeError:
        return g


def fix_columns(taxa_df: pd.DataFrame) -> pd.DataFrame:
    out_copy = taxa_df.copy(deep=True)
    # If taxanomic_status is accepted, make accepted_name and accepted_kew_id to taxon_name and kew_id

    out_copy['accepted_kew_id'] = np.where(out_copy['taxonomic_status'] == 'Accepted', taxa_df['kew_id'],
                                           taxa_df['accepted_kew_id'])
    out_copy['accepted_name'] = np.where(out_copy['taxonomic_status'] == 'Accepted', taxa_df['taxon_name'],
                                         taxa_df['accepted_name'])

    return out_copy


def get_all_taxa(families_of_interest: List[str] = None, ranks: List[str] = None,
                 accepted: bool = False, version: str = None, output_csv: str = None) -> pd.DataFrame:
    if output_csv is not None:
        if not os.path.isdir(os.path.dirname(output_csv)):
            os.mkdir(os.path.dirname(output_csv))

    if version is None:
        version = 'wcvp_v7_dec_2021'

    input_file = os.path.join(_inputs_path, version + '.txt')

    wcvp_link = 'http://sftp.kew.org/pub/data-repositories/WCVP/' + version + '.zip'

    # Download if doesn't exist
    if not os.path.exists(input_file):
        r = requests.get(wcvp_link, stream=True)
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall(_inputs_path)

    wcvp_data = pd.read_csv(input_file, sep='|')
    if families_of_interest is not None:
        wcvp_data = wcvp_data.loc[wcvp_data['family'].isin(families_of_interest)]

    if accepted:
        wcvp_data = wcvp_data[wcvp_data['taxonomic_status'] == 'Accepted']

    wcvp_data['rank'] = wcvp_data['rank'].apply(capitalize_first_letter_of_rank)

    if ranks is not None:
        wcvp_data = wcvp_data[wcvp_data['rank'].isin(ranks)]

    # Remove unplaced taxa
    wcvp_data = wcvp_data[wcvp_data['taxonomic_status'] != 'Unplaced']

    wcvp_data = fix_columns(wcvp_data)

    if output_csv is not None:
        wcvp_data.to_csv(output_csv)

    return wcvp_data


def main():
    get_all_taxa(families_of_interest=['Apocynaceae', 'Rubiaceae'], accepted=False,
                 output_csv=os.path.join(_outputs_path, 'wcvp_taxa_apocynaceae_rubiaceae.csv'))
    get_all_taxa(families_of_interest=['Apocynaceae', 'Rubiaceae'], accepted=True,
                 output_csv=os.path.join(_outputs_path, 'wcvp_accepted_taxa_apocynaceae_rubiaceae.csv'))
    get_all_taxa(families_of_interest=['Apocynaceae', 'Rubiaceae'], accepted=True, ranks=['Species'],
                 output_csv=os.path.join(_outputs_path, 'wcvp_accepted_species_apocynaceae_rubiaceae.csv'))


if __name__ == '__main__':
    main()
