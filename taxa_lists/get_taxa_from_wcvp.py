import os
import zipfile
import io

import pandas as pd
import requests
from pkg_resources import resource_filename

_inputs_path = resource_filename(__name__, 'inputs')
_outputs_path = resource_filename(__name__, 'outputs')


# Standardise rank names
def capitalize_first_letter_of_rank(g: str):
    try:
        l = g.lower()

        return l.capitalize()
    except AttributeError:
        return g


def fix_columns(taxa_df: pd.DataFrame):
    # If taxanomic_status is accepted, make accepted_name and accepted_kew_id to taxon_name and kew_id
    taxa_df[taxa_df['taxonomic_status'] == 'Accepted']['accepted_name'] = \
        taxa_df[taxa_df['taxonomic_status'] == 'Accepted']['taxon_name']

    taxa_df[taxa_df['taxonomic_status'] == 'Accepted']['accepted_kew_id'] = \
        taxa_df[taxa_df['taxonomic_status'] == 'Accepted']['kew_id']


def get_all_taxa(families_of_interest=None,
                 accepted=False, version=None, output_csv=None) -> pd.DataFrame:
    # TODO: Create/fill in various columns on loading, will make parsing later easier

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
    # Remove unplaced taxa
    wcvp_data = wcvp_data[wcvp_data['taxonomic_status'] != 'Unplaced']

    fix_columns(wcvp_data)

    if output_csv is not None:
        wcvp_data.to_csv(output_csv)

    return wcvp_data


def main():
    # get_accepted_taxa(output_csv=os.path.join(outputs_path, 'wcvp_accepted_taxa.csv'))
    get_all_taxa(families_of_interest=['Apocynaceae', 'Rubiaceae'], accepted=True,
                 output_csv=os.path.join(_outputs_path, 'wcvp_accepted_taxa_apocynaceae_rubiaceae.csv'))


if __name__ == '__main__':
    main()
