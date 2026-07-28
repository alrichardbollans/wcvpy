import os.path
import pathlib

import pandas as pd

from wcvpy.wcvp_download import get_all_taxa


def get_family_counts(wcvp_version):
    if wcvp_version is None:
        tag = 'latest'
        new_taxa = get_all_taxa(accepted=True, version=wcvp_version, get_new_version=True)
    else:
        tag = f'v{wcvp_version}'
        new_taxa = get_all_taxa(accepted=True, version=wcvp_version)

    out_dir = tag
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    species_df = new_taxa[new_taxa['accepted_rank'] == 'Species']
    genus_df = new_taxa[new_taxa['accepted_rank'] == 'Genus']

    # count the number of species in each family and output a dataframe of counts per family
    species_counts = species_df['accepted_family'].value_counts()
    species_counts.name = 'Species Count'
    genus_counts = genus_df['accepted_family'].value_counts()
    genus_counts.name = 'Genus Count'

    family_counts = pd.concat([species_counts, genus_counts], axis=1)
    family_counts.to_csv(os.path.join(out_dir, 'family_counts.csv'))


def main():
    get_family_counts(None)
    get_family_counts('15')


if __name__ == '__main__':
    main()
