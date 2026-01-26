import os.path

import matplotlib.pyplot as plt
import pandas as pd

from wcvpy.wcvp_download import get_all_taxa, wcvp_accepted_columns, plot_native_number_accepted_taxa_in_regions, plot_number_accepted_taxa_in_regions
import seaborn as sns

if __name__ == '__main__':
    wcvp_version = None
    if wcvp_version is None:
        tag = 'latest'
        new_taxa = get_all_taxa(accepted=True, version=wcvp_version, get_new_version=True)
    else:
        tag = f'v{wcvp_version}'
        new_taxa = get_all_taxa(accepted=True, version=wcvp_version)
    to_plot = new_taxa[~new_taxa[wcvp_accepted_columns['species']].isna()]
    plot_native_number_accepted_taxa_in_regions(to_plot, wcvp_accepted_columns['species'],
                                                '', f'all_species_native_distribution_{tag}.jpg', wcvp_version=wcvp_version,
                                                include_extinct=True)

    plot_number_accepted_taxa_in_regions(to_plot, wcvp_accepted_columns['species'],
                                         '', f'all_species_distribution_{tag}.jpg', wcvp_version=wcvp_version, include_introduced=True,
                                         include_extinct=True)

    native_df = pd.read_csv(f'all_species_native_distribution_{tag}.jpg_regions.csv', index_col=0)
    native_df.rename(columns={'Number of Taxa': 'Number of Native Taxa'}, inplace=True)
    intro_df = pd.read_csv(f'all_species_distribution_{tag}.jpg_regions.csv', index_col=0)
    all_df = pd.merge(native_df, intro_df, on=['Region'])
    all_df['Number Introduced Taxa'] = all_df['Number of Taxa'] - all_df['Number of Native Taxa']

    all_df['Proportion Introduced Taxa'] = all_df['Number Introduced Taxa'] / all_df['Number of Taxa'] * 100
    assert all(all_df['Number Introduced Taxa']>=0)
    all_df.to_csv(f'all_species_distribution_{tag}.csv')

    sns.displot(all_df['Proportion Introduced Taxa'])
    plt.savefig(f'proportions_of_introduced_taxa_per_region_{tag}.jpg')

