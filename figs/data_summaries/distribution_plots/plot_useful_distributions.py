import os.path
import pathlib

import matplotlib.pyplot as plt
import pandas as pd

from wcvpy.wcvp_download import get_all_taxa, wcvp_accepted_columns, plot_native_number_accepted_taxa_in_regions, plot_number_accepted_taxa_in_regions
import seaborn as sns


def do_plots(wcvp_version):
    if wcvp_version is None:
        tag = 'latest'
        new_taxa = get_all_taxa(accepted=True, version=wcvp_version, get_new_version=True)
    else:
        tag = f'v{wcvp_version}'
        new_taxa = get_all_taxa(accepted=True, version=wcvp_version)

    out_dir = tag
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    to_plot = new_taxa[~new_taxa[wcvp_accepted_columns['species']].isna()]

    if wcvp_version is None:
        to_plot = to_plot[~to_plot['accepted_species'].isin(['Lachnagrostis adamsonii', 'Hibiscus brackenridgei'])]

    if wcvp_version == '11':
        to_plot = to_plot[~to_plot['accepted_species'].isin(['Staurogyne polybotrya', 'Perityle bisetosa', 'Vernonia bainesii'])]

    plot_native_number_accepted_taxa_in_regions(to_plot, wcvp_accepted_columns['species'],
                                                '', os.path.join(out_dir, 'all_species_native_distribution.jpg'), wcvp_version=wcvp_version,
                                                include_extinct=True)

    plot_number_accepted_taxa_in_regions(to_plot, wcvp_accepted_columns['species'],
                                         '', os.path.join(out_dir, 'all_species_distribution.jpg'), wcvp_version=wcvp_version,
                                         include_introduced=True,
                                         include_extinct=True)

    native_df = pd.read_csv(os.path.join(out_dir, 'all_species_native_distribution.jpg_regions.csv'), index_col=0)
    native_df.rename(columns={'Number of Taxa': 'Number of Native Taxa'}, inplace=True)
    intro_df = pd.read_csv(os.path.join(out_dir, 'all_species_distribution.jpg_regions.csv'), index_col=0)
    all_df = pd.merge(native_df, intro_df, on=['Region'])
    all_df['Number Introduced Taxa'] = all_df['Number of Taxa'] - all_df['Number of Native Taxa']

    all_df['Proportion Introduced Taxa'] = all_df['Number Introduced Taxa'] / all_df['Number of Taxa'] * 100
    assert all(all_df['Number Introduced Taxa'] >= 0)
    all_df.to_csv(os.path.join(out_dir, 'all_species_distribution.csv'))

    sns.displot(all_df['Proportion Introduced Taxa'])
    plt.savefig(os.path.join(out_dir, 'proportions_of_introduced_taxa_per_region.jpg'))


def main():
    # do_plots(None)
    do_plots('15')


if __name__ == '__main__':
    main()
