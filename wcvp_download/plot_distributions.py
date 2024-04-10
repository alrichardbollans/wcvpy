import os

import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from wcvp_download import get_distributions_for_accepted_taxa, native_code_column, wcvp_accepted_columns

_inputs_path = resource_filename(__name__, 'inputs')


def _OHE_native_dists(df: pd.DataFrame) -> pd.DataFrame:
    df = df.reset_index(drop=True).dropna(subset=[native_code_column])

    def reformat_dist_col(given_val):
        # Only needs literal eval when using dist df read from csv
        if given_val == given_val:
            out = list(given_val)
        else:
            out = np.nan
        return out

    df[native_code_column] = df[native_code_column].apply(reformat_dist_col)

    multilabels = df[native_code_column].str.join('|').str.get_dummies()

    clean_df = df.join(multilabels)

    return clean_df


def get_native_region_distribution_dataframe_for_accepted_taxa(df: pd.DataFrame, acc_name_col: str, output_path: str = None,
                                                               include_doubtful: bool = False,
                                                               include_extinct: bool = False, wcvp_version: str = None):
    '''
    Gets the number of native species in df found in each region.
    :param output_path:
    :param df:
    :param acc_name_col:
    :param outpath:
    :param include_doubtful:
    :param include_extinct:
    :param wcvp_version:
    :return:
    '''
    df_with_dists = get_distributions_for_accepted_taxa(df.drop_duplicates(subset=[acc_name_col], keep='first'), acc_name_col,
                                                        include_doubtful=include_doubtful, include_extinct=include_extinct,
                                                        wcvp_version=wcvp_version)
    ohe = _OHE_native_dists(df_with_dists)

    region_sums = {}
    for c in ohe.columns:
        if c not in df_with_dists.columns:
            region_sums[c] = ohe[c].sum()

    dict_for_pandas = {'Region': region_sums.keys(), 'Number of Taxa': region_sums.values()}
    out_df = pd.DataFrame(dict_for_pandas)

    if output_path is not None:
        out_df.to_csv(output_path)

    return out_df


def plot_native_number_accepted_taxa_in_regions(df: pd.DataFrame, acc_name_col: str, output_dir: str, output_file_name: str,
                                                include_doubtful: bool = False,
                                                include_extinct: bool = False, wcvp_version: str = None, colormap: str = 'viridis'):
    '''
    An example of how to do some basic plotting using the output of get_native_region_distribution_dataframe_for_accepted_taxa
    :param df:
    :param acc_name_col:
    :param output_dir:
    :param output_file_name:
    :param title:
    :param include_doubtful:
    :param include_extinct:
    :param wcvp_version:
    :param colormap: matplotlib colormap name
    :return:
    '''
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cartopy.io.shapereader as shpreader
    import fiona  # Need fiona to read shapefiles

    df_with_region_data = get_native_region_distribution_dataframe_for_accepted_taxa(df, acc_name_col,
                                                                                     output_path=os.path.join(output_dir,
                                                                                                              output_file_name + '_regions.csv'),
                                                                                     include_doubtful=include_doubtful,
                                                                                     include_extinct=include_extinct,
                                                                                     wcvp_version=wcvp_version)
    #
    tdwg3_shp = shpreader.Reader(
        os.path.join(_inputs_path, 'wgsrpd-master', 'level3', 'level3.shp'))
    tdwg3_region_codes = df_with_region_data['Region'].values
    # min_val = df_with_region_data['Number of Taxa'].min()
    max_val = df_with_region_data['Number of Taxa'].max()
    min_val = 1
    if max_val == 1:
        min_val = 0
    norm = plt.Normalize(min_val, max_val)
    print('plotting countries')

    plt.figure(figsize=(15, 9.375))
    # if title is not None:
    #     plt.title(title, fontsize=40)
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, linewidth=2)

    cmap = mpl.colormaps[colormap]
    for country in tdwg3_shp.records():

        tdwg_code = country.attributes['LEVEL3_COD']
        if tdwg_code in tdwg3_region_codes:

            # print(country.attributes['name_long'], next(earth_colors))
            ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                              facecolor=cmap(
                                  norm(df_with_region_data.loc[df_with_region_data['Region'] == tdwg_code, 'Number of Taxa'].iloc[
                                           0])),
                              label=tdwg_code)

            x = country.geometry.centroid.x
            y = country.geometry.centroid.y

            # ax.text(x, y, tdwg_code, color='black', size=10, ha='center', va='center',
            #         transform=ccrs.PlateCarree())
        else:
            # print(f"code not in given malarial isocodes: {tdwg_code}")
            ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                              facecolor='white',
                              label=tdwg_code)

    all_map_isos = [country.attributes['LEVEL3_COD'] for country in tdwg3_shp.records()]
    missed_names = [x for x in tdwg3_region_codes if x not in all_map_isos]
    print(f'iso codes not plotted on map: {missed_names}')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    plt.tight_layout()
    fig = plt.gcf()
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.175, 0.02, 0.65])
    cbar1 = fig.colorbar(sm, cax=cbar_ax)
    cbar1.ax.tick_params(labelsize=30)

    plt.savefig(os.path.join(output_dir, output_file_name), dpi=400, bbox_inches='tight')
    plt.close()
    plt.cla()
    plt.clf()
