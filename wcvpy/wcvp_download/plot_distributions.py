import os

import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from wcvpy.wcvp_download import get_distributions_for_accepted_taxa, native_code_column, introduced_code_column

_inputs_path = resource_filename(__name__, 'inputs')


def _reformat_dist_col(given_val):
    """
    Reformats the input column value by checking its validity and converting it
    to a list. If the input value is invalid (e.g., NaN), the function assigns
    a NaN value instead.

    :param given_val: Input value that needs to be reformatted.
    :type given_val: Any
    :return: A list converted from the given input value if valid, or NaN if the
        input is invalid.
    :rtype: list or numpy.nan
    """
    # Only needs literal eval when using dist df read from csv
    if given_val == given_val:
        out = list(given_val)
    else:
        out = np.nan
    return out


def _OHE_native_dists(df: pd.DataFrame) -> pd.DataFrame:
    """
    One-hot encodes the data present in the native distributions column of the given
    DataFrame. This method processes the data by reformatting the column values, applying
    one-hot encoding to them, and then appending the resulting encoded columns to the
    original DataFrame.

    :param df: The input DataFrame containing the column to be reformatted
        and one-hot encoded.
    :type df: pd.DataFrame
    :return: A DataFrame that includes the original columns along with the
        one-hot encoded columns derived from the specified column.
    :rtype: pd.DataFrame
    """
    df = df.reset_index(drop=True).dropna(subset=[native_code_column])

    df[native_code_column] = df[native_code_column].apply(_reformat_dist_col)

    multilabels = df[native_code_column].str.join('|').str.get_dummies()

    clean_df = df.join(multilabels)

    return clean_df


def _OHE_introduced_dists(df: pd.DataFrame) -> pd.DataFrame:
    """
    One-hot encodes the data present in the introduced distributions column of the given
    DataFrame. This method processes the data by reformatting the column values, applying
    one-hot encoding to them, and then appending the resulting encoded columns to the
    original DataFrame.

    :param df: The input DataFrame containing the column to be reformatted
        and one-hot encoded.
    :type df: pd.DataFrame
    :return: A DataFrame that includes the original columns along with the
        one-hot encoded columns derived from the specified column.
    :rtype: pd.DataFrame
    """
    df = df.reset_index(drop=True).dropna(subset=[introduced_code_column])

    df[introduced_code_column] = df[introduced_code_column].apply(_reformat_dist_col)

    multilabels = df[introduced_code_column].str.join('|').str.get_dummies()

    clean_df = df.join(multilabels)

    return clean_df


def get_region_distribution_dataframe_for_accepted_taxa(df: pd.DataFrame, acc_name_col: str, output_path: str = None,
                                                        include_doubtful: bool = False,
                                                        include_extinct: bool = False, include_introduced: bool = False, wcvp_version: str = None):
    '''
    Gets the number of species in df found in each region.
    :param output_path:
    :param df:
    :param acc_name_col:
    :param outpath:
    :param include_doubtful:
    :param include_extinct:
    :param include_introduced:
    :param wcvp_version:
    :return:
    '''
    df_with_dists = get_distributions_for_accepted_taxa(df.drop_duplicates(subset=[acc_name_col], keep='first'), acc_name_col,
                                                        include_doubtful=include_doubtful, include_extinct=include_extinct,
                                                        wcvp_version=wcvp_version)
    ohe_native = _OHE_native_dists(df_with_dists)

    region_sums = {}
    for c in ohe_native.columns:
        if c not in df_with_dists.columns:
            region_sums[c] = ohe_native[c].sum()

    if include_introduced:
        ohe_introduced = _OHE_introduced_dists(df_with_dists)
        for c in ohe_introduced.columns:
            if c not in df_with_dists.columns:
                if c in region_sums:
                    region_sums[c] += ohe_introduced[c].sum()
                else:
                    region_sums[c] = ohe_introduced[c].sum()

    dict_for_pandas = {'Region': region_sums.keys(), 'Number of Taxa': region_sums.values()}
    out_df = pd.DataFrame(dict_for_pandas)

    if output_path is not None:
        out_df.to_csv(output_path)

    return out_df


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

    return get_region_distribution_dataframe_for_accepted_taxa(df, acc_name_col, output_path=output_path,
                                                               include_doubtful=include_doubtful,
                                                               include_extinct=include_extinct, include_introduced=False, wcvp_version=wcvp_version)


def plot_number_accepted_taxa_in_regions(df: pd.DataFrame, acc_name_col: str, output_dir: str, output_file_name: str,
                                         include_doubtful: bool = False,
                                         include_extinct: bool = False, include_introduced: bool = False, wcvp_version: str = None,
                                         colormap: str = 'viridis'):
    '''
        An example of how to do some basic plotting using the output of get_native_region_distribution_dataframe_for_accepted_taxa
        :param df:
        :param acc_name_col:
        :param output_dir:
        :param output_file_name:
        :param include_doubtful:
        :param include_extinct:
        :param include_introduced:
        :param wcvp_version:
        :param colormap: matplotlib colormap name
        :return:
        '''
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cartopy.io.shapereader as shpreader

    df_with_region_data = get_region_distribution_dataframe_for_accepted_taxa(df, acc_name_col,
                                                                              output_path=os.path.join(output_dir,
                                                                                                       output_file_name + '_regions.csv'),
                                                                              include_doubtful=include_doubtful,
                                                                              include_extinct=include_extinct,
                                                                              include_introduced=include_introduced,
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
    plot_number_accepted_taxa_in_regions(df, acc_name_col, output_dir, output_file_name, include_doubtful=include_doubtful,
                                         include_extinct=include_extinct, include_introduced=False, wcvp_version=wcvp_version, colormap=colormap)
