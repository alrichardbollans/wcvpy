import unittest

import numpy as np
import pandas as pd

from wcvp_download import add_distribution_list_to_wcvp, get_distributions_for_accepted_taxa
from wcvp_download import wcvp_columns, wcvp_accepted_columns, native_code_column, \
    introduced_code_column

with_dist = add_distribution_list_to_wcvp()

_output_path = 'test_outputs'

native_test_dict = {'Campomanesia thea': ('BZS',), 'Delonix regia': ('MDG',),
                    'Digitalis purpurea': (
                        'BGM', 'COR', 'CZE', 'FRA', 'GER', 'GRB', 'IRE', 'MOR', 'POR',
                        'SAR', 'SPA', 'SWE'), 'Digitalis purpurea var. purpurea': (
        'BGM', 'COR', 'CZE', 'FRA', 'GER', 'GRB', 'IRE', 'MOR', 'POR',
        'SAR', 'SPA', 'SWE'), 'Lebeckia': ('CPP',), 'Airyantha borneensis': ('BOR', 'PHI')}
intro_test_dict = {'Campomanesia thea': np.nan, 'Adenanthera microsperma': ('JAM', 'MAU', 'TAN'),
                   'Libidibia ferrea': ('IND', 'MLY', 'NWG', 'PAK', 'ZIM')}


def get_wcvp_id_from_ipni_id(all_taxa, ipni_id):
    return all_taxa[all_taxa['ipni_id'] == ipni_id]['plant_name_id'].values[0]


def get_acc_name_from_ipni_id(all_taxa, ipni_id):
    return all_taxa[all_taxa['ipni_id'] == ipni_id][wcvp_accepted_columns['name']].values[0]


class MyTestCase(unittest.TestCase):

    def test_getting_data_from_names(self):

        native_test_df = pd.DataFrame({'names': list(native_test_dict.keys())})
        out = get_distributions_for_accepted_taxa(native_test_df, 'names')

        for name in native_test_dict:
            row = out[out['names'] == name]
            print(row)
            print(name)
            if native_test_dict[name] == native_test_dict[name]:
                self.assertEqual(row[native_code_column].iloc[0], native_test_dict[name])
            else:
                self.assertTrue(np.isnan(row[native_code_column].iloc[0]))

        intro_test_df = pd.DataFrame({'names': list(intro_test_dict.keys())})
        out = get_distributions_for_accepted_taxa(intro_test_df, 'names')

        for name in intro_test_dict:
            row = out[out['names'] == name]
            if intro_test_dict[name] == intro_test_dict[name]:
                self.assertEqual(row[introduced_code_column].iloc[0], intro_test_dict[name])
            else:
                self.assertTrue(np.isnan(row[introduced_code_column].iloc[0]))

    def test_check_nonaccepted_have_dist(self):
        # Early check to verify wcvp distributions
        # This will break when dist fixed
        # Dist data is in general not given for non accepted taxa, except
        # where those taxa have no given accepted taxa or if taxa is local biotype

        subset_with_dist = with_dist[with_dist[wcvp_columns['status']] != 'Accepted']
        subset_with_dist = subset_with_dist[~subset_with_dist[wcvp_accepted_columns['name']].isna()]
        subset_with_dist = subset_with_dist.dropna(subset=[native_code_column, introduced_code_column],
                                                   how='all')
        print(subset_with_dist[wcvp_columns['status']].unique())
        print(subset_with_dist[native_code_column].unique())
        print(subset_with_dist[introduced_code_column].unique())
        self.assertNotEqual(len(subset_with_dist.index), 0)

    def test_check_accepted_and_synonyms_have_same_dist(self):
        # check distributions are the same for ids and associated accepted ids

        subset_with_dist = with_dist.dropna(subset=[wcvp_accepted_columns['ipni_id']])
        should_be_unique = subset_with_dist.groupby(wcvp_accepted_columns['ipni_id'])[native_code_column].apply(
            set).reset_index(
            name='dists')
        problem_df = should_be_unique[should_be_unique['dists'].map(len) > 1]
        self.assertEqual(len(problem_df.index), 0)

    def test_different_with_without_doubtful(self):

        with_doubtful = add_distribution_list_to_wcvp(include_doubtful=True)

        self.assertEqual(with_dist.columns.tolist(), with_doubtful.columns.tolist())
        df_diff = pd.concat([with_dist, with_doubtful]).drop_duplicates(keep=False)

        self.assertNotEqual(len(df_diff.index), 0)

    def test_instances(self):

        for t in native_test_dict:

            if native_test_dict[t] == native_test_dict[t]:
                row = with_dist[with_dist[wcvp_accepted_columns['name']] == t]
                self.assertEqual(row[native_code_column].iloc[0], native_test_dict[t])
            else:
                self.assertNotIn(t, with_dist[wcvp_accepted_columns['name']].values)

        for t in intro_test_dict:

            if intro_test_dict[t] == intro_test_dict[t]:
                row = with_dist[with_dist[wcvp_accepted_columns['name']] == t]
                self.assertEqual(row[introduced_code_column].iloc[0], intro_test_dict[t])


if __name__ == '__main__':
    unittest.main()
