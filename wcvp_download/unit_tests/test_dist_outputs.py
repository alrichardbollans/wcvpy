import os.path
import unittest

import numpy as np
import pandas as pd
from wcvp_download import get_all_taxa, wcvp_columns, wcvp_accepted_columns, native_code_column, \
    introduced_code_column
from wcvp_download import add_distribution_list_to_wcvp

with_dist = add_distribution_list_to_wcvp()

_output_path = 'test_outputs'


class MyTestCase(unittest.TestCase):
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

        subset_with_dist = with_dist.dropna(subset=[wcvp_accepted_columns['id']])
        should_be_unique = subset_with_dist.groupby(wcvp_accepted_columns['id'])[native_code_column].apply(
            set).reset_index(
            name='dists')
        problem_df = should_be_unique[should_be_unique['dists'].map(len) > 1]
        self.assertEqual(len(problem_df.index), 0)

    def test_different_with_without_doubtful(self):

        with_doubtful = add_distribution_list_to_wcvp(include_doubtful=True)

        pd.testing.assert_index_equal(with_dist.index, with_doubtful.index)
        self.assertEqual(with_dist.columns.tolist(), with_doubtful.columns.tolist())
        df_diff = pd.concat([with_dist, with_doubtful]).drop_duplicates(keep=False)

        self.assertNotEqual(len(df_diff.index), 0)

    def test_instances(self):
        native_test_dict = {'582307-1': np.nan, '491231-1': ('MDG',), '482884-1': ('MDG',),
                            '802077-1': (
                                'BGM', 'COR', 'CZE', 'FRA', 'GER', 'GRB', 'IRE', 'MOR', 'POR',
                                'SAR', 'SPA', 'SWE'), '77229223-1': (
                'BGM', 'COR', 'CZE', 'FRA', 'GER', 'GRB', 'IRE', 'MOR', 'POR',
                'SAR', 'SPA', 'SWE'), '801970-1': (
                'BGM', 'COR', 'CZE', 'FRA', 'GER', 'GRB', 'IRE', 'MOR', 'POR',
                'SAR', 'SPA', 'SWE'), '22751-1': ('CPP',), '480932-1': ('BOR', 'PHI')}
        for t in native_test_dict:
            row = with_dist[with_dist[wcvp_columns['id']] == t]
            if native_test_dict[t] == native_test_dict[t]:
                self.assertEqual(row[native_code_column].iloc[0], native_test_dict[t])
            else:
                self.assertTrue(np.isnan(row[native_code_column].iloc[0]))

        intro_test_dict = {'472162-1': ('JAM', 'MAU', 'TAN'),
                           '884567-1': ('JAM', 'MAU', 'TAN'),
                           '472169-1': ('JAM', 'MAU', 'TAN'),
                           '77100755-1': ('IND', 'MLY', 'NWG', 'PAK', 'ZIM')}
        for t in intro_test_dict:
            row = with_dist[with_dist[wcvp_columns['id']] == t]
            if intro_test_dict[t] == intro_test_dict[t]:
                self.assertEqual(row[introduced_code_column].iloc[0], intro_test_dict[t])
            else:
                self.assertTrue(np.isnan(row[introduced_code_column].iloc[0]))


if __name__ == '__main__':
    unittest.main()
