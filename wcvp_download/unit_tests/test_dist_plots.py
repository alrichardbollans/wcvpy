import os
import unittest

import numpy as np
import pandas as pd
import pandas.testing

from wcvp_download import plot_native_number_accepted_taxa_in_regions, get_native_region_distribution_dataframe_for_accepted_taxa, get_all_taxa

_output_path = 'test_outputs'

_v11_test_dict = {'AFG': 73,
                  'AGE': 249,
                  'AGW': 162,
                  'ALB': 63,
                  'ALD': 16,
                  'ALG': 61,
                  'AMU': 13,
                  'AND': 178,
                  'ANG': 529}
trait_df = pd.read_csv(os.path.join('test_inputs', 'all_traits.csv'))

new_taxa = get_all_taxa(accepted=True)
new_trait_df = trait_df[trait_df['accepted_name'].isin(new_taxa['accepted_name'].values)]


class MyTestCase(unittest.TestCase):

    def test_example(self):
        plot_native_number_accepted_taxa_in_regions(
            pd.DataFrame(['Campomanesia thea', 'Campomanesia thea', 'Campomanesia thea', 'Campomanesia thea'], columns=['name']), 'name',
            'test_outputs', 'Campomanesia.jpg')

        plot_native_number_accepted_taxa_in_regions(trait_df, 'accepted_name', 'test_outputs', 'testv11.jpg', wcvp_version='11')

        plot_native_number_accepted_taxa_in_regions(new_trait_df, 'accepted_name', 'test_outputs', 'test.jpg')

    def test_synonymy(self):
        try:
            cts_region_df = get_native_region_distribution_dataframe_for_accepted_taxa(
                pd.DataFrame(['Campomanesia thea', 'Citrosma thea', 'Siparuna thea', 'Citrosma lindenii'], columns=['name']), 'name',
                output_path=os.path.join('test_outputs', 'fails.csv'))
        except ValueError as e:
            print(e)
        else:
            raise ValueError

        # ct_region_df.to_csv(os.path.join('test_outputs', 'ct_region_df.csv'))
        #
        # cts_region_df.to_csv(os.path.join('test_outputs', 'cts_region_df.csv'))
        # pandas.testing.assert_frame_equal(ct_region_df, cts_region_df)

    def test_region_example(self):
        ct_region_df = get_native_region_distribution_dataframe_for_accepted_taxa(pd.DataFrame(['Campomanesia thea', 'Campomanesia thea'], columns=['name']), 'name',
                                                                                  output_path=os.path.join('test_outputs', 'Campomanesia_regions.csv'))
        pandas.testing.assert_frame_equal(ct_region_df, pd.DataFrame([['BZS', 1]], columns=['Region', 'Number of Native Species']))

        v11_region_df = get_native_region_distribution_dataframe_for_accepted_taxa(trait_df, 'accepted_name', wcvp_version='11',
                                                                                   output_path=os.path.join('test_outputs', 'v11_traits_regions.csv'))
        # Region examples
        for ex in _v11_test_dict:
            self.assertEqual(v11_region_df[v11_region_df['Region'] == ex]['Number of Native Species'].iloc[0], _v11_test_dict[ex])

        # Outputs differ between versions
        new_region_df = get_native_region_distribution_dataframe_for_accepted_taxa(new_trait_df, 'accepted_name',
                                                                                   output_path=os.path.join('test_outputs', 'new_traits_regions.csv'))
        pandas.testing.assert_frame_equal(new_region_df, new_region_df)
        try:
            pandas.testing.assert_frame_equal(new_region_df, v11_region_df)
        except AssertionError as e:
            print(e)
        else:
            raise ValueError



if __name__ == '__main__':
    unittest.main()
