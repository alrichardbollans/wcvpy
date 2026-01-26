import os
import unittest

import pandas as pd
import pandas.testing

from wcvpy.wcvp_download import plot_native_number_accepted_taxa_in_regions, get_native_region_distribution_dataframe_for_accepted_taxa, get_all_taxa, \
    wcvp_accepted_columns, wcvp_columns, get_region_distribution_dataframe_for_accepted_taxa, plot_number_accepted_taxa_in_regions

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
        plot_native_number_accepted_taxa_in_regions(trait_df, 'accepted_name', 'test_outputs', 'testv11.jpg', wcvp_version='11')

        # https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:44106-2
        plot_native_number_accepted_taxa_in_regions(
            pd.DataFrame(['Campomanesia thea', 'Campomanesia thea', 'Campomanesia thea', 'Campomanesia thea'], columns=['name']), 'name',
            'test_outputs', 'Campomanesia.jpg')

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
        ct_region_df = get_native_region_distribution_dataframe_for_accepted_taxa(
            pd.DataFrame(['Campomanesia thea', 'Campomanesia thea'], columns=['name']), 'name',
            output_path=os.path.join('test_outputs', 'Campomanesia_regions.csv'))
        pandas.testing.assert_frame_equal(ct_region_df, pd.DataFrame([['BZS', 1]], columns=['Region', 'Number of Taxa']))

        v11_region_df = get_native_region_distribution_dataframe_for_accepted_taxa(trait_df, 'accepted_name', wcvp_version='11',
                                                                                   output_path=os.path.join('test_outputs', 'v11_traits_regions.csv'))
        # Region examples
        for ex in _v11_test_dict:
            self.assertEqual(v11_region_df[v11_region_df['Region'] == ex]['Number of Taxa'].iloc[0], _v11_test_dict[ex])

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

    def test_plot_all(self):
        to_plot = new_taxa[new_taxa[wcvp_columns['rank']] == 'Species']
        to_plot = to_plot[~to_plot[wcvp_accepted_columns['name']].isna()]
        plot_native_number_accepted_taxa_in_regions(to_plot, wcvp_accepted_columns['name'],
                                                    'test_outputs', 'all_species_native_distribution.jpg',
                                                    include_extinct=True)

    def test_native_intro_example(self):

        test_df = pd.DataFrame(['Cinchona officinalis'], columns=['name'])
        # https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:58723-2
        plot_native_number_accepted_taxa_in_regions(test_df, 'name',
                                                    'test_outputs', 'Cinchona officinalis_native.jpg')
        test_input = pd.read_csv(os.path.join('test_inputs', 'Cinchona officinalis_native.jpg_regions.csv'),index_col=0)
        test_output = pd.read_csv(os.path.join('test_outputs', 'Cinchona officinalis_native.jpg_regions.csv'),index_col=0)
        pd.testing.assert_frame_equal(test_input, test_output)

        plot_number_accepted_taxa_in_regions(test_df, 'name',
                                             'test_outputs', 'Cinchona officinalis_all.jpg', include_introduced=True)
        test_input = pd.read_csv(os.path.join('test_inputs', 'Cinchona officinalis_all.jpg_regions.csv'), index_col=0)
        test_output = pd.read_csv(os.path.join('test_outputs', 'Cinchona officinalis_all.jpg_regions.csv'), index_col=0)
        pd.testing.assert_frame_equal(test_input, test_output)

    def test_genus_example(self):

        test_df = pd.DataFrame(['Campomanesia thea', 'Campomanesia thea', 'Campomanesia thea', 'Campomanesia thea'], columns=['name'])
        test_df['Genus'] = 'Campomanesia'
        # https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:30000055-2
        plot_native_number_accepted_taxa_in_regions(test_df, 'Genus',
                                                    'test_outputs', 'Campomanesia genus.jpg')

    def test_small_example(self):

        # https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:44106-2
        plot_native_number_accepted_taxa_in_regions(
            pd.DataFrame(['Campomanesia thea', 'Campomanesia adamantium'], columns=['name']), 'name',
            'test_outputs', 'Campomanesia_pairs.jpg')

    def test_native_vs_introduced(self):
        native1_df = get_native_region_distribution_dataframe_for_accepted_taxa(new_trait_df, 'accepted_name')

        native2_df = get_region_distribution_dataframe_for_accepted_taxa(new_trait_df, 'accepted_name', include_introduced=False)
        pd.testing.assert_frame_equal(native1_df, native2_df)
        native1_df.rename(columns={'Number of Taxa': 'Number of Native Taxa'}, inplace=True)

        introduced_df = get_region_distribution_dataframe_for_accepted_taxa(new_trait_df, 'accepted_name', include_introduced=True)

        all_df = pd.merge(native1_df, introduced_df, on=['Region'])

        self.assertTrue(all(all_df['Number of Taxa'] >= all_df['Number of Native Taxa']))
        self.assertFalse(all(all_df['Number of Taxa'] == all_df['Number of Native Taxa']))


if __name__ == '__main__':
    unittest.main()
