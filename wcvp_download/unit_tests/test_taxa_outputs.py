import os
import unittest

import pandas as pd

from wcvp_download import get_all_taxa, wcvp_columns, wcvp_accepted_columns

wcvp_data = get_all_taxa()
_output_path = 'test_outputs'


class MyTestCase(unittest.TestCase):
    def test_families(self):
        loganiaceae = get_all_taxa(families_of_interest=['Loganiaceae'])
        print(loganiaceae.drop_duplicates(subset=['family'], keep='first')[
                  ['taxon_name', 'family', wcvp_accepted_columns['family']]])
        print(loganiaceae['family'].unique())

        self.assertEqual(list(loganiaceae['family'].unique()),
                         ['Loganiaceae', 'Vitaceae', 'Gentianaceae', 'Acanthaceae', 'Rhamnaceae',
                          'Fabaceae', 'Sabiaceae', 'Boraginaceae', 'Plantaginaceae', 'Caprifoliaceae',
                          'Rutaceae', 'Primulaceae', ])

    def test_acc_cases(self):
        all_taxa_acc = wcvp_data[wcvp_data[wcvp_columns['status']] == 'Accepted']
        pd.testing.assert_series_equal(all_taxa_acc[wcvp_columns['id']],
                                       all_taxa_acc[wcvp_accepted_columns['id']],
                                       check_names=False)
        pd.testing.assert_series_equal(all_taxa_acc['plant_name_id'],
                                       all_taxa_acc['accepted_plant_name_id'],
                                       check_names=False)
        pd.testing.assert_series_equal(all_taxa_acc[wcvp_columns['name']],
                                       all_taxa_acc['accepted_name'],
                                       check_names=False)
        pd.testing.assert_series_equal(all_taxa_acc[wcvp_columns['family']],
                                       all_taxa_acc[wcvp_accepted_columns['family']],
                                       check_names=False)

        pd.testing.assert_series_equal(all_taxa_acc[wcvp_columns['rank']],
                                       all_taxa_acc['accepted_rank'],
                                       check_names=False)

        acc_species = all_taxa_acc[all_taxa_acc['taxon_rank'] == 'Species']
        pd.testing.assert_series_equal(acc_species['accepted_species'],
                                       acc_species['accepted_name'],
                                       check_names=False)
        pd.testing.assert_series_equal(acc_species['accepted_species_ipni_id'],
                                       acc_species[wcvp_accepted_columns['id']],
                                       check_names=False)

        # Currently breaks as accepted parents for some are hybrids, where genus is not given as hybrid
        difference = acc_species['accepted_parent'].compare(acc_species['genus'])
        pd.testing.assert_series_equal(acc_species['accepted_parent'], acc_species['genus'],
                                       check_names=False)

    def test_generic_cases(self):
        genus_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Genus']
        self.assertTrue(genus_accepted_rank['accepted_parent'].isnull().all())
        self.assertTrue(genus_accepted_rank['accepted_parent_ipni_id'].isnull().all())

        self.assertTrue(genus_accepted_rank['accepted_species'].isnull().all())
        self.assertTrue(genus_accepted_rank['accepted_species_ipni_id'].isnull().all())

    def test_species_cases(self):
        species_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Species']
        pd.testing.assert_series_equal(species_accepted_rank['accepted_species'],
                                       species_accepted_rank['accepted_name'],
                                       check_names=False)
        pd.testing.assert_series_equal(species_accepted_rank['accepted_species_ipni_id'],
                                       species_accepted_rank[wcvp_accepted_columns['id']],
                                       check_names=False)

    def test_subspecies_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Subspecies']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

    def test_Convariety_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Convariety']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

    def test_Form_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Form']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

    def test_proles_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'proles']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

    def test_Subform_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Subform']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

    def test_Variety_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Variety']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

    def test_Subvariety_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Subvariety']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

    def test_some_cases(self):
        bry_df = wcvp_data[wcvp_data[wcvp_columns['id']] == '545603-1']
        self.assertListEqual(bry_df[wcvp_accepted_columns['id']].values.tolist(), ['745133-1'])

        genus_df = wcvp_data[wcvp_data[wcvp_columns['id']] == '34250-1']
        self.assertListEqual(genus_df[wcvp_accepted_columns['id']].values.tolist(), ['34250-1'])

    def test_everything_has_ipni_id(self):
        # May fail: Not all plants have an ipni id
        without_ipni_id = wcvp_data[wcvp_data[wcvp_columns['id']].isna()]
        # without_ipni_id.to_csv(os.path.join(_output_path, 'without_id.csv'))
        self.assertEqual(len(without_ipni_id.index), 0)

    def test_everything_with_accepted_name_has_accepted_id(self):
        # May fail: Not all accepted plants have an ipni id
        with_accepted_name = wcvp_data[~wcvp_data[wcvp_accepted_columns['name']].isna()]
        without_acc_id = with_accepted_name[with_accepted_name[wcvp_accepted_columns['id']].isna()]
        without_acc_id.to_csv(os.path.join(_output_path, 'without_acc_id.csv'))
        self.assertEqual(len(without_acc_id.index), 0)

    def test_ranks(self):
        weird_df = wcvp_data[wcvp_data[wcvp_accepted_columns['rank']] == 'nothof.']
        weird_df.to_csv(os.path.join(_output_path, 'nothofs.csv'))
        self.assertEqual(len(weird_df.index), 0)

    def test_unusual_genera(self):
        logan_df = get_all_taxa(families_of_interest=['Loganiaceae'])
        logan_df[logan_df[wcvp_columns['genus']]=='Anthocleista'].to_csv(os.path.join(_output_path,'Anthocleista_in_logania.csv'))
        self.assertIn('Anthocleista', logan_df[wcvp_columns['genus']].values)
if __name__ == '__main__':
    unittest.main()
