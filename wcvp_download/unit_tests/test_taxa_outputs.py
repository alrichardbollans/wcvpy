import unittest

import pandas as pd

from wcvp_download import get_all_taxa, wcvp_columns

wcvp_data = get_all_taxa()


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
                                       all_taxa_acc['accepted_ipni_id'],
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
                                       acc_species['accepted_ipni_id'],
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
                                       species_accepted_rank['accepted_ipni_id'],
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
        self.assertListEqual(bry_df['accepted_ipni_id'].values.tolist(), ['745133-1'])

        genus_df = wcvp_data[wcvp_data[wcvp_columns['id']] == '34250-1']
        self.assertListEqual(genus_df['accepted_ipni_id'].values.tolist(), ['34250-1'])


if __name__ == '__main__':
    unittest.main()
