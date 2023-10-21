import os
import unittest

import pandas as pd
import pandas.testing

from wcvp_download import get_all_taxa, wcvp_columns, wcvp_accepted_columns, \
    wcvp_columns_used_in_direct_matching

wcvp_data = get_all_taxa()
_output_path = 'test_outputs'


class MyTestCase(unittest.TestCase):

    def test_check_all_have_unique_plant_name_ids(self):
        without_id = wcvp_data[wcvp_data[wcvp_columns['wcvp_id']].isna()]

        self.assertEqual(len(without_id.index), 0)

        duplicates = wcvp_data[wcvp_data.duplicated(subset=[wcvp_columns['wcvp_id']], keep=False)]
        duplicates.to_csv('dups.csv')
        self.assertEqual(len(duplicates.index), 0)

    def test_families(self):
        loganiaceae = get_all_taxa(families_of_interest=['Loganiaceae'])
        print(loganiaceae.drop_duplicates(subset=['family'], keep='first')[
                  ['taxon_name', 'family', wcvp_accepted_columns['family']]])
        # print(loganiaceae['family'].unique())

        self.assertEqual(sorted(list(loganiaceae['family'].unique())),
                         sorted(['Loganiaceae', 'Vitaceae', 'Gentianaceae', 'Acanthaceae', 'Rhamnaceae',
                                 'Fabaceae', 'Sabiaceae', 'Boraginaceae', 'Plantaginaceae', 'Caprifoliaceae',
                                 'Rutaceae', 'Primulaceae', 'Rubiaceae', 'Apocynaceae', 'Penaeaceae']))

    def test_acc_cases(self):
        all_taxa_acc = wcvp_data[wcvp_data[wcvp_columns['status']] == 'Accepted']
        pd.testing.assert_series_equal(all_taxa_acc[wcvp_columns['ipni_id']],
                                       all_taxa_acc[wcvp_accepted_columns['ipni_id']],
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
                                       acc_species[wcvp_accepted_columns['ipni_id']],
                                       check_names=False)

        pd.testing.assert_series_equal(acc_species['accepted_species_id'],
                                       acc_species[wcvp_accepted_columns['wcvp_id']],
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
                                       species_accepted_rank[wcvp_accepted_columns['ipni_id']],
                                       check_names=False)

        pd.testing.assert_series_equal(species_accepted_rank['accepted_species_id'],
                                       species_accepted_rank[wcvp_accepted_columns['wcvp_id']],
                                       check_names=False)

    def test_subspecies_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Subspecies']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_id'],
                                       subsp_accepted_rank['accepted_parent_id'],
                                       check_names=False)

    def test_Convariety_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Convariety']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_id'],
                                       subsp_accepted_rank['accepted_parent_id'],
                                       check_names=False)

    def test_Form_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Form']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_id'],
                                       subsp_accepted_rank['accepted_parent_id'],
                                       check_names=False)

    def test_proles_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'proles']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_id'],
                                       subsp_accepted_rank['accepted_parent_id'],
                                       check_names=False)

    def test_Subform_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Subform']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_id'],
                                       subsp_accepted_rank['accepted_parent_id'],
                                       check_names=False)

    def test_Variety_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Variety']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_id'],
                                       subsp_accepted_rank['accepted_parent_id'],
                                       check_names=False)

    def test_Subvariety_cases(self):
        subsp_accepted_rank = wcvp_data[wcvp_data['accepted_rank'] == 'Subvariety']

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species'],
                                       subsp_accepted_rank['accepted_parent'],
                                       check_names=False)
        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_ipni_id'],
                                       subsp_accepted_rank['accepted_parent_ipni_id'],
                                       check_names=False)

        pd.testing.assert_series_equal(subsp_accepted_rank['accepted_species_id'],
                                       subsp_accepted_rank['accepted_parent_id'],
                                       check_names=False)

    def test_some_cases(self):
        bry_df = wcvp_data[wcvp_data[wcvp_columns['ipni_id']] == '545603-1']
        self.assertListEqual(bry_df[wcvp_accepted_columns['ipni_id']].values.tolist(), ['745133-1'])

        genus_df = wcvp_data[wcvp_data[wcvp_columns['ipni_id']] == '34250-1']
        self.assertListEqual(genus_df[wcvp_accepted_columns['ipni_id']].values.tolist(), ['34250-1'])

    def test_unusual_genera(self):
        logan_df = get_all_taxa(families_of_interest=['Loganiaceae'])
        funny_genera = ['Anthocleista', 'Cissus', 'Ziziphus', 'Rhamnus']
        logan_df[logan_df[wcvp_columns['genus']].isin(['Anthocleista', 'Narcissus', 'Cissus'])].to_csv(
            os.path.join(_output_path, 'Ant_in_logania.csv'))
        for g in funny_genera:
            self.assertIn(g, logan_df[wcvp_columns['genus']].values)

    def test_author_information(self):
        # accepted cases
        # acc name with author = acc name + authors
        accepted = wcvp_data[wcvp_data[wcvp_columns['status']].isin(['Accepted', 'Artificial Hybrid'])]
        accepted['auth_check'] = accepted[wcvp_columns['name']].str.cat(
            [accepted[wcvp_columns['authors']].fillna('')],
            sep=' ').str.strip()
        problems1 = accepted[
            accepted[wcvp_accepted_columns['name_w_author']] != accepted['auth_check']]
        if len(problems1.index) > 0:
            problems1.to_csv('problems1.csv')
            self.assertEqual(len(problems1.index), 0)

        # acc species with author = acc name with author of accepted  species
        sp_check = accepted[['plant_name_id', wcvp_accepted_columns['name_w_author']]]
        sp_check = sp_check.rename(columns={'plant_name_id': 'sp_id',
                                            wcvp_accepted_columns['name_w_author']: 'sp_name_w_author'})
        merged_with_acc_species = pd.merge(accepted, sp_check, left_on=wcvp_accepted_columns['species'],
                                           right_on='sp_id')

        problems2 = merged_with_acc_species[
            merged_with_acc_species[wcvp_accepted_columns['species_w_author']] != merged_with_acc_species[
                'sp_name_w_author']]
        if len(problems2.index) > 0:
            problems2.to_csv('problems2.csv')
            self.assertEqual(len(problems2.index), 0)

        # synonyms
        # acc name with author  = acc name with author of accepted  name
        acc_check = accepted[['plant_name_id', wcvp_accepted_columns['name_w_author']]]
        acc_check = acc_check.rename(columns={'plant_name_id': 'acc_id',
                                              wcvp_accepted_columns['name_w_author']: 'acc_name_w_author'})
        merged_with_acc = pd.merge(wcvp_data, acc_check, left_on=wcvp_accepted_columns['name'],
                                   right_on='acc_id')

        problems3 = merged_with_acc[
            merged_with_acc[wcvp_accepted_columns['name_w_author']] != merged_with_acc[
                'acc_name_w_author']]
        if len(problems3.index) > 0:
            problems3.to_csv('problems3.csv')
            self.assertEqual(len(problems3.index), 0)

        # subspecies and varieties
        # acc parent w author = accepted species with author
        subsp_vars = wcvp_data[
            wcvp_data[wcvp_accepted_columns['rank']].isin(["Form", "Subspecies", "Subvariety", "Variety"])]

        problems4 = subsp_vars[
            subsp_vars[wcvp_accepted_columns['species_w_author']] != subsp_vars['accepted_parent_w_author']]

        # will catch some cases with no accepted parent info
        problems4 = problems4[~(problems4[wcvp_accepted_columns['species_w_author']].isna() & problems4[
            'accepted_parent_w_author'].isna())]

        if len(problems4.index) > 0:
            problems4.to_csv('problems4.csv')
            self.assertEqual(len(problems4.index), 0)
        # species
        # acc name with author  = accepted species with author
        species = wcvp_data[
            wcvp_data[wcvp_accepted_columns['rank']].isin(["Species"])]

        problems5 = species[
            species[wcvp_accepted_columns['species_w_author']] != species[
                wcvp_accepted_columns['name_w_author']]]

        if len(problems5.index) > 0:
            problems5.to_csv('problems5.csv')
            self.assertEqual(len(problems5.index), 0)

    def test_consistency_of_author_info(self):
        # Check that when some paranthet or prim author info is given this is reflected in taxon authors,
        # this is important for given acc names with authors
        without_authors = wcvp_data[
            wcvp_data[wcvp_columns['authors']].isna() & ~(wcvp_data[
                                                              wcvp_columns[
                                                                  'paranthet_author']].isna() &
                                                          wcvp_data[
                                                              wcvp_columns[
                                                                  'primary_author']].isna())]

        self.assertEqual(len(without_authors.index), 0)

    def test_getting_ranks(self):
        sp = get_all_taxa(ranks=['Species'])
        all_sp = wcvp_data[(wcvp_data[wcvp_columns['rank']].isin(['Species'])) | (
            wcvp_data[wcvp_accepted_columns['rank']].isin(['Species']))]

        pandas.testing.assert_frame_equal(sp,all_sp)


if __name__ == '__main__':
    unittest.main()
