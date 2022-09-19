import os
import unittest
import time
import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from automatchnames import id_lookup_wcvp, get_accepted_info_from_ids_in_column, \
    get_accepted_info_from_names_in_column, COL_NAMES

from automatchnames.get_accepted_info import _get_knms_matches_and_accepted_info_from_names_in_column, \
    _find_best_matches_from_multiples, _autoresolve_missing_matches
from taxa_lists import get_all_taxa

wcvp_taxa = get_all_taxa()
unittest_inputs = resource_filename(__name__, 'test_inputs')
unittest_outputs = resource_filename(__name__, 'test_outputs')


class MyTestCase(unittest.TestCase):
    @staticmethod
    def compare_series(s1: pd.Series, s2: pd.Series, **kwargs):
        try:
            pd.testing.assert_series_equal(s1, s2, **kwargs)
        except AssertionError as e:
            print(e)
            print('Problem taxa:')
            for taxa in s1.values:
                if taxa not in s2.values:
                    print(taxa)
            for taxa in s2.values:
                if taxa not in s1.values:
                    print(taxa)
            raise AssertionError


    @staticmethod
    def method_test_on_csv(method_function, input_csv_name, name_col: str, known_acc_name_col: str, fams=None):
        start = time.time()
        test_list = pd.read_csv(os.path.join(unittest_inputs, input_csv_name))
        if fams is not None:
            s = method_function(test_list, name_col, families_of_interest=fams)
        else:
            s = method_function(test_list, name_col)
        s.to_csv(os.path.join(unittest_outputs, input_csv_name))
        MyTestCase.compare_series(s['Accepted_Name'], test_list[known_acc_name_col], check_names=False)
        MyTestCase.compare_series(s[known_acc_name_col], s['Accepted_Name'], check_names=False)
        end = time.time()
        print(f'Time elapsed for method test: {end - start}s')

    @staticmethod
    def all_info_test(input_csv_name: str, name_col):
        '''
        Test all info is correct where input csv has appropriate column names
        :param input_csv_name:
        :param name_col:
        :return:
        '''
        start = time.time()

        test_df = pd.read_csv(os.path.join(unittest_inputs, input_csv_name))
        response = get_accepted_info_from_names_in_column(test_df, name_col)
        response.to_csv(os.path.join(unittest_outputs, input_csv_name))

        for k in COL_NAMES:

            MyTestCase.compare_series(test_df[k], response[COL_NAMES[k]], check_names=False)

        end = time.time()
        print(f'Time elapsed for all info test: {end - start}s')
    def test_id_lookup_wcvp(self):
        cap_dict = id_lookup_wcvp(wcvp_taxa, '44583-2')
        self.assertEqual(cap_dict['Accepted_Name'], 'Capirona macrophylla')
        self.assertEqual(cap_dict['Accepted_ID'], '77210192-1')
        self.assertEqual(cap_dict['Accepted_Rank'], 'Species')
        self.assertEqual(cap_dict['Accepted_Species'], 'Capirona macrophylla')
        self.assertEqual(cap_dict['Accepted_Species_ID'], '77210192-1')

        cap_dict = id_lookup_wcvp(wcvp_taxa, '77210192-1')
        self.assertEqual(cap_dict['Accepted_Name'], 'Capirona macrophylla')
        self.assertEqual(cap_dict['Accepted_ID'], '77210192-1')
        self.assertEqual(cap_dict['Accepted_Rank'], 'Species')
        self.assertEqual(cap_dict['Accepted_Species'], 'Capirona macrophylla')
        self.assertEqual(cap_dict['Accepted_Species_ID'], '77210192-1')

        cap_dict = id_lookup_wcvp(wcvp_taxa, '2217-1')
        self.assertEqual(cap_dict['Accepted_Name'], 'Aspidosperma')
        self.assertEqual(cap_dict['Accepted_ID'], '2217-1')
        self.assertEqual(cap_dict['Accepted_Rank'], 'Genus')
        self.assertTrue(np.isnan(cap_dict['Accepted_Species']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species_ID']))

        cap_dict = id_lookup_wcvp(wcvp_taxa, '41511-1')
        self.assertEqual(cap_dict['Accepted_Name'], 'Aspidosperma')
        self.assertEqual(cap_dict['Accepted_ID'], '2217-1')
        self.assertEqual(cap_dict['Accepted_Rank'], 'Genus')
        self.assertTrue(np.isnan(cap_dict['Accepted_Species']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species_ID']))

        cap_dict = id_lookup_wcvp(wcvp_taxa, '35260-1')
        self.assertEqual(cap_dict['Accepted_Name'], 'Richardia')
        self.assertEqual(cap_dict['Accepted_ID'], '35260-1')
        self.assertEqual(cap_dict['Accepted_Rank'], 'Genus')
        self.assertTrue(np.isnan(cap_dict['Accepted_Species']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species_ID']))

        cap_dict = id_lookup_wcvp(wcvp_taxa, '30000008-2')
        self.assertTrue(np.isnan(cap_dict['Accepted_Name']))
        self.assertTrue(np.isnan(cap_dict['Accepted_ID']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Rank']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species_ID']))

        cap_dict = id_lookup_wcvp(wcvp_taxa, '30000422-2')
        self.assertTrue(np.isnan(cap_dict['Accepted_Name']))
        self.assertTrue(np.isnan(cap_dict['Accepted_ID']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Rank']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species_ID']))

        cap_dict = id_lookup_wcvp(wcvp_taxa, 'abc')
        self.assertTrue(np.isnan(cap_dict['Accepted_Name']))
        self.assertTrue(np.isnan(cap_dict['Accepted_ID']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Rank']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species_ID']))

        cap_dict = id_lookup_wcvp(wcvp_taxa, '')
        self.assertTrue(np.isnan(cap_dict['Accepted_Name']))
        self.assertTrue(np.isnan(cap_dict['Accepted_ID']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Rank']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species_ID']))

        cap_dict = id_lookup_wcvp(wcvp_taxa, ' ')
        self.assertTrue(np.isnan(cap_dict['Accepted_Name']))
        self.assertTrue(np.isnan(cap_dict['Accepted_ID']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Rank']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species_ID']))

        cap_dict = id_lookup_wcvp(wcvp_taxa, np.nan)
        self.assertTrue(np.isnan(cap_dict['Accepted_Name']))
        self.assertTrue(np.isnan(cap_dict['Accepted_ID']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Rank']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species']))
        self.assertTrue(np.isnan(cap_dict['Accepted_Species_ID']))

    def test_get_accepted_info_from_ids_in_column(self):
        unstandardised = pd.read_csv(os.path.join(unittest_inputs, 'powo_medicinal.csv'), index_col=False)
        standardised = pd.read_csv(os.path.join(unittest_inputs, 'powo_medicinal_cleaned.csv'), index_col=False)

        unstandardised.rename(columns={'Unnamed: 0': 'X'}, inplace=True)
        start = time.time()
        x = get_accepted_info_from_ids_in_column(unstandardised, 'fqId')
        end = time.time()
        print(f'Time elapsed for test: {end - start}s')
        print(standardised.columns)
        print(x.columns)
        pd.testing.assert_frame_equal(standardised, x)
        start = time.time()
        garbage = get_accepted_info_from_ids_in_column(unstandardised, 'powo_Snippet')
        end = time.time()
        print(f'Time elapsed for test: {end - start}s')
        unstandardised['Accepted_Name'] = np.nan
        unstandardised['Accepted_ID'] = np.nan
        unstandardised['Accepted_Rank'] = np.nan
        unstandardised['Accepted_Species'] = np.nan
        unstandardised['Accepted_Species_ID'] = np.nan

        pd.testing.assert_frame_equal(unstandardised, garbage)

    def test_find_best_matches(self):
        multiple_names = {
            'submitted': ['Asclepias curassavica', 'Asclepias curassavica', 'Condylocarpon', 'Condylocarpon',
                          'Condylocarpon'],
            'match_state': ['multiple_matches', 'multiple_matches', 'multiple_matches', 'multiple_matches',
                            'multiple_matches'],
            'ipni_id': ['urn:lsid:ipni.org:names:94213-1', 'urn:lsid:ipni.org:names:94212-1',
                        'urn:lsid:ipni.org:names:39836-1', 'urn:lsid:ipni.org:names:328988-2',
                        'urn:lsid:ipni.org:names:11637-1'],
            'matched_name': ['Asclepias curassavica Sp. Pl.: 215 (1753) L. 1753',
                             'Asclepias curassavica Abh. Königl. Ges. Wiss. Göttingen 19: 159 (1874) Griseb. 1874',
                             'Condylocarpus Gen. Pl. Umbell., ed. 2: 202 (1816) Hoffm. 1816',
                             'Condylocarpon Mém. Mus. Hist. Nat. 8: 119 (1822) Desf. 1822',
                             'Condylocarpus Descr. Pinus, ed. 3, 1: 120 (1832) Salisb. ex Lamb. 1832']}

        multiple_names_df = pd.DataFrame(multiple_names)
        multiple_match_records = _find_best_matches_from_multiples(multiple_names_df)

        self.assertEqual(
            len(multiple_match_records[multiple_match_records['submitted'] == 'Asclepias curassavica'].index), 1)
        self.assertEqual(len(multiple_match_records[multiple_match_records['submitted'] == 'Condylocarpon'].index), 1)

        self.assertEqual(
            multiple_match_records.loc[multiple_match_records['submitted'] == 'Condylocarpon', 'Accepted_ID'].iloc[0],
            '328988-2')
        self.assertEqual(
            multiple_match_records.loc[
                multiple_match_records['submitted'] == 'Asclepias curassavica', 'Accepted_ID'].iloc[0],
            '94213-1')

        self.assertEqual(
            multiple_match_records.loc[multiple_match_records['submitted'] == 'Condylocarpon', 'Accepted_Name'].iloc[0],
            'Condylocarpon')
        self.assertEqual(
            multiple_match_records.loc[
                multiple_match_records['submitted'] == 'Asclepias curassavica', 'Accepted_Name'].iloc[0],
            'Asclepias curassavica')

    def test_get_knms_names_and_accepted_info_from_names_in_column(self):
        self.method_test_on_csv(_get_knms_matches_and_accepted_info_from_names_in_column, 'genera_list.csv',
                                'Unlabelled', 'Unlabelled')

        self.method_test_on_csv(_get_knms_matches_and_accepted_info_from_names_in_column, 'species_list.csv',
                                'Labelled', 'Labelled')

    def test_hybrids(self):
        self.method_test_on_csv(get_accepted_info_from_names_in_column, 'hybrid_list.csv', 'name', 'Know_acc_name')

    def test_simple_genera(self):
        self.method_test_on_csv(get_accepted_info_from_names_in_column, 'genera_list.csv', 'Unlabelled', 'Unlabelled',
                                fams=['Rubiaceae',
                                      'Apocynaceae'])

    def test_synonyms(self):
        self.method_test_on_csv(get_accepted_info_from_names_in_column, 'synonym_list.csv', 'syn', 'Know_acc_name')

    def test_known_errors(self):
        self.method_test_on_csv(get_accepted_info_from_names_in_column, 'examples_to_fix.csv', 'name', 'acc_name')

    def test_known_unresolved(self):
        self.method_test_on_csv(get_accepted_info_from_names_in_column, 'unaccepted_examples.csv', 'name', 'acc_name')


    def test_hard_genera(self):
        self.method_test_on_csv(get_accepted_info_from_names_in_column, 'hard_genera_list.csv', 'genera', 'acc_name',
                                fams=['Rubiaceae',
                                      'Apocynaceae'])

    def test_simple_species(self):
        self.method_test_on_csv(get_accepted_info_from_names_in_column, 'species_list.csv', 'Labelled', 'Labelled')

    def test_misspellings(self):
        self.method_test_on_csv(get_accepted_info_from_names_in_column, 'misspellings.csv', 'name', 'acc_name')


    def test_hard_cases_all_info(self):
        self.all_info_test('hard_cases.csv', 'name')

    def test_capitals(self):
        self.all_info_test('test_capitals_db.csv', 'name')

    def test_varieties(self):
        self.all_info_test('test_variety_db.csv', 'name')

    def test_our_data(self):
        self.all_info_test('our_data_test.csv', 'Name')

    def test_subspecies(self):
        self.all_info_test('subspecies.csv', 'name')

    def test_unmatched_resolutions(self):
        self.method_test_on_csv(get_accepted_info_from_names_in_column, 'unmatched.csv', 'submitted', 'acc_name',
                                fams=['Rubiaceae',
                                      'Apocynaceae'])


if __name__ == '__main__':

    unittest.main()
