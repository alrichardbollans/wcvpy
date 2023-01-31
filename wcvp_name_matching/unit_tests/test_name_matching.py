import os
import unittest
import time
from typing import List

import numpy as np
import pandas as pd
import pandas.testing
from pkg_resources import resource_filename

from wcvp_name_matching import id_lookup_wcvp, get_accepted_wcvp_info_from_ids_in_column, \
    get_accepted_info_from_names_in_column, acc_info_col_names, clean_urn_ids
from wcvp_name_matching.get_accepted_info import _get_knms_matches_and_accepted_info_from_names_in_column, \
    _find_best_matches_from_multiple_knms_matches

from wcvp_download import get_all_taxa, wcvp_accepted_columns

wcvp_taxa = get_all_taxa()

unittest_inputs = resource_filename(__name__, 'test_inputs')
unittest_outputs = resource_filename(__name__, 'test_outputs')

# columns used in testing csvs
test_columns = {'acc_id': wcvp_accepted_columns['id'],
                'acc_name': 'accepted_name',
                'acc_fam': wcvp_accepted_columns['family'],
                'acc_rank': 'accepted_rank',
                'acc_parent': 'accepted_parent',
                'acc_species': 'accepted_species',
                'acc_species_id': 'accepted_species_ipni_id',
                'tax_stat': 'taxon_status'}


# # Check adding familiy doesn't ruin resolutions outside family
class MyTestCase(unittest.TestCase):

    def compare_series(self, s1: pd.Series, s2: pd.Series):
        try:
            if s1.isnull().all():
                s1 = s1.astype(object)
            if s2.isnull().all():
                s2 = s2.astype(object)
            self.assertTrue(s1.equals(s2))
        except AssertionError as e:
            print(e)
            print(s1.name)
            print(s2.name)

            problem_taxa = [p for p in s1.values if p not in s2.values] + [p for p in s2.values if
                                                                           p not in s1.values]
            if len(problem_taxa) > 0:
                print(f'Problem taxa: {problem_taxa}')
            else:
                print('Likely incorrect ordering')
            raise AssertionError

    def _test_get_acc_info_names_on_csv(self, input_csv_name, name_col: str, known_acc_name_col: str,
                                        fams: List[str] = None, **kwargs):
        start = time.time()
        test_list = pd.read_csv(os.path.join(unittest_inputs, input_csv_name))

        if fams is None:
            s = get_accepted_info_from_names_in_column(test_list, name_col, **kwargs)
        else:
            s = get_accepted_info_from_names_in_column(test_list, name_col, families_of_interest=fams,
                                                       **kwargs)
        s.to_csv(os.path.join(unittest_outputs, input_csv_name))
        self.compare_series(s['accepted_name'], test_list[known_acc_name_col])
        self.compare_series(s[known_acc_name_col], s['accepted_name'])
        end = time.time()
        print(f'Time elapsed for method test: {end - start}s')

    def _test_get_knms_names_on_csv(self, input_csv_name, name_col: str, known_acc_name_col: str,
                                    taxa_df: pd.DataFrame = None):
        start = time.time()
        test_list = pd.read_csv(os.path.join(unittest_inputs, input_csv_name))

        if taxa_df is None:
            s = _get_knms_matches_and_accepted_info_from_names_in_column(test_list, name_col, name_col,
                                                                         wcvp_taxa)
        else:
            s = _get_knms_matches_and_accepted_info_from_names_in_column(test_list, name_col, name_col,
                                                                         taxa_df)
        s.to_csv(os.path.join(unittest_outputs, input_csv_name))
        # Doesn't preserve order
        self.assertEqual(sorted(list(s['accepted_name'])), sorted(list(test_list[known_acc_name_col])))
        end = time.time()
        print(f'Time elapsed for method test: {end - start}s')

    def all_info_test(self, input_csv_name: str, name_col: str, **kwargs):
        '''
        Test all info is correct where input csv has appropriate column names
        :param input_csv_name:
        :param name_col:
        :return:
        '''
        start = time.time()

        test_df = pd.read_csv(os.path.join(unittest_inputs, input_csv_name))
        response = get_accepted_info_from_names_in_column(test_df, name_col, **kwargs)
        response.to_csv(os.path.join(unittest_outputs, input_csv_name))

        for k in test_columns:
            self.compare_series(test_df[k], response[test_columns[k]])

        end = time.time()
        print(f'Time elapsed for all info test: {end - start}s')

        return test_df, response

    def all_info_id_test(self, input_csv_name: str, id_col: str, taxa_df: pd.DataFrame):
        '''
        Test all info is correct where input csv has appropriate column names
        :param taxa_df:
        :param id_col:
        :param input_csv_name:
        :return:
        '''
        start = time.time()

        test_df = pd.read_csv(os.path.join(unittest_inputs, input_csv_name))
        response = get_accepted_wcvp_info_from_ids_in_column(test_df, id_col, taxa_df)
        response.to_csv(os.path.join(unittest_outputs, input_csv_name))

        for k in test_columns:
            self.compare_series(test_df[k], response[test_columns[k]])

        end = time.time()
        print(f'Time elapsed for all info test: {end - start}s')

    def test_id_lookup_wcvp(self):
        cap_record = id_lookup_wcvp(wcvp_taxa, '44583-2')
        self.assertEqual(cap_record['accepted_name'].iloc[0], 'Capirona macrophylla')
        self.assertEqual(cap_record[wcvp_accepted_columns['id']].iloc[0], '77210192-1')
        self.assertEqual(cap_record['accepted_rank'].iloc[0], 'Species')
        self.assertEqual(cap_record['accepted_species'].iloc[0], 'Capirona macrophylla')
        self.assertEqual(cap_record['accepted_species_ipni_id'].iloc[0], '77210192-1')

        cap_record = id_lookup_wcvp(wcvp_taxa, '77210192-1')
        self.assertEqual(cap_record['accepted_name'].iloc[0], 'Capirona macrophylla')
        self.assertEqual(cap_record[wcvp_accepted_columns['id']].iloc[0], '77210192-1')
        self.assertEqual(cap_record['accepted_rank'].iloc[0], 'Species')
        self.assertEqual(cap_record['accepted_species'].iloc[0], 'Capirona macrophylla')
        self.assertEqual(cap_record['accepted_species_ipni_id'].iloc[0], '77210192-1')

        cap_record = id_lookup_wcvp(wcvp_taxa, '2217-1')
        self.assertEqual(cap_record['accepted_name'].iloc[0], 'Aspidosperma')
        self.assertEqual(cap_record[wcvp_accepted_columns['id']].iloc[0], '2217-1')
        self.assertEqual(cap_record['accepted_rank'].iloc[0], 'Genus')
        self.assertTrue(np.isnan(cap_record['accepted_species'].iloc[0]))
        self.assertTrue(np.isnan(cap_record['accepted_species_ipni_id'].iloc[0]))

        cap_record = id_lookup_wcvp(wcvp_taxa, '41511-1')
        self.assertEqual(cap_record['accepted_name'].iloc[0], 'Aspidosperma')
        self.assertEqual(cap_record[wcvp_accepted_columns['id']].iloc[0], '2217-1')
        self.assertEqual(cap_record['accepted_rank'].iloc[0], 'Genus')
        self.assertTrue(np.isnan(cap_record['accepted_species'].iloc[0]))
        self.assertTrue(np.isnan(cap_record['accepted_species_ipni_id'].iloc[0]))

        cap_record = id_lookup_wcvp(wcvp_taxa, '35260-1')
        self.assertEqual(cap_record['accepted_name'].iloc[0], 'Richardia')
        self.assertEqual(cap_record[wcvp_accepted_columns['id']].iloc[0], '35260-1')
        self.assertEqual(cap_record['accepted_rank'].iloc[0], 'Genus')
        self.assertTrue(np.isnan(cap_record['accepted_species'].iloc[0]))
        self.assertTrue(np.isnan(cap_record['accepted_species_ipni_id'].iloc[0]))

        cap_record = id_lookup_wcvp(wcvp_taxa, '30000008-2')
        self.assertEqual(len(cap_record.index), 0)

        cap_record = id_lookup_wcvp(wcvp_taxa, '30000422-2')
        self.assertEqual(len(cap_record.index), 0)

        cap_record = id_lookup_wcvp(wcvp_taxa, 'abc')
        self.assertEqual(len(cap_record.index), 0)

        cap_record = id_lookup_wcvp(wcvp_taxa, '')
        self.assertEqual(len(cap_record.index), 0)

        cap_record = id_lookup_wcvp(wcvp_taxa, ' ')
        self.assertEqual(len(cap_record.index), 0)

        cap_record = id_lookup_wcvp(wcvp_taxa, np.nan)
        self.assertEqual(len(cap_record.index), 0)

    # @unittest.skip
    def test_get_accepted_info_from_ids_in_column(self):
        self.all_info_id_test('powo_medicinal_cleaned.csv', 'fqId', wcvp_taxa)

        standardised = pd.read_csv(os.path.join(unittest_inputs, 'powo_medicinal_cleaned.csv'),
                                   index_col=False)

        start = time.time()
        garbage = get_accepted_wcvp_info_from_ids_in_column(standardised, 'powo_Snippet', wcvp_taxa)
        end = time.time()
        print(f'Time elapsed for test: {end - start}s')

        for k in test_columns:
            self.assertTrue(garbage[test_columns[k]].isnull().all())

    # @unittest.skip
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
        multiple_names_df['ipni_id'] = multiple_names_df['ipni_id'].apply(clean_urn_ids)
        multiple_names_df = get_accepted_wcvp_info_from_ids_in_column(multiple_names_df, 'ipni_id', wcvp_taxa)
        multiple_match_records = _find_best_matches_from_multiple_knms_matches(multiple_names_df, 'submitted')

        self.assertEqual(
            len(multiple_match_records[multiple_match_records['submitted'] == 'Asclepias curassavica'].index),
            1)
        self.assertEqual(
            len(multiple_match_records[multiple_match_records['submitted'] == 'Condylocarpon'].index), 1)

        self.assertEqual(
            multiple_match_records.loc[
                multiple_match_records['submitted'] == 'Condylocarpon', wcvp_accepted_columns['id']].iloc[
                0],
            '328988-2')
        self.assertEqual(
            multiple_match_records.loc[
                multiple_match_records['submitted'] == 'Asclepias curassavica', wcvp_accepted_columns[
                    'id']].iloc[0],
            '94213-1')

        self.assertEqual(
            multiple_match_records.loc[
                multiple_match_records['submitted'] == 'Condylocarpon', 'accepted_name'].iloc[0],
            'Condylocarpon')
        self.assertEqual(
            multiple_match_records.loc[
                multiple_match_records['submitted'] == 'Asclepias curassavica', 'accepted_name'].iloc[0],
            'Asclepias curassavica')

    # @unittest.skip
    def test_get_knms_names_and_accepted_info_from_names_in_column(self):
        self._test_get_knms_names_on_csv(
            'genera_list.csv',
            'Unlabelled', 'Unlabelled')

        self._test_get_knms_names_on_csv(
            'species_list.csv',
            'Labelled', 'Labelled')

    def test_hybrids(self):
        self._test_get_acc_info_names_on_csv('hybrid_list.csv', 'Name',
                                             'Know_acc_name')

    def test_simple_genera(self):

        self._test_get_acc_info_names_on_csv('genera_list.csv',
                                             'Unlabelled',
                                             'Unlabelled', fams=['Rubiaceae',
                                                                 'Apocynaceae'])

    def test_synonyms(self):
        self._test_get_acc_info_names_on_csv('synonym_list.csv', 'syn',
                                             'Know_acc_name')

    def test_examples_to_fix(self):
        self._test_get_acc_info_names_on_csv('examples_to_fix.csv',
                                             'Name',
                                             'acc_name')

    def test_fam_examples_to_fix(self):
        self._test_get_acc_info_names_on_csv('fam_examples_to_fix.csv',
                                             'Name',
                                             'acc_name', family_column='Family')

    def test_spacelike(self):
        self._test_get_acc_info_names_on_csv('spacelike_cases.csv',
                                             'Name',
                                             'acc_name', family_column='Family')

    def test_hard_genera(self):
        self._test_get_acc_info_names_on_csv('hard_genera_list.csv',
                                             'genera',
                                             'acc_name')

    def test_simple_species(self):
        self._test_get_acc_info_names_on_csv('species_list.csv',
                                             'Labelled',
                                             'Labelled')

    def test_misspellings(self):
        self._test_get_acc_info_names_on_csv('misspellings.csv',
                                             'Name',
                                             'acc_name')

    def test_hard_cases_all_info(self):
        self.all_info_test('hard_cases.csv', 'Name')

    # @unittest.skip
    def test_capitals(self):
        self.all_info_test('test_capitals_db.csv', 'Name')

    def test_varieties(self):
        self.all_info_test('test_variety_db.csv', 'Name')

    def test_our_data(self):
        self.all_info_test('our_data_example.csv', 'Name')

    def test_our_data_with_fams(self):
        self.all_info_test('our_data_example_with_fams.csv', 'Name', family_column='Family')

    def test_subspecies(self):
        self.all_info_test('subspecies.csv', 'Name')

    def test_unmatched_resolutions(self):
        self._test_get_acc_info_names_on_csv('unmatched.csv',
                                             'submitted_to_knms',
                                             'acc_name', fams=['Rubiaceae',
                                                               'Apocynaceae'])

    def test_some_cases(self):
        test_df = pd.DataFrame(['Anthocleista brieyi'], columns=['Name'])
        acc_df = get_accepted_info_from_names_in_column(test_df, name_col='Name')
        self.assertListEqual(acc_df[wcvp_accepted_columns['family']].values.tolist(), ['Rubiaceae'])

    def test_direct_match(self):
        self.all_info_test('wcvp_level.csv', 'Name', match_level='direct')
        try:
            get_accepted_info_from_names_in_column(pd.DataFrame(), 'Name', match_level='garbage')
        except ValueError:
            pass
        else:
            raise ValueError

    # @unittest.skip
    def test_fam_testing(self):
        self.all_info_test('family_test.csv', 'Name', family_column='Family')

    def test_powo_example_wtih_fams(self):
        self.all_info_test('powo_medicinal_cleaned.csv', 'Name', family_column='Family')

    def test_manual_additions(self):
        manual_template = os.path.join('..', 'matching data', 'manual_match_template.csv')
        test_df, result = self.all_info_test('manual_test.csv', 'Name', manual_resolution_csv=manual_template)

        self.assertTrue((result['matched_by'] == 'manual').all())

    def test_other_column_values_stay_the_same(self):
        test_df, result = self.all_info_test('powo_medicinal_cleaned.csv', 'Name')
        result = result.drop(columns=acc_info_col_names + ['matched_by'])
        pandas.testing.assert_frame_equal(test_df, result)

    def test_knms_synonyms(self):
        self.all_info_test('knms_synonyms.csv', 'Name')

    def test_ecbot(self):

        test_df = pd.read_excel(os.path.join(unittest_inputs, 'Adam-ex-data-ZI.xlsx'))
        response = get_accepted_info_from_names_in_column(test_df, 'Taxon ', family_column='Family')
        response.to_csv(os.path.join(unittest_outputs, 'ecbot.csv'))

    def test_autoresolutions(self):
        self.all_info_test('autoresolutions.csv', 'Name')

        self.all_info_test('autoresoltuions_with_fams.csv', 'Name', family_column='Family')

    def test_wcvp_families(self):
        # Currently fails
        # See readme on this issue
        # One must be careful with usage of families_of_interest --- there are some genera that are often considered to be in
        # different families. For example, 'Anthocleista' is an accepted genus in Gentianaceae but is often considered
        # to be in Loganiaceae. If you try to match e.g. 'Anthocleista procera' within Loganiaceae it will be unresolved.
        # Note however that examples like 'Anthocleista brieyi' are synonyms within Loganiaceae whose accepted family
        # is rubiaceae and in this case the program will find the match
        self.all_info_test('outside_family_test.csv', 'Name', families_of_interest=['Loganiaceae'])

    def test_matched_to_example(self):
        test_df = pd.read_csv(os.path.join(unittest_inputs, 'matched_name_example.csv'))
        response = get_accepted_info_from_names_in_column(test_df, 'Name')
        response.to_csv(os.path.join(unittest_outputs, 'matched_name_example.csv'))


        self.compare_series(test_df['m_name'], response['matched_name'])


if __name__ == '__main__':
    unittest.main()
