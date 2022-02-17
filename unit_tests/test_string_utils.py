import unittest

import numpy as np

from name_matching_cleaning import get_genus_from_full_name, clean_urn_ids, get_species_from_full_name


class MyTestCase(unittest.TestCase):

    def test_genus_name(self):
        self.assertIsNone(get_genus_from_full_name(None))
        self.assertTrue(np.isnan(get_genus_from_full_name(np.nan)))
        correct_dict = {'': '', 'x': 'x', 'x ': 'x', ' x': 'x', ' x x': 'x',
                        'Hoodia Sweet ex Decne.': 'Hoodia', '× Sarcorhiza Anon.': '× Sarcorhiza',
                        'Medinilla sarcorhiza Cogn.': 'Medinilla', 'Clematis × pinnata': 'Clematis'}

        for k in correct_dict:
            print(k)
            self.assertEqual(correct_dict[k], get_genus_from_full_name(k))
    def test_species_name(self):
        self.assertIsNone(get_species_from_full_name(None))
        self.assertTrue(np.isnan(get_species_from_full_name(np.nan)))
        correct_dict = {'': '', 'x': '', 'x ': '', ' x': '', ' x y': 'y',
                        'Hoodia Sweet ex Decne': 'Sweet', '× Sarcorhiza Anon.': 'Anon.',
                        'Medinilla sarcorhiza Cogn.': 'sarcorhiza', 'Clematis × pinnata': '× pinnata'}

        for k in correct_dict:
            print(k)
            self.assertEqual(correct_dict[k], get_species_from_full_name(k))

    def test_clean_urn_ids(self):
        self.assertIsNone(clean_urn_ids(None))
        self.assertTrue(np.isnan(clean_urn_ids(np.nan)))
        correct_dict = {'urn:lsid:ipni.org:names:1': '1', 'urn:lsid:ipni.org:names:': '',
                        ' urn:lsid:ipni.org:names:2': '2', ' http://ipni.org/urn:lsid:ipni.org:names:2171-1': '2171-1'}

        for k in correct_dict:
            print(k)
            self.assertEqual(correct_dict[k], clean_urn_ids(k))


if __name__ == '__main__':
    unittest.main()
