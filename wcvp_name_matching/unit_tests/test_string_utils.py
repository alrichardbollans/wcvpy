import unittest

import numpy as np

from wcvp_name_matching import get_genus_from_full_name, clean_urn_ids, get_species_from_full_name
from wcvp_name_matching.string_utils import _capitalize_first_letter_of_taxon


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

    def test_capitalising(self):
        test_dict = {'X': 'X', '3': '3', np.nan: np.nan,
                     'Xa': 'Xa', 'Xan and': 'Xan and',
                     'XAn XAn': 'Xan xan', "× genus": '× Genus',
                     "× genus species": '× Genus species',
                     "areaceae": 'Areaceae',
                     "× genus SPecies": '× Genus species',
                     'genus Species auth.': 'Genus species Auth.',
                     'genus auth1. Species auth2.': 'Genus Auth1. species Auth2.',
                     'genus auth1. Genus auth2. Pub.': 'Genus Auth1. genus Auth2. Pub.',
                     'Pub. Genus auth1. Genus auth2. Pub.': 'Pub. genus Auth1. genus Auth2. Pub.',
                     'Pub. Genus auth1. var. Genus auth2. Pub.': 'Pub. genus Auth1. var. genus Auth2. Pub.',
                     'Genus sp1 subsp. sp2': 'Genus sp1 subsp. sp2',
                     'Genus sp1 subs. sp2': 'Genus sp1 Subs. sp2',
                     }
        for t in test_dict:
            out = _capitalize_first_letter_of_taxon(t)
            out_two = _capitalize_first_letter_of_taxon(out)
            if out != out:
                self.assertTrue(np.isnan(test_dict[t]))
            else:
                self.assertEqual(out, test_dict[t])
                self.assertEqual(out_two, test_dict[t])

    def test_clean_urn_ids(self):
        self.assertIsNone(clean_urn_ids(None))
        self.assertTrue(np.isnan(clean_urn_ids(np.nan)))
        correct_dict = {'urn:lsid:ipni.org:names:1': '1', 'urn:lsid:ipni.org:names:': '', '3-1': '3-1',
                        ' urn:lsid:ipni.org:names:2': '2',
                        ' http://ipni.org/urn:lsid:ipni.org:names:2171-1': '2171-1'}

        for k in correct_dict:
            print(k)
            self.assertEqual(correct_dict[k], clean_urn_ids(k))

        self.assertEqual(clean_urn_ids('urn:lsid:ipni.org:names:30479151-2'), '30479151-2')
        self.assertIs(clean_urn_ids(':lsid:ipni.org:names:30479151-2'), ':lsid:ipni.org:names:30479151-2')
        self.assertIs(clean_urn_ids('a'), 'a')
        self.assertIs(clean_urn_ids(''), '')
        self.assertIsInstance(clean_urn_ids('urn:lsid:ipni.org:names:30479151-2'), str)
        self.assertIsInstance(clean_urn_ids(np.NAN), type(np.NAN))


if __name__ == '__main__':
    unittest.main()
