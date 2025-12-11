import os
import unittest

import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from wcvpy.wcvp_download import clean_whitespaces_in_names
from wcvpy.wcvp_name_matching import get_genus_from_full_name, clean_urn_ids, get_species_epithet_from_full_name
from wcvpy.wcvp_name_matching.string_utils import _capitalize_first_letter_of_taxon, tidy_authors, \
    get_word_combinations, remove_spacelike_chars, add_space_around_hybrid_chars_and_infraspecific_epithets, \
    get_species_binomial_from_full_name

unittest_inputs = resource_filename(__name__, 'test_inputs')


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
        self.assertIsNone(get_species_epithet_from_full_name(None))
        self.assertTrue(np.isnan(get_species_epithet_from_full_name(np.nan)))
        correct_dict = {'': '', 'x': '', 'x ': '', ' x': '', ' x y': 'y',
                        'Hoodia Sweet ex Decne': 'Sweet', '× Sarcorhiza Anon.': 'Anon.',
                        'Medinilla sarcorhiza Cogn.': 'sarcorhiza', 'Clematis × pinnata': '× pinnata'}

        for k in correct_dict:
            print(k)
            self.assertEqual(correct_dict[k], get_species_epithet_from_full_name(k))
    def test_species_binomial_name(self):
        self.assertIsNone(get_species_binomial_from_full_name(None))
        self.assertTrue(np.isnan(get_species_binomial_from_full_name(np.nan)))
        correct_dict = {'': '', 'g': 'g', 'g ': 'g', ' g': 'g', ' g y': 'g y',
                        'Hoodia Sweet ex Decne': 'Hoodia Sweet', '× Sarcorhiza Anon.': '× Sarcorhiza Anon.',
                        'Medinilla sarcorhiza Cogn.': 'Medinilla sarcorhiza', 'Clematis × pinnata': 'Clematis × pinnata'}

        for k in correct_dict:
            print(k)
            self.assertEqual(correct_dict[k], get_species_binomial_from_full_name(k))

    def test_capitalising(self):
        test_dict = {'': '', ' ': ' ', 'Abies .abies. (L. ) Druce': 'Abies .Abies. (L. ) druce',
                     'Abies abies (L. ) Druce': 'Abies abies (L. ) druce',
                     'X': 'X', '3': '3', np.nan: np.nan,
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
                     'Genus sp1 subs. sp2': 'Genus sp1 Subs. sp2'
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

    def test_tidying_authors(self):
        test_dict = {'Abies abies (L. ) druce': 'Abies abies (L.) druce','Abies abies (L. A druce': 'Abies abies (L.A druce',
                     'Strychnos axillaris': 'Strychnos axillaris', np.nan: np.nan, None: None,
                     'Strychnos axillaris Dalzell. & A.Gibson.': 'Strychnos axillaris Dalzell. & A.Gibson.',
                     'Amsonia tabernaemontana Walter var. gattingeri Woodson': 'Amsonia tabernaemontana Walter var. gattingeri Woodson',
                     'Palicourea gracilenta (Müll. Arg.) Delprete & J. H. Kirkbr.': 'Palicourea gracilenta (Müll.Arg.) Delprete & J.H.Kirkbr.',
                     'ROTHMANNIA ENGLERIANA (K. SCHUM.) KEAV': 'ROTHMANNIA ENGLERIANA (K.SCHUM.) KEAV'}

        for d in test_dict:
            if test_dict[d] == test_dict[d]:
                self.assertEqual(test_dict[d], tidy_authors(d))
            else:
                self.assertTrue(np.isnan(tidy_authors(d)))

    def test_get_word_combintations(self):
        test_dict = {'first second third': ['first', 'first second', 'first second third'],
                     'A   B B  C': ['A', 'A B', 'A B B', 'A B B C']}
        for t in test_dict:
            self.assertEqual(get_word_combinations(t), test_dict[t])

    def test_spacelike(self):
        df = pd.read_csv(os.path.join(unittest_inputs, 'spacelike_cases.csv'))
        for x in df['Name'].values:
            clean = remove_spacelike_chars(x)
            pass

    def test_adding_space(self):
        examples = {'Acokanthera deflersii Schweinf. ex Lewin': 'Acokanthera deflersii Schweinf. ex Lewin',
                    'Asubsp.B': 'A subsp. B'}
        for t in examples:
            self.assertEqual(add_space_around_hybrid_chars_and_infraspecific_epithets(t), examples[t])

    def test_whitespace_removal(self):
        test_dict = {'first second third': 'first second third',
                     'A   B B  C ': 'A B B C', 2: 2, '  ': ''}
        for t in test_dict:
            self.assertEqual(clean_whitespaces_in_names(t), test_dict[t])


if __name__ == '__main__':
    unittest.main()
