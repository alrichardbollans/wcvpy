import os
import time
import unittest

from pkg_resources import resource_filename

from OpenRefineMatching import *

unittest_inputs = resource_filename(__name__, 'test_inputs')
unittest_outputs = resource_filename(__name__, 'test_outputs')


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

    def _test_get_names_on_csv(self, input_csv_name, name_col: str):
        in_csv = input_csv_name + '_in.csv'
        correct_csv = input_csv_name + '_correct.csv'
        # assert frames are equal
        start = time.time()
        test_list = pd.read_csv(os.path.join(unittest_inputs, in_csv))

        s = openrefine_match_full_names(test_list, name_col)

        s.to_csv(os.path.join(unittest_outputs, input_csv_name + '.csv'))

        correct_df = pd.read_csv(os.path.join(unittest_inputs, correct_csv))
        pd.testing.assert_frame_equal(s.reset_index(drop=True), correct_df.reset_index(drop=True),
                                      check_dtype=False)
        end = time.time()
        print(f'Time elapsed for method test: {end - start}s')

    def test_examples(self):
        self._test_get_names_on_csv('simple_examples', 'Name')

    def test_gui_diffs(self):
        self._test_get_names_on_csv('gui_diffs', 'Name')

    def test_fuzzy(self):
        self._test_get_names_on_csv('fuzzy_matches', 'Name')

    def test_unmatched(self):
        # Cases where openrefine doesn't find matches
        self._test_get_names_on_csv('unmatched', 'Name')


if __name__ == '__main__':
    unittest.main()
