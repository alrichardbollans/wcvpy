import os
import time
import unittest

import pandas as pd
from pkg_resources import resource_filename

from wcvpy.OpenRefineMatching import openrefine_match_full_names

from wcvpy.wcvp_name_matching import resolve_openrefine_to_best_matches

unittest_inputs = resource_filename(__name__, 'test_inputs')
unittest_outputs = resource_filename(__name__, 'test_outputs')
from wcvpy.wcvp_download import get_all_taxa

_all_taxa = get_all_taxa()


class MyTestCase(unittest.TestCase):

    def _test_get_names_on_csv(self, input_csv_name, name_col: str):
        in_csv = input_csv_name + '_in.csv'
        correct_csv = input_csv_name + '_correct.csv'
        # assert frames are equal
        start = time.time()
        test_list = pd.read_csv(os.path.join(unittest_inputs, in_csv))

        s = openrefine_match_full_names(test_list, name_col)

        out_df = resolve_openrefine_to_best_matches(s, _all_taxa)

        out_df.to_csv(os.path.join(unittest_outputs, input_csv_name + '.csv'))

        correct_df = pd.read_csv(os.path.join(unittest_inputs, correct_csv), dtype={'plant_name_id': object})

        out_df = out_df.sort_values('Name').reset_index(drop=True)
        correct_df = correct_df.sort_values('Name').reset_index(drop=True)
        pd.testing.assert_frame_equal(out_df, correct_df,
                                      check_dtype=False)
        end = time.time()
        print(f'Time elapsed for method test: {end - start}s')

    def test_example(self):
        self._test_get_names_on_csv('openrefine_examples', 'Name')


if __name__ == '__main__':
    unittest.main()
