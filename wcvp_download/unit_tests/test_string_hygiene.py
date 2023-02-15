import os
import re
import unittest

import pandas as pd

from wcvp_download import get_all_taxa, wcvp_columns_used_in_direct_matching, infraspecific_chars, \
    hybrid_characters, wcvp_columns

wcvp_data = get_all_taxa()
_output_path = 'test_outputs'


class MyTestCase(unittest.TestCase):
    def dataset_tests(self, df, output_dir):
        things_not_in_checklist_taxon_names = ['  ', ' .']
        notin_problem_dfs = []
        for col in wcvp_columns_used_in_direct_matching:
            for unused_string in things_not_in_checklist_taxon_names:
                problem_df = df[df[col].str.contains(unused_string, na=False, regex=False)]

                if len(problem_df.index) > 0:
                    problem_df[['accepted_ipni_id', 'accepted_name', col]].to_csv(
                        os.path.join(output_dir,
                                     '(' + unused_string + ')' + col + 'wcvp_problems.csv'))
                    notin_problem_dfs.append(problem_df)

        things_that_should_be_followed_by_spaces_in_names = [c for c in infraspecific_chars if
                                                             '.' in c] + hybrid_characters + ['.']
        problem_dfs = []
        for col in wcvp_columns_used_in_direct_matching:
            if col not in [wcvp_columns['authors'], wcvp_columns['paranthet_author'],
                           wcvp_columns['primary_author']]:
                for test_string in things_that_should_be_followed_by_spaces_in_names:
                    problem_df = df[df[col].str.contains(re.escape(test_string) + r'(?!$|\s)', na=False)]

                    if len(problem_df.index) > 0:
                        problem_df[['accepted_ipni_id', 'accepted_name', col]].to_csv(
                            os.path.join(output_dir, test_string + col + 'wcvp_problems.csv'))
                        problem_dfs.append(problem_df)
        self.assertEqual(len(problem_dfs), 0)
        self.assertEqual(len(notin_problem_dfs), 0)

    def test_things_that_shouldnt_be_in_strings(self):
        self.dataset_tests(wcvp_data, os.path.join(_output_path, 'strings'))

    def test_given_WCVP(self):
        wcvp_given_data = get_all_taxa(clean_strings=False)
        self.dataset_tests(wcvp_given_data, os.path.join(_output_path, 'strings', 'wcvp_problems'))


if __name__ == '__main__':
    unittest.main()
