import os
import re
import unittest

from wcvp_download import get_all_taxa, wcvp_columns_used_in_direct_matching, infraspecific_chars, \
    hybrid_characters, wcvp_columns, wcvp_accepted_columns

wcvp_data = get_all_taxa()
_output_path = 'test_outputs'


def string_hygeine_tests(df, output_dir):
    things_not_in_checklist = ['  ', ' .',  '\t']
    notin_problem_dfs = []
    for col in wcvp_columns_used_in_direct_matching:
        for unused_string in things_not_in_checklist:
            problem_df = df[df[col].str.contains(unused_string, na=False, regex=False)]

            if len(problem_df.index) > 0:
                problem_df[['accepted_ipni_id', 'accepted_name', col]].to_csv(
                    os.path.join(output_dir,
                                 '(' + unused_string + ')' + col + 'wcvp_problems.csv'))
                notin_problem_dfs.append(problem_df)

    things_that_should_be_followed_by_spaces_in_names = [c for c in infraspecific_chars if
                                                         '.' in c] + hybrid_characters + ['.']
    things_that_should_be_preceded_by_spaces_in_names = [c for c in infraspecific_chars if
                                                         '.' in c] + hybrid_characters
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
            for test_string in things_that_should_be_preceded_by_spaces_in_names:
                problem_df2 = df[
                    df[col].str.contains(r'(?<!\s)(?<!^)\.' + re.escape(test_string), na=False)]
                if len(problem_df2.index) > 0:
                    problem_df2[['accepted_ipni_id', 'accepted_name', col]].to_csv(
                        os.path.join(output_dir, test_string + col + '_prec_wcvp_problems.csv'))
                    problem_dfs.append(problem_df2)
    return problem_dfs, notin_problem_dfs


class MyTestCase(unittest.TestCase):

    def test_things_that_shouldnt_be_in_strings(self):
        p1_df, p2_df = string_hygeine_tests(wcvp_data, _output_path)
        self.assertEqual(len(p1_df), 0)
        self.assertEqual(len(p2_df), 0)

    def test_whitespace(self):
        for c in (wcvp_columns_used_in_direct_matching + list(wcvp_accepted_columns.values())):
            to_check = wcvp_data[~wcvp_data[c].isna()]
            problems = to_check[to_check[c].str.endswith(' ')]
            if len(problems.index) > 0:
                problems.to_csv(os.path.join(_output_path, str(c) + '.csv'))
                raise ValueError


if __name__ == '__main__':
    unittest.main()
