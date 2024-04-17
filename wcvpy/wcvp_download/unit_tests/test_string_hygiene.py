import os
import re
import unittest

import pandas.testing

from wcvpy.wcvp_download import get_all_taxa, wcvp_columns_used_in_direct_matching, infraspecific_chars, \
    hybrid_characters, wcvp_columns, wcvp_accepted_columns, clean_whitespaces_in_names

wcvp_data = get_all_taxa()
_output_path = 'test_outputs'

things_not_in_checklist = ['  ', ' .', '. )', ' )', '\t']


def string_hygeine_tests(df, output_dir):
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

    def test_authors(self):
        test_columns = [wcvp_columns['paranthet_author'], wcvp_columns['primary_author'], wcvp_columns['authors']]
        taxa_to_test = wcvp_data.dropna(subset=test_columns, how='any').copy()
        taxa_to_test['()'] = '(' + taxa_to_test[wcvp_columns['paranthet_author']] + ')'

        for c in test_columns:
            taxa_to_test[c] = taxa_to_test[c].apply(clean_whitespaces_in_names)
        taxa_to_test['test_col'] = taxa_to_test['()'].str.cat([taxa_to_test[wcvp_columns['primary_author']]],
                                                              sep=' ')

        pandas.testing.assert_series_equal(taxa_to_test['test_col'], taxa_to_test[wcvp_columns['authors']], check_names=False)

    # @unittest.skip('Known to fail')
    def test_publication_years(self):
        # This currently fails
        def parse_publication_year(given_string: str):
            if given_string == '(1981 publ. 1082)':  # An exception that returns a valid date
                return 'fails'
            elif given_string in ['(19166)',
                                  '(19543)',
                                  '(19553)',
                                  '(19667)',
                                  '(19983)',
                                  ]:  # Some obivous errors that return valid dates
                return 'fails'

            try:
                return given_string[-5:-1]
            except TypeError:
                return given_string

        from datetime import datetime
        wcvp_data['pyear'] = wcvp_data['first_published'].apply(parse_publication_year)
        format = "%Y"

        failed = []

        def test_year_parsing(given_str):
            if given_str == '' or given_str is None or given_str != given_str:
                pass
            else:
                try:
                    datetime.strptime(given_str, format)
                except ValueError:
                    failed.append(given_str)

        wcvp_data['pyear'].apply(test_year_parsing)

        failed_df = wcvp_data[wcvp_data['pyear'].isin(failed)]
        if len(failed_df.index) > 0:
            failed_df.to_csv(os.path.join(_output_path, 'publication_year_issues.csv'))
            raise Exception


if __name__ == '__main__':
    unittest.main()
