import re
import string

import pandas as pd

from wcvp_download import wcvp_accepted_columns, wcvp_columns, hybrid_characters, infraspecific_chars

acc_info_col_names = [wcvp_accepted_columns['ipni_id'],
                      wcvp_accepted_columns['name'],
                      wcvp_accepted_columns['family'],
                      wcvp_accepted_columns['rank'],
                      wcvp_accepted_columns['species'],
                      wcvp_accepted_columns['species_ipni_id'],
                      wcvp_accepted_columns['parent_name'],
                      'accepted_parent_ipni_id']
output_record_col_names = acc_info_col_names + [wcvp_columns['status']]

submitted_name_col_id = 'submitted_name_col_id'
recapitalised_name_col = 'recap_name_col'
lowercase_name_col = 'lower_case_name_col'
tidied_taxon_authors_col = 'tidied_taxon_authors'
submitted_family_name_col_id = 'submitted_family_name_col_id'
unique_submission_index_col = 'unique_submission_index_col'


def get_genus_from_full_name(full_name_beginning_with_genus: str) -> str:
    try:
        genus_plus = remove_whitespace_at_beginning_and_end(full_name_beginning_with_genus)
        y = genus_plus[:1]
        if genus_plus[:1] in hybrid_characters:
            hybrid_char = genus_plus[:1]
            if genus_plus[:2] == hybrid_char + ' ':
                g = genus_plus[2:].partition(' ')[0]
                return hybrid_char + ' ' + g
            else:
                raise ValueError(genus_plus)
        else:
            return genus_plus.partition(' ')[0]

    except (TypeError, AttributeError):
        return full_name_beginning_with_genus


def get_species_from_full_name(full_name_beginning_with_genus: str) -> str:
    genus = get_genus_from_full_name(full_name_beginning_with_genus)

    try:
        try:

            species_plus = full_name_beginning_with_genus.partition(genus)[2]
        except ValueError:
            species_plus = ''
        return get_genus_from_full_name(species_plus)

    except (TypeError, AttributeError):
        return full_name_beginning_with_genus


def remove_whitespace_at_beginning_and_end(value: str) -> str:
    try:
        return value.strip()
    except AttributeError:
        return value


def remove_double_spaces(given_str: str) -> str:
    if pd.isnull(given_str):
        return given_str
    else:
        return " ".join(given_str.split())


def add_space_around_hybrid_chars_and_infraspecific_epithets(value: str):
    try:
        out = value
        for h_char in [c for c in infraspecific_chars if '.' in c and c != 'f.'] + hybrid_characters:
            # These things should be preceded and followed by space
            if h_char in value:
                out = re.sub(r'(?<={})(?=[^\s])'.format(re.escape(h_char)), r' ', out)
                out = re.sub(r'(?={})(?<=[^\s])'.format(re.escape(h_char)), r' ', out)

        return out
    except AttributeError:
        return value


def _capitalize_first_letter_of_taxon(g: str) -> str:
    '''
    Aims to captilise first letter of genus and lower case everything else
    :param g:
    :return:
    '''
    try:

        append_to_beginning = ''
        if any(g.startswith(x + ' ') for x in hybrid_characters):
            hybrid_char = g[:1]
            append_to_beginning = hybrid_char + ' '
            g = g[2:]
        l = g.lower()

        words = l.split()
        capitalised_words = [w.capitalize() if w.endswith('.') and w not in infraspecific_chars else w for w
                             in words]
        for i in range(len(capitalised_words)):
            c = capitalised_words[i]
            if any(c.startswith(p) for p in string.punctuation):
                try:
                    c_list = list(c)
                    c_list[1] = c_list[1].capitalize()
                    capitalised_words[i] = ''.join(c_list)
                except IndexError:
                    pass
        capitalised_words[0] = capitalised_words[0].capitalize()
        return append_to_beginning + ' '.join(capitalised_words)
    except AttributeError:
        return g


def remove_spacelike_chars(given_name: str) -> str:
    try:

        # new_name = given_name.replace('\xa0', ' ')
        # new_name = new_name.replace('\t', ' ')
        new_name = ' '.join(given_name.split())
        return new_name
    except AttributeError:
        return given_name


def remove_fullstop(given_name: str) -> str:
    try:
        return given_name.replace('.', '')
    except AttributeError:
        return given_name


def tidy_families_in_column(df: pd.DataFrame, fam_column: str):
    df[submitted_family_name_col_id] = df[fam_column]
    df[fam_column] = df[fam_column].apply(remove_spacelike_chars)
    df[fam_column] = df[fam_column].apply(remove_whitespace_at_beginning_and_end)
    df[fam_column] = df[fam_column].apply(_capitalize_first_letter_of_taxon)


def tidy_names_in_column(df: pd.DataFrame, name_col: str):
    df[submitted_name_col_id] = df[name_col]
    df[name_col] = df[name_col].apply(remove_spacelike_chars)
    df[name_col] = df[name_col].apply(add_space_around_hybrid_chars_and_infraspecific_epithets)
    df[name_col] = df[name_col].apply(remove_double_spaces)
    df[name_col] = df[name_col].apply(remove_whitespace_at_beginning_and_end)
    df[recapitalised_name_col] = df[name_col].apply(_capitalize_first_letter_of_taxon)


def clean_urn_ids(given_value: str) -> str:
    '''
    Strips urn:lsid:ipni.org:names: from id
    '''

    try:
        if re.search('urn:lsid:ipni.org:names:', given_value):
            pos = re.search('urn:lsid:ipni.org:names:', given_value).end()
            return given_value[pos:]
        else:
            return given_value
    except TypeError:
        return given_value


def tidy_authors(given_string: str):
    # Remove spaces after full stops if full stop isn't part of infraspecific epithet
    # and after the space is a letter
    my_regex = ""
    for ch in infraspecific_chars:
        my_regex += "(?<!" + re.escape(ch.replace('.', '')) + ")"
    my_regex += "\.\s(?=[a-z]|[A-Z]|[)])"
    try:
        return re.sub(my_regex, ".", given_string)
    except TypeError:
        return given_string


def get_word_combinations(given_string: str):
    splitted = given_string.split()
    combinations = []
    for i in range(1, len(splitted) + 1):
        to_append = splitted[0]
        for j in range(1, i):
            to_append = ' '.join([to_append, splitted[j]])
        combinations.append(to_append)
    return combinations
