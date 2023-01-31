import re

import pandas as pd

from wcvp_download import wcvp_accepted_columns, wcvp_columns

acc_info_col_names = [wcvp_accepted_columns['id'],
                      wcvp_accepted_columns['name'],
                      wcvp_accepted_columns['family'],
                      wcvp_accepted_columns['rank'],
                      wcvp_accepted_columns['species'],
                      wcvp_accepted_columns['species_id'],
                      wcvp_accepted_columns['parent_name'],
                      'accepted_parent_ipni_id',
                      wcvp_columns['status']]

hybrid_characters = ["×", "+"]

infraspecific_chars = ['agamosp.', 'convar.', 'ecas.', 'f.', 'grex', 'group', 'lusus', 'microf.', 'microgene',
                       'micromorphe', 'modif.', 'monstr.', 'mut.', 'nid', 'nothof.', 'nothosubsp.',
                       'nothovar.', 'positio', 'proles', 'provar.', 'psp.', 'stirps', 'subf.', 'sublusus',
                       'subproles', 'subsp.', 'subspecioid', 'subvar.', 'unterrasse', 'var.']

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


def remove_whitespace_at_beginning_and_end(value):
    try:
        v = value.rstrip()
        out = v.lstrip()
        return out
    except AttributeError:
        return value


def add_space_after_hybrid_char(value: str):
    try:
        out = value
        for h_char in hybrid_characters:
            if h_char in value:
                out = re.sub(r'(?<=[{}])(?=[^\s])'.format(h_char), r' ', out)
        return out
    except AttributeError:
        return value


def _capitalize_first_letter_of_taxon(g: str):
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

        words = l.split(' ')
        capitalised_words = [w.capitalize() if (w.endswith('.') and w not in infraspecific_chars) else w for w
                             in words]
        capitalised_words[0] = capitalised_words[0].capitalize()
        return append_to_beginning + ' '.join(capitalised_words)
    except AttributeError:
        return g


def remove_spacelike_chars(given_name: str):
    try:
        new_name = given_name.replace('\xa0', ' ')
        new_name = new_name.replace('\t', ' ')

        return new_name
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
    df[name_col] = df[name_col].apply(add_space_after_hybrid_char)
    df[name_col] = df[name_col].apply(remove_whitespace_at_beginning_and_end)
    df[recapitalised_name_col] = df[name_col].apply(_capitalize_first_letter_of_taxon)
    df[lowercase_name_col] = df[name_col].str.lower()


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
    my_regex = "\.\s(?=[a-z]|[A-Z])"
    for ch in infraspecific_chars:
        ch_reg_form = re.escape(ch)
        my_regex += "(?<!" + ch_reg_form + "\s)"
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
