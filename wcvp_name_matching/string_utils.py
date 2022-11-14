import re

import pandas as pd

acc_info_col_names = ['accepted_ipni_id',
                      'accepted_name',
                      'accepted_family',
                      'accepted_rank',
                      'accepted_species',
                      'accepted_species_ipni_id',
                      'accepted_parent',
                      'accepted_parent_ipni_id',
                      'taxon_status']

hybrid_characters = ["×", "+"]

infraspecific_chars = ['agamosp.', 'convar.', 'ecas.', 'f.', 'grex', 'group', 'lusus', 'microf.', 'microgene',
                       'micromorphe', 'modif.', 'monstr.', 'mut.', 'nid', 'nothof.', 'nothosubsp.',
                       'nothovar.', 'positio', 'proles', 'provar.', 'psp.', 'stirps', 'subf.', 'sublusus',
                       'subproles', 'subsp.', 'subspecioid', 'subvar.', 'unterrasse', 'var.']

submitted_name_col_id = 'submitted_name_col_id'
recapitalised_name_col = 'recap_name_col'


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


def tidy_names_in_column(df: pd.DataFrame, name_col: str):
    df[submitted_name_col_id] = df[name_col]
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
