import re

import numpy as np
import pandas as pd

COL_NAMES = {'acc_name': 'Accepted_Name',
             'acc_species': 'Accepted_Species',
             'acc_species_id': 'Accepted_Species_ID',
             'acc_id': 'Accepted_ID',
             'acc_rank': 'Accepted_Rank',
             'single_source': 'Source',
             'sources': 'Sources'}


def get_genus_from_full_name(full_name_beginning_with_genus: str) -> str:
    try:
        genus_plus = remove_whitespace_at_beginning_and_end(full_name_beginning_with_genus)
        y = genus_plus[:1]
        if genus_plus[:2] == '× ':
            g = genus_plus[2:].partition(' ')[0]
            return '× ' + g
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

    v = value.rstrip()
    out = v.lstrip()
    return out

def _capitalize_first_letter_of_taxon(g: str, check_string_is_uppercase=False):
    try:
        if check_string_is_uppercase:
            if not g.isupper():
                return g

        append_to_beginning = ''
        if g.startswith('× '):
            append_to_beginning = '× '
            g = g[2:]
        l = g.lower()

        if len([x for x in g if x == " "]) > 1:
            # KNMS does not return matches where authors are incorrectly capitalised
            return g

        return append_to_beginning + l.capitalize()
    except AttributeError:
        return g

def tidy_names_in_column(df:pd.DataFrame,col:str):
    df[col] = df[col].apply(_capitalize_first_letter_of_taxon, check_string_is_uppercase=True)
    df[col] = df[col].apply(remove_whitespace_at_beginning_and_end)

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
