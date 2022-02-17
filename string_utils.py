import re

import numpy as np


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
    # try:

    v = value.rstrip()
    out = v.lstrip()
    return out


def clean_urn_ids(given_value: str) -> str:
    '''
    Strips urn:lsid:ipni.org:names: from id
    '''
    try:
        if re.search('urn:lsid:ipni.org:names:', given_value):
            pos = re.search('urn:lsid:ipni.org:names:', given_value).end()
            return given_value[pos:]
        else:
            return np.nan
    except TypeError:
        return given_value
