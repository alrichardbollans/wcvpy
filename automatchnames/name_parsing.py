import json

import pandas as pd
import requests


def get_epithets_for_name(name: str):
    epithet1 = ''
    epithet2 = ''
    epithet3 = ''
    try:
        server_url = 'https://parser.globalnames.org/api/v1/'
        name_request = {"names": [name],
                        "withDetails": True}
        r = requests.post(server_url, json=name_request)

        content = json.loads(r.content.decode('utf-8'))

        for d in content[0]['words']:
            if d['wordType'] == "GENUS":
                epithet1 = d['normalized']
            elif d['wordType'] == "SPECIES":
                epithet2 = d['normalized']
            elif d['wordType'] == "INFRASPECIES":
                epithet3 = d['normalized']
    except:
        pass

    return [epithet1, epithet2, epithet3]


def get_epithets_for_names_in_df(in_df: pd.DataFrame, name_col: str):
    df_copy = in_df.copy()
    df_copy['gnparse_results'] = df_copy.apply(lambda row: get_epithets_for_name(row[name_col]),
                                               axis=1)
    mask = (df_copy.gnparse_results.notnull())
    df_copy.loc[mask, 'gn_genus'] = df_copy[mask].gnparse_results.apply(lambda x: x[0])
    df_copy.loc[mask, 'gn_species'] = df_copy[mask].gnparse_results.apply(lambda x: x[1])
    df_copy.loc[mask, 'gn_infraspecies'] = df_copy[mask].gnparse_results.apply(lambda x: x[2])

    return df_copy

if __name__ == '__main__':

    from io import StringIO

    data = """
    id|fullname_w_auth|genus|species|infra|bas_auth|pub_auth
    0|Hedera helix|Hedera|helix
    1|Quercus robur L.|Quercus|robur|||L.
    2|Medinilla sarcorhiza var. sarcorhiza|Ilex|aquifolia
    3|Vaccinium vitis-idaea L.|||
    4|Clematis × pinnata|||
    """
    df = pd.read_csv(StringIO(data), sep="|")
    df = get_epithets_for_names_in_df(df,'fullname_w_auth')
    print(df)
