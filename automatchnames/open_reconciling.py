# Taken from https://gist.github.com/nickynicolson/11fe9e57a198d31fa010fb3feaa65d94
import json

import pandas as pd
import requests

from automatchnames import get_epithets_for_names_in_df

_KEW_RECONCILE_SERVICE_URL = 'http://data1.kew.org/reconciliation/reconcile/IpniName'


def _buildQuery(row, col_name, col_prop_mapper):
    query = {'query': row[col_name]}
    # Add properties
    properties = []
    for p_cname, p_pname in col_prop_mapper.items():
        if p_cname in row:
            if pd.notnull(row[p_cname]):
                property = {"p": p_pname, "pid": p_pname, "v": row[p_cname]}
                properties.append(property)
    query['properties'] = properties
    return query


def _reconcile(row, col, props):
    id = None
    query = _buildQuery(row, col_name=col, col_prop_mapper=props)
    query_json = json.dumps(query)
    # print(query_json)
    # pass to reconciliation service
    url = _KEW_RECONCILE_SERVICE_URL + '?query=' + query_json
    r = requests.get(url=url)
    reco_results = []
    try:
        res = r.json()['result']
        if res is not None:
            for result in res:
                id = result['id']
                name = result['name']
                score = result['score']
                reco_results.append([id, name, score])
    except:
        pass
    return reco_results


def get_reconciliations(in_df: pd.DataFrame, full_name_col: str, genus_col: str = None, species_col: str = None,
                        infraspecies_col: str = None, keep_all=False):
    r = requests.get(_KEW_RECONCILE_SERVICE_URL)
    service_metadata = r.json()

    df_copy = in_df.copy()

    if genus_col is None or species_col is None or infraspecies_col is None:
        df_copy = get_epithets_for_names_in_df(df_copy, name_col=full_name_col)

    if genus_col is None:
        genus_col = 'gn_genus'

    if species_col is None:
        species_col = 'gn_species'

    if infraspecies_col is None:
        infraspecies_col = 'gn_infraspecies'

    col_prop_mapper = dict()
    col_prop_mapper[full_name_col] = 'full_name'
    col_prop_mapper[genus_col] = 'epithet_1'
    col_prop_mapper[species_col] = 'epithet_2'
    col_prop_mapper[infraspecies_col] = 'epithet_3'

    # Reconcile
    df_copy['reco_results'] = df_copy.apply(lambda row: _reconcile(row, col=full_name_col, props=col_prop_mapper),
                                            axis=1)
    # Explode reconciliation results so that each in own row
    df_copy = df_copy.explode('reco_results')
    # Extract ID and name from exploded reconciliation results
    mask = (df_copy.reco_results.notnull())
    df_copy.loc[mask, 'reco_id'] = df_copy[mask].reco_results.apply(lambda x: x[0])
    df_copy.loc[mask, 'reco_name'] = df_copy[mask].reco_results.apply(lambda x: x[1])
    df_copy.loc[mask, 'reco_score'] = df_copy[mask].reco_results.apply(lambda x: x[2])
    # Drop source column as it is no longer needed
    df_copy.drop(columns=['reco_results'], inplace=True)
    # Add a link to the reconciled entity
    mask = (df_copy.reco_id.notnull())
    df_copy.loc[mask, 'reco_link'] = df_copy[mask].reco_id.apply(
        lambda reco_id: service_metadata["view"]["url"].replace('{{id}}', reco_id))

    if not keep_all:
        df_copy.drop_duplicates(keep=False, subset=[full_name_col], inplace=True)
        df_copy.dropna(subset=['reco_id'], inplace=True)
    return df_copy


if __name__ == '__main__':
    from io import StringIO

    data = """
    id|fullname_w_auth|genus|species|infra|bas_auth|pub_auth
    0|Hedera helix|Hedera|helix
    1|Quercus robur L.|Quercus|robur|||L.
    2|Ilex aquifolia|Ilex|aquifolia
    3|Vaccinium vitis-idaea L.|||
    4|Vaccinium L.|||
    """
    df = pd.read_csv(StringIO(data), sep="|")
    acc_df = get_reconciliations(df, 'fullname_w_auth', keep_all=True)
    acc_df.to_csv('test1.csv')

    # df = pd.read_csv('unit_tests/test_inputs/genera_list.csv')
    # acc_df = get_reconciliations(df, 'Unlabelled')
    # acc_df.to_csv('test3.csv')

    # df = pd.read_csv('unit_tests/test_inputs/standardised_order.csv')
    # acc_df = get_reconciliations(df, 'Name')
    # acc_df.to_csv('test2.csv')
