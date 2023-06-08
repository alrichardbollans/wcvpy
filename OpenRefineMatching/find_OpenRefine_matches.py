import json

import pandas as pd
import requests

_openrefine_ipni_service_url = 'http://data1.kew.org/reconciliation/reconcile/IpniName'
_openrefine_ipni_request = requests.get(_openrefine_ipni_service_url)
_openrefine_ipni_service_metadata = _openrefine_ipni_request.json()

reco_submitted_name_col_id = 'reco_submitted_name_col_id'


def _buildQuery(row: pd.Series, full_name_col: str):
    query = {'query': row[full_name_col]}
    # # Add properties - not used for full names
    # properties = []
    # for p_cname, p_pname in col_prop_mapper.items():
    #     if pd.notnull(row[p_cname]):
    #         property = {"p": p_pname, "pid": p_pname, "v": row[p_cname]}
    #         properties.append(property)
    # query['properties'] = properties
    return query


def _reconcile(row, full_name_col):
    query = _buildQuery(row, full_name_col=full_name_col)
    query_json = json.dumps(query)
    # pass to reconciliation service
    url = _openrefine_ipni_service_url + '?query=' + query_json
    r = requests.get(url=url)
    reco_results = []
    try:
        res = r.json()['result']
        if res is not None:
            for result in res:
                id = result['id']
                name = result['name']
                reco_results.append([id, name, result['score']])
    except json.decoder.JSONDecodeError:
        pass
    return reco_results


def openrefine_match_full_names(df: pd.DataFrame, full_name_col: str,
                                output_csv: str = None) -> pd.DataFrame:
    out_df = df.copy(deep=True)
    out_df[reco_submitted_name_col_id] = df[full_name_col]
    out_df = out_df.drop_duplicates(subset=[reco_submitted_name_col_id])

    # Reconcile
    out_df['reco_results'] = out_df.apply(
        lambda row: _reconcile(row, full_name_col=full_name_col),
        axis=1)
    # Explode reconciliation results so that each in own row
    out_df = out_df.explode('reco_results')
    # Extract ID, name and score from exploded reconciliation results
    mask = (out_df.reco_results.notnull())
    out_df.loc[mask, 'reco_id'] = out_df[mask].reco_results.apply(lambda x: x[0])
    out_df.loc[mask, 'reco_name'] = out_df[mask].reco_results.apply(lambda x: x[1])
    out_df.loc[mask, 'reco_score'] = out_df[mask].reco_results.apply(lambda x: x[2])
    # Drop source column as it is no longer needed
    out_df.drop(columns=['reco_results'], inplace=True)
    # Add a link to the reconciled entity
    # mask = (out_df.reco_id.notnull())
    # out_df.loc[mask, 'reco_link'] = out_df[mask].reco_id.apply(
    #     lambda id: _openrefine_ipni_service_metadata["view"]["url"].replace('{{id}}', id))

    # out_df = get_accepted_wcvp_info_from_ipni_ids_in_column(out_df, 'reco_id', _all_taxa)
    if output_csv is not None:
        out_df.to_csv(output_csv)

    return out_df
