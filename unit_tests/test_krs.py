import unittest
from io import StringIO

import numpy as np
import pandas as pd

from name_matching_cleaning import get_genus_from_full_name, clean_urn_ids, get_species_from_full_name, \
    get_reconciliations


class MyTestCase(unittest.TestCase):

    def test_unmatched_resolutions(self):
        data = """
            id|fullname_w_auth|genus|species|infra|bas_auth|known_full_name
            0|Hedera helix|Hedera|helix
            1|Quercus robur L.|Quercus|robur|||L.
            2|Ilex aquifolia|Ilex|aquifolia
            3|Vaccinium vitis-idaea L.|||
            4|Vaccinium L.|||
            """
        df = pd.read_csv(StringIO(data), sep="|")
        acc_df = get_reconciliations(df, 'fullname_w_auth', keep_all=True)
        pd.testing.assert_series_equal(acc_df['reco_name'], acc_df['known_full_name'])


if __name__ == '__main__':
    unittest.main()
