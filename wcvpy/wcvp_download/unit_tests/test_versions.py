import unittest

from wcvpy.wcvp_download import get_all_taxa


class MyTestCase(unittest.TestCase):

    def test_downloads(self):
        self.assertRaises(ValueError, get_all_taxa,get_new_version=True, version='11')

        v11 = get_all_taxa(version='11')
        v10 = get_all_taxa(version='10')
        v12 = get_all_taxa(version='12')
        assert v10['accepted_name'].unique().size < v12['accepted_name'].unique().size
        assert v11['accepted_name'].unique().size < v12['accepted_name'].unique().size
        assert v10['accepted_name'].unique().size < v11['accepted_name'].unique().size


if __name__ == '__main__':
    unittest.main()
