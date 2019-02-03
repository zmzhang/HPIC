"""
python -m unittest -v hpic.tests.hpic_test
"""

import unittest
from hpic.hpic import hpic
class hpic_test(unittest.TestCase):

    def test_hpic(self):
        peak_list = hpic('D:/data/Mspos_MM48_20uM_1-B,1_01_14614.mzxml','D:/data/result',250,3)
        peak_list_sorted = peak_list.sort_values(by='intensity', ascending=False)
        self.assertAlmostEqual(peak_list_sorted.iloc[0]['intensity'], 654353)

if __name__ == '__main__':
    unittest.main()