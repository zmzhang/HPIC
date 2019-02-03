"""
python -m unittest -v hpic.tests.hpic_test
"""

import unittest
from hpic.hpic import hpic
class hpic_test(unittest.TestCase):

    def test_hpic(self):
        hpic(r'D:/data/Mspos_MM48_20uM_1-B,1_01_14614.mzxml','D:/data/result',250,3)
        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()