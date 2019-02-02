"""
python -m unittest -v hpic.tests.hpic_test
"""

import unittest
class hpic_test(unittest.TestCase):

    def test_hpic(self):
        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()