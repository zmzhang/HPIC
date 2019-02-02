"""
python -m unittest -v hpic.tests.test_mspd
"""
import unittest

class mspd_test(unittest.TestCase):

    def test_mspd(self):
        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()