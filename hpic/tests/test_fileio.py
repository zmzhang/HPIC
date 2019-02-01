import unittest

"""
python -m unittest -v hpic.tests.test_fileio
"""
class fileio_test(unittest.TestCase):

    def test_mzxml(self):
        self.assertEqual(1, 1)

    def test_mzml(self):
        self.assertTrue(True)

    def test_mzdata(self):
        self.assertEqual(5, 5)

if __name__ == '__main__':
    unittest.main()