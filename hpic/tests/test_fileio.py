"""
python -m unittest -v hpic.tests.test_fileio
"""
import unittest
from hpic.fileio import readms
import os

class fileio_test(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        
    def test_mzxml(self):
        filename = os.path.join(self.dirname, "tiny1.mzXML2.0.mzXML")
        ms,intensity,rt,rt_mean_interval = readms(filename)
        self.assertEqual(ms[0].shape[0],1313)
        self.assertEqual(intensity[0].shape[0],1313)
        self.assertEqual(len(rt),1)
        self.assertEqual(rt_mean_interval,0.0)

    def test_mzml(self):
        filename = os.path.join(self.dirname, "tiny.pwiz.1.1.mzML")
        ms,intensity,rt,rt_mean_interval = readms(filename)
        self.assertEqual(ms[0].shape[0],15)
        self.assertEqual(intensity[0].shape[0],15)
        self.assertEqual(len(rt),3)
        self.assertTrue(rt_mean_interval<0.0)

    def test_mzdata(self):
        filename = os.path.join(self.dirname, "tiny1.mzData")
        ms,intensity,rt,rt_mean_interval = readms(filename)

if __name__ == '__main__':
    unittest.main()