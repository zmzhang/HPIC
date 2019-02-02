"""
python -m unittest -v hpic.tests.test_fileio
"""
import unittest
from hpic.fileio import readms
import os

class fileio_test(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        
    def test_mzxml2(self):
        filename = os.path.join(self.dirname, "tiny1.mzXML2.0.mzXML")
        print(filename)
        ms,intensity,rt,rt_mean_interval = readms(filename)
        print(ms)
        
    def test_mzxml3(self):
        filename = os.path.join(self.dirname, "tiny1.mzXML3.0.mzXML")
        ms,intensity,rt,rt_mean_interval = readms(filename)
        print(ms)

    def test_mzml(self):
        filename = os.path.join(self.dirname, "tiny.pwiz.1.1.mzML")
        ms,intensity,rt,rt_mean_interval = readms(filename)
        print(ms)

    def test_mzdata(self):
        filename = os.path.join(self.dirname, "tiny1.mzData")
        ms,intensity,rt,rt_mean_interval = readms(filename)
        print(ms)

if __name__ == '__main__':
    unittest.main()