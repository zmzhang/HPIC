"""
python -m unittest -v hpic.tests.test_mspd
"""
import unittest
import numpy as np
from scipy.stats import norm
from hpic.mspd import peaks_detection

class mspd_test(unittest.TestCase):

    def setUp(self):
        size = 200
        rt = np.arange(0, size, 1)
        mz = np.full((size,), 500.00) + np.random.normal(0,0.01,size)
        noise = np.random.normal(0,1,size)
        baseline = rt*0.1
        rv1 = norm(loc = 50, scale = 3)
        rv2 = norm(loc = 100, scale = 5)
        rv3 = norm(loc = 150, scale = 10)
        ints = (rv1.pdf(rt)+rv2.pdf(rt)+rv3.pdf(rt))*500 + noise + baseline
        self.data = np.column_stack((rt,mz,ints))


    def test_mspd(self):
        peaks = peaks_detection(self.data, np.arange(1, 30), 2, 20)
        self.assertAlmostEqual(int(peaks[0][1]), 50, delta=3)
        self.assertAlmostEqual(int(peaks[1][1]), 100, delta=3)
        self.assertAlmostEqual(int(peaks[2][1]), 150, delta=3)

if __name__ == '__main__':
    unittest.main()