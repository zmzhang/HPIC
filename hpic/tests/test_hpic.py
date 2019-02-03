import unittest, os, wget, tarfile
from hpic.hpic import hpic

def makedir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

class hpic_test(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        dir_data = os.path.join(dirname, "tmp/data")
        dir_result = os.path.join(dirname, "tmp/result")
        tar_path = os.path.join(dir_data, "Sample-1.tar.bz2")
        self.data_path = os.path.join(dir_data, "MM14_20um.mzxml")
        self.dir_result = dir_result
        if not(os.path.isfile(self.data_path)):
            makedir(dir_data)
            makedir(dir_result)
            wget.download(url="http://msbi.ipb-halle.de/download/Sample-1.tar.bz2", out=dir_data)
            tar = tarfile.open(tar_path, "r:bz2")
            tar.extractall(dir_data)
            tar.close()

    def test_hpic(self):
        peak_list = hpic(self.data_path, self.dir_result, 250, 3)
        peak_list_sorted = peak_list.sort_values(by='intensity', ascending=False)
        self.assertAlmostEqual(peak_list_sorted.iloc[0]['intensity'], 208942.0)

if __name__ == '__main__':
    unittest.main()