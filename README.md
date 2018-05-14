HPIC
===========================
A new method based on HDBSCAN has been developed for PIC extraction from raw LC-MS data, which is a reliable and sensitive feature detection algorithm implemented in Python
# Install
## Required Dependencies
* Python2.7(https://www.python.org/ftp/python/2.7.11/python-2.7.11.amd64.msi)
* install requirements(numpy, scipy, scikit-learn, pyOpenMS and hdbscan) using pip
	```shell
	pip install numpy
	pip install scipy
	pip install scikit-learn
	pip install pyopenms
	pip install hdbscan
	```
# Download
* Download [HPIC](https://codeload.github.com/wangronggit/HPIC/zip/master)
* Unzip it into HPIC directory
# Usage
* Go to HPIC directory
* Download MM14 dataset from this [url](https://msbi.ipb-halle.de/download/Sample-1.tar.bz2) and unzip it
* Run following Python code to extract PICs from it
	```python
	import os  
	os.chdir("D:/HPIC")
	import hpic
	lcfile = "MM14_20um.mzxml"
	feature_out = "D:/feature" 
	feature = hpic.hpic(lcfile,feature_out,500,3,1,15)
	```

