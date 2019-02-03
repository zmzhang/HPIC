HPIC
===========================
HPIC is a reliable and sensitive method to extract PIC from raw LC-MS dataset directly. Based on the concept that meaningful metabolites tend to generate ions with high density in the m/z and retention space, HDBSCAN is used to find these high density regions and distinguish real signals from noises. It can provide a reasonable and flexible way to determine the m/z tolerance range rather than to give or estimate a fixed m/z tolerance value. 

![Architecture of HPIC method](https://user-images.githubusercontent.com/6937141/52176768-4bdf3900-27f2-11e9-95c6-be94b717fb93.png)


# Install

## [Anaconda Python (Python version 3.7.1)](https://repo.continuum.io/archive/Anaconda3-2018.12-Windows-x86_64.exe)
## Run following command in Anaconda Prompt

		```shell
		pip install git+git://github.com/zmzhang/HPIC@master
		```

# Unittests and Docstrings
Run the unittests in ![tests folder](https://github.com/zmzhang/HPIC/tree/master/hpic/tests)
Usages of major functions are provoided as their ![Docstrings](https://github.com/zmzhang/HPIC/blob/master/hpic/hpic.py#L195)

# Contact

For any questions, please contact:
[zmzhang@csu.edu.cn](mailto:zmzhang@csu.edu.cn)