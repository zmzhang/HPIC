"""
Module to parse mzXML, mzML and mzData files based on pyopenms
"""
# pylint: disable=no-name-in-module
import os
import numpy as np
from pyopenms import MSExperiment, MzXMLFile, MzMLFile, MzDataFile


def readms(file_path):
    """
    Read mzXML, mzML and mzData files.

    Arguments:
        file_path: path to the dataset locally

    Returns:
        Tuple of Numpy arrays: (m/z, intensity, retention time, mean interval of retention time).
    
    Examples:
        >>> from hpic.fileio import readms
        >>> ms,intensity,rt,rt_mean_interval = readms("MM14_20um.mzxml")
    """
    ms_format = os.path.splitext(file_path)[1]
    ms_format = ms_format.lower()
    msdata = MSExperiment()
    if ms_format == '.mzxml':
        file = MzXMLFile()
    elif ms_format == '.mzml':
        file = MzMLFile()
    elif ms_format == '.mzdata':
        file = MzDataFile()
    else:
        raise Exception('ERROR: %s is wrong format' % file_path)
    file.load(r'%s' % file_path, msdata)
    m_s = []
    intensity = []
    r_t = []
    for spectrum in msdata:
        if spectrum.getMSLevel() == 1:
            r_t.append(spectrum.getRT())
            p_ms = []
            p_intensity = []
            for peak in spectrum:
                if peak.getIntensity() != 0:
                    p_ms.append(peak.getMZ())
                    p_intensity.append(peak.getIntensity())
            ms_index = np.argsort(np.negative(p_intensity))
            m_s.append(np.array(p_ms)[ms_index])
            intensity.append(np.array(p_intensity)[ms_index])
    rt1 = np.array(r_t)
    if rt1.shape[0] > 1:
        rt_mean_interval = np.mean(np.diff(rt1))
    else:
        rt_mean_interval = 0.0
    return m_s, intensity, r_t, rt_mean_interval
