import os
import numpy as np
import pyopenms

def readms(input_file):#only 'mzml,mzdata or mzxml' format
    #ms_format = re.search('\.\w+',input_file)
    ms_format= os.path.splitext(input_file)[1]
    #ms_format = ms_format.group()
    ms_format = ms_format.lower()
    msdata = pyopenms.MSExperiment()
    if ms_format == '.mzxml':
        file =  pyopenms.MzXMLFile()
    elif ms_format == '.mzml':
        file =  pyopenms.MzMLFile()
    elif ms_format =='.mzdata':
        file =  pyopenms.MzDataFile()
    else:
        raise Exception('ERROR: %s is wrong format' % input_file)
    file.load(r'%s' % input_file,msdata)
    ms = []
    intensity = []
    rt = []
    for spectrum in msdata:
        if spectrum.getMSLevel() == 1:
            rt.append(spectrum.getRT())
            p_ms = []
            p_intensity = []
            for peak in spectrum:
                if peak.getIntensity()!= 0:
                    p_ms.append(peak.getMZ())
                    p_intensity.append(peak.getIntensity())
            #print len(p_intensity)
            ms_index = np.argsort(-np.array(p_intensity))
            ms.append(np.array(p_ms)[ms_index])
            intensity.append(np.array(p_intensity)[ms_index])
            #scan+=1        
    rt1 = np.array(rt)
    if rt1.shape[0] > 1:
        rt_mean_interval = np.mean(np.diff(rt1))
    else:
        rt_mean_interval = 0.0
    #print  rt_mean_interval
    #rt_mean_interval = np.mean(rt1[1:]-rt1[:-1])
    #return ms,intensity,rt,scan,rt_max_interval
    return ms,intensity,rt,rt_mean_interval
#返回质谱，强度，保留时间，扫描次数，以及保留时间的最大间隔，和平均的时间间隔