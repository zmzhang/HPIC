# -*- coding: utf-8 -*-
"""
Created on Wed Jun 07 08:40:23 2017

@author: Wangrong
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 11 08:18:18 2017

@author: Wangrong
"""

import time,os
import numpy as np
#import pyopenms
import hdbscan
from sklearn.preprocessing import minmax_scale
from retention_time import to_deque
import matplotlib.pyplot as plt
import seaborn as sns

def readms(input_file):#only 'mzml,mzdata or mzxml' format
    #ms_format = re.search('\.\w+',input_file)
    import pyopenms
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
    rt_mean_interval = np.mean(np.diff(rt1))
    #print  rt_mean_interval
    #rt_mean_interval = np.mean(rt1[1:]-rt1[:-1])
    #return ms,intensity,rt,scan,rt_max_interval
    return ms,intensity,rt,rt_mean_interval
#返回质谱，强度，保留时间，扫描次数，以及保留时间的最大间隔，和平均的时间间隔
def maxI(intensity):
    #print intensity[:2]
    max_i = np.array([intensity_i[0] if len(intensity_i)>0 else 0 for intensity_i in intensity])
    
    max_intensity_rt = np.argmax(max_i)
    max_intensity_intensity = max_i[max_intensity_rt]
    return max_intensity_rt,max_intensity_intensity
#返回强度最大值得保留时间，以及质谱索引，和强度值
#输入质谱，强度，保留时间，质谱宽度，保留时间宽度，以及最高强度的保留时间，和质谱索引，和扫描次数
def choosedata(ms,intensity,rt,inv_ms,inv_rt,max_intensity_rt,scan):
    max_int_ms =ms[max_intensity_rt][0]
    #max_int_ms_intensity = max(max_i)
    #mz_max = max_int_ms+inv_ms
    #mz_min = max_int_ms-inv_ms
    #choose_ms = []
    #choose_intensity = []
    start = max_intensity_rt-inv_rt
    end = max_intensity_rt+inv_rt
    if start <0:
        start = 0
    if end>scan:
        end = scan
    choose_spec = []
    for index in range(start,end):
        ms_index = np.where(np.abs(ms[index]-max_int_ms)<inv_ms)[0]
        if ms_index.shape[0]>0:
            rt_1 = np.full(ms_index.shape[0],rt[index], dtype=np.double)
            choose_spec.extend(zip(rt_1,ms[index][ms_index],intensity[index][ms_index]))
        #choose_ms.append(ms[index][ms_index])
        #choose_intensity.append(intensity[index][ms_index])
    #print h_intensity_rt
    choose_spec = np.array(choose_spec)
    choose_rt = rt[start:end]
    h_intensity_rt = rt[max_intensity_rt]
    return choose_spec,choose_rt,max_int_ms,h_intensity_rt
def hdbscan_lc(choose_spec,choose_rt,h_intensity_rt,max_int_ms,rt_inv,mis_gap):
    #import hdbscan
    #from sklearn.preprocessing import minmax_scale
    #print 'what wrong'
    choose_spec_trans = np.zeros_like(choose_spec[:,:2])
    choose_spec_trans[:,1] = np.abs(choose_spec[:,1]-max_int_ms)
    choose_spec_trans[:,0] = choose_spec[:,0]-h_intensity_rt
    data = minmax_scale(choose_spec_trans)
    data[:,1] = data[:,1]*50
    clusterer = hdbscan.HDBSCAN(min_cluster_size=5)
    clusterer.fit(data)
    labels = clusterer.labels_
    label_index = np.argmax(choose_spec[:,2])
    label = labels[label_index]
    if label == -1:
        choose_spec_1 = choose_spec[label_index].reshape(1,3)        
    else:
        choose_spec_1 = choose_spec[labels==label,:]
        #plt.scatter(choose_spec_1[:,0],choose_spec_1[:,1],marker = 'o', color = 'r')
        #plt.scatter(choose_spec[labels!=label,:][:,0],choose_spec[labels!=label,:][:,1])
        #plt.savefig(r'D:/new direction/look/%s_scatter.png' % max_int_ms)
        #plt.close()
        choose_spec_1_rt = choose_spec_1[:,0]
        choose_rt_1 =choose_rt[choose_rt.index(choose_spec_1_rt[0]):choose_rt.index(choose_spec_1_rt[-1])+1]
        no_common = np.setdiff1d(np.array(choose_rt_1),choose_spec_1_rt)
        #print no_common.shape[0]
        if no_common.shape[0]>mis_gap:
            left = no_common[no_common<h_intensity_rt]
            if left.shape[0]>mis_gap:
                pre_rt = np.searchsorted(choose_spec_1_rt,left[-mis_gap-1])
            else:
                pre_rt = 0
            right = no_common[no_common>h_intensity_rt]
            if right.shape[0]>mis_gap:
                dif_rt = np.searchsorted(choose_spec_1_rt,right[mis_gap])
            else:
                dif_rt = len(choose_spec_1_rt)
            choose_spec_1 = choose_spec_1[pre_rt:dif_rt,:]
            choose_spec_1_rt = choose_spec_1[:,0]
        choose_1_rt = np.diff(choose_spec_1_rt)
        #delete the data while there are more than one point at the rentention time
        if not np.all(choose_1_rt):
            choose_1_rt_index = range(choose_spec_1_rt.shape[0])
            choose_1_rt = np.where(choose_1_rt ==0)[0]
            choose_spec_1_rt_1 = choose_spec_1_rt[choose_1_rt]
            k= 0
            while k<choose_spec_1_rt_1.shape[0]:
                p_multiple = choose_spec_1[choose_spec_1_rt==choose_spec_1_rt_1[k]]
                p_index = range(choose_1_rt[k],choose_1_rt[k]+ p_multiple.shape[0])
                p_index.remove(np.argmin(np.abs(p_multiple[:,1]-max_int_ms))+choose_1_rt[k])
                for index in p_index:
                    choose_1_rt_index.remove(index)
                    
                k +=(p_multiple.shape[0]-1)
            choose_spec_1= choose_spec_1[choose_1_rt_index,:]
            #plt.scatter(choose_spec_1[:,0],choose_spec_1[:,2],marker = 'o', color = 'r')
            #plt.savefig(r'D:/new direction/look/%s_peak.png' % max_int_ms)
            #plt.close()	
		
    return choose_spec_1
def PIC(inputfile,file_t,min_intensity=200,mass_inv=1,rt_inv=15,mis_gap=3):
    #start = time.time() 
    ms,intensity,rt,rt_mean_interval=readms(inputfile)
    rt_inv = int(rt_inv/rt_mean_interval)
    scan = len(rt)
    #print scan
    while True:
        max_intensity_rt,max_intensity_intensity = maxI(intensity)
        if max_intensity_intensity<min_intensity:
            break
        else:
            choose_spec,choose_rt,max_int_ms,h_intensity_rt = choosedata(ms,intensity,rt,mass_inv,rt_inv,max_intensity_rt,scan)
            if choose_spec.shape[0]<5:
                choose_spec_2= choose_spec[np.argmax(choose_spec[:,2])].reshape(1,3)
            else:
                choose_spec_2 = hdbscan_lc(choose_spec,choose_rt,h_intensity_rt,max_int_ms,rt_inv,mis_gap)
                #np.savetxt('%s/%s_%s_%s_1.txt' % (file_t,max_int_ms,h_intensity_rt,max_intensity_intensity),choose_spec)
            np.savetxt('%s/%s_%s_%s_%s_%s_%s.txt' % (file_t,max_int_ms,h_intensity_rt,max_intensity_intensity,choose_spec_2[0,0],choose_spec_2[-1,0],choose_spec_2.shape[0]),choose_spec_2)                
            del_1_index = rt.index(choose_spec_2[0,0])
            del_2_index = rt.index(choose_spec_2[-1,0])+1
            del_not = np.setdiff1d(np.array(rt[del_1_index:del_2_index]),choose_spec_2[:,0])
            del_all_ms = choose_spec_2[:,1]
            del_not_count = 0
            for index,del_ms_index in enumerate(range(del_1_index,del_2_index)):
                if rt[del_ms_index] in del_not:
                    del_not_count+=1
                else:
                    index = index-del_not_count
                    intensity[del_ms_index] = intensity[del_ms_index][ms[del_ms_index]!=del_all_ms[index]]
                    ms[del_ms_index] = ms[del_ms_index][ms[del_ms_index]!=del_all_ms[index]]
    #end = time.time()
    #print (end-start)/60
    return rt_mean_interval
#file_ms = 'D:/1/plot_0.2'
#file_dir = 'MM48_MSS.mzML'
#file_result = '1/test'
#file_all = os.listdir('%s' % file_dir)
#print file_all
'''for i in file_all[6:8]:
    #print i
    PIC('%s/%s' % (file_dir,i),os.path.splitext(os.path.basename(i))[0],500)'''
#EIC(r'D:\data\pos\Mspos_MM48_20uM_1-B,1_01_14615.mzdata',1,50,2,'D:/1/10_50_5_2')
#PIC('D:/HDBSCAN_code/MM48_MSS.mzML','D:/HDBSCAN_code/2',200)
def hpic(file_in,file_out,min_intensity=250,min_snr=3,mass_inv=1,rt_inv=15):
    file_t = os.path.splitext(os.path.basename(file_in))[0]
    file_t = '%s/%s' % (file_out,file_t)
    try:
        os.makedirs('%s' % file_t )
    except:
        pass
    interval = PIC(file_in,file_t,min_intensity)*1.5
    return to_deque(file_t,file_out,min_snr,interval,min_intensity)
#hpic(r'D:/data/pos/Mspos_MM48_20uM_1-B,1_01_14616.mzdata','D:/5/Mspos_MM48_20uM_1-B,1_01_14616',250,3)
    
      

    
     

        
            

    
    
    
    




