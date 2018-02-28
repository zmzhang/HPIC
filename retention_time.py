# -*- coding: utf-8 -*-
"""
Created on Sat Jun 03 23:50:05 2017

@author: Wangrong
"""
import numpy as np
import pandas as pd
import os
import shutil
import time
from collections import deque
from lc_ms_peak import lc_ms_peak 
def to_deque(file_in,file_out,min_snr,rt_v,intensity):
    number = 10
    start = time.time()
    width = np.arange(1,60)
    os.chdir(file_in)
    files = os.listdir(file_in)
    peak_list = []
    alls = []
    for each in files:
        each = map(float,each[:-4].split('_'))
        alls.append(each)
    frame = np.array(alls)
    #print frame.shape
    files = np.array(files)[np.argsort(frame[:,0])]
    frame = frame[np.argsort(frame[:,0])]
    frame_mz = frame[:,0]
    l_frame = frame.shape[0]
    i = 0
    count = []
    while i<l_frame:
        if i not in count:
            mz = frame_mz[i]
            mz_range = mz +0.1
            p_frame =  frame[i:,:][frame_mz[i:]<mz_range]
            length = p_frame.shape[0]
            if length >1:
                pre = p_frame[0,3]
                dif = p_frame[0,4]
                condition = np.array(list(p_frame[:,3]-dif) + list(pre-p_frame[:,4]))
                if condition[(condition>0)&(condition<rt_v)].shape[0] >0 :
                    file_sep = deque()
                    file_sep.append(0)
                    a = np.where(((pre-p_frame[:,4])>0)&((pre-p_frame[:,4])<rt_v))[0]
                    while a.shape[0]>0 and p_frame[p_frame[:,3]==pre].shape[0]<2:
                        file_sep.appendleft(a[0])
                        pre = p_frame[a[0],3]
                        a = np.where(((pre-p_frame[:,4])>0)&((pre-p_frame[:,4])<rt_v))[0]
                    b = np.where(((p_frame[:,3]-dif)>0)&((p_frame[:,3]-dif)<rt_v))[0]
                    while b.shape[0]>0 and p_frame[p_frame[:,4]==dif].shape[0]<2 :
                        file_sep.append(b[0])
                        dif = p_frame[b[0],4]
                        b = np.where(((p_frame[:,3]-dif)>0)&((p_frame[:,3]-dif)<rt_v))[0]
                    new = np.zeros(6)
                    new[:3] = p_frame[file_sep[np.argmax(p_frame[file_sep,2])],:3]
                    index = [m+i for m in file_sep]
                    mass = np.vstack([np.loadtxt(file_i) for file_i in list(files[index])])
                    new[3] = mass[0,0]
                    new[4] = mass[-1,0]
                    new[5] = mass.shape[0]
                    if new[5] >=number:
                        #np.savetxt('%s/%s.txt' % (new_dir,'_'.join(map(str,list(new)))),mass)
                        #peak_list.extend(lc_ms_peak(new,np.arange(1,100),3,12,True,mass))
                        peak_list.extend(lc_ms_peak(mass,width,min_snr,True,intensity))
                        #peak_list.extend(lc_ms_peak(mass,width,min_snr,True))
                        #print peak_list
                        
                    count.extend(index)
                else:
                    #print frame[i]
                    #print files[i]
                    if frame[i,5]>=number:
                        
                        #peak_list.extend(lc_ms_peak(frame[i],np.arange(1,100),3,12,False,old_dir+'/'+files[i]))
                        peak_list.extend(lc_ms_peak(files[i],width,min_snr,False,intensity))
                        #peak_list.extend(lc_ms_peak(old_dir+'/'+files[i],width,min_snr,False))
                        #shutil.copyfile('%s' % old_dir+'/'+files[i],'%s' % new_dir+'/'+files[i])   
            else:
                #print frame[i]
                #print files[i]
                if frame[i,5]>=number:
                    
                    peak_list.extend(lc_ms_peak(files[i],width,min_snr,False,intensity))
                    #peak_list.extend(lc_ms_peak(old_dir+'/'+files[i],width,min_snr,False))
                    #peak_list.extend(lc_ms_peak(frame[i],np.arange(1,100),3,12,False,old_dir+'/'+files[i]))
                    #Yshutil.copyfile('%s' % old_dir+'/'+files[i],'%s' % new_dir+'/'+files[i])
                
        i+=1
    #print len(peak_list)
    #print peak_list[0]
    #shutil.rmtree(file_in)
    peak_list = pd.DataFrame(peak_list,columns = ['mz','rt','intensity','rt1','rt2','signal','snr','ind','shape'])
    os.chdir(file_out)
    peak_list.to_csv(r'%s/%s.csv' % (file_out,os.path.splitext(os.path.basename(file_in))[0]))
    #os.chdir(file_in)
    shutil.rmtree(file_in)
    #print time.time()-start
    return peak_list
'''file_out = 'D:/HDBSCAN_code/test'
files = 'D:/HDBSCAN_code/test'
file_all = os.listdir('%s' % files)
min_snr = 3
rt_v = 0.3
intensity = 200
intensity = intensity
print file_all'''
'''for i in file_all[1:]:
    #print i
    to_deque('%s/%s' % (files,i),'%s/%s' % (file_out,i),min_snr,rt_v,intensity)
    #to_deque('%s/%s' % (files,i),'%s/%s' % (file_out,i),min_snr,rt_v)'''
#to_deque('D:/HDBSCAN_code/2','D:/HDBSCAN_code',min_snr,rt_v,intensity)
                
                


