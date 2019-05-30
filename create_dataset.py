#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#function to shift data by 6hrs and zero out any non zero value between 8pm till 4am
#and interpolate any negative values
"""
Created on Mon May  6 17:03:24 2019

@author: tarun
"""





import numpy as np
import h5py
from pathlib import Path

#file paths
base_file_path = '/Users/tarun/Desktop/india_tmy.h5'
new_file_path = 'data.h5'
log_file = 'log_h5py.txt'

#opening base file
h5f = h5py.File(base_file_path,'r')

#variable
total_column = 100  #np.shape(h5f['dhi'])[1]  #total column
column_chunk_size=100

#constants
keys={0:'dni',1:'dhi',2:'ghi',3:'surface_temperature',4:'wspd'}
day=365

def rearrangeArray(data_array,days,cols):
    no_of_partion = days*4   #a day is partitioned in 4 parts
    temp=np.array_split(data_array,no_of_partion,axis=0)
    mat_1_4=temp[0:no_of_partion:4] # first 6 hours of each day
    mat_2_4=temp[1:no_of_partion:4] # 7am till 12 noon of each day
    mat_3_4=temp[2:no_of_partion:4] # 12noon till 6pm of each day
    mat_4_4=temp[3:no_of_partion:4] # 7pm 12 midnight of each day
    return np.array(list(zip(mat_4_4,mat_1_4,mat_2_4,mat_3_4))).reshape(days*24,cols)


def processArray(h5f,hf,col,key):
    print('key %s\n'%(key),file=log)
    for i in range(col,total_column,column_chunk_size):
        chunk_size = column_chunk_size
        if(total_column-i < column_chunk_size):
            chunk_size = total_column-i
        hf[key].resize(i+chunk_size,axis=1)
        old_data=h5f[key][:,i:i+chunk_size]
        if (key != 'wspd' and key != 'surface_temperature'):
            old_data[id_8pm_4am==1] = 0
        new_data=rearrangeArray(old_data,day,chunk_size)
        if (key != 'surface_temperature'):
            count_negValue = (new_data<0).sum()
            if (count_negValue > 0) :
                idx = np.where(new_data<0)
                interp_value(new_data,idx)
        hf[key][:,-chunk_size:] = new_data.reshape(24*day,chunk_size)
        print('col %d,chunk size %d\n'%(i,chunk_size),file=log)




def interp_value(arr,idx):
    row_index=idx[0]  #row index of neg values
    col_index=idx[1]  #column index of neg values
    n=len(row_index)  #count of neg values discovered
    x=list()          #empty list which will be updated with index of column of non negative values
    for i in range(n):
        for j in range(-5,5): # using 10 values to estimate the 11th value
            if (row_index[i]+j > np.shape(arr)[0]-1) or (row_index[i]+j < 0):
                continue #checking if the index is not the initial or last row
            if arr[row_index[i]+j,col_index[i]] >= 0:
                x.append(row_index[i]+j) #neglecting neg values and using only non negative value loation for interpolation
        arr[row_index[i],col_index[i]] = np.interp(row_index[i],x,arr[x,col_index[i]])
        x=list()





#creating new data file if it doenst exist with 1 empty column and 8760 rows
file_h5 = Path(new_file_path)
if not file_h5.is_file():
    print('creating new file')
    hf = h5py.File('data.h5', 'w')
    hf.create_dataset('dni',(24*day,1), maxshape=(24*day, total_column))
    hf.create_dataset('dhi', (24*day,1),maxshape=(24*day, total_column))
    hf.create_dataset('ghi', (24*day,1), maxshape=(24*day, total_column))
    hf.create_dataset('surface_temperature', (24*day,1), maxshape=(24*day, total_column))
    hf.create_dataset('wspd', (24*day,1), maxshape=(24*day, total_column))
    hf.close()


#index calculation from 8pm till 4am for old data
# 14 corresponds to 8pm, 23 corresponds to 4am in old data
id_8pm_4am = (np.arange(day*24)%24)//14
index_5am = np.arange(23,day*24,24)
id_8pm_4am[index_5am] = 0



log = open(log_file,'w')
with h5py.File(new_file_path, 'a') as hf:
    for i in keys:
        key=keys[i]
        dataset=hf[key]
        col=np.shape(dataset)[1]
        print(col)
        if col<total_column:
            print('appending')
            processArray(h5f,hf,col-1,key)
    hf.create_dataset('meta', data=h5f['meta'])

log.close()
h5f.close()
