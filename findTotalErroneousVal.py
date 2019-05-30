#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#functions to total non-zero values in dni,dhi,ghi between 8pm till 4am
#find negative values if any in dni,dhi,ghi,wspd and interpolate it
"""
Created on Fri May  3 11:35:53 2019

@author: tarun
"""
import numpy as np
import pandas as pd
import h5py




def rearrangeArray(data_array,days,cols):
    no_of_partion = days*4   #a day is partitioned in 4 parts
    temp=np.array_split(data_array,no_of_partion,axis=0)
    mat_1_4=temp[0:no_of_partion:4] # first 6 hours of each day
    mat_2_4=temp[1:no_of_partion:4] # 7am till 12 noon of each day
    mat_3_4=temp[2:no_of_partion:4] # 12noon till 6pm of each day
    mat_4_4=temp[3:no_of_partion:4] # 7pm 12 midnight of each day
    return np.array(list(zip(mat_4_4,mat_1_4,mat_2_4,mat_3_4))).reshape(days*24,cols)
 


def findIndexofErrorValue(h5f,key,id_8pm_4am,array0,array1):
    non_zero_val = 0
    neg_val = 0
    for i in range(starting_column,last_column,column_chunk_size):
        chunk_size = column_chunk_size
        if(last_column-i<column_chunk_size):
            chunk_size = last_column-i 
        data = h5f[key][:,i:i+chunk_size]
        if key != 'wspd':
            val = np.sum(data[np.where(id_8pm_4am==1)]!=0)
            if val > 0:
                non_zero_val = non_zero_val + val
                tu = np.where(data[np.where(id_8pm_4am==1)]!=0)
                col = tu[1] + i
                array0[col] = 1
            
        val = np.sum(data < 0)
        if val > 0:
            neg_val = neg_val + val
            tu = np.where(data < 0)
            col = tu[1] + i
            array1[col] = 1
    return [non_zero_val,neg_val,array0,array1]



def findNonZeroIndexWrapper():
    keys={0:'dni',1:'dhi',2:'ghi',3:'wspd'}
    for i in keys:
        key=keys[i]
        z_key = key+'_NonZeroErr'
        Err_LatLong[z_key] = np.zeros((np.shape(h5f[key])[1]))
        n_key = key+'_NegValueErr'
        Err_LatLong[n_key] = np.zeros((np.shape(h5f[key])[1]))
        temp0 = np.zeros((np.shape(h5f[key])[1]))
        temp1 = np.zeros((np.shape(h5f[key])[1]))
        result = findIndexofErrorValue(h5f,key,id_8pm_4am,temp0,temp1)
        #total_nonzero = findIndexofNegValue(h5f,key,id_8pm_4am,temp1)
        Err_LatLong[z_key]=result[2]
        Err_LatLong[n_key]=result[3]
        
        if key != 'wspd':
            print('total non zero values %d, in %s' %(result[0], key))
        print('total negative values %d, in %s' %(result[1], key))
        if key == 'wspd':
            Err_LatLong.drop(z_key,axis=1,inplace=True)
        print('\n')

h5f = h5py.File('/Users/tarun/Desktop/india_tmy.h5','r')
day=365
starting_column =0
last_column=np.shape(h5f['dni'])[1]
column_chunk_size=200

id_8pm_4am = (np.arange(day*24)%24)//14  #index_8pm = np.arange(14,day*24,24)
index_5am = np.arange(23,day*24,24)
id_8pm_4am[index_5am] = 0

Err_LatLong = pd.DataFrame()

#t1=time.time()

findNonZeroIndexWrapper()
h5f.close()  
