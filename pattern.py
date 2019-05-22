#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 16:21:48 2019

@author: tarun
"""

import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d


base_file_path = '/Users/tarun/Desktop/india_tmy.h5'
h5f = h5py.File(base_file_path,'r')

#variable
total_column = np.shape(h5f['ghi'])[1]  #total column
column_chunk_size=330


day=365
id_8pm_4am = (np.arange(day*24)%24)//14
index_5am = np.arange(23,day*24,24)
id_8pm_4am[index_5am] = 0

mean_data = np.zeros((330,310))
col = 0
key = 'dhi'
start_day = 50
end_day = 60


def processArray(h5f,mean_data,key,x-1,y-1):
    column_chunk_size = 330
    col = 0
    j=0
    x=x*24
    y=y*24
    for i in range(col,total_column,column_chunk_size):
        data=h5f[key][x:y,i:i+column_chunk_size]
        avg = np.flip(data.mean(axis=0))
        mean_data.T[j]=avg
        j=j+1


processArray(h5f,mean_data,key,start_day,end_day)
#print(mean_GHI)
y = np.arange(37.95,5.05,-0.1)
x = np.arange(67.05, 97.95, 0.1)

X, Y = np.meshgrid(x, y)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.plot_wireframe(X,Y,mean_GHI, rstride=2, cstride=2)
#ax.plot_wireframe(Y,X,mean_GHI, rstride=2, cstride=2)
#colors = matplotlib.cm.jet(np.hypot(x,y))
#ax.plot_surface(X,Y,mean_GHI, rstride=2, cstride=2)
ax.contourf(X, Y, mean_data)

#ax.scatter3D(X, Y, mean_GHI, c=mean_GHI, cmap='Greens');
plt.show()
#def rearrangeArray(data_array,days,cols):
#    no_of_partion = days*4   #a day is partitioned in 4 parts
#    temp=np.array_split(data_array,no_of_partion,axis=0)
#    mat_1_4=temp[0:no_of_partion:4] # first 6 hours of each day
#    mat_2_4=temp[1:no_of_partion:4] # 7am till 12 noon of each day
#    mat_3_4=temp[2:no_of_partion:4] # 12noon till 6pm of each day
#    mat_4_4=temp[3:no_of_partion:4] # 7pm 12 midnight of each day
#    return np.array(list(zip(mat_4_4,mat_1_4,mat_2_4,mat_3_4))).reshape(days*24,cols)
#
#
#def processArray(h5f,col,key):
#    print('key %s\n'%(key),file=log)
#    for i in range(col,total_column,column_chunk_size):
#        chunk_size = column_chunk_size
#        if(total_column-i < column_chunk_size):
#            chunk_size = total_column-i
#        hf[key].resize(i+chunk_size,axis=1)
#        old_data=h5f[key][:,i:i+chunk_size]
#        if (key != 'wspd' and key != 'surface_temperature'):
#            old_data[id_8pm_4am==1] = 0
#        new_data=rearrangeArray(old_data,day,chunk_size)
#
#
