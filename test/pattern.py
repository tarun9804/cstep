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
from mpl_toolkits.basemap import Basemap


base_file_path = '/Users/tarun/Desktop/india_tmy.h5'
h5f = h5py.File(base_file_path,'r')

#variable
total_column = np.shape(h5f['ghi'])[1]  #total column
column_chunk_size=330
key='tmy_year'
yr=h5f[key][2,:]
idx=np.where((yr==2014)|(yr==2011)|(yr==2009)|(yr==2008)|(yr==2006)|\
(yr==2003)|(yr==2000))
mask=np.full(total_column,0)
mask[idx]=1

day=365
id_8pm_4am = (np.arange(day*24)%24)//14
index_5am = np.arange(23,day*24,24)
id_8pm_4am[index_5am] = 0


start_day = 60
end_day = 61


def processArray(h5f,mean_data,key,x,y,var):
    x=x-1
    y=y-1
    column_chunk_size = 330
    col = 0
    j=0
    x=x*24
    y=y*24
    for i in range(col,total_column,column_chunk_size):
        data=h5f[key][x:y,i:i+column_chunk_size]
        ms=mask[i:i+column_chunk_size]
        data=data*ms
        avg = np.flip(data.mean(axis=0))
        var_t = np.flip(data.var(axis=0))
        var.T[j]=var_t
        mean_data.T[j]=avg
        j=j+1

mean_dhi = np.zeros((330,310))
var_dhi = np.zeros((330,310))
key = 'dhi'
processArray(h5f,mean_dhi,key,start_day,end_day,var_dhi)

mean_dni = np.zeros((330,310))
var_dni = np.zeros((330,310))
key = 'dni'
processArray(h5f,mean_dni,key,start_day,end_day,var_dni)

mean_temp = np.zeros((330,310))
var_temp = np.zeros((330,310))
key = 'surface_temperature'#'dni'
processArray(h5f,mean_temp,key,start_day,end_day,var_temp)

mean_dew = np.zeros((330,310))
var_dew = np.zeros((330,310))
key = 'dew_point'#'dni'
processArray(h5f,mean_dew,key,start_day,end_day,var_dew)

#out = ((var_temp>30)*(mean_dhi>100)).astype(int)
out = ((mean_dhi-mean_dni))
#out = var_temp
#out = (((var_dni/1000)<100)*(mean_dhi>100)).astype(int)
#out = (mean_dni/mean_dhi)
#out = (mean_dhi>100).astype(int)
#out = ((mean_dhi>100)*(mean_temp>20)).astype(int)
#out = ((mean_dni<600)*(mean_temp>30)).astype(int)
#out =  (mean_dni/mean_dhi)*((mean_temp>30)).astype(int)#(mean_dni/mean_dhi)
#print(mean_GHI)
y = np.arange(37.95,5.05,-0.1)
x = np.arange(67.05, 97.95, 0.1)

X, Y = np.meshgrid(x, y)


#ax.plot_wireframe(X,Y,mean_GHI, rstride=2, cstride=2)
#ax.plot_wireframe(Y,X,mean_GHI, rstride=2, cstride=2)
#colors = matplotlib.cm.jet(np.hypot(x,y))
#ax.plot_surface(X,Y,mean_GHI, rstride=2, cstride=2)
fig, ax = plt.subplots()
earth = Basemap(ax=ax,llcrnrlat=5.05,urcrnrlat=37.95,llcrnrlon=67.05,urcrnrlon=97.95)
earth.drawcoastlines( linewidth=1)
earth.drawcountries( linewidth=1)
im=ax.contourf(X, Y, out)
fig.colorbar(im,ax=ax)

ax.set_title("High DHI Area")
ax.plot()
f_name = 'file_'+str(start_day)+'.png'
#plt.savefig(f_name)
plt.show()





#gdf.plot(ax=ax, color='red')

#plt.show()

#ax.scatter3D(X, Y, mean_GHI, c=mean_GHI, cmap='Greens');

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
