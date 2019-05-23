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


day=365
id_8pm_4am = (np.arange(day*24)%24)//14
index_5am = np.arange(23,day*24,24)
id_8pm_4am[index_5am] = 0



y = np.arange(37.95,5.05,-0.1)
x = np.arange(67.05, 97.95, 0.1)

X, Y = np.meshgrid(x, y)

def processArray(h5f,mean_data,key,x,y):
    x=x-1
    y=y-1
    column_chunk_size = 330
    col = 0
    j=0
    x=x*24
    y=y*24
    for i in range(col,total_column,column_chunk_size):
        data=h5f[key][x+6:y+6:24,i:i+column_chunk_size]
        avg = np.flip(data.mean(axis=0))
        mean_data.T[j]=avg
        j=j+1



start_day = 60
end_day = 70
fig, ax = plt.subplots()
earth = Basemap(ax=ax,llcrnrlat=5.05,urcrnrlat=37.95,llcrnrlon=67.05,urcrnrlon=97.95)
earth.drawcoastlines( linewidth=1)
earth.drawcountries( linewidth=1)
ax.set_title("High DHI Area")
mean_dhi = np.zeros((330,310))

for i in range(5):
    start_day = start_day + i
    end_day = start_day + 1
    key = 'dhi'
    if i>1:
        key = 'dni'
    mean_dhi = np.random.randint(0,10,(330,310))
    #processArray(h5f,mean_dhi,key,start_day,end_day)
    out = mean_dhi

    fig, ax = plt.subplots()
    earth = Basemap(ax=ax,llcrnrlat=5.05,urcrnrlat=37.95,llcrnrlon=67.05,urcrnrlon=97.95)
    earth.drawcoastlines( linewidth=1)
    earth.drawcountries( linewidth=1)
    ax.set_title("High DHI Area")
    im=ax.contourf(X, Y, out)
    #fig.colorbar(im,ax=ax)
    ax.plot()
    plt.show()
    plt.pause(0.001)


'''
mean_dhi = np.zeros((330,310))
key = 'dhi'
processArray(h5f,mean_dhi,key,start_day,end_day)
mean_dni = np.zeros((330,310))
key = 'dni'
processArray(h5f,mean_dni,key,start_day,end_day)
mean_temp = np.zeros((330,310))
key = 'surface_temperature'#'dni'
processArray(h5f,mean_temp,key,start_day,end_day)
out = mean_dhi
'''
#out = ((mean_dni<600)*(mean_temp>30)).astype(int)
#out =  (mean_dni/mean_dhi)*((mean_temp>30)).astype(int)#(mean_dni/mean_dhi)
#print(mean_GHI)



#ax.plot_wireframe(X,Y,mean_GHI, rstride=2, cstride=2)
#ax.plot_wireframe(Y,X,mean_GHI, rstride=2, cstride=2)
#colors = matplotlib.cm.jet(np.hypot(x,y))
#ax.plot_surface(X,Y,mean_GHI, rstride=2, cstride=2)
'''
fig, ax = plt.subplots()
earth = Basemap(ax=ax,llcrnrlat=5.05,urcrnrlat=37.95,llcrnrlon=67.05,urcrnrlon=97.95)
earth.drawcoastlines( linewidth=1)
earth.drawcountries( linewidth=1)
im=ax.contourf(X, Y, out)
fig.colorbar(im,ax=ax)

ax.set_title("High DHI Area")
ax.plot()
plt.show()
'''
