#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:03:19 2019

@author: tarun
"""




import pandas as pd
import numpy as np
import time
import user_input as ui

row_count = 8760

def con2time(x):
    temp = np.modf(x)
    hr=temp[1]
    temp1=np.modf(temp[0]*60)
    mi = temp1[1]
    sec=np.modf(temp1[0]*60)
    return hr,mi,sec[1]

def getEOT(x):
    return (229.2 * (0.000075 + 0.001868 * np.cos(x) - 0.032077 * np.sin(x) - 0.014615 * np.cos(2*x) - 0.040849 * np.sin(2*x)))

def getAngDeclination(x):
    return ((57.29577951308232)*(0.006918 - 0.399912 * np.cos(x) + 0.070257 * np.sin(x) - 0.006758 * np.cos(2*x) + 0.000907 * np.sin(2*x) - 0.002697 * np.cos(3*x) + 0.00148 * np.sin(3*x)))

def monthorder(x):
    month = [31,59,90,120,151,181,212,243,273,304,334]
    temp=np.zeros((12,31))
    temp[0]=x[0:month[0]]
    temp[1,0:28]=x[month[0]:month[1]]
    temp[2]=x[month[1]:month[2]]
    temp[3,0:30]=x[month[2]:month[3]]
    temp[4]=x[month[3]:month[4]]
    temp[5,0:30]=x[month[4]:month[5]]
    temp[6]=x[month[5]:month[6]]
    temp[7]=x[month[6]:month[7]]
    temp[8,0:30]=x[month[7]:month[8]]
    temp[9]=x[month[8]:month[9]]
    temp[10,0:30]=x[month[9]:month[10]]
    temp[11]=x[month[10]:]
    return temp



def ZoneTime_SolarAngles(Location_Latitude, Location_Longitude, Module_Tilt, \
                          Module_Surface_Azimuth, Ref_Latitude, Ref_Longitude):

        Location_Latitude = Location_Latitude#x[0]#6.05
        Location_Longitude = Location_Longitude #x[1]#-67.05
        #print('processing',x[1])
        #Module_Tilt = 12.966
        #Module_Surface_Azimuth = 0
        row_count = 8760
        #Reference_Latitude = 25.15
        Reference_Longitude = Ref_Longitude #x[5]#-82.58

        EOLD = 4 * (Reference_Longitude - Location_Longitude)
        Location_Latitude_r = np.deg2rad(Location_Latitude)
        #Hour = pd.DataFrame(np.arange(row_count),columns=['hr'])
        Hour=np.arange(row_count)
        #Hour['ZT_Day_Hour']=Hour['Hour']%24
        ZT_Day_Hour = Hour%24
        #Hour['ZT_day']=(Hour['hr']//24)+1
        ZT_day = (Hour//24)+1
        temp=360/365
        #Hour['ZT_B']=np.deg2rad((Hour['ZT_day']-1)*temp)
        ZT_B = np.deg2rad((ZT_day-1)*temp)
        #Hour['Zone_Time']=pd.date_range(start='2015-1-1 00:00:00',freq='1H',periods=row_count)
        Zone_Time=pd.date_range(start='2015-1-1 00:00:00',freq='1H',periods=row_count)

        #Hour['EOT'] = Hour['ZT_B'].apply(getEOT)   #(229.2 * (0.000075 + 0.001868 * np.cos(Hour['ZT_B_rad']) - 0.032077 * np.sin(Hour['ZT_B_rad']) - 0.014615 * np.cos(2*Hour['ZT_B_rad']) - 0.040849 * np.sin(2*Hour['ZT_B_rad'])))
        #Hour['EOT'] = getEOT(Hour['ZT_B'])
        EOT = getEOT(ZT_B)
        #Hour['Correction'] = Hour['EOT'] + EOLD
        Correction = EOT + EOLD
        #Hour['Correction/60']=Hour['Correction']/60
        Correction_60 = Correction/60
        #day_12['corr_hr'],day_12['corr_min'],day_12['corr_sec']=con2time(day_12['Correction'])
        Correction_Hour,Correction_Min,Correction_Sec = con2time(Correction_60)

        #Hour['Solar_Time'] =  Hour['Zone_Time'] + pd.to_timedelta(Correction_60,unit='h')
        Solar_Time =  Zone_Time + pd.to_timedelta(Correction_60,unit='h')

        #Hour['ST_day_hr'] = pd.to_datetime(Hour['Solar_Time']).dt.hour + pd.to_datetime(Hour['Solar_Time']).dt.minute/60 + pd.to_datetime(Hour['Solar_Time']).dt.second/3600
        #ST_day_hr = pd.to_datetime(Hour['Solar_Time']).dt.hour + pd.to_datetime(Hour['Solar_Time']).dt.minute/60 + pd.to_datetime(Hour['Solar_Time']).dt.second/3600
        ST_day_hr = Solar_Time.hour + Solar_Time.minute/60 + Solar_Time.second/3600

        #Hour['ST_day'] = Hour.Solar_Time.dt.dayofyear
        ST_day = Solar_Time.dayofyear

        temp=360/365
        #Hour['ST_B'] = np.deg2rad((Hour.ST_day - 1)*temp)
        ST_B = np.deg2rad((ST_day - 1)*temp)
        #Hour['ST_Ang_Declination'] = getAngDeclination(Hour['ST_B'])
        ST_Ang_Declination = getAngDeclination(ST_B)
        #Hour['ST_Ang_Hour_f'] = (-15*(12-Hour['ST_day_hr']))
        ST_Ang_Hour_f = (-15*(12-ST_day_hr))
        #ST_Ang_Declination_r = np.deg2rad(Hour['ST_Ang_Declination'])
        ST_Ang_Declination_r = np.deg2rad(ST_Ang_Declination)

        #ST_Ang_Hour_f_r = np.deg2rad(Hour['ST_Ang_Hour_f'])
        ST_Ang_Hour_f_r = np.deg2rad(ST_Ang_Hour_f)
        #Hour['ST_Ang_Sol_Zenith'] = np.rad2deg(np.arccos(np.cos(Location_Latitude_r)*np.cos(ST_Ang_Hour_f_r)*np.cos(ST_Ang_Declination_r) + np.sin(Location_Latitude_r)*np.sin(ST_Ang_Declination_r)))
        ST_Ang_Sol_Zenith = np.rad2deg(np.arccos(np.cos(Location_Latitude_r)*np.cos(ST_Ang_Hour_f_r)*np.cos(ST_Ang_Declination_r) + np.sin(Location_Latitude_r)*np.sin(ST_Ang_Declination_r)))

        #Hour['ST_Ang_Sol_Altitude'] = 90 - Hour['ST_Ang_Sol_Zenith']
        ST_Ang_Sol_Altitude = 90 - ST_Ang_Sol_Zenith

        #Hour['Sun_Hour_Status'] = 0
        #Hour.loc[Hour['ST_Ang_Sol_Altitude'] >= 0 , 'Sun_Hour_Status']=1

        Sun_Hour_Status = (ST_Ang_Sol_Altitude>=0).astype(int)
        Annual_Sun_Hours = np.sum(ST_Ang_Sol_Altitude>=0)

        #Hour['ST_Ang_Sol_Azimuth'] = 0
        ST_Ang_Sol_Azimuth=np.zeros(row_count)
        #ST_Ang_Sol_Zenith_r = np.deg2rad(Hour['ST_Ang_Sol_Zenith'])
        ST_Ang_Sol_Zenith_r = np.deg2rad(ST_Ang_Sol_Zenith)

        #Hour['ST_Ag_Sol_Az_temp'] = (np.cos(ST_Ang_Sol_Zenith_r)*np.sin(Location_Latitude_r) - np.sin(ST_Ang_Declination_r))/(np.sin(ST_Ang_Sol_Zenith_r)*np.cos(Location_Latitude_r))
        ST_Ag_Sol_Az_temp = (np.cos(ST_Ang_Sol_Zenith_r)*np.sin(Location_Latitude_r) - np.sin(ST_Ang_Declination_r))/(np.sin(ST_Ang_Sol_Zenith_r)*np.cos(Location_Latitude_r))

        #Hour.loc[Hour['ST_Ag_Sol_Az_temp'] < -1,'ST_Ag_Sol_Az_temp'] = -1
        #Hour.loc[Hour['ST_Ag_Sol_Az_temp'] > 1,'ST_Ag_Sol_Az_temp'] = 1
        ST_Ag_Sol_Az_temp = np.array(ST_Ag_Sol_Az_temp)
        ST_Ag_Sol_Az_temp[ST_Ag_Sol_Az_temp<-1]=-1
        ST_Ag_Sol_Az_temp[ST_Ag_Sol_Az_temp>1]=1
        #b = np.where(a<3,0,1)
        #ST_Ag_Sol_Az_temp = np.where(ST_Ag_Sol_Az_temp<-1,-1)
        #ST_Ag_Sol_Az_temp = np.where(ST_Ag_Sol_Az_temp>1,1)
        #Hour['ST_Ang_Sol_Azimuth'] = np.sign(Hour['ST_Ang_Hour_f'])*np.abs(np.rad2deg(np.arccos(Hour['ST_Ag_Sol_Az_temp'])))
        ST_Ang_Sol_Azimuth = np.sign(ST_Ang_Hour_f)*np.abs(np.rad2deg(np.arccos(ST_Ag_Sol_Az_temp)))


        #Hour.loc[Hour['ST_Ang_Sol_Azimuth'] == np.inf,'ST_Ang_Sol_Azimuth'] = 0
        ST_Ang_Sol_Azimuth = np.array(ST_Ang_Sol_Azimuth)
        ST_Ang_Sol_Azimuth[ST_Ang_Sol_Azimuth==np.inf] = 0

        Module_Tilt_r =  np.deg2rad(Module_Tilt)
        Module_Surface_Azimuth_r =  np.deg2rad(Module_Surface_Azimuth)
        ST_Ang_Incidence =np.arccos(np.sin(ST_Ang_Declination_r)*np.sin(Location_Latitude_r)*np.cos(Module_Tilt_r) -\
                                    np.sin(ST_Ang_Declination_r)*np.cos(Location_Latitude_r)*np.sin(Module_Tilt_r)*np.cos(Module_Surface_Azimuth_r)+\
                                    np.cos(ST_Ang_Declination_r)*np.cos(Location_Latitude_r)*np.cos(Module_Tilt_r)*np.cos(ST_Ang_Hour_f_r)+\
                                    np.cos(ST_Ang_Declination_r)*np.sin(Location_Latitude_r)*np.sin(Module_Tilt_r)*np.cos(Module_Surface_Azimuth_r)*np.cos(ST_Ang_Hour_f_r)+\
                                    np.cos(ST_Ang_Declination_r)*np.sin(Module_Tilt_r)*np.sin(Module_Surface_Azimuth_r)*np.sin(ST_Ang_Hour_f_r))
        ST_Ang_Incidence = np.rad2deg(ST_Ang_Incidence)
        '''
        ST_Ang_Incidence.append(acosd(sind(ST_Ang_Declination[i])*sind(Location_Latitude)*cosd(Module_Tilt) -\
                                      sind(ST_Ang_Declination[i])*cosd(Location_Latitude)*sind(Module_Tilt)*cosd(Module_Surface_Azimuth)+\
                                      cosd(ST_Ang_Declination[i])*cosd(Location_Latitude)*cosd(Module_Tilt)*cosd(ST_Ang_Hour_f[i])+\
                                      cosd(ST_Ang_Declination[i])*sind(Location_Latitude)*sind(Module_Tilt)*cosd(Module_Surface_Azimuth)*cosd(ST_Ang_Hour_f[i])+\
                                      cosd(ST_Ang_Declination[i])*sind(Module_Tilt)*sind(Module_Surface_Azimuth)*sind(ST_Ang_Hour_f[i])))
        '''

        index_12 = np.arange(12,row_count,24)
        #day_12 = pd.DataFrame()
        Correction_60_12 = Correction_60[index_12]
        Correction_HHMMSS = (Correction_Hour[index_12],Correction_Min[index_12],Correction_Sec[index_12])

        ST_Wset = np.rad2deg(np.arccos(-np.tan(Location_Latitude_r)*np.tan(ST_Ang_Declination_r[index_12])))

        #day_12['ST_Wrise'] = -day_12['ST_Wset']
        ST_Wrise = -ST_Wset

        #day_12['T_Set'] = 12 + day_12['ST_Wset']/15
        T_Set = 12 + ST_Wset/15

        #day_12['Lis_Tset'] = day_12['T_Set'] - Hour['Correction'][index_12]/60
        #day_12['Lis_Tset'] = day_12['T_Set'] - day_12['Correction']
        Lis_Tset = T_Set - Correction_60_12

        ST_Tset = Solar_Time[index_12].normalize() + pd.to_timedelta(T_Set,unit='h')

        #day_12['T_Rise'] = 12 + day_12['ST_Wrise']/15
        T_Rise = 12 + ST_Wrise/15

        #day_12['Lis_Trise'] = day_12['T_Rise'] - day_12['Correction']
        Lis_Trise = T_Rise - Correction_60_12

        #day_12['ST_Trise'] = pd.to_datetime(Hour['Solar_Time'][index_12]).dt.normalize() + pd.to_timedelta(day_12['T_Rise'],unit='h')
        ST_Trise = Solar_Time[index_12].normalize() + pd.to_timedelta(T_Rise,unit='h')


        #day_12['ZT_Trise'] = day_12['ST_Trise'] - pd.to_timedelta(day_12['Correction'],unit='h')
        ZT_Trise = ST_Trise - pd.to_timedelta(Correction_60_12,unit='h')


        #day_12['ZT_Tset'] = day_12['ST_Tset'] - pd.to_timedelta(day_12['Correction'],unit='h')
        ZT_Tset = ST_Tset - pd.to_timedelta(Correction_60_12,unit='h')

        #day_12['ZT_Date'] = pd.to_datetime(Hour['Solar_Time'][index_12]).dt.normalize()
        ZT_Date = Solar_Time[index_12].normalize()

        #day_12['Day_length'] = day_12['ST_Wset']*2/15
        Day_length = ST_Wset*0.133333#2/15

        #day_12['ST_Day_Length'] = pd.to_datetime(Hour['Solar_Time'][index_12]).dt.normalize() + pd.to_timedelta(day_12['Day_length'],unit='h')
        ST_Day_Length = Solar_Time[index_12].normalize() + pd.to_timedelta(Day_length,unit='h')

        #day_12['Day_Length_Tset_Trise'] = day_12['Lis_Tset'] - day_12['Lis_Trise']
        Day_Length_Tset_Trise = Lis_Tset - Lis_Trise

        #day_12.reset_index(inplace=True,drop=True)
        #indexs=['MaxCorrection','MinCorrection','MaxTrise','MaxTset','MinTrise','MinTset','MaxDayLength','MinDayLength','Sun_Window_Length']

        #i=day_12['Correction'].idxmax(axis=0)
        i=np.argmax(Correction_60_12)
        #dt = day_12['ZT_Date'][i]
        #dt = ZT_Date[i]
        Max_Correction = ZT_Date[i] + pd.to_timedelta(Correction_60_12[i],unit='h')
        #result.loc['MaxCorrection','DateTime'] = dt + pd.to_timedelta(day_12['Correction'][i],unit='h') #+ pd.to_timedelta(day_12['corr_min'],unit='m') + pd.to_timedelta(day_12['corr_sec'],unit='s')

        #i=day_12['Correction'].idxmin(axis=0)
        #dt = day_12['ZT_Date'][i]
        #result.loc['MinCorrection','DateTime'] = dt + pd.to_timedelta(day_12['Correction'][i],unit='h') #+ pd.to_timedelta(day_12['corr_min'],unit='m') + pd.to_timedelta(day_12['corr_sec'],unit='s')
        i=np.argmin(Correction_60_12)
        Min_Correction = ZT_Date[i] + pd.to_timedelta(Correction_60_12[i],unit='h')

        #i=day_12['Lis_Trise'].idxmax(axis=0)
        #dt = day_12['ZT_Date'][i]
        #result.loc['MaxTrise','DateTime'] = dt + pd.to_timedelta(day_12['Lis_Trise'][i],unit='h')
        i=np.argmax(Lis_Trise)
        Max_Trise = ZT_Date[i] + pd.to_timedelta(Lis_Trise[i],unit='h')

        #i=day_12['Lis_Tset'].idxmax(axis=0)
        #dt = day_12['ZT_Date'][i]
        #result.loc['MaxTset','DateTime'] = dt + pd.to_timedelta(day_12['Lis_Tset'][i],unit='h')
        i=np.argmax(Lis_Tset)
        Max_Tset = ZT_Date[i] + pd.to_timedelta(Lis_Tset[i],unit='h')


        #i=day_12['Lis_Trise'].idxmin(axis=0)
        #dt = day_12['ZT_Date'][i]
        #result.loc['MinTrise','DateTime'] = dt + pd.to_timedelta(day_12['Lis_Trise'][i],unit='h')
        i=np.argmin(Lis_Trise)
        Min_Trise = ZT_Date[i] + pd.to_timedelta(Lis_Trise[i],unit='h')


        #i=day_12['Lis_Tset'].idxmin(axis=0)
        #dt = day_12['ZT_Date'][i]
        #result.loc['MinTset','DateTime'] = dt + pd.to_timedelta(day_12['Lis_Tset'][i],unit='h')
        i=np.argmin(Lis_Tset)
        Min_Tset = ZT_Date[i] + pd.to_timedelta(Lis_Tset[i],unit='h')


        #i=day_12['Day_Length_Tset_Trise'].idxmax(axis=0)
        #dt = day_12['ZT_Date'][i]
        #result.loc['MaxDayLength','DateTime'] = dt + pd.to_timedelta(day_12['Day_Length_Tset_Trise'][i],unit='h')
        i=np.argmax(Day_Length_Tset_Trise)
        Max_DayLength = ZT_Date[i] + pd.to_timedelta(Day_Length_Tset_Trise[i],unit='h')



        Max_DL_SRSS_Trise = ZT_Trise.time[i]
        Max_DL_SRSS_Tset = ZT_Tset.time[i]
        Max_DL_SRSS = (Max_DL_SRSS_Trise,Max_DL_SRSS_Tset)

        i=np.argmin(Day_Length_Tset_Trise)
        Min_DayLength = ZT_Date[i] + pd.to_timedelta(Day_Length_Tset_Trise[i],unit='h')


        Min_DL_SRSS_Trise = ZT_Trise.time[i]
        Min_DL_SRSS_Tset = ZT_Tset.time[i]
        Min_DL_SRSS = (Min_DL_SRSS_Trise,Min_DL_SRSS_Tset)


        Lis_Trise=np.array(Lis_Trise)
        Lis_Tset=np.array(Lis_Tset)
        Sun_Window = (np.floor(np.min(Lis_Trise)),np.ceil(np.max(Lis_Tset)))
        l=Sun_Window[1]-Sun_Window[0]
        Sun_Window_Length = pd.to_timedelta(l,unit='h')


        return (Zone_Time, ZT_Day_Hour,Solar_Time, ST_Ang_Hour_f, ST_Ang_Sol_Altitude, \
                ST_Ang_Sol_Azimuth, ST_Ang_Incidence, Sun_Hour_Status,\
                Annual_Sun_Hours,Correction_HHMMSS, ZT_Trise, ZT_Tset, ST_Day_Length,\
                Max_Correction, Min_Correction, Min_Trise, Max_Trise, Min_Tset, Max_Tset, \
                Max_DayLength, Min_DayLength, Sun_Window, Sun_Window_Length, \
                Max_DL_SRSS, Min_DL_SRSS)






def Resource_Estimation(ghi,dhi,dni,ws,tamb,shs,sw,zdh):
    #Resource_param.loc[shs==1,'Sun_Hour_Status'] =1
    #print(shs[0:20])
    #shs=shs.values
    ghi= ghi.values
    dhi = dhi.values
    dni = dni.values
    tamb = tamb.values
    ws = ws.values

    gh=shs*ghi
    dh=shs*dhi
    dn=shs*dni
    ta=shs*tamb
    wss=shs*ws


    GHI_SH=gh#ghi[idx[0]],ghi[shs==1]
    DHI_SH=dh
    DNI_SH=dn
    Tamb_SH=ta
    WS_SH=wss

    Daily_Sum_Tamb = np.sum(np.reshape(ta,(365,24)),axis=1)
    Day_Ag_GHI = np.sum(np.reshape(gh*0.001,(365,24)),axis=1)
    Day_Ag_DHI = np.sum(np.reshape(dh*0.001,(365,24)),axis=1)
    Day_Ag_DNI = np.sum(np.reshape(dn*0.001,(365,24)),axis=1)
    Daily_Sum_WS = np.sum(np.reshape(wss,(365,24)),axis=1)
    Day_Sunshine_Hours = np.sum(np.reshape(shs,(365,24)),axis=1)



    Day_Max_WS = np.max(np.reshape(ws,(365,24)),axis=1)
    Day_Min_WS = np.min(np.reshape(ws,(365,24)),axis=1)
    Day_Max_Tamb = np.max(np.reshape(ta,(365,24)),axis=1)
    Day_Min_Tamb = np.min(np.reshape(ta,(365,24)),axis=1)



    Mon_Ag_GHI=np.sum(monthorder(Day_Ag_GHI),axis=1)
    Mon_Ag_DNI=np.sum(monthorder(Day_Ag_DNI),axis=1)
    Mon_Ag_DHI=np.sum(monthorder(Day_Ag_DHI),axis=1)
    Mon_Ag_Sun_Hours=np.sum(monthorder(Day_Sunshine_Hours),axis=1)
    Mon_Sum_Tamb=np.sum(monthorder(Daily_Sum_Tamb),axis=1)
    Mon_Sum_WS=np.sum(monthorder(Daily_Sum_WS),axis=1)

    Mon_Max_Tamb = np.max(monthorder(Day_Max_Tamb),axis=1)
    Mon_Min_Tamb = np.min(monthorder(Day_Min_Tamb),axis=1)
    Mon_Max_WS = np.max(monthorder(Day_Max_WS),axis=1)
    Mon_Min_WS = np.min(monthorder(Day_Min_WS),axis=1)

    Day_Av_Tamb = Daily_Sum_Tamb/Day_Sunshine_Hours
    Day_Av_WS = Daily_Sum_WS/Day_Sunshine_Hours

    Mon_Av_Tamb = Mon_Sum_Tamb/Mon_Ag_Sun_Hours
    Mon_Av_WS = Mon_Sum_WS/Mon_Ag_Sun_Hours


    An_Ag_GHI_SH = np.sum(Day_Ag_GHI)/1000
    An_Ag_DHI_SH = np.sum(Day_Ag_DHI)/1000
    An_Ag_DNI_SH = np.sum(Day_Ag_DNI)/1000
    An_Ag_SunHours = np.sum(Day_Sunshine_Hours)

    temp = np.nonzero(Tamb_SH)
    Min_Tamb_SH = np.min(Tamb_SH[temp[0]])
    Max_Tamb_SH = np.max(Tamb_SH[temp[0]])
    Ave_Tamb_SH = np.mean(Tamb_SH[temp[0]])

    Min_WS_SH = np.min(WS_SH)
    Max_WS_SH = np.max(WS_SH)
    Ave_WS_SH = np.sum(Daily_Sum_WS)/np.sum(Day_Sunshine_Hours)

    Av_Ag_GHI_Per_Day = np.mean(Day_Ag_GHI)
    Av_Ag_DHI_Per_Day = np.mean(Day_Ag_DHI)
    Av_Ag_DNI_Per_Day = np.mean(Day_Ag_DNI)
    Av_Sunshine_Hours_Per_Day = np.sum(Day_Sunshine_Hours)/365


    #print(Resource_minmax['param'])
    #print('Sun Window',sw,len(shs))
    return (An_Ag_GHI_SH, An_Ag_DHI_SH, An_Ag_DNI_SH, An_Ag_SunHours, \
            Max_Tamb_SH, Min_Tamb_SH, Ave_Tamb_SH, Max_WS_SH, Min_WS_SH, \
            Ave_WS_SH, Av_Ag_GHI_Per_Day, Av_Ag_DHI_Per_Day, Av_Ag_DNI_Per_Day,\
            Av_Sunshine_Hours_Per_Day, Mon_Ag_GHI, Mon_Ag_DHI, Mon_Ag_DNI, \
            Mon_Ag_Sun_Hours, Mon_Max_Tamb, Mon_Min_Tamb, Mon_Av_Tamb, \
            Mon_Max_WS, Mon_Min_WS, Mon_Av_WS, Day_Ag_GHI, Day_Ag_DHI, \
            Day_Ag_DNI, Day_Sunshine_Hours, Day_Max_WS, Day_Min_WS, \
            Day_Av_WS, Day_Max_Tamb, Day_Min_Tamb, Day_Av_Tamb, GHI_SH, \
            DHI_SH, DNI_SH, Tamb_SH, WS_SH)


def Spacing_Factor(Module_Tilt, Ang_Sol_Azimuth, Ang_Sol_Altitude, ZT_Day_Hour):

    mask = (Ang_Sol_Altitude>=1)
    Pgen_Hour_Status = mask.astype(int)
    Pgen_Day_Hour = ZT_Day_Hour*mask
    Module_Tilt_r = np.deg2rad(Module_Tilt)
    Ang_Sol_Azimuth_r = np.deg2rad(Ang_Sol_Azimuth)
    Ang_Sol_Altitude_r = np.deg2rad(Ang_Sol_Altitude.values)
    test = (np.sin(Module_Tilt_r)*np.cos(Ang_Sol_Azimuth_r))/np.tan(Ang_Sol_Altitude_r)
    Lrow_Factor = test*mask
    test = (np.sin(Module_Tilt_r)*np.sin(Ang_Sol_Azimuth_r))/np.tan(Ang_Sol_Altitude_r)
    Lcol_Factor = test*mask
    Pgen_Window = (np.min(Pgen_Day_Hour[np.nonzero(Pgen_Day_Hour)]), \
                            np.max(Pgen_Day_Hour[np.nonzero(Pgen_Day_Hour)]))
    Annual_Pgen_Hours = np.sum(Pgen_Hour_Status)

    return(Pgen_Hour_Status, Pgen_Day_Hour, Lrow_Factor, Lcol_Factor, \
           Pgen_Window, Annual_Pgen_Hours)


#------------------------------------------------------------------------------
# Determination of net effective radiation on the tilted panel
#------------------------------------------------------------------------------
#
# This module is only for fixed tilt configuration
# Seperate functions need to be developed for module tracking
#
def Net_Effective_Radiation (Module_Tilt, Array_Height, STAngIncidence, \
                             GHI_SH, DHI_SH, DNI_SH, Albedo, Lmod, \
                             Ground_Clearance, Pgen_Hour_Status):

    #Rb = np.zeros(len(Pgen_Hour_Status))
    #Gt = np.zeros(len(Pgen_Hour_Status))
    Module_Tilt_r = np.deg2rad(Module_Tilt)

    #Rd = 0.5 * (1 + cosd(Module_Tilt))
    #Rg = 0.5 * (1 - cosd(Module_Tilt))
    Rd = 0.5 * (1 + np.cos(Module_Tilt_r))
    Rg = 0.5 * (1 - np.cos(Module_Tilt_r))

    #n = math.floor((Array_Height)/(Lmod * sind(Module_Tilt)))
    n = np.floor((Array_Height)/(Lmod * np.sin(Module_Tilt_r)))

    status = (Pgen_Hour_Status==1)
    Rb = np.cos(np.deg2rad(STAngIncidence.values))
    Gt = (DNI_SH * Rb + DHI_SH * Rd + GHI_SH * Rg*Albedo)*status
    Annual_Gt = np.sum(Gt)*1e-06#/1000000

    return (Gt, n, Annual_Gt)

#-------------------------- End of Net Effective Radiation function -----------

#------------------------------------------------------------------------------
# Determination of cell temperature and resource to module power factor
#------------------------------------------------------------------------------
#
def Tcell_ResToModPower (Pgen_Hour_Status, Gt, Tamb_SH, WS_SH, a_CT, b_CT,
                         DelT, Pmod, Kt_Pmod, Kt_Isc, Kt_Voc, Isc, Voc):

    Gref = 1000
    Tref = 25


    status = (Pgen_Hour_Status==1)
    temp = Gt/Gref
    Tcell = (Gt*np.exp(a_CT + b_CT * WS_SH) + Tamb_SH + DelT * temp)*status
    R_Pmod = temp * (1+(Kt_Pmod*0.01)*(Tcell-Tref))*status
    Agg_R_Pmod = np.sum(R_Pmod)
    Mod_Isc = Isc*temp*(1+(Kt_Isc*0.01) *(Tcell-Tref))*status
    Mod_Voc = (Voc + (Kt_Voc*0.01) * (Tcell - Tref))*status

    return (Tcell, Mod_Isc, Mod_Voc, R_Pmod, Agg_R_Pmod)

#----- End of cell temperature and resource to module power function ----------
#----- End of cell temperature and resource to module power function ----------

#------------------------------------------------------------------------------
# Determination of plant sizing parameters
#------------------------------------------------------------------------------
#
def Plant_Sizing_Parameters (V_PCU_Choice, Vuser, n, R_Pmod, Mod_Isc, Mod_Voc, \
                             Pplant_Target, P_nomDC, P_PCU_AC, Vmp, Imp, Pmod, \
                             V_maxDC, I_maxDC, SL_PC, EL_PC, Lmod, Bmod, \
                             Vmppmax, Vmppmin):

    if V_PCU_Choice == 1:
        V_PCU_Ref = 0.5 * (Vmppmax + Vmppmin)
    elif V_PCU_Choice == 2:
        V_PCU_Ref = Vuser

    Mod_Pmax = np.max(R_Pmod) * Pmod
    Mod_Voc_min = np.min(Mod_Voc[np.nonzero(Mod_Voc)])
    Mod_Voc_max = np.max(Mod_Voc[np.nonzero(Mod_Voc)])
    Mod_Isc_max = np.max(Mod_Isc[np.nonzero(Mod_Voc)])
    N_PCU = Pplant_Target *1000 // P_nomDC
    m_max = V_maxDC // Mod_Voc_max
    m = np.ceil(V_PCU_Ref / Vmp )
    I_PCU_Ref = P_nomDC * 1000 / V_PCU_Ref
    y_max = I_maxDC // (n * Mod_Isc_max)
    y = np.ceil (I_PCU_Ref / (n * Imp))
    N_mod_PCU = n * m * y
    m_change = 0
    DC_Losses = (1 - SL_PC)
    P_PCU_DC_Max = N_mod_PCU * Mod_Pmax * DC_Losses * 0.001# 1000

#    m_Pmod_max = Mod_Pmax * m/1000
#    P_PCU_DC_Max = 230

    Original_P_plant = N_mod_PCU * N_PCU * Pmod * 1e-06 #/ 1000000

    if (m > m_max) or (y > y_max):
        print('Invalid choice for PCU reference voltage')

    if (P_PCU_DC_Max > P_nomDC):
        while(P_PCU_DC_Max > P_nomDC):
            N_mod_PCU_new = N_mod_PCU - m
            m_change = m_change - 1
            P_PCU_DC_Max = N_mod_PCU_new * Mod_Pmax * DC_Losses*0.001 #/1000
            N_mod_PCU = N_mod_PCU_new
            y = N_mod_PCU / (n * m)

    else:
        while (P_PCU_DC_Max <= P_nomDC):
            N_mod_PCU_new =  N_mod_PCU + m
            m_change = m_change + 1
            P_PCU_DC_Max = N_mod_PCU_new * Mod_Pmax * DC_Losses *0.001 #/1000
            N_mod_PCU = N_mod_PCU_new
            y = N_mod_PCU / (n * m)

    Pmax_Deviation_NewCap =  (P_PCU_DC_Max/P_nomDC) - 1
    y_area = np.ceil (y)

    N_Plant = N_mod_PCU * N_PCU
    P_plant = N_Plant * Pmod * 1e-06#/ 1000000
    Pure_mod_area = N_Plant * Lmod * Bmod
    Pure_mod_area_acres = Pure_mod_area * 0.00024710
    DC_Scaling_Factor = N_mod_PCU * Pmod / (1000 * P_PCU_AC)
    P_plant_PCU = N_PCU * P_PCU_AC * 0.001#/ 1000


    return (Mod_Pmax, Mod_Voc_min, Mod_Voc_max, Mod_Isc_max, N_PCU, m, \
            m_change, N_mod_PCU, y, y_area, N_Plant, Pure_mod_area, \
            Pure_mod_area_acres, DC_Scaling_Factor, P_plant, P_plant_PCU, \
            Original_P_plant, Pmax_Deviation_NewCap, P_PCU_DC_Max)

#-------- End of plant sizing parameters function -----------------------------
#------------------------------------------------------------------------------
# Estimation of inter-row and inter-column spacing
#------------------------------------------------------------------------------
#
def Inter_Row_Col_Spacing (Pgen_Window, Sun_Hour_Status, Pgen_Hour_Status, \
                           ZT_Day_Hour, L_row_Factor, L_col_Factor, n, Lmod, \
                           R_Pmod, Gt):

    Pgen_2hr_Window = (Pgen_Window [0] + 1, Pgen_Window [1] - 1)
    Pgen_4hr_Window = (Pgen_Window [0] + 2, Pgen_Window [1] - 2)
    Agg_R_Pmod = np.zeros (3)
    Annual_GenHr = np.zeros (3)
   # index = 0 --> Pgen hrs, 1 --> Pgen - 2 hrs,  2 --> Pgen - 4hrs

    Drow = np.zeros(3)
    Dcol = np.zeros(3)



    temp = (Sun_Hour_Status==1)
    Lrow = n*Lmod*L_row_Factor*temp
    Lcol = n*Lmod*L_col_Factor*temp

    temp = (Pgen_Hour_Status==1)
    Agg_R_Pmod[0] = np.sum(R_Pmod*temp)
    Annual_GenHr[0] = np.sum(temp)
    Gt_Pgen = Gt*temp
    Drow[0] = np.max(Lrow*temp)
    Dcol [0] = np.max(Lcol*temp)

    temp = ((ZT_Day_Hour >= Pgen_2hr_Window [0])*(ZT_Day_Hour <= Pgen_2hr_Window [1]))
    Agg_R_Pmod[1] = np.sum(R_Pmod*temp)
    Annual_GenHr[1] = np.sum(temp)
    Gt_Pgen_2hr = Gt*temp
    Drow[1] = np.max(Lrow*temp)
    Dcol [1] = np.max(Lcol*temp)

    temp = ((ZT_Day_Hour >= Pgen_4hr_Window [0])*(ZT_Day_Hour <= Pgen_4hr_Window [1]))
    Agg_R_Pmod [2] = np.sum(R_Pmod*temp)
    Annual_GenHr[2] = np.sum(temp)
    Gt_Pgen_4hr = Gt*temp
    Drow[2] = np.max(Lrow*temp)
    Dcol [2] = np.max(Lcol*temp)

    # Annaul radiation on the tilted panel in MW / sq.m
    An_Gt = np.zeros (3)
    An_Gt [0] = np.sum(Gt_Pgen)*1e-06#/1000000
    An_Gt [1] = np.sum(Gt_Pgen_2hr)*1e-06#/1000000
    An_Gt [2] = np.sum(Gt_Pgen_4hr)*1e-06#/1000000

    return (Drow, Dcol, Annual_GenHr, Agg_R_Pmod, Pgen_2hr_Window, \
            Pgen_4hr_Window, Gt_Pgen, Gt_Pgen_2hr, Gt_Pgen_4hr, An_Gt)

#-------- End of Inter-row and Inter-col spacing  function -------------------

#------------------------------------------------------------------------------
# Estimation of Plant area (includes spiral subfunction)
#------------------------------------------------------------------------------
#
# Defining the rectangular spiral
#
def spiral (X):
    #
    # Basic initialisations
    # .....................
    # Coefficients for tracing area components of the block
    A_LB = 0
    A_LDc = 0
    A_BDr = 0
    #
    # Coefficients for tracing transition set for LDc and BDr components
    A_LDc_temp = 0
    A_BDr_temp = 0
    #
    # Counter to trace the growth of units along length (y- axis) and
    # rectangular matrix formation (Default assumption rectangular matrix)
    R = 1
    #
    # Counter to trace the growth of units along breadth (x - axis) and
    # square matrix formation
    S = 0
    #
    # Counter to trace the second entity of the transition set
    count = 0
    #
    No = 0      # Number of outlying units
    NLe = 0     # No. of units along length contributing to enclosed set
    NLo = 0     # No. of units along length contributing to outlying set
    NBe = 0     # No. of units along breadth contributing to enclosed set
    NBo = 0     # No. of units along breadth contributing to outlying set
    #
    # No.of inter-row spacing sec. along length contributing to enclosed set
    Nre = 0
    # No.of inter-row spacing sec. along length contributing to outlying set
    Nro = 0
    # No.of inter-col spacing sec. along breadth contributing to enclosed set
    Nce = 0
    # No.of inter-col spacing sec. along breadth contributing to outlying set
    Nco = 0
    #
    # Flag to trace outlying set additon along length
    L_flag = 0
    # Flag to trace outlying set addition along breadth
    B_flag = 0
    #
    # Tracing the growth along length and breadth to determine R, S
    #
    ii = 1
    while ii<=X:
        # Factoring area contribution by the unit itself
        A_LB += 1
        # Check if the no of units is a perfect square --> square matrix check
        if np.sqrt(A_LB)% 1 == 0:
            S += 1

        # Underpinning condition to set the growth of matrix along length
        # (y- axis) at X = 2 as a rectangular matrix and set the growth of
        # inter-row spacing
        if ii == 2:
            A_BDr = 1
            R += 1

        # Accounting for inter-row & inter column spacing for successive block
        if ii >= 3:
            A_BDr_temp += 1
            A_BDr = A_BDr_temp
            A_LDc_temp += 1
            A_LDc = A_LDc_temp
            #
            # At the end of every rectangular matrix, the area component LDc
            # would be a perfect square, the next element would share the same
            # LDc area component. These elements form the 'transition set(TS)'
            # Listed transition sets X = 6,7; 12,13,; 21,20; ... The first
            # number of the transition set forms a rectangular matrix. Hence R
            # incremented by 1. Since the spacing is shared, the LDc and BDr
            # component addition is reduced by 1 for the first element of the
            # transition set and it is incremented by 1 for the second element
            # of the transition set and the count proceeds.
            #
            if A_LDc > 1:
                if np.sqrt(A_LDc_temp)%1.0 == 0 and count < 1:
                    A_LDc_temp -= 1
                    A_BDr_temp -= 1
                    count += 1 # First element of TS traced
                    R += 1

                # Resetting the TS counter when crossing the 2nd element of TS
                if np.sqrt(A_LDc_temp)%1.0 == 0 and count == 1:
                    count = 0

        ii += 1

    # Estimating the number of outlying units
    No = X - R*S
    #
    # If the enclosed set forms a perfect square: (R = S)
    # The outlying units are added along the breadth till the TS is reached.
    #
    if R == S:
        # Accounting for enclosed set
        NLe = R
        Nre = R - 1
        NBe = S
        Nce = S - 1
        # Accounting for outlying set
        if No > 0:
            NLo = 1     # Fixing the dimension of units along length
            NBo = No    # Tracing unit growth along breadth
            Nro = NLo
            Nco = NBo - 1
            B_flag = 1


    # If the enclosed set forms a perfect rectangle: (R > S)
    # The outlying units are added along length till perfect square is reached.
    #
    if R > S:
        NLe = R
        Nre = R - 1
        NBe = S
        Nce = S - 1
        #
        # Accounting for outlying set
        #
        if No > 0:
            NLo = No
            NBo = 1
            Nro = NLo - 1
            Nco = NBo
            L_flag = 1


    # Consolidated coefficients
    # Contribution due to units
    LB = NLe * NBe + NLo * NBo
    # Contribution due to unit breadth and inter-row spacing
    BDrow = NBe * Nre + NBo * Nro
    # Contribution due to unit length and inter-col spacing
    LDcol = NLe * Nce + NLo * Nco

    # Contribution due to inter-row and inter-col spacing.
    # Initialising inter-row inter-col area component at X = 4 units onwards
    if X < 4:
        DrowDcol = 0

    if X == 4:
        DrowDcol = 1

    if X > 4:
        DrowDcol = Nre * Nce + Nro * Nco

    return (NLe, NBe, NLo, NBo, Nre, Nce, Nro, Nco, L_flag, B_flag)

#----------------End of spiral sub-function definition-------------------------

#------------------------------------------------------------------------------
#   Estimating Plant area
#------------------------------------------------------------------------------

def Plant_Area_Estimation (n, Lmod, Bmod, Module_Tilt, D_row, D_col, m, y_area, \
                           N_PCU, BS_a, BS_b, P_plant, Pure_mod_area_acres, \
                           Pgen_Window, Pgen_2hr_Window, Pgen_4hr_Window, \
                           Gt_Pgen, Gt_Pgen_2hr, Gt_Pgen_4hr):

    N_Time_Windows = len(D_row)

    LArray_e = np.zeros (N_Time_Windows)
    BArray_e = np.zeros (N_Time_Windows)
    LArray_o = np.zeros (N_Time_Windows)
    BArray_o = np.zeros (N_Time_Windows)
    PCU_area = np.zeros (N_Time_Windows)
    PCU_area_acres = np.zeros (N_Time_Windows)
    L_PCU = np.zeros (N_Time_Windows)
    B_PCU = np.zeros (N_Time_Windows)
    Gross_PCU_area = np.zeros (N_Time_Windows)
    Gross_PCU_area_acres = np.zeros (N_Time_Windows)

    LPCU_e = np.zeros (N_Time_Windows)
    BPCU_e = np.zeros (N_Time_Windows)
    LPCU_o = np.zeros (N_Time_Windows)
    BPCU_o = np.zeros (N_Time_Windows)
    P_area = np.zeros (N_Time_Windows)
    P_area_acres = np.zeros (N_Time_Windows)
    L_P = np.zeros (N_Time_Windows)
    B_P = np.zeros (N_Time_Windows)
    Gross_P_area = np.zeros (N_Time_Windows)
    Gross_P_area_acres = np.zeros (N_Time_Windows)

    B_L_P = np.zeros (N_Time_Windows)
    B_B_P = np.zeros (N_Time_Windows)
    Total_Plant_area = np.zeros (N_Time_Windows)
    Total_Plant_area_acres = np.zeros (N_Time_Windows)
    Aux_P = np.zeros (N_Time_Windows)
    Total_Plant_area_with_Aux = np.zeros (N_Time_Windows)
    Total_Plant_area_with_Aux_acres = np.zeros (N_Time_Windows)
    Del_Total_P_gross_P = np.zeros (N_Time_Windows)
    Del_Total_P_Aux_Total_P = np.zeros (N_Time_Windows)
    Aspect_ratio = np.zeros (N_Time_Windows)
    Packing_density_net = np.zeros (N_Time_Windows)
    Packing_density_gross  = np.zeros (N_Time_Windows)
    Packing_density_Total_aux = np.zeros (N_Time_Windows)
    Deviation_Factor_Total_aux = np.zeros (N_Time_Windows)

    Abs_Deviation = np.zeros(N_Time_Windows)
    Gt_TW = np.zeros(row_count)
    CERC_Benchmark = 5

    # Application of ulam spiral for placement of arrays for a PCU block
    NLe_PCU, NBe_PCU, NLo_PCU, NBo_PCU, Nre_PCU, Nce_PCU, Nro_PCU, Nco_PCU, \
    L_flag_PCU, B_flag_PCU =  spiral (y_area)

    # Calculating the base length and breadth dimensions of an array
    L = n * Lmod * np.cos(np.deg2rad(Module_Tilt))
    B = m * Bmod

    ii = 0
    # Calculating the PCU block dimensions and area
    while ii < N_Time_Windows:
        LArray_e [ii] = NLe_PCU * L + Nre_PCU * D_row [ii]
        BArray_e [ii] = NBe_PCU * B + Nce_PCU * D_col [ii]
        LArray_o [ii] = NLo_PCU * L + Nro_PCU * D_row [ii]
        BArray_o [ii] = NBo_PCU * B + Nco_PCU * D_col [ii]
        PCU_area [ii] = LArray_e [ii] * BArray_e [ii] + LArray_o [ii] * BArray_o [ii]
        PCU_area_acres [ii] = PCU_area [ii] * 0.00024710
        L_PCU [ii] = LArray_e [ii] + B_flag_PCU * (L + D_row [ii])
        B_PCU [ii] = BArray_e [ii] + L_flag_PCU * (B + D_col [ii])
        Gross_PCU_area [ii] = L_PCU [ii] * B_PCU [ii]
        Gross_PCU_area_acres [ii] = Gross_PCU_area [ii] * 0.00024710
        ii += 1

    # Application of Ulam Spiral for placement of PCU block for the plant
    NLe_P, NBe_P, NLo_P, NBo_P, Nre_P, Nce_P, Nro_P, Nco_P, L_flag_P, B_flag_P = spiral (N_PCU)
    ii = 0
    #
    # Calculating plant dimensions and area
    while ii < N_Time_Windows:
        LPCU_e [ii] = NLe_P * L_PCU [ii] + Nre_P * D_row [ii]
        BPCU_e [ii] = NBe_P * B_PCU [ii] + Nce_P * D_col [ii]
        LPCU_o [ii] = NLo_P * L_PCU [ii] + Nro_P * D_row [ii]
        BPCU_o [ii] = NBo_P * B_PCU [ii] + Nco_P * D_col [ii]
        P_area [ii] = LPCU_e [ii] * BPCU_e [ii] + LPCU_o [ii] * BPCU_o [ii]     # Net Plant Area in sq.m
        P_area_acres [ii] = P_area [ii]* 0.00024710                             # Net Plant Area in acres
        L_P [ii] = LPCU_e [ii] + B_flag_P * (L_PCU [ii] + D_row [ii])
        B_P [ii] = BPCU_e [ii] + L_flag_P * (B_PCU [ii] + D_col [ii])
        Gross_P_area [ii] = L_P [ii] * B_P [ii]                                 # Effective or Gross Plant Area in sq.m
        Gross_P_area_acres [ii] = Gross_P_area [ii] * 0.00024710                # Effective or Gross Plant Area in acres
        ii += 1

    ii = 0
    # Application of boundary area and auxiliary area of plant
    while ii < N_Time_Windows:
        B_L_P [ii] = 2 * BS_b * (L_P[ii] + 2 * BS_a)
        B_B_P [ii] = 2 * BS_a * B_P[ii]
        Total_Plant_area [ii] = B_L_P[ii] + B_B_P[ii] + Gross_P_area[ii]
        Total_Plant_area_acres [ii] = Total_Plant_area [ii] * 0.00024710
        if P_plant <= 1:
            Aux_P [ii] = 0.165 * Total_Plant_area [ii]
        elif P_plant > 1 and P_plant <= 100:
            Aux_P [ii] = 0.16273 * np.exp(-0.027 * P_plant) * Total_Plant_area[ii]
        else:
            Aux_P [ii] = 0.01 * Total_Plant_area [ii]

        # Estimating total plant area with auxiliary area requirements
        Total_Plant_area_with_Aux [ii] = Total_Plant_area [ii] + Aux_P [ii]
        Total_Plant_area_with_Aux_acres [ii] = Total_Plant_area_with_Aux [ii] * 0.00024710
        #
        # Accounting for excess auxiliary area
        Del_Total_P_gross_P [ii] = Total_Plant_area_acres [ii] - Gross_P_area_acres [ii]
        Del_Total_P_Aux_Total_P [ii] = Total_Plant_area_with_Aux_acres [ii] - Total_Plant_area_acres [ii]
        #
        if Del_Total_P_gross_P [ii] >= Del_Total_P_Aux_Total_P[ii]:
            Total_Plant_area_with_Aux_acres [ii] = Total_Plant_area_acres [ii]
        elif Del_Total_P_gross_P [ii] < Del_Total_P_Aux_Total_P[ii]:
            Total_Plant_area_with_Aux_acres [ii] = Gross_P_area_acres [ii] + Del_Total_P_Aux_Total_P[ii]

        ii += 1

    # Area related metrics
    ii = 0
    while ii < N_Time_Windows:
        #
        # Spiral related metrics
        Aspect_ratio [ii] = L_P[ii]/B_P[ii]
        Packing_density_net [ii] = Pure_mod_area_acres/P_area_acres[ii]
        Packing_density_gross [ii] = Pure_mod_area_acres/Gross_P_area_acres[ii]
        Packing_density_Total_aux [ii] = Pure_mod_area_acres/Total_Plant_area_with_Aux_acres [ii]
        Deviation_Factor_Total_aux [ii] = (Total_Plant_area_with_Aux_acres [ii] - P_plant*CERC_Benchmark)/(P_plant*CERC_Benchmark)
        ii += 1

    # Selection of time window of interest
    Abs_Deviation = np.absolute(Deviation_Factor_Total_aux)
    TW = 0
    Del_BA_TA_TW = 999

    for i in range(N_Time_Windows):
        if Abs_Deviation [i] < Del_BA_TA_TW :
            Del_BA_TA_TW = Abs_Deviation [i]
            TW = i
        if TW == 0:
            TW_DayHrs = (Pgen_Window [0], Pgen_Window [1])
            Gt_TW = Gt_Pgen
        elif TW == 1:
            TW_DayHrs = (Pgen_2hr_Window [0], Pgen_2hr_Window [1])
            Gt_TW = Gt_Pgen_2hr
        elif TW == 2:
            TW_DayHrs = (Pgen_4hr_Window [0], Pgen_4hr_Window [1])
            Gt_TW = Gt_Pgen_4hr


    return (Total_Plant_area_with_Aux_acres, Packing_density_Total_aux, \
            Deviation_Factor_Total_aux, TW, TW_DayHrs, Gt_TW )

#------------------------ End of Plant area estimation function ---------------

#------------------------------------------------------------------------------
# Estimation of Plant output
#------------------------------------------------------------------------------
#
def Plant_output_year_0 (TW, An_Gt, Pure_Mod_Area, SL_PC, EL_PC, TW_DayHrs, \
                         N_Plant, ZT_Day_Hour, R_Pmod, Pmod, Eta_PCU, \
                         P_plant, Gt_TW):

    # Data initialisation
    #TW_Status = np.zeros(len(Gt_TW))

    #R_Pplant_TW = np.zeros (len(Gt_TW))
    #Pmod_out = np.zeros (len(Gt_TW))
    #P_plant_PV_DC = np.zeros (len(Gt_TW))
    #P_plant_PV_DC_inc_SL = np.zeros (len(Gt_TW))
    #P_plant_PV_AC = np.zeros (row_count)
    Agg_R_Pplant_AC_TW = 0.0
    #Hourly_Gt_all_pmod = np.zeros(len(Gt_TW))
    #Hourly_Gen_PC = np.zeros(len(Gt_TW))
    #Hourly_SEE_PV = np.zeros (len(Gt_TW))

    # Aggregate solar radiation incident on solar panels in MW/sq.m
    Agg_Gt_all_modules_TW = An_Gt [TW] * Pure_Mod_Area

    Soiling_Loss_Factor = 1 - SL_PC
    Electrical_Loss_Factor = 1 - EL_PC

    hour = (ZT_Day_Hour >= TW_DayHrs[0])*(ZT_Day_Hour <=TW_DayHrs [1])
    R_Pplant_TW = N_Plant * R_Pmod * hour * 1e-06
    Pmod_out = Pmod * R_Pmod * hour
    P_plant_PV_DC = Pmod * R_Pplant_TW * hour
    P_plant_PV_DC_inc_SL = P_plant_PV_DC * hour * Soiling_Loss_Factor
    P_plant_PV_AC = P_plant_PV_DC_inc_SL * hour * Electrical_Loss_Factor * Eta_PCU
    Hourly_Gen_PC = (P_plant_PV_AC * hour)/P_plant
    Hourly_Gt_all_pmod = Gt_TW * hour * Pure_Mod_Area * 1e-06
    temp = (Gt_TW*hour!=0)
    temp1=Hourly_Gt_all_pmod
    Hourly_Gt_all_pmod[temp==False] = 1
    Hourly_SEE_PV = (temp*P_plant_PV_AC)/Hourly_Gt_all_pmod
    Hourly_Gt_all_pmod[temp==False] = 0

    Agg_R_Pplant_AC_TW = np.sum(R_Pplant_TW * Soiling_Loss_Factor * \
                            Electrical_Loss_Factor * Eta_PCU)
    # Annual parameters in MWh
    An_Pplant_DC = np.sum(P_plant_PV_DC)
    An_Pplant_DC_inc_SL = np.sum(P_plant_PV_DC_inc_SL)
    An_Pplant_AC = np.sum(P_plant_PV_AC)
    An_soiling_loss = An_Pplant_DC - An_Pplant_DC_inc_SL
    An_electrical_loss = An_Pplant_DC_inc_SL - An_Pplant_AC

    An_PR_PV = (An_Pplant_AC * 1000)/ (An_Gt [TW] * N_Plant * Pmod )
    An_CUF_PV = An_Pplant_AC / (8760 * P_plant)
    An_SEE_PV = An_Pplant_AC / (Agg_Gt_all_modules_TW)

    # Adding An_Pplant_DC_inc_SL in the return list, 8 May 2019
    return (Hourly_Gen_PC, Hourly_SEE_PV, Hourly_Gt_all_pmod, \
            P_plant_PV_AC, Agg_Gt_all_modules_TW, An_Pplant_DC, \
            An_Pplant_DC_inc_SL, An_Pplant_AC, An_soiling_loss, \
            An_electrical_loss, An_PR_PV, An_CUF_PV, An_SEE_PV)

#------------------------ End of Plant output for year 0 function -------------

#------------------------------------------------------------------------------
# Estimation of Plant output factoring module degradation
#------------------------------------------------------------------------------
#
def Degradation_Results (BaseDegRate, Plant_Life, Pmod, Mod_Rating_EoY1, \
                        N_Plant, Agg_R_Pmod, Eta_PCU, SL_PC, EL_PC, P_plant, \
                        Pure_Mod_Area, An_Gt):
    Degradation_Rates = np.array([-0.005, BaseDegRate, -0.01, -0.015, -0.02, \
                                  -0.025, -0.03, -0.035, -0.04])

    Degradation_Rates = sorted (Degradation_Rates, reverse = True)
    Base_DegCase_Index = Degradation_Rates.index(BaseDegRate)

    Agg_Gt_all_Modules = np.zeros(3)

    Net_module_rating = np.zeros((len(Degradation_Rates), Plant_Life + 1))
    Module_Rating_PC = np.zeros((len(Degradation_Rates), Plant_Life + 1))
    Annual_PV_Gen_Deg = np.zeros((3,len(Degradation_Rates), Plant_Life + 1))
    Annual_CUF_PV = np.zeros((3,len(Degradation_Rates), Plant_Life + 1))
    Annual_PR_PV = np.zeros((3,len(Degradation_Rates), Plant_Life + 1))
    Annual_SEE_PV = np.zeros((3,len(Degradation_Rates), Plant_Life + 1))

    Agg_RPF = np.zeros(3)

    for i in range(3):
        Agg_Gt_all_Modules [i] = An_Gt[i] * Pure_Mod_Area
        Agg_RPF [i] = N_Plant * Agg_R_Pmod [i] * (1 - SL_PC) * (1 - EL_PC) * \
                   Eta_PCU*1e-06 #1000000


    for Deg_ind in range (len(Degradation_Rates)):
        for yr in range (0, Plant_Life + 1):
            if yr == 0:
                Module_Rating_PC [Deg_ind][yr] = 1
                Net_module_rating [Deg_ind][yr] = Pmod

            elif yr == 1:
                Module_Rating_PC [Deg_ind][yr] = Mod_Rating_EoY1
                Net_module_rating [Deg_ind][yr] = Pmod * \
                                                Module_Rating_PC [Deg_ind][yr]

            elif yr > 1:
                Module_Rating_PC [Deg_ind][yr] = \
                Module_Rating_PC [Deg_ind][yr - 1] + Degradation_Rates [Deg_ind]
                Net_module_rating [Deg_ind][yr] = Pmod * \
                                                Module_Rating_PC [Deg_ind][yr]

    for TimeWindow in range(3):
        for Deg_ind in range(0, len(Degradation_Rates)):
            for yr in range(0, Plant_Life + 1):
                Annual_PV_Gen_Deg [TimeWindow][Deg_ind][yr] = \
                    Net_module_rating [Deg_ind][yr] * Agg_RPF [TimeWindow]

                Annual_CUF_PV [TimeWindow][Deg_ind][yr] = \
                    Annual_PV_Gen_Deg [TimeWindow][Deg_ind][yr] /(8760 * Pplant)

                Annual_PR_PV [TimeWindow][Deg_ind][yr] = \
                    Annual_PV_Gen_Deg [TimeWindow][Deg_ind][yr]*1000 / \
                    (An_Gt [TimeWindow] \
                            * N_Plant * Net_module_rating [Deg_ind][yr] )

                Annual_SEE_PV [TimeWindow][Deg_ind][yr] = \
                    Annual_PV_Gen_Deg [TimeWindow][Deg_ind][yr] / \
                                    (Agg_Gt_all_Modules [TimeWindow])

    return (Degradation_Rates, Base_DegCase_Index, Agg_RPF, Annual_PV_Gen_Deg,\
            Net_module_rating, Annual_CUF_PV, Module_Rating_PC, Annual_PR_PV, \
            Annual_SEE_PV)

# the financial model uses only generation (Annual_PV_Gen_Deg)
# for time window = TW, and base degradation case index = Base_DegCase_Index

# -------------------------- end of degradation results function --------------

#------------------------------------------------------------------------------
# Estimating monthly, daily aggregate gen, Histogram for frequency of gen
#------------------------------------------------------------------------------
#
def Monthly_Daily_Gen_Hist1 (P_plant_PV_AC, ZT_Day_Hour, P_plant_PCU, \
                             ZoneTime):

    PV_Gen_Yr0_Daily = np.sum (P_plant_PV_AC.reshape(365,-1), axis = 1)

    Gen_Hrs_Yr0_TW_PC = np.zeros(11)

    PV_Gen_Yr0_Monthly = np.sum(monthorder(PV_Gen_Yr0_Daily),axis=1)

    Gen_Hrs_Yr0_TW_PC [0] = np.sum((P_plant_PV_AC >= 0.01 * P_plant_PCU) * (P_plant_PV_AC <= 0.1 * P_plant_PCU))
    Gen_Hrs_Yr0_TW_PC [1] = np.sum((P_plant_PV_AC > 0.1 * P_plant_PCU) * (P_plant_PV_AC <= 0.2 * P_plant_PCU))
    Gen_Hrs_Yr0_TW_PC [2] = np.sum((P_plant_PV_AC > 0.2 * P_plant_PCU) * (P_plant_PV_AC <= 0.3 * P_plant_PCU))
    Gen_Hrs_Yr0_TW_PC [3] = np.sum((P_plant_PV_AC > 0.3 * P_plant_PCU) * (P_plant_PV_AC <= 0.4 * P_plant_PCU))
    Gen_Hrs_Yr0_TW_PC [4] = np.sum((P_plant_PV_AC > 0.4 * P_plant_PCU) * (P_plant_PV_AC <= 0.5 * P_plant_PCU))
    Gen_Hrs_Yr0_TW_PC [5] = np.sum((P_plant_PV_AC > 0.5 * P_plant_PCU) * (P_plant_PV_AC <= 0.6 * P_plant_PCU))
    Gen_Hrs_Yr0_TW_PC [6] = np.sum((P_plant_PV_AC > 0.6 * P_plant_PCU) * (P_plant_PV_AC <= 0.7 * P_plant_PCU))
    Gen_Hrs_Yr0_TW_PC [7] = np.sum((P_plant_PV_AC > 0.7 * P_plant_PCU) * (P_plant_PV_AC <= 0.8 * P_plant_PCU))
    Gen_Hrs_Yr0_TW_PC [8] = np.sum((P_plant_PV_AC > 0.8 * P_plant_PCU) * (P_plant_PV_AC <= 0.9 * P_plant_PCU))
    Gen_Hrs_Yr0_TW_PC [9] = np.sum((P_plant_PV_AC > 0.9 * P_plant_PCU) * (P_plant_PV_AC <= 1.0 * P_plant_PCU))
    Gen_Hrs_Yr0_TW_PC [10] = np.sum(P_plant_PV_AC > P_plant_PCU)
    Ag_Gen_Hrs_1to100PC_TW = np.sum(Gen_Hrs_Yr0_TW_PC)

    return (PV_Gen_Yr0_Daily, PV_Gen_Yr0_Monthly, Gen_Hrs_Yr0_TW_PC, \
            Ag_Gen_Hrs_1to100PC_TW)

#---End of function which estimates monthly and daily gen and hist 1-----------


#------------------------------------------------------------------------------
# Histogram of power generation vs hour of the day
#------------------------------------------------------------------------------
#
def Gen_vs_DayHr_Hist2 (P_plant_PV_AC):

    Gen_PV_AC_TW_24x365 = P_plant_PV_AC.reshape(365, -1)
    Gen_PV_AC_TW_24x365_T = Gen_PV_AC_TW_24x365.transpose()
    Min_Gen_PV_AC_TW_24x365 = \
            np.min(Gen_PV_AC_TW_24x365_T.reshape(24,-1), axis = 1)
    Max_Gen_PV_AC_TW_24x365 = \
            np.max(Gen_PV_AC_TW_24x365_T.reshape(24,-1), axis = 1)
    Mean_Gen_PV_AC_TW_24x365 = \
            np.mean(Gen_PV_AC_TW_24x365_T.reshape(24,-1), axis = 1)

    # Added sum and percentage share of hourly generation on 6 May 2019
    Sum_Gen_PV_AC_TW_24x365 = \
        np.sum(Gen_PV_AC_TW_24x365_T.reshape(24, -1), axis = 1)

    PC_Share_Gen_AC_TW_24x365 = \
        Sum_Gen_PV_AC_TW_24x365/sum(Sum_Gen_PV_AC_TW_24x365)


    return(Gen_PV_AC_TW_24x365, Gen_PV_AC_TW_24x365_T, Min_Gen_PV_AC_TW_24x365,\
            Max_Gen_PV_AC_TW_24x365, Mean_Gen_PV_AC_TW_24x365, \
            Sum_Gen_PV_AC_TW_24x365, PC_Share_Gen_AC_TW_24x365 )



#----------------End of function Gen_vs_DayHr_Hist2 ---------------------------

def Misc_Tech_Parameters (Plant_Area, An_PV_Gen_Deg, TW, \
                          Base_Deg_Ind, An_Gt, An_Gen_Hr, \
                          Hr_Gen_PC, Hr_Gtallmod, Hr_SEE_PV, Gt_TW):

    An_PV_Gen_per_acre = np.zeros(len(Plant_Area))

    # Creating 24x365 matrix for hourly generation percentage, solar to electric
    # efficiency and incident radiation on titled panel, added on 7 May 2019
    Gen_PC_TW_24x365 = Hr_Gen_PC.reshape(365, -1)
    Gen_PC_TW_24x365_T = Gen_PC_TW_24x365.transpose()

    Gtallmod_TW_24x365 = Hr_Gtallmod.reshape(365, -1)
    Gtallmod_TW_24x365_T = Gtallmod_TW_24x365.transpose()

    SEE_PV_TW_24x365 = Hr_SEE_PV.reshape(365, -1)
    SEE_PV_TW_24x365_T = SEE_PV_TW_24x365.transpose()


    for i in range(len(Plant_Area)):
        An_PV_Gen_per_acre [i] = \
            An_PV_Gen_Deg [i][Base_Deg_Ind][0] / Plant_Area [i]

    Area_Pgen_Pgen2hr_Factor = Plant_Area [0] / Plant_Area [1]
    Area_Pgen2hr_Pgen4hr_Factor = Plant_Area [1] / Plant_Area [2]
    Area_Pgen_Pgen4hr_Factor = Plant_Area [0] / Plant_Area [2]

    Gt_Pgen_Pgen2hr_Factor = An_Gt [0] / An_Gt [1]
    Gt_Pgen2hr_Pgen4hr_Factor = An_Gt [1] / An_Gt [2]
    Gt_Pgen_Pgen4hr_Factor = An_Gt [0] / An_Gt [2]

    Gen_Pgen_Pgen2hr_Factor = \
        An_PV_Gen_Deg [0][Base_Deg_Ind][0] / An_PV_Gen_Deg [1][Base_Deg_Ind][0]
    Gen_Pgen2hr_Pgen4hr_Factor = \
        An_PV_Gen_Deg [1][Base_Deg_Ind][0] / An_PV_Gen_Deg [2][Base_Deg_Ind][0]
    Gen_Pgen_Pgen4hr_Factor = \
        An_PV_Gen_Deg [0][Base_Deg_Ind][0] / An_PV_Gen_Deg [2][Base_Deg_Ind][0]

    Gen_PC_Pgen2hr_Pgen = \
        An_PV_Gen_Deg [1][Base_Deg_Ind][0] / An_PV_Gen_Deg [0][Base_Deg_Ind][0]

    Gen_PC_Pgen4hr_Pgen = \
        An_PV_Gen_Deg [2][Base_Deg_Ind][0] / An_PV_Gen_Deg [0][Base_Deg_Ind][0]

    # Generation hour factor, added on 7 May 2019
    GenHr_Pgen_Pgen2hr_Factor = An_Gen_Hr [0] / An_Gen_Hr [1]
    GenHr_Pgen2hr_Pgen4hr_Factor = An_Gen_Hr [1] / An_Gen_Hr [2]
    GenHr_Pgen_Pgen4hr_Factor = An_Gen_Hr [0] / An_Gen_Hr [2]

    return (An_PV_Gen_per_acre, Area_Pgen_Pgen2hr_Factor, \
            Area_Pgen2hr_Pgen4hr_Factor, Area_Pgen_Pgen4hr_Factor, \
            Gt_Pgen_Pgen2hr_Factor, Gt_Pgen2hr_Pgen4hr_Factor, \
            Gt_Pgen_Pgen4hr_Factor, Gen_Pgen_Pgen2hr_Factor, \
            Gen_Pgen2hr_Pgen4hr_Factor, Gen_Pgen_Pgen4hr_Factor, \
            Gen_PC_Pgen2hr_Pgen, Gen_PC_Pgen4hr_Pgen, GenHr_Pgen_Pgen2hr_Factor, \
            GenHr_Pgen2hr_Pgen4hr_Factor, GenHr_Pgen_Pgen4hr_Factor, \
            Gen_PC_TW_24x365_T, Gtallmod_TW_24x365_T, SEE_PV_TW_24x365_T)

#--------- End of function which calculates misc parameters for tech model-----




ui.User_Assumed_Inputs()
ts1=time.time_ns()
ZoneTime, ZTDayHour,SolarTime, STAng_Hour_f, STAngSolAltitude, \
            STAngSolAzimuth, STAngIncidence, SunHourStatus,\
            AnnualSunHours,CorrectionHHMMSS, ZTTrise, ZTTset, STDayLength,\
            MaxCorrection, MinCorrection, MinTrise, MaxTrise, MinTset, MaxTset, \
            MaxDayLength, MinDayLength, SunWindow, SunWindowLength, MaxDLSRSS, \
            MinDLSRSS \
            = ZoneTime_SolarAngles(ui.User_Assumed_Inputs.Location_Latitude,\
                                    ui.User_Assumed_Inputs.Location_Longitude,\
                                    ui.User_Assumed_Inputs.Module_Tilt, \
                                    ui.User_Assumed_Inputs.Module_Surface_Azimuth,\
                                    ui.User_Assumed_Inputs.Reference_Latitude, \
                                    ui.User_Assumed_Inputs.Reference_Longitude)
ts2=time.time_ns()
print((ts2-ts1)/1e6,"- ZoneTime_SolarAngles(ms)")


ts1=time.time_ns()
AnAgGHI, AnAgDHI, AnAgDNI, AnAgSunHours, MaxTambSH, MinTambSH, AveTambSH, \
 MaxWSSH, MinWSSH, AveWSSH, AvAgGHIPerDay, AvAgDHIPerDay, AvAgDNIPerDay,\
 AvSunshineHoursPerDay, MonAgGHI, MonAgDHI, MonAgDNI,MonAgSunHours, MonMaxTamb,\
 MonMinTamb, MonAvTamb, MonMaxWS, MonMinWS, MonAvWS, DayAgGHI, DayAgDHI, \
 DayAgDNI, DaySunshineHours, DayMaxWS, DayMinWS, DayAvWS, DayMaxTamb, \
 DayMinTamb, DayAvTamb, GHI_SH, DHI_SH, DNI_SH, Tamb_SH, WS_SH \
            = Resource_Estimation(ui.User_Assumed_Inputs.GHI, \
                                   ui.User_Assumed_Inputs.DHI, \
                                   ui.User_Assumed_Inputs.DNI, \
                                   ui.User_Assumed_Inputs.WS, \
                                   ui.User_Assumed_Inputs.Tamb, \
                                   SunHourStatus, SunWindow, \
                                   ZTDayHour)
ts2=time.time_ns()
print((ts2-ts1)/1e6,"- Resource_Estimation(ms)")

# Calling Function that estimates spacing factor and generation hours
ts1=time.time_ns()
PgenHourStatus, PgenDayHour, LrowFactor, LcolFactor, PgenWindow, \
    AnnualPgenHours = Spacing_Factor(ui.User_Assumed_Inputs.Module_Tilt, \
                                       STAngSolAzimuth, STAngSolAltitude, \
                                       ZTDayHour)
ts2=time.time_ns()
print((ts2-ts1)/1e6,"- Spacing_Factor(ms)")

# Calling Function that estimates net effective radiation
ts1=time.time_ns()
Gt, n, AnnualGt = Net_Effective_Radiation (ui.User_Assumed_Inputs.Module_Tilt, \
                         ui.User_Assumed_Inputs.Misc_Array_height, \
                         STAngIncidence, GHI_SH, DHI_SH, DNI_SH, \
                         ui.User_Assumed_Inputs.Misc_Albedo, \
                         ui.User_Assumed_Inputs.Module_Lmod, \
                         ui.User_Assumed_Inputs.Misc_Ground_Clearance, \
                         PgenHourStatus)
ts2=time.time_ns()
print((ts2-ts1)/1e6,"- Net_Effective_Radiation(ms)")

# Calling Function that estimates cell temperature and resource to
# module power function
ts1=time.time_ns()
ModTCell, ModIsc, ModVoc, RPmod, AggRPmod = \
    Tcell_ResToModPower (PgenHourStatus, Gt, Tamb_SH, WS_SH, \
                         ui.User_Assumed_Inputs.Mount_a_CT, \
                         ui.User_Assumed_Inputs.Mount_b_CT, \
                         ui.User_Assumed_Inputs.Mount_DelT, \
                         ui.User_Assumed_Inputs.Module_Pmod, \
                         ui.User_Assumed_Inputs.Module_Kt_Pmax, \
                         ui.User_Assumed_Inputs.Module_Kt_Isc, \
                         ui.User_Assumed_Inputs.Module_Kt_Voc, \
                         ui.User_Assumed_Inputs.Module_Isc, \
                         ui.User_Assumed_Inputs.Module_Voc)
ts2=time.time_ns()
print((ts2-ts1)/1e6,"- Tcell_ResToModPower(ms)")

# Calling function that estimates the plant sizing parameters
# This function optimises the plant design to generate maximum power for the
# available resource
ts1=time.time_ns()
ModPmax, ModVocMin, ModVocMax, ModIscMax, N_PCU, m, m_change, NmodPCU, y, \
y_area, NPlant, PureModArea, PureModAreaAcres, DCScalingFactor, Pplant, \
PplantPCU, OriginalPplant, PmaxDeviationNewCap, P_PCU_DCMax = \
    Plant_Sizing_Parameters (ui.User_Assumed_Inputs.V_PCU_Choice, \
                         ui.User_Assumed_Inputs.Vuser, n, RPmod, ModIsc, ModVoc, \
                         ui.User_Assumed_Inputs.Pplant_Target, \
                         ui.User_Assumed_Inputs.PCU_PnomDC, \
                         ui.User_Assumed_Inputs.PCU_P_PCU_AC, \
                         ui.User_Assumed_Inputs.Module_Vmp, \
                         ui.User_Assumed_Inputs.Module_Imp, \
                         ui.User_Assumed_Inputs.Module_Pmod, \
                         ui.User_Assumed_Inputs.PCU_VmaxDC, \
                         ui.User_Assumed_Inputs.PCU_ImaxDC, \
                         ui.User_Assumed_Inputs.Misc_SL_PC, \
                         ui.User_Assumed_Inputs.Misc_EL_PC, \
                         ui.User_Assumed_Inputs.Module_Lmod, \
                         ui.User_Assumed_Inputs.Module_Bmod, \
                         ui.User_Assumed_Inputs.PCU_Vmppmax, \
                         ui.User_Assumed_Inputs.PCU_Vmppmin)
ts2=time.time_ns()
print((ts2-ts1)/1e3,"- Plant_Sizing_Parameters(microseconds)")

# Calling the function that estimates inter-row, inter-column spacing;
# aggregate resource to module power factor and annual generation hours
# for all 3 time windows
ts1=time.time_ns()
Drow, Dcol, AnnualGenHr, AggRPmod, Pgen2hrWindow, Pgen4hrWindow, GtPgen, \
GtPgen2hr, GtPgen4hr, AnGt = \
        Inter_Row_Col_Spacing (PgenWindow, SunHourStatus, PgenHourStatus, \
                           ZTDayHour, LrowFactor, LcolFactor, n, \
                           ui.User_Assumed_Inputs.Module_Lmod, RPmod, Gt)

ts2=time.time_ns()
print((ts2-ts1)/1e6,"- Inter_Row_Col_Spacing(ms)")

# Calling the function that estimates plant area, packing density and suitable
# time window
ts1=time.time_ns()
PlantArea, PackingDensity, DeviationFactor, TW, TW_DayHrs, Gt_TW = \
                Plant_Area_Estimation (n, ui.User_Assumed_Inputs.Module_Lmod, \
                       ui.User_Assumed_Inputs.Module_Bmod, \
                       ui.User_Assumed_Inputs.Module_Tilt, \
                       Drow, Dcol, m, y_area, N_PCU, \
                       ui.User_Assumed_Inputs.Misc_BS_AlongL_a, \
                       ui.User_Assumed_Inputs.Misc_BS_AlongB_b, \
                       Pplant, PureModAreaAcres, PgenWindow, Pgen2hrWindow,\
                       Pgen4hrWindow, GtPgen, GtPgen2hr, GtPgen4hr )

ts2=time.time_ns()
print((ts2-ts1)/1e3,"- Plant_Area_Estimation(microseconds)")

# Calling the function that estimates plant generation in year 0 and other
# performance metrics
ts1=time.time_ns()
HourlyGenPC, HourlySEEPV, HourlyGtallpmod, PplantPVAC, AggGtallmodulesTW, \
AnPplantDC, AnPplantDCincSL, AnPplantAC, AnSoilingLoss, AnElectricalLoss, \
AnPRPV, AnCUFPV,  AnSEEPV = Plant_output_year_0 (TW, AnGt, PureModArea, \
                                 ui.User_Assumed_Inputs.Misc_SL_PC, \
                                 ui.User_Assumed_Inputs.Misc_EL_PC, \
                                 TW_DayHrs, NPlant, ZTDayHour, RPmod, \
                                 ui.User_Assumed_Inputs.Module_Pmod, \
                                 ui.User_Assumed_Inputs.PCU_Efficiency,\
                                 Pplant, Gt_TW)

ts2=time.time_ns()
print((ts2-ts1)/1e6,"- Plant_output_year_0(ms)")

# Calling the function that estimates the PV plant output factoring
# Module degradation
ts1=time.time_ns()
DegradationRates, BaseDegCaseIndex, AggRPF, AnnualPVGenDeg, NetModuleRating, \
AnnualCUFPV, ModuleRatingPC, AnnualPRPV, AnnualSEEPV  = Degradation_Results\
                     (ui.User_Assumed_Inputs.Module_YoY_Degrdn_rate, \
                      ui.User_Assumed_Inputs.Plant_Life, \
                      ui.User_Assumed_Inputs.Module_Pmod, \
                      ui.User_Assumed_Inputs.Module_Rating_EoYr1, \
                      NPlant, AggRPmod,\
                      ui.User_Assumed_Inputs.PCU_Efficiency, \
                      ui.User_Assumed_Inputs.Misc_SL_PC, \
                      ui.User_Assumed_Inputs.Misc_EL_PC, Pplant, PureModArea, AnGt)
ts2=time.time_ns()
print((ts2-ts1)/1e6,"- Degradation_Results(ms)")

# Calling function which estimates the monthly and daily generation and
# generates histogram to measure frequency of generation
ts1=time.time_ns()
PVGenYr0Daily, PVGenYr0Monthly, GenHrsYr0TWPC, AgGenHrs1to100PCTW \
     = Monthly_Daily_Gen_Hist1 (PplantPVAC, ZTDayHour, PplantPCU, ZoneTime)
ts2=time.time_ns()
print((ts2-ts1)/1e6,"- Monthly_Daily_Gen_Hist1(ms)")

ts1=time.time_ns()
GenPVACTW24x365, GenPVACTW24x365T, MinGenPVACTW24x365, MaxGenPVACTW24x365, \
    MeanGenPVACTW24x365, SumGenPVACTW24x365, PCShareGenACTW24x365 \
    = Gen_vs_DayHr_Hist2 (PplantPVAC)
ts2=time.time_ns()
print((ts2-ts1)/1e3,"- Gen_vs_DayHr_Hist2(microseconds)")

# Calling the function estimating miscellaneous parameters in the tech model
ts1=time.time_ns()
AnPVGenPerAcre, AreaPgenPgen2hrFactor, AreaPgen2hrPgen4hrFactor, \
 AreaPgenPgen4hrFactor, GtPgenPgen2hrFactor, GtPgen2hrPgen4hrFactor, \
 GtPgenPgen4hrFactor, GenPgenPgen2hrFactor, GenPgen2hrPgen4hrFactor, \
 GenPgenPgen4hrFactor, GenPCPgen2hrPgen, GenPCPgen4hrPgen, \
 GenHrPgenPgen2hrFactor, GenHrPgen2hrPgen4hrFactor, GenHrPgenPgen4hrFactor, \
 GenPCTW24x365T, GtallmodTW24x365T, SEEPVTW24x365T \
 = Misc_Tech_Parameters (PlantArea, AnnualPVGenDeg, TW, \
                          BaseDegCaseIndex, AnGt, AnnualGenHr,  HourlyGenPC, \
                          HourlyGtallpmod, HourlySEEPV, Gt_TW)
ts2=time.time_ns()
print((ts2-ts1)/1e3,"- Misc_Tech_Parameters(microseconds)")
