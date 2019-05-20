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
print(ts2-ts1)


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
print(ts2-ts1)
ts1=time.time_ns()
# Calling Function that estimates spacing factor and generation hours
PgenHourStatus, PgenDayHour, LrowFactor, LcolFactor, PgenWindow, \
    AnnualPgenHours = Spacing_Factor(ui.User_Assumed_Inputs.Module_Tilt, \
                                       STAngSolAzimuth, STAngSolAltitude, \
                                       ZTDayHour)
ts2=time.time_ns()
print(ts2-ts1)
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
print(ts2-ts1)

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
print(ts2-ts1)

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
print(ts2-ts1)
