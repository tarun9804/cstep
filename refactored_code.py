#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:03:19 2019

@author: tarun
"""




import pandas as pd
import numpy as np
import time

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




def ZoneTime_SolarAngles_(Location_Latitude, Location_Longitude, Module_Tilt, \
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

def getRadiationAvgData(radiationRawData, sunHourStatus):
    radiationRawData = radiationRawData.values
    radiationDataSunHours = sunHourStatus * radiationRawData
    Day_Ag_radiation = np.sum(np.reshape(radiationDataSunHours*0.001,(365,24)),axis=1)
    Mon_Ag_radiation=np.sum(monthorder(Day_Ag_radiation),axis=1)
    An_Ag_radiation_SH = np.sum(Day_Ag_radiation)/1000
    Av_Ag_radiation_Per_Day = np.mean(Day_Ag_radiation)
    return (Day_Ag_radiation,Mon_Ag_radiation,An_Ag_radiation_SH,Av_Ag_radiation_Per_Day)


def getClimatDataReport(rawData,sunHourStatus):
    rawData = rawData
    rawDataSunHoursMasked = rawData * sunHourStatus

    Day_Sunshine_Hours = np.sum(np.reshape(sunHourStatus, (365, 24)), axis=1)
    Mon_Ag_Sun_Hours = np.sum(monthorder(Day_Sunshine_Hours), axis=1)

    Daily_Sum = np.sum(np.reshape(rawDataSunHoursMasked, (365, 24)), axis=1)
    Day_Max = np.max(np.reshape(rawDataSunHoursMasked,(365,24)),axis=1)
    Day_Min = np.min(np.reshape(rawDataSunHoursMasked,(365,24)),axis=1)
    Mon_Sum=np.sum(monthorder(Daily_Sum),axis=1)
    Mon_Max = np.max(monthorder(Day_Max), axis=1)
    Mon_Min = np.min(monthorder(Day_Min), axis=1)
    Day_Av = Daily_Sum/Day_Sunshine_Hours
    Mon_Av = Mon_Sum / Mon_Ag_Sun_Hours

    temp = np.nonzero(rawDataSunHoursMasked)
    Min = np.min(rawDataSunHoursMasked[temp[0]])
    Max = np.max(rawDataSunHoursMasked[temp[0]])
    Ave = np.mean(rawDataSunHoursMasked[temp[0]])

    return (Max, Min, Ave, Mon_Max, Mon_Min, Mon_Av, Day_Max, Day_Min,
            Day_Av, rawDataSunHoursMasked)

def getSunHourData(sunHourStatus):
    Day_Sunshine_Hours = np.sum(np.reshape(sunHourStatus, (365, 24)), axis=1)
    Mon_Ag_Sun_Hours = np.sum(monthorder(Day_Sunshine_Hours), axis=1)
    An_Ag_SunHours = np.sum(Day_Sunshine_Hours)
    Av_Sunshine_Hours_Per_Day = np.sum(Day_Sunshine_Hours) / 365
    return (Day_Sunshine_Hours,Mon_Ag_Sun_Hours,An_Ag_SunHours,Av_Sunshine_Hours_Per_Day)


# ui.User_Assumed_Inputs()
#
# ts1=time.time_ns()
#
# ZoneTime, ZTDayHour,SolarTime, STAng_Hour_f, STAngSolAltitude, \
#             STAngSolAzimuth, STAngIncidence, SunHourStatus,\
#             AnnualSunHours,CorrectionHHMMSS, ZTTrise, ZTTset, STDayLength,\
#             MaxCorrection, MinCorrection, MinTrise, MaxTrise, MinTset, MaxTset, \
#             MaxDayLength, MinDayLength, SunWindow, SunWindowLength, MaxDLSRSS, \
#             MinDLSRSS \
#             = ZoneTime_SolarAngles_(ui.User_Assumed_Inputs.Location_Latitude,\
#                                     ui.User_Assumed_Inputs.Location_Longitude,\
#                                     ui.User_Assumed_Inputs.Module_Tilt, \
#                                     ui.User_Assumed_Inputs.Module_Surface_Azimuth,\
#                                     ui.User_Assumed_Inputs.Reference_Latitude, \
#                                     ui.User_Assumed_Inputs.Reference_Longitude)
# ts2=time.time_ns()
# print(ts2-ts1)
#
#
# ts1=time.time_ns()
# AnAgGHI, AnAgDHI, AnAgDNI, AnAgSunHours, MaxTambSH, MinTambSH, AveTambSH, \
#  MaxWSSH, MinWSSH, AveWSSH, AvAgGHIPerDay, AvAgDHIPerDay, AvAgDNIPerDay,\
#  AvSunshineHoursPerDay, MonAgGHI, MonAgDHI, MonAgDNI,MonAgSunHours, MonMaxTamb,\
#  MonMinTamb, MonAvTamb, MonMaxWS, MonMinWS, MonAvWS, DayAgGHI, DayAgDHI, \
#  DayAgDNI, DaySunshineHours, DayMaxWS, DayMinWS, DayAvWS, DayMaxTamb, \
#  DayMinTamb, DayAvTamb, GHI_SH, DHI_SH, DNI_SH, Tamb_SH, WS_SH \
#             = Resource_Estimation_ (ui.User_Assumed_Inputs.GHI, \
#                                    ui.User_Assumed_Inputs.DHI, \
#                                    ui.User_Assumed_Inputs.DNI, \
#                                    ui.User_Assumed_Inputs.WS, \
#                                    ui.User_Assumed_Inputs.Tamb, \
#                                    SunHourStatus, SunWindow, \
#                                    ZTDayHour)
# ts2=time.time_ns()
# print(ts2-ts1)
