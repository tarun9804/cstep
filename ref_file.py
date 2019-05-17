#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CSTEP's Solar Techno-Economic Model for Photovoltaics (CSTEM PV)
Technical model
Created on Wed Apr 17 22:04 2019
Updated on Tue Apr 30 10:02 2019
@author: Harshid Sridhar
"""
import numpy as np
import time as ti

# for importing data from excel
import pandas as pd

# for solar angles and time related calculations:
import math
from datetime import datetime
from datetime import timedelta
from datetime import date
from datetime import time
# Function definition
def cosd(x):return math.cos(math.radians(x))
def sind(x):return math.sin(math.radians(x))
def tand(x):return math.tan(math.radians(x))
def acosd(x):return math.degrees(math.acos(x))
def asind(x):return math.degrees(math.asin(x))
def atand(x):return math.degrees(math.atan(x))

#------------------------------------------------------------------------------
# Reading the Input data
#------------------------------------------------------------------------------

def User_Assumed_Inputs():
    
    # Importing data from base input file
    User_Inp = pd.read_excel ('Model_Inputs.xlsx', 'Base', index_col = None)
    Solar_Data = pd.read_excel('Model_Inputs.xlsx','Solar Data', index_col = 0)
    
    
    # Mapping inputs
    # Base Data Inputs
    User_Assumed_Inputs.Location_Name = User_Inp ['Data'][0]
    User_Assumed_Inputs.Location_Latitude = User_Inp ['Data'][1]
    User_Assumed_Inputs.Location_Longitude = User_Inp ['Data'][2]
    User_Assumed_Inputs.Reference_Latitude = 25.15      # Assumed
    User_Assumed_Inputs.Reference_Longitude = -82.58    # Assumed
    User_Assumed_Inputs.Pplant_Target = User_Inp ['Data'][3]
    User_Assumed_Inputs.Plant_Life = User_Inp ['Data'][4]
    User_Assumed_Inputs.Module_Tilt = User_Inp ['Data'][5]
    User_Assumed_Inputs.Module_Surface_Azimuth = User_Inp ['Data'][6]
    
    # Module Data
    User_Assumed_Inputs.Module_Manufacturer = User_Inp ['Data'][7]
    User_Assumed_Inputs.Module_Technology = User_Inp ['Data'][8]
    User_Assumed_Inputs.Module_Model = User_Inp ['Data'][9]
    User_Assumed_Inputs.Module_Pmod = User_Inp ['Data'][10]
    User_Assumed_Inputs.Module_Voc = User_Inp ['Data'][11]
    User_Assumed_Inputs.Module_Isc = User_Inp ['Data'][12]
    User_Assumed_Inputs.Module_Vmp = User_Inp ['Data'][13]
    User_Assumed_Inputs.Module_Imp = User_Inp ['Data'][14]
    User_Assumed_Inputs.Module_Lmod = User_Inp ['Data'][15]
    User_Assumed_Inputs.Module_Bmod = User_Inp ['Data'][16]
    User_Assumed_Inputs.Module_Kt_Pmax = User_Inp ['Data'][17]
    User_Assumed_Inputs.Module_Kt_Voc = User_Inp ['Data'][18]
    User_Assumed_Inputs.Module_Kt_Isc = User_Inp ['Data'][19]
    User_Assumed_Inputs.Module_Rating_EoYr1 = User_Inp ['Data'][20]
    User_Assumed_Inputs.Module_YoY_Degrdn_rate = User_Inp ['Data'][21]
    
    # PCU Data
    User_Assumed_Inputs.PCU_Manufacturer = User_Inp ['Data'][22]
    User_Assumed_Inputs.PCU_Model = User_Inp ['Data'][23]
    User_Assumed_Inputs.PCU_Efficiency = User_Inp ['Data'][24]
    User_Assumed_Inputs.PCU_P_PCU_AC = User_Inp ['Data'][25]
    User_Assumed_Inputs.PCU_PnomDC = User_Inp ['Data'][26]
    User_Assumed_Inputs.PCU_Vmppmin = User_Inp ['Data'][27]
    User_Assumed_Inputs.PCU_Vmppmax = User_Inp ['Data'][28]
    User_Assumed_Inputs.PCU_Vstart = User_Inp ['Data'][29]
    User_Assumed_Inputs.PCU_VmaxDC = User_Inp ['Data'][30]
    User_Assumed_Inputs.PCU_InomDC = User_Inp ['Data'][31]
    User_Assumed_Inputs.PCU_ImaxDC = User_Inp ['Data'][32]
    
    # Mount Data
    User_Assumed_Inputs.Mount_Module_Type = User_Inp ['Data'][33]
    User_Assumed_Inputs.Mount_Style = User_Inp ['Data'][34]
    User_Assumed_Inputs.Mount_a_CT = User_Inp ['Data'][35]
    User_Assumed_Inputs.Mount_b_CT = User_Inp ['Data'][36]
    User_Assumed_Inputs.Mount_DelT = User_Inp ['Data'][37]
    
    # Miscellaneous Data
    User_Assumed_Inputs.Misc_Albedo = User_Inp ['Data'][38]
    User_Assumed_Inputs.Misc_Array_height = User_Inp ['Data'][39]
    User_Assumed_Inputs.Misc_Ground_Clearance = User_Inp ['Data'][40]
    User_Assumed_Inputs.Misc_BS_AlongL_a = User_Inp ['Data'][41]
    User_Assumed_Inputs.Misc_BS_AlongB_b = User_Inp ['Data'][42]
    User_Assumed_Inputs.Misc_SL_PC = User_Inp ['Data'][43]
    User_Assumed_Inputs.Misc_EL_PC = User_Inp ['Data'][44]
    User_Assumed_Inputs.Misc_Aux_MWhR = User_Inp ['Data'][45]
    
    # Solar Data
    User_Assumed_Inputs.GHI = Solar_Data ['GHI']
    User_Assumed_Inputs.DHI = Solar_Data ['DHI']
    User_Assumed_Inputs.DNI = Solar_Data ['DNI']
    User_Assumed_Inputs.Tamb = Solar_Data ['Tamb']
    User_Assumed_Inputs.WS = Solar_Data ['WS']
    
    # User inputs
    User_Assumed_Inputs.V_PCU_Choice = 1
    User_Assumed_Inputs.Vuser = 0
    
#-----------------------------------------------------------------------------
# Estimation of Solar angles and zone time related parameters
#------------------------------------------------------------------------------

"""
Solar angles, Sun hours
--------------------------------
The following code considers the zone time (ZT) reference data & estimates the
solar time (ST) equivalent for a location with a given latitude and longitude

Further it estimates the:
    Declination angle (ST)
    Hour angle (ST)
    Sun rise and sun set hour angles (ST)
    Sun rise and sun set time (ST and ZT)
    Day length in hours (ST)
    Solar Zenith angle (ST)
    Solar Altitude angle (ST)
    Solar Azimuth angle (ST)
    Sun hour status (ST)
    Annual Sun hours
    
Considering the Module tilt angle and surface azimuth angle of the array,
it also estimates the 'Incidence angle (ST)' of the sun rays on the module
    
"""
                      
# Defining time conversion function 
def Conv2Time (Var):
    Var_Hour = int(Var)
    Var_Min = int((Var - Var_Hour)*60)
    Var_Sec = int(((Var - Var_Hour)*60 - Var_Min)*60)
    return(Var_Hour, Var_Min, Var_Sec)

def ZoneTime_SolarAngles (Location_Latitude, Location_Longitude, Module_Tilt, \
                          Module_Surface_Azimuth, Ref_Latitude, Ref_Longitude):
    
    # Declarations:
    # =============
    Hour = list()
    ZT_Day_Hour = list()
    ZT_Day = list()
    ZT_B = list()
    Zone_Time = list()
    EOT = list()
    Correction = list()
    Correction_HHMMSS = list()
    Solar_Time = list()
    Sun_Hour_Status = list()
    Annual_Sun_Hours = 0
    ST_Day = list()
    ST_B = list()
    ST_Ang_Declination = list()
    ST_Day_Hour_f = list()
    ST_Ang_Hour_f = list()
    ST_Ang_Sol_Zenith = list()
    ST_Ang_Sol_Altitude = list()
    ST_Ang_Sol_Azimuth = list()
    ST_Ang_Incidence = list()
    ST_Wset = list()
    ST_Wrise = list()
    Lis_Trise = list()
    Lis_Tset = list()
    Day_Length_Tset_Trise = list()
    ST_Tset = list()
    ST_Trise = list()
    ZT_Tset = list()
    ZT_Trise = list()
    ZT_Date = list()
    ST_Day_Length = list()
    
    EOLD = 4 * (Ref_Longitude - Location_Longitude)
    ii = 0
    # Defining Zone Time and Relevant angles
    for i in range(8760):
        Hour.append(i)
        ZT_Day_Hour.append(Hour[i] % 24)
        ZT_Day.append((Hour[i]//24)+1)
        ZT_B.append((ZT_Day[i]-1)*360/365)    
        if i == 0: 
            # Syntax: datetime.datetime(year, month, day, hour=0, minute=0, 
            # second=0, microsecond=0, tzinfo=None, *, fold=0)
            Zone_Time.append(datetime (2015,1,1,00,00,00))        
        else:
            # Syntax: datetime.timedelta(days=0, seconds=0, microseconds=0, 
            # milliseconds=0, minutes=0, hours=0, weeks=0)
            Zone_Time.append(Zone_Time[i-1] + timedelta(0,0,0,0,0,1,0))
        EOT.append(229.2 * (0.000075 + 0.001868 * cosd(ZT_B[i]) - 0.032077 * sind(ZT_B[i]) - 0.014615 * cosd(2*ZT_B[i]) - 0.040849 * sind(2*ZT_B[i])))
        Correction.append(EOLD + EOT[i])
        Correction_Hour, Correction_Min, Correction_Sec = Conv2Time (Correction[i]/60)
        Solar_Time.append(Zone_Time[i] + timedelta(0,Correction_Sec,0,0, Correction_Min,Correction_Hour,0))
        ST_Day_Hour_f.append(Solar_Time[i].hour + Solar_Time[i].minute/60 + Solar_Time[i].second/3600)
        
        # timetuple().tm_yday --> extracts the julian day 1 to 365 (or 366) from datetime object
        ST_Day.append(Solar_Time[i].timetuple().tm_yday)
        ST_B.append((ST_Day[i]-1)*360/365) 
        ST_Ang_Declination.append((180/math.pi)*(0.006918 - 0.399912 * cosd(ST_B[i]) + 0.070257 * sind(ST_B[i]) - 0.006758 * cosd(2*ST_B[i]) + 0.000907 * sind(2*ST_B[i]) - 0.002697 * cosd(3*ST_B[i]) + 0.00148 * sind(3 * ST_B[i]))  )
        ST_Ang_Hour_f.append(-15*(12-ST_Day_Hour_f[i]))
        
        ST_Ang_Sol_Zenith.append(acosd(cosd(Location_Latitude)*cosd(ST_Ang_Hour_f[i])*cosd(ST_Ang_Declination[i]) + sind(Location_Latitude)*sind(ST_Ang_Declination[i])))
        ST_Ang_Sol_Altitude.append(90 - ST_Ang_Sol_Zenith[i])
        if ST_Ang_Sol_Altitude[i] >= 0:
            Sun_Hour_Status.append(1)
            Annual_Sun_Hours = Annual_Sun_Hours + 1
        else:
            Sun_Hour_Status.append(0)
        try:
            ST_Ang_Sol_Azimuth.append((np.sign(ST_Ang_Hour_f[i]))*abs(acosd((cosd(ST_Ang_Sol_Zenith[i])*sind(Location_Latitude) - sind(ST_Ang_Declination[i]))/(sind(ST_Ang_Sol_Zenith[i])*cosd(Location_Latitude)))))
        except:
            ST_Ang_Sol_Azimuth.append(0)
            
        ST_Ang_Incidence.append(acosd(sind(ST_Ang_Declination[i])*sind(Location_Latitude)*cosd(Module_Tilt) -\
                                      sind(ST_Ang_Declination[i])*cosd(Location_Latitude)*sind(Module_Tilt)*cosd(Module_Surface_Azimuth)+\
                                      cosd(ST_Ang_Declination[i])*cosd(Location_Latitude)*cosd(Module_Tilt)*cosd(ST_Ang_Hour_f[i])+\
                                      cosd(ST_Ang_Declination[i])*sind(Location_Latitude)*sind(Module_Tilt)*cosd(Module_Surface_Azimuth)*cosd(ST_Ang_Hour_f[i])+\
                                      cosd(ST_Ang_Declination[i])*sind(Module_Tilt)*sind(Module_Surface_Azimuth)*sind(ST_Ang_Hour_f[i])))
        if ZT_Day_Hour[i] == 12:
            Correction_HHMMSS.append((Correction_Hour, Correction_Min, Correction_Sec))
            ST_Wset.append(acosd(-tand(Location_Latitude)*tand(ST_Ang_Declination[i])))
            ST_Wrise.append(-ST_Wset[ii])
            T_Set = 12 + ST_Wset[ii]/15
            Lis_Tset.append(T_Set - Correction[i]/60)
            T_Set_Hour, T_Set_Min, T_Set_Sec = Conv2Time(T_Set)
            ST_Tset.append(datetime(Solar_Time[i].year, Solar_Time[i].month, Solar_Time[i].day, T_Set_Hour, T_Set_Min, T_Set_Sec))
            T_Rise = 12 + ST_Wrise[ii]/15
            Lis_Trise.append(T_Rise - Correction[i]/60)
            T_Rise_Hour, T_Rise_Min, T_Rise_Sec = Conv2Time(T_Rise)
            ST_Trise.append(datetime(Solar_Time[i].year, Solar_Time[i].month, Solar_Time[i].day, T_Rise_Hour, T_Rise_Min, T_Rise_Sec))
            ZT_Trise.append(ST_Trise[ii] - timedelta(0,Correction_Sec,0,0,Correction_Min,Correction_Hour,0))  
            ZT_Tset.append(ST_Tset[ii] - timedelta(0,Correction_Sec,0,0,Correction_Min,Correction_Hour,0))  
            ZT_Date.append(date(Solar_Time[i].year, Solar_Time[i].month, Solar_Time[i].day))
            Day_Length = ST_Wset[ii]*2/15
            Day_Length_Hour, Day_Length_Min, Day_Length_Sec = Conv2Time (Day_Length)
            ST_Day_Length.append(datetime(Solar_Time[i].year, Solar_Time[i].month, Solar_Time[i].day, Day_Length_Hour,Day_Length_Min,Day_Length_Sec,0))
            Day_Length_Tset_Trise.append(Lis_Tset[ii] - Lis_Trise[ii])
            
            ii = ii + 1
   

    Max_Correction = (ZT_Date[Correction_HHMMSS.index(max(Correction_HHMMSS))],max(Correction_HHMMSS))
    Min_Correction = (ZT_Date[Correction_HHMMSS.index(min(Correction_HHMMSS))],min(Correction_HHMMSS))
    
    Buff_Hour, Buff_Min, Buff_Sec = Conv2Time (max(Lis_Trise))
    Max_Trise = (ZT_Date[Lis_Trise.index(max(Lis_Trise))],time(Buff_Hour, Buff_Min, Buff_Sec))
    
    Buff_Hour, Buff_Min, Buff_Sec = Conv2Time (max(Lis_Tset))
    Max_Tset = (ZT_Date[Lis_Tset.index(max(Lis_Tset))],time(Buff_Hour, Buff_Min, Buff_Sec))
    
    Buff_Hour, Buff_Min, Buff_Sec = Conv2Time (min(Lis_Trise))
    Min_Trise = (ZT_Date[Lis_Trise.index(min(Lis_Trise))],time(Buff_Hour, Buff_Min, Buff_Sec))
    
    Buff_Hour, Buff_Min, Buff_Sec = Conv2Time (min(Lis_Tset))
    Min_Tset = (ZT_Date[Lis_Tset.index(min(Lis_Tset))],time(Buff_Hour, Buff_Min, Buff_Sec))
    
    Buff_Hour, Buff_Min, Buff_Sec = Conv2Time (max(Day_Length_Tset_Trise))
    Max_DayLength = \
        (ZT_Date[Day_Length_Tset_Trise.index(max(Day_Length_Tset_Trise))], \
                 time(Buff_Hour, Buff_Min, Buff_Sec))
    
#--------------------------
# Sun rise and sun set time for max and min daylength added on 30 Apr 2019 
    Buff_Ind = Day_Length_Tset_Trise.index(max(Day_Length_Tset_Trise))
    
    Max_DL_SRSS  = (time(ZT_Trise[Buff_Ind].hour, ZT_Trise[Buff_Ind].minute, 0), \
                    time(ZT_Tset[Buff_Ind].hour, ZT_Tset[Buff_Ind].minute, 0))
    
    
    Buff_Hour, Buff_Min, Buff_Sec = Conv2Time (min(Day_Length_Tset_Trise))
    
    Min_DayLength = \
        (ZT_Date[Day_Length_Tset_Trise.index(min(Day_Length_Tset_Trise))],\
                 time(Buff_Hour, Buff_Min, Buff_Sec))
        
    Buff_Ind = Day_Length_Tset_Trise.index(min(Day_Length_Tset_Trise))
    
    Min_DL_SRSS = (time(ZT_Trise[Buff_Ind].hour, ZT_Trise[Buff_Ind].minute, ZT_Trise[Buff_Ind].second), \
                time(ZT_Tset[Buff_Ind].hour, ZT_Tset[Buff_Ind].minute, ZT_Tset[Buff_Ind].second))
        

#------------------------    
    Sun_Window = (math.floor(min(Lis_Trise)), math.ceil(max(Lis_Tset)))
    
    Buff_Hour, Buff_Min, Buff_Sec = Conv2Time (Sun_Window[1] - Sun_Window[0])
    Sun_Window_Length = (time(Buff_Hour, Buff_Min, Buff_Sec))

    return (Zone_Time, ZT_Day_Hour,Solar_Time, ST_Ang_Hour_f, ST_Ang_Sol_Altitude, \
            ST_Ang_Sol_Azimuth, ST_Ang_Incidence, Sun_Hour_Status,\
            Annual_Sun_Hours,Correction_HHMMSS, ZT_Trise, ZT_Tset, ST_Day_Length,\
            Max_Correction, Min_Correction, Min_Trise, Max_Trise, Min_Tset, Max_Tset, \
            Max_DayLength, Min_DayLength, Sun_Window, Sun_Window_Length, \
            Max_DL_SRSS, Min_DL_SRSS)

# -------------------- End of functions related to time and solar angles ------
    
#Obtaining user inputs
'''
User_Assumed_Inputs()

ts1=ti.time_ns()
# Calling Function related to ZoneTime and Solar Angles
#for i in range(50):
ZoneTime, ZTDayHour,SolarTime, STAng_Hour_f, STAngSolAltitude, \
            STAngSolAzimuth, STAngIncidence, SunHourStatus,\
            AnnualSunHours,CorrectionHHMMSS, ZTTrise, ZTTset, STDayLength,\
            MaxCorrection, MinCorrection, MinTrise, MaxTrise, MinTset, MaxTset, \
            MaxDayLength, MinDayLength, SunWindow, SunWindowLength, MaxDLSRSS, \
            MinDLSRSS \
            = ZoneTime_SolarAngles (User_Assumed_Inputs.Location_Latitude,\
                                    User_Assumed_Inputs.Location_Longitude,\
                                    User_Assumed_Inputs.Module_Tilt, \
                                    User_Assumed_Inputs.Module_Surface_Azimuth,\
                                    User_Assumed_Inputs.Reference_Latitude, \
                                    User_Assumed_Inputs.Reference_Longitude)
ts2=ti.time_ns()
print(ts2-ts1)
'''
'''

ts1=ti.time_ns()

# Calling Function related to ZoneTime and Solar Angles
for i in range(10):
    ZoneTime, ZTDayHour,SolarTime, STAng_Hour_f, STAngSolAltitude, \
                STAngSolAzimuth, STAngIncidence, SunHourStatus,\
                AnnualSunHours,CorrectionHHMMSS, ZTTrise, ZTTset, STDayLength,\
                MaxCorrection, MinCorrection, MinTrise, MaxTrise, MinTset, MaxTset, \
                MaxDayLength, MinDayLength, SunWindow, SunWindowLength, MaxDLSRSS, \
                MinDLSRSS \
                = ZoneTime_SolarAngles (User_Assumed_Inputs.Location_Latitude,\
                                        User_Assumed_Inputs.Location_Longitude,\
                                        User_Assumed_Inputs.Module_Tilt, \
                                        User_Assumed_Inputs.Module_Surface_Azimuth,\
                                        User_Assumed_Inputs.Reference_Latitude, \
                                        User_Assumed_Inputs.Reference_Longitude)
ts2=ti.time_ns()
print(ts2-ts1)
'''