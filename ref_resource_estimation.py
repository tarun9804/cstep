#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 10:42:56 2019

@author: tarun
"""
import numpy as np

import user_input as ui
import time
import old_code as rf
# for importing data from excel


# for solar angles and time related calculations:
import math

# Function definition
def cosd(x):return math.cos(math.radians(x))
def sind(x):return math.sin(math.radians(x))
def tand(x):return math.tan(math.radians(x))
def acosd(x):return math.degrees(math.acos(x))
def asind(x):return math.degrees(math.asin(x))
def atand(x):return math.degrees(math.atan(x))


#------------------------------------------------------------------------------
# Determination of spacing factor
#------------------------------------------------------------------------------
#
# This module is only for fixed tilt configuration
# Seperate functions need to be developed for module tracking
#
def Spacing_Factor(Module_Tilt, Ang_Sol_Azimuth, Ang_Sol_Altitude, ZT_Day_Hour):
    Lrow_Factor = np.zeros(len(Ang_Sol_Azimuth))
    Lcol_Factor = np.zeros(len(Ang_Sol_Azimuth))
    Pgen_Hour_Status = np.zeros(len(Ang_Sol_Azimuth))
    Pgen_Day_Hour = np.zeros(len(Ang_Sol_Azimuth))

    for hour in range(len(Ang_Sol_Azimuth)):
        if Ang_Sol_Altitude [hour] >= 1:
            Pgen_Hour_Status [hour] = 1
            Pgen_Day_Hour [hour] = ZT_Day_Hour [hour]
            Lrow_Factor [hour] = sind (Module_Tilt) * \
                cosd (Ang_Sol_Azimuth[hour]) / tand (Ang_Sol_Altitude [hour])
            Lcol_Factor [hour] = sind (Module_Tilt) * \
                sind (Ang_Sol_Azimuth[hour]) / tand (Ang_Sol_Altitude [hour])


    Pgen_Window = (Pgen_Day_Hour[np.nonzero(Pgen_Day_Hour)].min(), \
                            Pgen_Day_Hour[np.nonzero(Pgen_Day_Hour)].max())
    Annual_Pgen_Hours = Pgen_Hour_Status.sum()

    return(Pgen_Hour_Status, Pgen_Day_Hour, Lrow_Factor, Lcol_Factor, \
           Pgen_Window, Annual_Pgen_Hours)

#-------------------------- End of Spacing factor function --------------------






ui.User_Assumed_Inputs()
ts1=time.time_ns()
ZoneTime, ZTDayHour,SolarTime, STAng_Hour_f, STAngSolAltitude, \
            STAngSolAzimuth, STAngIncidence, SunHourStatus,\
            AnnualSunHours,CorrectionHHMMSS, ZTTrise, ZTTset, STDayLength,\
            MaxCorrection, MinCorrection, MinTrise, MaxTrise, MinTset, MaxTset, \
            MaxDayLength, MinDayLength, SunWindow, SunWindowLength, MaxDLSRSS, \
            MinDLSRSS \
            = rf.ZoneTime_SolarAngles(ui.User_Assumed_Inputs.Location_Latitude,\
                                    ui.User_Assumed_Inputs.Location_Longitude,\
                                    ui.User_Assumed_Inputs.Module_Tilt, \
                                    ui.User_Assumed_Inputs.Module_Surface_Azimuth,\
                                    ui.User_Assumed_Inputs.Reference_Latitude, \
                                    ui.User_Assumed_Inputs.Reference_Longitude)


ts2=time.time_ns()
print(ts2-ts1)
ts1=time.time_ns()
# Calling Function that estimates solar resource at different temporal resolution
AnAgGHI, AnAgDHI, AnAgDNI, AnAgSunHours, MaxTambSH, MinTambSH, AveTambSH, \
 MaxWSSH, MinWSSH, AveWSSH, AvAgGHIPerDay, AvAgDHIPerDay, AvAgDNIPerDay,\
 AvSunshineHoursPerDay, MonAgGHI, MonAgDHI, MonAgDNI,MonAgSunHours, MonMaxTamb,\
 MonMinTamb, MonAvTamb, MonMaxWS, MonMinWS, MonAvWS, DayAgGHI, DayAgDHI, \
 DayAgDNI, DaySunshineHours, DayMaxWS, DayMinWS, DayAvWS, DayMaxTamb, \
 DayMinTamb, DayAvTamb, GHI_SH, DHI_SH, DNI_SH, Tamb_SH, WS_SH \
            = rf.Resource_Estimation(ui.User_Assumed_Inputs.GHI, \
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
