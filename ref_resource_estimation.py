#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 10:42:56 2019

@author: tarun
"""
import numpy as np

import user_input as ui
import time
import ref_file as rf
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
# Estimating solar resource
#------------------------------------------------------------------------------

def Resource_Estimation (GHI, DHI, DNI, WS, Tamb, Sun_Hour_Status, Sun_Window, \
                        ZT_Day_Hour):

    GHI = GHI.values
    DHI = DHI.values
    DNI = DNI.values
    Tamb = Tamb.values
    WS = WS.values

    GHI_SH = np.zeros(len(GHI))
    DHI_SH = np.zeros(len(DHI))
    DNI_SH = np.zeros(len(DNI))
    Tamb_SH = np.zeros(len(Tamb))
    WS_SH = np.zeros(len(WS))

    # Annual Aggregate Parameters for base dataset in MW/sq.m
    An_Ag_GHI_base = (np.sum(GHI)/1000000)
    An_Ag_DHI_base = (np.sum(DHI)/1000000)
    An_Ag_DNI_base = (np.sum(DNI)/1000000)

    # Defining daily parameters during sunshine hours
    Day_Ag_GHI = np.zeros(365)
    Day_Ag_DHI = np.zeros(365)
    Day_Ag_DNI = np.zeros(365)
    Daily_Sum_Tamb = np.zeros(365)
    Daily_Sum_WS = np.zeros(365)
    Day_Av_WS = np.zeros(365)
    Day_Av_Tamb = np.zeros(365)
    Day_Sunshine_Hours = np.zeros(365)
    Day_Min_Tamb = np.full(365, 999.0)
    Day_Max_Tamb = np.zeros(365)
    Day_Min_WS = np.full(365, 999.0)
    Day_Max_WS = np.zeros(365)

    # Defining monthly parameters during sunshine hours
    Mon_Ag_GHI = np.zeros(12)
    Mon_Ag_DHI = np.zeros(12)
    Mon_Ag_DNI = np.zeros(12)
    Mon_Ag_Sun_Hours = np.zeros(12)
    Mon_Sum_Tamb = np.zeros(12)
    Mon_Sum_WS = np.zeros(12)
    Mon_Av_WS = np.zeros(12)
    Mon_Av_Tamb = np.zeros(12)
    Mon_Max_Tamb = np.zeros(12)
    Mon_Min_Tamb = np.full(12, 999.0)
    Mon_Max_WS = np.zeros(12)
    Mon_Min_WS = np.full(12, 999.0)

    Day = 0
    Month = 0
    for i in range(len(Sun_Hour_Status)):
        Day = (i//24)
        if Sun_Hour_Status [i] == 1:

            # Slicing solar data for sun-hours only
            GHI_SH [i] = GHI [i]
            DHI_SH [i] = DHI [i]
            DNI_SH [i] = DNI [i]
            Tamb_SH [i] = Tamb [i]
            WS_SH [i] = WS [i]

            # Estimating daily parameters (to get in kW/sq.m, sol data x 0.001)
            Daily_Sum_Tamb [Day] = Daily_Sum_Tamb [Day] + Tamb [i]
            Day_Ag_GHI [Day] = Day_Ag_GHI [Day] + GHI [i] * 0.001
            Day_Ag_DHI [Day] = Day_Ag_DHI [Day] + DHI [i] * 0.001
            Day_Ag_DNI [Day] = Day_Ag_DNI [Day] + DNI [i] * 0.001
            Daily_Sum_WS [Day] = Daily_Sum_WS [Day] + WS [i]
            Day_Sunshine_Hours [Day] = Day_Sunshine_Hours [Day] + 1

            if Day_Max_WS [Day] < WS [i]:
                Day_Max_WS [Day] = WS [i]
            if Day_Min_WS [Day] > WS [i]:
                Day_Min_WS [Day] = WS [i]
            if Day_Max_Tamb [Day]< Tamb [i]:
                Day_Max_Tamb [Day] = Tamb [i]
            if Day_Min_Tamb [Day] > Tamb [i]:
                Day_Min_Tamb [Day] = Tamb [i]

            #Estimating monthly parameters
            if Day <= 30:
                Month = 0
            if Day > 30 and Day <= 58:
                Month = 1
            if Day > 58 and Day <= 89:
                Month = 2
            if Day > 89 and Day <= 119:
                Month = 3
            if Day > 119 and Day <= 150:
                Month = 4
            if Day > 150 and Day <= 180:
                Month = 5
            if Day > 180 and Day <= 211:
                Month = 6
            if Day > 211 and Day <= 242:
                Month = 7
            if Day > 242 and Day <= 272:
                Month = 8
            if Day > 272 and Day <= 303:
                Month = 9
            if Day > 303 and Day <= 333:
                Month = 10
            if Day > 333 and Day <= 364:
                Month = 11

            Mon_Ag_GHI [Month] = Mon_Ag_GHI [Month] + GHI [i] * 0.001
            Mon_Ag_DNI [Month] = Mon_Ag_DNI [Month] + DNI [i] * 0.001
            Mon_Ag_DHI [Month] = Mon_Ag_DHI [Month] + DHI [i] * 0.001
            Mon_Ag_Sun_Hours [Month] = Mon_Ag_Sun_Hours [Month] + 1
            Mon_Sum_Tamb [Month] = Mon_Sum_Tamb [Month] + Tamb [i]
            Mon_Sum_WS [Month] = Mon_Sum_WS [Month] + WS [i]
            if Mon_Max_Tamb [Month] < Tamb [i]:
                    Mon_Max_Tamb [Month] = Tamb [i]
            if Mon_Min_Tamb [Month] > Tamb [i]:
                    Mon_Min_Tamb [Month] = Tamb [i]
            if Mon_Max_WS [Month] < WS [i]:
                    Mon_Max_WS [Month] = WS [i]
            if Mon_Min_WS [Month] > WS [i]:
                    Mon_Min_WS [Month] = WS [i]

    for i in range(len(Day_Sunshine_Hours)):
        Day_Av_Tamb [i] = Daily_Sum_Tamb [i]/Day_Sunshine_Hours[i]
        Day_Av_WS [i]= Daily_Sum_WS [i]/Day_Sunshine_Hours [i]

    for i in range(len(Mon_Ag_Sun_Hours)):
        Mon_Av_Tamb [i] = Mon_Sum_Tamb [i]/Mon_Ag_Sun_Hours[i]
        Mon_Av_WS [i] = Mon_Sum_WS [i]/Mon_Ag_Sun_Hours [i]

    # Annual Aggregate radiation parameters for sunshine hours in MW/sq.m
    An_Ag_GHI_SH = (np.sum(Day_Ag_GHI)/1000)
    An_Ag_DHI_SH = (np.sum(Day_Ag_DHI)/1000)
    An_Ag_DNI_SH = (np.sum(Day_Ag_DNI)/1000)
    An_Ag_SunHours = Day_Sunshine_Hours.sum()

    # Annual Min and Max parameters
    Min_Tamb_SH = Tamb_SH[np.nonzero(Tamb_SH)].min()
    Max_Tamb_SH = Tamb_SH[np.nonzero(Tamb_SH)].max()
    Ave_Tamb_SH = Tamb_SH[np.nonzero(Tamb_SH)].mean()

    Min_WS_SH = WS_SH.min()
    Max_WS_SH = WS_SH.max()
    Ave_WS_SH = Daily_Sum_WS.sum()/Day_Sunshine_Hours.sum()

    # Per Day average
    Av_Ag_GHI_Per_Day = Day_Ag_GHI.mean()
    Av_Ag_DHI_Per_Day = Day_Ag_DHI.mean()
    Av_Ag_DNI_Per_Day = Day_Ag_DNI.mean()
    Av_Sunshine_Hours_Per_Day = Day_Sunshine_Hours.sum()/365

    return (An_Ag_GHI_SH, An_Ag_DHI_SH, An_Ag_DNI_SH, An_Ag_SunHours, \
            Max_Tamb_SH, Min_Tamb_SH, Ave_Tamb_SH, Max_WS_SH, Min_WS_SH, \
            Ave_WS_SH, Av_Ag_GHI_Per_Day, Av_Ag_DHI_Per_Day, Av_Ag_DNI_Per_Day,\
            Av_Sunshine_Hours_Per_Day, Mon_Ag_GHI, Mon_Ag_DHI, Mon_Ag_DNI, \
            Mon_Ag_Sun_Hours, Mon_Max_Tamb, Mon_Min_Tamb, Mon_Av_Tamb, \
            Mon_Max_WS, Mon_Min_WS, Mon_Av_WS, Day_Ag_GHI, Day_Ag_DHI, \
            Day_Ag_DNI, Day_Sunshine_Hours, Day_Max_WS, Day_Min_WS, \
            Day_Av_WS, Day_Max_Tamb, Day_Min_Tamb, Day_Av_Tamb, GHI_SH, \
            DHI_SH, DNI_SH, Tamb_SH, WS_SH)

# -------------------- End of resource estimation function --------------------

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
            = Resource_Estimation (ui.User_Assumed_Inputs.GHI, \
                                   ui.User_Assumed_Inputs.DHI, \
                                   ui.User_Assumed_Inputs.DNI, \
                                   ui.User_Assumed_Inputs.WS, \
                                   ui.User_Assumed_Inputs.Tamb, \
                                   SunHourStatus, SunWindow, \
                                   ZTDayHour)

ts2=time.time_ns()
print(ts2-ts1)
