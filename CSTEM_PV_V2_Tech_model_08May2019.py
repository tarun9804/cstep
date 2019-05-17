"""
CSTEP's Solar Techno-Economic Model for Photovoltaics (CSTEM PV)
Technical model
Created on Wed Apr 17 22:04 2019
Updated on Tue Apr 30 10:02 2019
@author: Harshid Sridhar
"""
import numpy as np

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
User_Assumed_Inputs()

# Calling Function related to ZoneTime and Solar Angles
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


# Calling Function that estimates solar resource at different temporal resolution
AnAgGHI, AnAgDHI, AnAgDNI, AnAgSunHours, MaxTambSH, MinTambSH, AveTambSH, \
 MaxWSSH, MinWSSH, AveWSSH, AvAgGHIPerDay, AvAgDHIPerDay, AvAgDNIPerDay,\
 AvSunshineHoursPerDay, MonAgGHI, MonAgDHI, MonAgDNI,MonAgSunHours, MonMaxTamb,\
 MonMinTamb, MonAvTamb, MonMaxWS, MonMinWS, MonAvWS, DayAgGHI, DayAgDHI, \
 DayAgDNI, DaySunshineHours, DayMaxWS, DayMinWS, DayAvWS, DayMaxTamb, \
 DayMinTamb, DayAvTamb, GHI_SH, DHI_SH, DNI_SH, Tamb_SH, WS_SH \
            = Resource_Estimation (User_Assumed_Inputs.GHI, \
                                   User_Assumed_Inputs.DHI, \
                                   User_Assumed_Inputs.DNI, \
                                   User_Assumed_Inputs.WS, \
                                   User_Assumed_Inputs.Tamb, \
                                   SunHourStatus, SunWindow, \
                                   ZTDayHour)
            

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

# Calling Function that estimates spacing factor and generation hours
PgenHourStatus, PgenDayHour, LrowFactor, LcolFactor, PgenWindow, \
    AnnualPgenHours = Spacing_Factor(User_Assumed_Inputs.Module_Tilt, \
                                       STAngSolAzimuth, STAngSolAltitude, \
                                       ZTDayHour)

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
    
    Rb = np.zeros(len(Pgen_Hour_Status))
    Gt = np.zeros(len(Pgen_Hour_Status))
    Annual_Gt = 0
    
    Rd = 0.5 * (1 + cosd(Module_Tilt))
    Rg = 0.5 * (1 - cosd(Module_Tilt))
    n = math.floor((Array_Height)/(Lmod * sind(Module_Tilt)))

    for hour in range(len(Gt)):
        if Pgen_Hour_Status [hour] == 1:
            Rb [hour] = cosd(STAngIncidence [hour])
            Gt [hour] = DNI_SH [hour] * Rb [hour] + DHI_SH [hour] * Rd + \
                    GHI_SH [hour] * Rg * Albedo
            # in MW/sq.m
            Annual_Gt = Annual_Gt + Gt [hour]/1000000
                    
    return (Gt, n, Annual_Gt)

#-------------------------- End of Net Effective Radiation function -----------

# Calling Function that estimates net effective radiation 
    
Gt, n, AnnualGt = Net_Effective_Radiation (User_Assumed_Inputs.Module_Tilt, \
                         User_Assumed_Inputs.Misc_Array_height, \
                         STAngIncidence, GHI_SH, DHI_SH, DNI_SH, \
                         User_Assumed_Inputs.Misc_Albedo, \
                         User_Assumed_Inputs.Module_Lmod, \
                         User_Assumed_Inputs.Misc_Ground_Clearance, \
                         PgenHourStatus)
    
#------------------------------------------------------------------------------
# Determination of cell temperature and resource to module power factor
#------------------------------------------------------------------------------
#
def Tcell_ResToModPower (Pgen_Hour_Status, Gt, Tamb_SH, WS_SH, a_CT, b_CT, 
                         DelT, Pmod, Kt_Pmod, Kt_Isc, Kt_Voc, Isc, Voc):
   
    Gref = 1000
    Tref = 25
    
    Tcell = np.zeros(len(Pgen_Hour_Status))
    Mod_Isc = np.zeros(len(Pgen_Hour_Status))
    Mod_Voc = np.zeros(len(Pgen_Hour_Status))
    R_Pmod = np.zeros(len(Pgen_Hour_Status))
    Agg_R_Pmod = 0
    
    for hour in range(len(Pgen_Hour_Status)):
        if Pgen_Hour_Status [hour] == 1:
            Tcell [hour] = Gt [hour] * math.exp(a_CT + b_CT * WS_SH [hour]) + \
                    Tamb_SH [hour] + DelT * Gt [hour]/Gref
            R_Pmod [hour] = (Gt [hour]/Gref) * (1 + (Kt_Pmod/100) * (Tcell [hour] -\
               Tref))
            Agg_R_Pmod = Agg_R_Pmod + R_Pmod [hour]
            Mod_Isc [hour] = Isc * (Gt [hour]/Gref) * (1 + (Kt_Isc/100) * \
                (Tcell [hour] - Tref))
            Mod_Voc [hour] = Voc + (Kt_Voc/100) * (Tcell [hour] - Tref)
    
    return (Tcell, Mod_Isc, Mod_Voc, R_Pmod, Agg_R_Pmod)

#----- End of cell temperature and resource to module power function ----------

# Calling Function that estimates cell temperature and resource to 
# module power function
ModTCell, ModIsc, ModVoc, RPmod, AggRPmod = \
    Tcell_ResToModPower (PgenHourStatus, Gt, Tamb_SH, WS_SH, \
                         User_Assumed_Inputs.Mount_a_CT, \
                         User_Assumed_Inputs.Mount_b_CT, \
                         User_Assumed_Inputs.Mount_DelT, \
                         User_Assumed_Inputs.Module_Pmod, \
                         User_Assumed_Inputs.Module_Kt_Pmax, \
                         User_Assumed_Inputs.Module_Kt_Isc, \
                         User_Assumed_Inputs.Module_Kt_Voc, \
                         User_Assumed_Inputs.Module_Isc, \
                         User_Assumed_Inputs.Module_Voc)
    
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
    if V_PCU_Choice == 2:
        V_PCU_Ref = Vuser
    
    Mod_Pmax = R_Pmod.max() * Pmod
    Mod_Voc_min = Mod_Voc[np.nonzero(Mod_Voc)].min()
    Mod_Voc_max = Mod_Voc[np.nonzero(Mod_Voc)].max()
    Mod_Isc_max = Mod_Isc[np.nonzero(Mod_Voc)].max()
    N_PCU = Pplant_Target *1000 // P_nomDC
    m_max = V_maxDC // Mod_Voc_max
    m = math.ceil(V_PCU_Ref / Vmp )
    I_PCU_Ref = P_nomDC * 1000 / V_PCU_Ref
    y_max = I_maxDC // (n * Mod_Isc_max)
    y = math.ceil (I_PCU_Ref / (n * Imp))
    N_mod_PCU = n * m * y
    m_change = 0
    DC_Losses = (1 - SL_PC)
    P_PCU_DC_Max = N_mod_PCU * Mod_Pmax * DC_Losses / 1000
    
#    m_Pmod_max = Mod_Pmax * m/1000
#    P_PCU_DC_Max = 230
    
    Original_P_plant = N_mod_PCU * N_PCU * Pmod / 1000000
    
    if (m > m_max) or (y > y_max):
        print('Invalid choice for PCU reference voltage')
        
    if (P_PCU_DC_Max > P_nomDC):
        while(P_PCU_DC_Max > P_nomDC):
            N_mod_PCU_new = N_mod_PCU - m
            m_change = m_change - 1
            P_PCU_DC_Max = N_mod_PCU_new * Mod_Pmax * DC_Losses /1000
            N_mod_PCU = N_mod_PCU_new
            y = N_mod_PCU / (n * m)
                        
    else:
        while (P_PCU_DC_Max <= P_nomDC):
            N_mod_PCU_new =  N_mod_PCU + m
            m_change = m_change + 1
            P_PCU_DC_Max = N_mod_PCU_new * Mod_Pmax * DC_Losses /1000
            N_mod_PCU = N_mod_PCU_new
            y = N_mod_PCU / (n * m)
            
    Pmax_Deviation_NewCap =  (P_PCU_DC_Max/P_nomDC) - 1
    y_area = math.ceil (y)

    N_Plant = N_mod_PCU * N_PCU 
    P_plant = N_Plant * Pmod / 1000000
    Pure_mod_area = N_Plant * Lmod * Bmod
    Pure_mod_area_acres = Pure_mod_area * 0.00024710
    DC_Scaling_Factor = N_mod_PCU * Pmod / (1000 * P_PCU_AC)
    P_plant_PCU = N_PCU * P_PCU_AC / 1000
    
    
    return (Mod_Pmax, Mod_Voc_min, Mod_Voc_max, Mod_Isc_max, N_PCU, m, \
            m_change, N_mod_PCU, y, y_area, N_Plant, Pure_mod_area, \
            Pure_mod_area_acres, DC_Scaling_Factor, P_plant, P_plant_PCU, \
            Original_P_plant, Pmax_Deviation_NewCap, P_PCU_DC_Max)

#-------- End of plant sizing parameters function -----------------------------

# Calling function that estimates the plant sizing parameters
# This function optimises the plant design to generate maximum power for the 
# available resource

ModPmax, ModVocMin, ModVocMax, ModIscMax, N_PCU, m, m_change, NmodPCU, y, \
y_area, NPlant, PureModArea, PureModAreaAcres, DCScalingFactor, Pplant, \
PplantPCU, OriginalPplant, PmaxDeviationNewCap, P_PCU_DCMax = \
    Plant_Sizing_Parameters (User_Assumed_Inputs.V_PCU_Choice, \
                         User_Assumed_Inputs.Vuser, n, RPmod, ModIsc, ModVoc, \
                         User_Assumed_Inputs.Pplant_Target, \
                         User_Assumed_Inputs.PCU_PnomDC, \
                         User_Assumed_Inputs.PCU_P_PCU_AC, \
                         User_Assumed_Inputs.Module_Vmp, \
                         User_Assumed_Inputs.Module_Imp, \
                         User_Assumed_Inputs.Module_Pmod, \
                         User_Assumed_Inputs.PCU_VmaxDC, \
                         User_Assumed_Inputs.PCU_ImaxDC, \
                         User_Assumed_Inputs.Misc_SL_PC, \
                         User_Assumed_Inputs.Misc_EL_PC, \
                         User_Assumed_Inputs.Module_Lmod, \
                         User_Assumed_Inputs.Module_Bmod, \
                         User_Assumed_Inputs.PCU_Vmppmax, \
                         User_Assumed_Inputs.PCU_Vmppmin)

#------------------------------------------------------------------------------
# Estimation of inter-row and inter-column spacing
#------------------------------------------------------------------------------
#
def Inter_Row_Col_Spacing (Pgen_Window, Sun_Hour_Status, Pgen_Hour_Status, \
                           ZT_Day_Hour, L_row_Factor, L_col_Factor, n, Lmod, \
                           R_Pmod, Gt):
    
    Pgen_2hr_Window = (Pgen_Window [0] + 1, Pgen_Window [1] - 1)
    Pgen_4hr_Window = (Pgen_Window [0] + 2, Pgen_Window [1] - 2)
    
   # index = 0 --> Pgen hrs, 1 --> Pgen - 2 hrs,  2 --> Pgen - 4hrs
    Agg_R_Pmod = np.zeros (3)
        
    Lrow = np.zeros(len(Sun_Hour_Status))
    Lcol = np.zeros(len(Sun_Hour_Status))
    Gt_Pgen = np.zeros(len(Sun_Hour_Status))
    Gt_Pgen_2hr = np.zeros (len(Sun_Hour_Status))
    Gt_Pgen_4hr = np.zeros (len(Sun_Hour_Status))
    
    Drow = np.zeros(3)
    Dcol = np.zeros(3)
    
    Annual_GenHr = np.zeros (3)
        
    for hour in range(len(Sun_Hour_Status)):
        if Sun_Hour_Status [hour] == 1:
            Lrow [hour] = n * Lmod * L_row_Factor [hour]
            Lcol [hour] = n * Lmod * L_col_Factor [hour]
        
        if Pgen_Hour_Status [hour] == 1:
            Agg_R_Pmod [0] = Agg_R_Pmod [0] + R_Pmod [hour]
            Annual_GenHr [0] = Annual_GenHr [0] + 1
            Gt_Pgen [hour] = Gt [hour]
            if Lrow [hour] > Drow [0]:
                Drow [0] = Lrow [hour]
            if Lcol [hour] > Dcol [0]:
                Dcol [0] = Lcol [hour]
        
        if (ZT_Day_Hour [hour] >= Pgen_2hr_Window [0]) and (ZT_Day_Hour [hour] <= \
           Pgen_2hr_Window [1]):
            Agg_R_Pmod [1] = Agg_R_Pmod [1] + R_Pmod [hour]
            Annual_GenHr [1] = Annual_GenHr [1] + 1
            Gt_Pgen_2hr [hour] = Gt [hour]
            if Lrow [hour] > Drow [1]:
                Drow [1] = Lrow [hour]
            if Lcol [hour] > Dcol [1]:
                Dcol [1] = Lcol [hour]
        
        if (ZT_Day_Hour [hour] >= Pgen_4hr_Window [0]) and (ZT_Day_Hour [hour] <= \
           Pgen_4hr_Window [1]):
            Agg_R_Pmod [2] = Agg_R_Pmod [2] + R_Pmod [hour]
            Annual_GenHr [2] = Annual_GenHr [2] + 1
            Gt_Pgen_4hr [hour] = Gt [hour]
            if Lrow [hour] > Drow [2]:
                Drow [2] = Lrow [hour]
            if Lcol [hour] > Dcol [2]:
                Dcol [2] = Lcol [hour]
    
    # Annaul radiation on the tilted panel in MW / sq.m  
    An_Gt = np.zeros (3)
    An_Gt [0] = Gt_Pgen.sum()/1000000
    An_Gt [1] = Gt_Pgen_2hr.sum()/1000000
    An_Gt [2] = Gt_Pgen_4hr.sum()/1000000
    
    return (Drow, Dcol, Annual_GenHr, Agg_R_Pmod, Pgen_2hr_Window, \
            Pgen_4hr_Window, Gt_Pgen, Gt_Pgen_2hr, Gt_Pgen_4hr, An_Gt)

#-------- End of Inter-row and Inter-col spacing  function -------------------
    
# Calling the function that estimates inter-row, inter-column spacing;
# aggregate resource to module power factor and annual generation hours 
# for all 3 time windows
    
Drow, Dcol, AnnualGenHr, AggRPmod, Pgen2hrWindow, Pgen4hrWindow, GtPgen, \
GtPgen2hr, GtPgen4hr, AnGt = \
        Inter_Row_Col_Spacing (PgenWindow, SunHourStatus, PgenHourStatus, \
                           ZTDayHour, LrowFactor, LcolFactor, n, \
                           User_Assumed_Inputs.Module_Lmod, RPmod, Gt)

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
        if math.sqrt(A_LB)% 1 == 0:
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
                if math.sqrt(A_LDc_temp)%1.0 == 0 and count < 1:
                    A_LDc_temp -= 1
                    A_BDr_temp -= 1
                    count += 1 # First element of TS traced
                    R += 1
                    
                # Resetting the TS counter when crossing the 2nd element of TS
                if math.sqrt(A_LDc_temp)%1.0 == 0 and count == 1:
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
    Gt_TW = np.zeros(len(Gt_Pgen))
    
    CERC_Benchmark = 5
    
    # Application of ulam spiral for placement of arrays for a PCU block
    NLe_PCU, NBe_PCU, NLo_PCU, NBo_PCU, Nre_PCU, Nce_PCU, Nro_PCU, Nco_PCU, \
    L_flag_PCU, B_flag_PCU =  spiral (y_area)
    
    # Calculating the base length and breadth dimensions of an array
    L = n * Lmod * cosd (Module_Tilt)
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
            Aux_P [ii] = 0.16273 * math.exp(-0.027 * P_plant) * Total_Plant_area[ii]
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


# Calling the function that estimates plant area, packing density and suitable
# time window
PlantArea, PackingDensity, DeviationFactor, TW, TW_DayHrs, Gt_TW = \
                Plant_Area_Estimation (n, User_Assumed_Inputs.Module_Lmod, \
                       User_Assumed_Inputs.Module_Bmod, \
                       User_Assumed_Inputs.Module_Tilt, \
                       Drow, Dcol, m, y_area, N_PCU, \
                       User_Assumed_Inputs.Misc_BS_AlongL_a, \
                       User_Assumed_Inputs.Misc_BS_AlongB_b, \
                       Pplant, PureModAreaAcres, PgenWindow, Pgen2hrWindow,\
                       Pgen4hrWindow, GtPgen, GtPgen2hr, GtPgen4hr )



#------------------------------------------------------------------------------
# Estimation of Plant output
#------------------------------------------------------------------------------
#
def Plant_output_year_0 (TW, An_Gt, Pure_Mod_Area, SL_PC, EL_PC, TW_DayHrs, \
                         N_Plant, ZT_Day_Hour, R_Pmod, Pmod, Eta_PCU, \
                         P_plant, Gt_TW):
    
    # Data initialisation
    #TW_Status = np.zeros(len(Gt_TW))
    R_Pplant_TW = np.zeros (len(Gt_TW))
    Pmod_out = np.zeros (len(Gt_TW))
    P_plant_PV_DC = np.zeros (len(Gt_TW))
    P_plant_PV_DC_inc_SL = np.zeros (len(Gt_TW))
    P_plant_PV_AC = np.zeros (len(Gt_TW))
    Agg_R_Pplant_AC_TW = 0.0
    Hourly_Gt_all_pmod = np.zeros(len(Gt_TW))
    Hourly_Gen_PC = np.zeros(len(Gt_TW))
    Hourly_SEE_PV = np.zeros (len(Gt_TW))
    
    # Aggregate solar radiation incident on solar panels in MW/sq.m
    Agg_Gt_all_modules_TW = An_Gt [TW] * Pure_Mod_Area
    
    Soiling_Loss_Factor = 1 - SL_PC
    Electrical_Loss_Factor = 1 - EL_PC
    
    for hour in range(len(Gt_TW)):
        if (ZT_Day_Hour [hour] >= TW_DayHrs [0]) and (ZT_Day_Hour [hour] <= \
        TW_DayHrs [1]):
            R_Pplant_TW [hour] = N_Plant * R_Pmod [hour] /1000000
            Pmod_out [hour] = Pmod * R_Pmod [hour]
            P_plant_PV_DC [hour] = Pmod * R_Pplant_TW [hour]
            P_plant_PV_DC_inc_SL [hour] = P_plant_PV_DC [hour] * \
                                            Soiling_Loss_Factor
            P_plant_PV_AC [hour] = P_plant_PV_DC_inc_SL [hour] * \
                                    Electrical_Loss_Factor * Eta_PCU
            Hourly_Gen_PC [hour] = P_plant_PV_AC [hour]/P_plant
#            Hourly_Gen_PC [hour] = 
            Hourly_Gt_all_pmod [hour] = Gt_TW [hour] * Pure_Mod_Area/1000000
            if Gt_TW [hour] == 0:
                Hourly_SEE_PV [hour] = 0
            else:
                Hourly_SEE_PV [hour] = \
                    P_plant_PV_AC [hour] / (Hourly_Gt_all_pmod [hour])
        
        else:
            P_plant_PV_AC [hour] = 0
            Hourly_Gen_PC [hour] = 0
            Hourly_SEE_PV [hour] = 0
            
    
    Agg_R_Pplant_AC_TW = (R_Pplant_TW * Soiling_Loss_Factor * \
                            Electrical_Loss_Factor * Eta_PCU).sum()
    print(Agg_R_Pplant_AC_TW)
    # Annual parameters in MWh
    An_Pplant_DC = P_plant_PV_DC.sum()
    An_Pplant_DC_inc_SL = P_plant_PV_DC_inc_SL.sum()
    An_Pplant_AC = P_plant_PV_AC.sum()
    An_soiling_loss = An_Pplant_DC - An_Pplant_DC_inc_SL
    An_electrical_loss = An_Pplant_DC_inc_SL - An_Pplant_AC
    
    An_PR_PV = An_Pplant_AC / (An_Gt [TW] * N_Plant * Pmod / 1000)
    An_CUF_PV = An_Pplant_AC / (8760 * P_plant)
    An_SEE_PV = An_Pplant_AC / (Agg_Gt_all_modules_TW)
    
    # Adding An_Pplant_DC_inc_SL in the return list, 8 May 2019
    return (Hourly_Gen_PC, Hourly_SEE_PV, Hourly_Gt_all_pmod, \
            P_plant_PV_AC, Agg_Gt_all_modules_TW, An_Pplant_DC, \
            An_Pplant_DC_inc_SL, An_Pplant_AC, An_soiling_loss, \
            An_electrical_loss, An_PR_PV, An_CUF_PV, An_SEE_PV)

#------------------------ End of Plant output for year 0 function -------------

# Calling the function that estimates plant generation in year 0 and other
# performance metrics

HourlyGenPC, HourlySEEPV, HourlyGtallpmod, PplantPVAC, AggGtallmodulesTW, \
AnPplantDC, AnPplantDCincSL, AnPplantAC, AnSoilingLoss, AnElectricalLoss, \
AnPRPV, AnCUFPV,  AnSEEPV = Plant_output_year_0 (TW, AnGt, PureModArea, \
                                 User_Assumed_Inputs.Misc_SL_PC, \
                                 User_Assumed_Inputs.Misc_EL_PC, \
                                 TW_DayHrs, NPlant, ZTDayHour, RPmod, \
                                 User_Assumed_Inputs.Module_Pmod, \
                                 User_Assumed_Inputs.PCU_Efficiency,\
                                 Pplant, Gt_TW) 
     

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
                   Eta_PCU/1000000
                   
    
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
                    Annual_PV_Gen_Deg [TimeWindow][Deg_ind][yr] / \
                    (An_Gt [TimeWindow] \
                            * N_Plant * Net_module_rating [Deg_ind][yr] / 1000)
    
                Annual_SEE_PV [TimeWindow][Deg_ind][yr] = \
                    Annual_PV_Gen_Deg [TimeWindow][Deg_ind][yr] / \
                                    (Agg_Gt_all_Modules [TimeWindow])
    
    return (Degradation_Rates, Base_DegCase_Index, Agg_RPF, Annual_PV_Gen_Deg,\
            Net_module_rating, Annual_CUF_PV, Module_Rating_PC, Annual_PR_PV, \
            Annual_SEE_PV)

# the financial model uses only generation (Annual_PV_Gen_Deg) 
# for time window = TW, and base degradation case index = Base_DegCase_Index

# -------------------------- end of degradation results function --------------
 
# Calling the function that estimates the PV plant output factoring 
# Module degradation 
           
DegradationRates, BaseDegCaseIndex, AggRPF, AnnualPVGenDeg, NetModuleRating, \
AnnualCUFPV, ModuleRatingPC, AnnualPRPV, AnnualSEEPV  = Degradation_Results\
                     (User_Assumed_Inputs.Module_YoY_Degrdn_rate, \
                      User_Assumed_Inputs.Plant_Life, \
                      User_Assumed_Inputs.Module_Pmod, \
                      User_Assumed_Inputs.Module_Rating_EoYr1, \
                      NPlant, AggRPmod,\
                      User_Assumed_Inputs.PCU_Efficiency, \
                      User_Assumed_Inputs.Misc_SL_PC, \
                      User_Assumed_Inputs.Misc_EL_PC, Pplant, PureModArea, AnGt)
            
            
#------------------------------------------------------------------------------
# Estimating monthly, daily aggregate gen, Histogram for frequency of gen
#------------------------------------------------------------------------------
#
def Monthly_Daily_Gen_Hist1 (P_plant_PV_AC, ZT_Day_Hour, P_plant_PCU, \
                             ZoneTime):
    
    PV_Gen_Yr0_Daily = np.sum (P_plant_PV_AC.reshape(365,-1), axis = 1)
    PV_Gen_Yr0_Monthly = np.zeros(12)
    mt = np.zeros(12)
    
    Ag_Gen_Hrs_1to100PC_TW = 0
    Gen_Hrs_Yr0_TW_PC = np.zeros(11)
    
    for hr in range(len(P_plant_PV_AC)):
        
        Gen = P_plant_PV_AC [hr]
        # modified monthly aggregate generation estimation on 6 May 2019
        mon = ZoneTime[hr].month - 1
        day = ZoneTime[hr].day
        mt[mon] = day
        PV_Gen_Yr0_Monthly[mon] = PV_Gen_Yr0_Monthly[mon] + Gen
        
        if (Gen >= 0.01 * P_plant_PCU) and (Gen <= 0.1 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [0] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
            
        if (Gen > 0.1 * P_plant_PCU) and (Gen <= 0.2 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [1] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
        
        if (Gen > 0.2 * P_plant_PCU) and (Gen <= 0.3 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [2] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
            
        if (Gen > 0.3 * P_plant_PCU) and (Gen <= 0.4 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [3] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
            
        if (Gen > 0.4 * P_plant_PCU) and (Gen <= 0.5 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [4] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
        
        if (Gen > 0.5 * P_plant_PCU) and (Gen <= 0.6 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [5] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
        
        if (Gen > 0.6 * P_plant_PCU) and (Gen <= 0.7 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [6] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
            
        if (Gen > 0.7 * P_plant_PCU) and (Gen <= 0.8 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [7] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
            
        if (Gen > 0.8 * P_plant_PCU) and (Gen <= 0.9 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [8] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
            
        if (Gen > 0.9 * P_plant_PCU) and (Gen <= 1 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [9] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
            
        if (Gen > 1 * P_plant_PCU):
            Gen_Hrs_Yr0_TW_PC [10] += 1
            Ag_Gen_Hrs_1to100PC_TW += 1
        
    return (PV_Gen_Yr0_Daily, PV_Gen_Yr0_Monthly, Gen_Hrs_Yr0_TW_PC, \
            Ag_Gen_Hrs_1to100PC_TW)

#---End of function which estimates monthly and daily gen and hist 1-----------

# Calling function which estimates the monthly and daily generation and 
# generates histogram to measure frequency of generation 
    
PVGenYr0Daily, PVGenYr0Monthly, GenHrsYr0TWPC, AgGenHrs1to100PCTW \
     = Monthly_Daily_Gen_Hist1 (PplantPVAC, ZTDayHour, PplantPCU, ZoneTime)
    
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
    
    
GenPVACTW24x365, GenPVACTW24x365T, MinGenPVACTW24x365, MaxGenPVACTW24x365, \
    MeanGenPVACTW24x365, SumGenPVACTW24x365, PCShareGenACTW24x365 \
    = Gen_vs_DayHr_Hist2 (PplantPVAC)

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
 
# Calling the function estimating miscellaneous parameters in the tech model 

AnPVGenPerAcre, AreaPgenPgen2hrFactor, AreaPgen2hrPgen4hrFactor, \
 AreaPgenPgen4hrFactor, GtPgenPgen2hrFactor, GtPgen2hrPgen4hrFactor, \
 GtPgenPgen4hrFactor, GenPgenPgen2hrFactor, GenPgen2hrPgen4hrFactor, \
 GenPgenPgen4hrFactor, GenPCPgen2hrPgen, GenPCPgen4hrPgen, \
 GenHrPgenPgen2hrFactor, GenHrPgen2hrPgen4hrFactor, GenHrPgenPgen4hrFactor, \
 GenPCTW24x365T, GtallmodTW24x365T, SEEPVTW24x365T \
 = Misc_Tech_Parameters (PlantArea, AnnualPVGenDeg, TW, \
                          BaseDegCaseIndex, AnGt, AnnualGenHr,  HourlyGenPC, \
                          HourlyGtallpmod, HourlySEEPV, Gt_TW)
 
