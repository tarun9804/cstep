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
