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
print((ts2-ts1)/1e6,"- ZoneTime_SolarAngles(ms)")

# Calling Function that estimates solar resource at different temporal resolution
ts1=time.time_ns()
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
print((ts2-ts1)/1e3,"- Plant_Sizing_Parameters(micro sec)")


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
print((ts2-ts1)/1e3,"- Plant_Area_Estimation(micro sec)")

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
