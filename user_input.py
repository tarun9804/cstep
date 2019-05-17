#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:56:29 2019

@author: tarun
"""
import pandas as pd


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