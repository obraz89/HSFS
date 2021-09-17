# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 14:37:29 2021

@author: vmdev
"""

import configparser

SECT_FLOW = 'Flow'
KEY_P0 = 'P0'
KEY_T0 = 'T0'
KEY_MACH = 'Mach'
SECT_GAS = 'Gas'
KEY_GAMMA = 'Gamma'
KEY_MOLWGT = 'MolWgt'
KEY_R = "R"
def calc_flow_params():
    
    config = configparser.ConfigParser()    
    config.read("main.ini")
    check_config(config)
       
    T0 = config.getfloat(SECT_FLOW, KEY_T0)
    P0 = config.getfloat(SECT_FLOW, KEY_P0)
    M  = config.getfloat(SECT_FLOW, KEY_MACH)
    gm = config.getfloat(SECT_GAS, KEY_GAMMA)
    R = config.getfloat(SECT_GAS, KEY_R)
    MolWgt = config.getfloat(SECT_GAS, KEY_MOLWGT)
    
    coef = 1.0 + 0.5*(gm - 1.0)*M*M
    T_inf = T0/coef
    print("T_inf[K]={:.6e}".format(T_inf))
    Mu_inf = calc_visc_suther_dim(T_inf)
    print("Visc_inf[kg/ms]={:.6e}".format(Mu_inf))
    c_inf = (gm*R*T_inf/MolWgt)**0.5
    print("C_inf[m/s]={:.6e}".format(c_inf))
    u_inf = c_inf*M
    print("u_inf[m/s]={:.6e}".format(u_inf))
    rho_inf = P0*MolWgt*coef**(-1.0/(gm-1.0))/(R*T0)
    print("rho_inf[kg/m^3]={:.6e}".format(rho_inf))
    Re1 = rho_inf*u_inf*1.0/Mu_inf
    print("Re1[1/m]={:.6e}".format(Re1))
    
    
# check config and write defaults if param is missing
def check_config(config):
    # flow
    if (not SECT_FLOW in config):
        config[SECT_FLOW] = {}
    if (not config[SECT_FLOW].get(KEY_P0)):
        config[SECT_FLOW][KEY_P0] = '1.0e5'
    if (not config[SECT_FLOW].get(KEY_T0)):
        config[SECT_FLOW][KEY_T0] = '300.0'
    if (not config[SECT_FLOW].get(KEY_MACH)):
        config[SECT_FLOW][KEY_MACH] = '6.0'
    # gas
    if (not SECT_GAS in config):
        config[SECT_GAS] = {}
    if (not config[SECT_GAS].get(KEY_GAMMA)):    
        config[SECT_GAS][KEY_GAMMA] = '1.4'
    if (not config[SECT_GAS].get(KEY_MOLWGT)):
        config[SECT_GAS][KEY_MOLWGT] = '0.029'
    if (not config[SECT_GAS].get(KEY_R)):
        config[SECT_GAS][KEY_R] = '8.31'
        
    with open("main.ini", 'w') as configfile:
        config.write(configfile)

def calc_visc_suther_dim(T):
    mu_ref = 1.716e-05
    T_ref = 273.15
    S = 110.4
    return mu_ref*(T/T_ref)**1.5*(T_ref+S)/(T+S)
    
calc_flow_params()