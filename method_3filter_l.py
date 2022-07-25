# -*- coding: utf-8 -*-
""" 
Author : Luis A. Gutiérrez-Soto
Script to recover the Ha emission by appliying the 3fiters method, based on Villela-Rojo et al. 2015 (VR+15)
"""
import splus_estimator_l as spe

# Read the file with the transmision curves
r =  spe.file_band("rSDSS.dat")
f660 = spe.file_band("F0660.dat")
i = spe.file_band("iSDSS.dat")

# Estmating alpha and beta parameters
#ALPHAx
alpha_r = spe.calc_alphax(r)     
alpha_f660 = spe.calc_alphax(f660)
alpha_i = spe.calc_alphax(i)       
#BETAxx
beta_r = spe.calc_betax(r)   
beta_f660 = spe.calc_betax(f660)

# for alpha
a = (alpha_r - alpha_i) / (alpha_f660 - alpha_i)

def method_3filter(Fr, Ff660, Fi):
    """ Estimating the Ha + [NII] emission"""
    halpha = ((Fr - Fi) - a * (Ff660 - Fi)) / (beta_r - beta_f660 * a)
    return halpha
    
    
    
