# -*- coding: utf-8 -*-
"""
Author : Luis A. Gutiérrez -Soto 
based on Carlos Eduardo Barbosa code
Methods to extract H-alpha and other line with SPLUS filters according to
methods presented in Villela-Rojo et al. 2015 (VL+2015).
"""
from __future__ import division, print_function
import os
import astropy.units as u
import numpy as np
from scipy.interpolate import interp1d

def calc_deltax(trans, wline = 6562.8):
    """ Numerical integration of equation (4) in VR+2015 for SPLUS
       system -- beta. """
    wave = trans["wl"]
    trans = trans["Flux"]
    curve = interp1d(wave, trans, kind="linear", fill_value=0.,
                     bounds_error=False)

    deltax = np.trapz(trans * wave, wave) / (curve(wline) * wline) 
    return deltax

def calc_alphax(trans):
    """ Numerical integration of equation (4) in VR+2015 for SPLUS
       system -- alpha """
    wave = trans["wl"]
    trans = trans["Flux"]
    term1 = np.trapz(trans * wave * wave, wave) 
    term2 = np.trapz(trans * wave, wave) 
    alphax = term1 / term2
    return alphax

def calc_betax(trans, wline = 6562.8):
    """ Numerical integration of equation (4) in VR+2015 for SPLUS
       system -- beta. """
    wave = trans["wl"]
    trans = trans["Flux"]
    curve = interp1d(wave, trans, kind="linear", fill_value=0.,
                     bounds_error=False)

    betax = (wline * curve(wline)) / np.trapz(trans * wave, wave)
    return betax

# Estimating the estimator for each filter
def file_band(file_, Path="filter_curves/"):
    """ Read the band-file"""
    dt = np.dtype([ ('wl', 'f8'), ('Flux', 'f8')])
    band = np.loadtxt(os.path.join(Path, file_), dtype = dt)
    return band

if __name__ == "__main__":
    r =  file_band("rSDSS.dat")
    f660 =  file_band("F0660.dat")
    i =  file_band("iSDSS.dat")

    # estimator
    print("Delta_f660:", calc_deltax(f660))
    print("----------------------------------------------------------------")
    print("Alpha_r:", calc_alphax(r))
    print("Alpha_f660:", calc_alphax(f660))
    print("Alpha_i:", calc_alphax(i))
    print("----------------------------------------------------------------")
    print("Beta_r:", calc_betax(r))
    print("Beta_f660:", calc_betax(f660))
    
    
    
