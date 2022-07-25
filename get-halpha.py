# -*- coding: utf-8 -*-
""" 
Author : Luis A. Gutiérrez-Soto
Script to recover the Ha emission by applying the 3-fiters method, based on Villela-Rojo et al. 2015 (VR+15).
"""
import splusdata
import getpass
import pandas as pd
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import aplpy
from scipy.ndimage import gaussian_filter
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy.wcs import WCS
import os
import argparse
import sys
from pathlib import Path
import astropy.constants as const
import cv2
import misc 
from method_3filter_l import method_3filter
ROOT_PATH = Path("..")

# Connecting with SPLUS database
username = str(input("Login: "))
password = getpass.getpass("Password: ")

conn = splusdata.connect(username, password)

parser = argparse.ArgumentParser(
    description="""Estomatting the Halpha + [N II] in SPLUS""")

parser.add_argument("table", type=str,
                    default="table with the sources",
                    help="Name of table, taken the prefix ")

parser.add_argument("--Object", type=str,
                    default=None,
                    help="Id object of the source under interest ")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info about each line in region file")

args = parser.parse_args()
file_ = args.table + ".ecsv"

try:
    data = Table.read(file_, format="ascii.ecsv")
except FileNotFoundError:
    file_ = args.table + ".dat"
    data = Table.read(file_, format="ascii")

if args.Object is not None:
    Object_ = str(args.Object)
    mask = np.array([source in Object_ for source in data["ID"]])
    data = data[mask]
else:
    data = data

# speed light A / s
c = 2.99792458e18
############################################################
# Definition to decompress the images ######################
############################################################
def fz2fits(image):
    """
    It converts SPLUS images
    from .fz to .fits
    """
    data = fits.open(image)[1].data
    head = fits.open(image)[1].header
    imageout = image[:-2] + 'fits'
    print ('Creating file: ')
    print (imageout)
    fits.writeto(imageout, data, head, overwrite=True)
############################################################
# wave_eff = {"F378": 3771.0, "F395": 3941.0, "F410": 4094.0,
#             "F430": 4292.0, "F515": 5133.0, "F660": 6614.0,
#             "F861": 8611.0, "G": 4756.0, "I": 7692.0, "R": 6260.0,
#             "U": 3574.0, "Z": 8733.0}
wave_eff = {"R": 6260.0,
            "F660": 6614.0,
            "I": 7692.0 }

for tab in data:
    ra = tab["RA"]
    dec = tab["DEC"]
    Field = tab["Field"]
    Name = tab["ID"].split("R3.")[-1].replace(".", "-")
    
    for b, v in wave_eff.items():
        filter_ = conn.get_cut(ra, dec, 150, b)
        filter_.writeto('{}_{}.fz'.format(Name, b), overwrite=True) # write to fits
        
        # Decompress
        fz2fits('{}_{}.fz'.format(Name, b))
    ##################################################################
    #Read the FITS file
    ##################################################################
    data = np.array([fits.open('{}_{}.fits'.format(Name, bandss))[0].data for bandss in wave_eff.keys()])

    # Getting the zero point
    #zp = 21
    zp = misc.get_zps_dr3(Field, [bands for bands in wave_eff.keys()])
    
    f0 = np.power(10, -0.4 * (48.6 - 20)) # fluxt const without zero point 
    fnu = f0 * data
    flm_r = np.power(10, -0.4*zp[0]) * fnu[0] * c / wave_eff["R"]**2
    flm_f660 = np.power(10, -0.4*zp[1]) * fnu[1] * c / wave_eff["F660"]**2
    flm_i = np.power(10, -0.4*zp[2]) * fnu[2] * c / wave_eff["I"]**2
    #idx = np.array([bands.index(band) for band in bands])
    
    #Getting the Ha + [NII]
    halpha = method_3filter(flm_r, flm_f660, flm_i)

    #Plotting
    fig, ax = plt.subplots()
    m, s = np.mean(halpha), np.std(halpha)
    inverted_image = cv2.bitwise_not(halpha)
    #plt.imshow(inverted_image, vmin=m-s, vmax=m+s, interpolation = 'nearest', origin='lower', filternorm=True, cmap='gray')
    plt.imshow(halpha, vmin=m-s, vmax=m+s, interpolation = 'nearest', origin='lower', filternorm=True, cmap='coolwarm')
    
    plt.savefig("halpha-image-" + Name + ".pdf")
    
    # Deleting files
    for b, v in wave_eff.items():
        os.remove('{}_{}.fz'.format(Name, b))
        print("FZ File Removed!")
    
        

        
        
    
