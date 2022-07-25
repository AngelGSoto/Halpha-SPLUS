# -*- coding: utf-8 -*-
""" 
Created on 26/06/19
"""
from __future__ import print_function, division

import os

import numpy as np
from astropy.table import Table
from astropy.io import fits

tables_dir = "zp"
file_ = "zp.ecsv"
def get_zps_dr3(Field, bands):
    """ Read the table containing the zero points for a given tile and given
    bands. """
    zpfile = os.path.join(tables_dir, file_)

    zpdata = Table.read(zpfile, format="ascii.ecsv")
    zpdic = {a: {"R": b,
                 "F660": c,
                 "I": d}
              for a, b, c, d in zip(zpdata["field"],
                                    zpdata["R"],
                                    zpdata["F660"],
                                    zpdata["I"])}

    zps = np.array([zpdic[Field][band] for band in bands])
    return zps

if __name__ == "__main__":
    bands = ["R", "F660", "I"]
    test_zp = get_zps_dr3('STRIPE82-0001', bands)
    print(test_zp)
  
