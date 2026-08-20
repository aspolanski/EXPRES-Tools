#!/usr/bin/env python3

"""
Helper script that will take a directory of EXPRES CCF files,
fit the CCFS, and return the recommended fitting window size.

Author: Alex Polanski
Date: August, 20 2026


"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys, argparse, glob

from astropy.io import fits
from astropy import units as u

from lmfit.models import GaussianModel,ConstantModel

ccf_model = GaussianModel() + ConstantModel()

parser = argparse.ArgumentParser()
parser.add_argument('--ccf_dir', dest='ccf_dir', type=str, help='Directory containing CCF fits files. (format = /*/*/)')
parser.add_argument('--fwhm_factor', dest='fwhm_factor', type=float, default=1.5, help='Fraction of the FWHM to use as the fitting window.')
parser.add_argument('--plot', dest='plot', action='store_true', default=False, help='Generate plots.')

args = parser.parse_args()

# First check if directory exists

if not os.path.isdir(args.ccf_dir):
    print(f"{args.ccf_dir} not found.")
    exit()
else:
    ccfs = glob.glob(f"{args.ccf_dir}*.fits")
    ccfs.sort()

    if len(ccfs) == 0:
        print(f"No CCF files found in {args.ccf_dir}")
        exit()

if args.plot:
    fig, ax = plt.subplots()

fwhms = []

for i, f in enumerate(ccfs):
    
    hdu = fits.open(f)
    
    v_grid, ccf = hdu[1].data['V_grid']/100000 , hdu[1].data['ccf']
    
    min_rv = v_grid[np.argmin(ccf)]
    
    amp_guess = min(ccf)-np.nanmedian(ccf)
    
    params = ccf_model.make_params(center=min_rv, amplitude=amp_guess,c=np.nanmedian(ccf))
    result = ccf_model.fit(ccf, params, x=v_grid)

    print(f"FWHM: {result.values['fwhm']:0.4} km/s | Center: {result.values['center']:0.4} km/s")
    
    fwhms.append(result.values['fwhm'])

    if args.plot:
    
        ax.plot(v_grid-result.values['center'],ccf)
        ax.plot(v_grid-result.values['center'],result.best_fit)
        ax.axvline(x=result.values['fwhm'])
        ax.axvline(x=-result.values['fwhm'])
        
        ax.axvline(x=(1.5*result.values['fwhm'])/2)
        ax.axvline(x=-(1.5*result.values['fwhm'])/2)
        
        ax.set_xlim(-3*result.values['fwhm'],3*result.values['fwhm'])
        
    
print(f"\nReccomended CCF Window: {(args.fwhm_factor*np.mean(fwhms))/2*100000:0.3}")

if args.plot:
    plt.show()