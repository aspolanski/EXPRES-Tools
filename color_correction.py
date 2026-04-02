#/usr/bin/env python3


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tqdm

try:
    with suppress_stdout_stderr(): #This does not work, need to figure that out
        import smplotlib #For fancy plots.
except:
    pass

from astropy.io import fits
import astropy.units as u
import astropy.constants as const

import os, sys, argparse, glob

"""
Script for correcting EXPRES spectra for variations in flux 
levels (e.g. resulting from airmass differences). Mixture of prescriptions
from the HARPS/ESPRESSO DRP and NEID DRP:

Creates a low-resolution spectrum from all the observations and compares it
to the expected spectral energy distribution for a star of the same effective
temperature. This is done by using a B-star spectrum multiplied with a 
blackbody. The differences between expected and observed flux levels are
fit with a 6th order polynomial and used to rescale the observed flux.

Color-corrected flux is stored as an additional fits rec entry: 'spectrum_ccor'

Author: Alex Polanski
2026-April-1
"""

parser = argparse.ArgumentParser()
parser.add_argument('--in_dir', dest='in_dir', type=str, help='The input FITS file.')
parser.add_argument('--out_dir',dest='out_dir', type=str, help='Directory of the registered spectra.')
parser.add_argument('--plot',dest='if_plot', type=bool, help='Whether to make shift plots.')
parser.add_argument('--target-temp',dest='target_temp', type=int, help='Effective temperature of the target star (K).')

args = parser.parse_args()


def planck(wav_aa, T):
    wav_cm = wav_aa * 1e-8
    h = 6.626e-27
    c = 2.998e10
    k = 1.381e-16
    return (2*h*c**2 / wav_cm**5) / (np.exp(h*c / (wav_cm*k*T)) - 1)


# Read in the EXPRES throughput file
expres_hdu = fits.open("./expres_throughput.fits")
inst_response = expres_hdu[1].data['throughput']


fits_files = glob.glob(f"{args.in_dir}/*.fits")

lowR_flux = np.zeros((86,len(fits_files)))
lowR_wave = expres_hdu[1].data['lowR_wave']
idx_low, idx_high = 1980, 5940  #use central third of order for integrating flux

target_bb = planck(lowR_wave,args.target_temp)
target_bb /= np.median(target_bb)

for f, fits_file in enumerate(tqdm.tqdm(fits_files)):
    hdu = fits.open(fits_file)

    for i,order in enumerate(range(hdu[1].data['spectrum'].shape[0])):

        flux = hdu[1].data['spectrum'][order,idx_low:idx_high] * hdu[1].data['blaze'][order,idx_low:idx_high]
        flux = flux / hdu[1].data['tellurics'][order,idx_low:idx_high]
        lowR_flux[i,f] = np.nansum(flux)

    lowR_flux[:,f] = lowR_flux[:,f]/np.median(lowR_flux[:,f])

    p = np.poly1d(np.polyfit(lowR_wave[1:-3],(inst_response*target_bb)[1:-3]/lowR_flux[1:-3,f],deg=6))
    
    scaled_flux = ((hdu[1].data['spectrum']*hdu[1].data['blaze'])*p(hdu[1].data['bary_wavelength']))/hdu[1].data['blaze']


    new_col = fits.Column(name='spectrum_ccor', format='7920D', array=scaled_flux)
    orig_cols = hdu[1].columns
    new_hdu = fits.BinTableHDU.from_columns(orig_cols + new_col, header=hdu[1].header)

    hdu[1] = new_hdu
    
    hdu.writeto(f'{args.out_dir}/{fits_file.split("/")[-1]}',overwrite=True)
    
    




