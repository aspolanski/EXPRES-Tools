#/usr/bin/env python3


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

try:
    with suppress_stdout_stderr(): #This does not work, need to figure that out
        import smplotlib #For fancy plots.
except:
    pass

from scipy.interpolate import interp1d
from scipy.signal import correlate

from astropy.io import fits
import astropy.units as u
import astropy.constants as const

from astroquery.simbad import Simbad
from lmfit.models import ConstantModel,VoigtModel
import os, sys, argparse


"""
Script for registering EXPRES spectra into the rest frame.
Cross-correlates with Solar spectrum to derive bulk RV.
Saves the FITS file with an additional extension 'REGISTERED'

Author: Alex Polanski
2025-May-21
"""

parser = argparse.ArgumentParser()
parser.add_argument('--in_file', dest='in_file', type=str, help='The input FITS file.')
parser.add_argument('--out_dir',dest='out_dir', type=str, help='Directory of the registered spectra.')
parser.add_argument('--plot',dest='if_plot', type=bool, help='Whether to make shift plots.')
parser.add_argument('--target',dest='target', type=str, help='FITS file to register to (defaults to Solar spectrum).')

args = parser.parse_args()

if args.in_file is None:
    print("Requires an input spectrum (--in_file)")
else:
    if not os.path.exists(args.in_file):
        sys.exit(f"{args.in_file} not found")
    else:
        hdu = fits.open(args.in_file)
        l=args.in_file.split('/')[-1].split('_')
        l.insert(1,"reg")
        new_file_name = '_'.join(l)   #Incredibly EXPRES-specific way to get a new file name.

if args.out_dir is None:
    out_dir = './fits_files/registered'
else:
    out_dir = args.out_dir

if args.if_plot is None:
    if_plot = True
else:
    if_plot = args.if_plot

if args.target is None:
    print("No target spectrum given. Defaulting to Solar spectrum.")
    target_hdu = fits.open("./Sun_210616.5095.fits") #The Solar spectrum currently be used, need to change to "idolic" Solar spectrum.
else:
    if not os.path.exists(args.target):
        sys.exit(f"{args.in_file} not found")
    else:
        target_hdu = fits.open(args.target)



 
def cross_correlate(x1,x2):

    """
    Perform cross-correlation of two arrays.
    """
    
    x1_sans_mean = x1 - np.nanmean(x1)
    x2_sans_mean = x2 - np.nanmean(x2)
    
    numerator = np.nansum(x1_sans_mean * x2_sans_mean)
    denominator = np.sqrt(np.nansum(x1_sans_mean**2)) * np.sqrt(np.nansum(x2_sans_mean**2))
    
    return numerator/denominator


def plot_ccf(rv,ccf,ccf_model,order,ax,offset):
    
    """
    Function to plot an individual CCF and its model.
    X and Y are flipped.
    """
    
    #get the order label position
    pos = int(len(rv)/6)
    
    ax.plot(ccf+offset,rv,c='k',linewidth=1.6)
    ax.plot(ccf_model+offset,rv,color='C10',linewidth=1.6,linestyle='dashed')
    ax.text(s=f'{order}',y=rv[pos],x=ccf_model[pos]+offset,ha='center', va='center',
        bbox=dict(facecolor='white', edgecolor='black',boxstyle='square,pad=0.2'),fontsize=10)
    
    
    
def dict_to_hdu(data_dict):

    """
    Takes a Python dictionary and returns a set of FITS columns.
    """
    #assert all_equal([len(data_dict[key]) for key in data_dict]), "Arrays must be the same length"
    cols = []
    for key in data_dict:
        arr = np.asarray(data_dict[key])
        #u = units[key] if key in units else None
        dim = None if arr.ndim < 3 else str(arr.shape[1:][::-1])
        cols.append(fits.Column(name=key, format='7920D',
                                array=arr, dim=dim))
        
    return fits.ColDefs(cols)

def get_radial_velocity(object_name):

    """
    Query Simbad to get the bulk RV of a star in order to center RV grid.
    """
    # Customize Simbad to include radial velocity
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('rv_value')  # Radial velocity field

    try:
        result = custom_simbad.query_object(object_name)
        if result is None:
            print(f"Object '{object_name}' not found in SIMBAD.")
            return None

        rv = result['RV_VALUE'][0]
        if rv is None:
            print(f"No radial velocity found for '{object_name}'.")
            return None

        return rv
    
    except Exception as e:
        print(f"Error querying SIMBAD: {e}")
        return None


#------------- Try to get a known RV for the star -------------------
#  
#  Assumes the star name is simbad resolvable. 
#  If not found for any reason, defaults to an RV grid +/- 30 km/s
#



rv_guess = get_radial_velocity(hdu[0].header['OBJECT']) 

if rv_guess is None:
    print("WARNING: radial velocity not found on Simbad. Default RV grid may not be ideal for this spectrum.")
    rv_grid = np.arange(-30,30,0.5)
else:
    rv_grid = np.arange(rv_guess-20,rv_guess+20,0.5)


# begin plotting object
if if_plot:
    fig = plt.figure(figsize=(13, 8))

    # Define GridSpec with height ratios
    gs = gridspec.GridSpec(3, 3, height_ratios=[1,1,1], hspace=0.1,wspace=0.1)

    ccf_ax_top = fig.add_subplot(gs[0,:])
    ccf_ax_bot = fig.add_subplot(gs[1,:])
    ccf_ax_center = fig.add_subplot(gs[2,0:2])
    ccf_ax_hist = fig.add_subplot(gs[2,2])

    ccf_ax_top.set_xticks([])
    ccf_ax_bot.set_xticks([])

    ccf_ax_top.set_ylabel("RV [km/s]",fontsize=15)
    ccf_ax_top.yaxis.set_label_coords(-0.05,0)

    ccf_ax_center.set_ylabel(r"V$_{rad}$ [km/s]",fontsize=15)
    ccf_ax_center.set_xlabel("Order",fontsize=15)

    ccf_ax_hist.set_xlabel(r"V$_{rad}$ [km/s]",fontsize=15)
    ccf_ax_hist.set_yticks([])


#------------- Find the order-by-order bulk RV -------------------
#  
# The first 10 orders we skip since they tend to be lower SNR
# Tellurics are corrected for unless they exceed 50 percent of continuum.
# We also build the fits header here for the new extension
#

vrads = []
offset = 0
skip_orders = [0,1,2,3,4,5,6,7,8,9]
new_hdr = fits.Header()

for order in range(0,hdu[1].data['bary_wavelength'].shape[0]):
    
    if order in skip_orders:
        
        #print(f"Skipping Order {order}")
        vrads.append(np.nan)
        new_hdr.set(f'ORDER {order} VRAD', 9999, f'VRAD for Order {order}. Skipped due to low SNR.')
        continue
    
    
    targ_wave = target_hdu[1].data['bary_wavelength'][order,:]
    targ_spec = target_hdu[1].data['spectrum'][order,:]
    targ_cont = target_hdu[1].data['continuum'][order,:]
    targ_tell = target_hdu[1].data['tellurics'][order,:]

    targ_spec = targ_spec/targ_cont
    
    
    wave = hdu[1].data['bary_wavelength'][order,1000:7290-1000]
    spec = hdu[1].data['spectrum'][order,1000:7290-1000]
    cont = hdu[1].data['continuum'][order,1000:7290-1000]
    tell = hdu[1].data['tellurics'][order,1000:7290-1000]
    
    spec = spec/cont/tell
    
    if sum(tell<0.5) > 0:
        
        #print(f"Skipping Order {order}")
        vrads.append(np.nan)
        new_hdr.set(f'ORDER {order} VRAD', 9999, f'VRAD for Order {order}. Skipped due to tellurics.')
        continue
        
    
    cc_val = np.ones(len(rv_grid))

    for i, rv in enumerate(rv_grid):
        
        new_wave = wave - (wave * (rv/const.c.to(u.km/u.s).value))

        new_spec = interp1d(new_wave, spec, bounds_error=False,fill_value=(np.nan,np.nan))(targ_wave)

        cc_val[i] = cross_correlate(new_spec,targ_spec)
    
    vrad = rv_grid[np.argmax(cc_val)]
    ccf_model =VoigtModel() + ConstantModel()
    params = ccf_model.make_params(center=vrad, amplitude=np.max(cc_val), sigma=1)
    params['gamma'].set(value=0.1, vary=True)
    
    try:
        result = ccf_model.fit(cc_val, params, x=rv_grid)
        
        if if_plot:

            if order < 43:
                
                plot_ccf(rv_grid,cc_val,result.best_fit,order,ccf_ax_top,offset)

            elif order == 43:
                offset=0 #reset
                plot_ccf(rv_grid,cc_val,result.best_fit,order,ccf_ax_bot,offset)
                
            else:
                plot_ccf(rv_grid,cc_val,result.best_fit,order,ccf_ax_bot,offset)
     
        
        vrads.append(result.values['center'])
        new_hdr.set(f'ORDER {order} VRAD', result.values['center'], f'VRAD for Order {order}.')
        offset+=1
        
        #print(order, sum(tell<0.99),sum(tell<0.5), result.values['center'])
        
    except ValueError:
        
        print(f"Value error on order {order}")
        vrads.append(np.nan)
        new_hdr.set(f'ORDER {order} VRAD', 9999, f'VRAD for Order {order}. CCF fitting error.')
        pass

vrads = np.array(vrads)
vrad_bulk = np.nanmedian(vrads)
vrad_errs = np.diff(np.nanpercentile(vrads,[16,50,84]))


if if_plot:
    ccf_ax_center.plot(vrads,marker='o',color='k')
    ccf_ax_hist.hist(np.array(vrads),histtype='step',color='k',hatch='//')
    ccf_ax_hist.text(s=rf"""V$_{{rad}}$: {vrad_bulk:0.2f}$^{{+{vrad_errs[0]:.2f}}}_{{-{vrad_errs[1]:.2f}}}$
    No. Orders: {len(vrads[~np.isnan(vrads)])}""",
                    x=0.05,
                    y=0.7,
                    transform = ccf_ax_hist.transAxes,
                    fontsize=12,
                    bbox=dict(facecolor='white', edgecolor='black',boxstyle='square,pad=0.2'))
    
    fig.savefig(f"{out_dir}/plots/{hdu[0].header['OBJECT']}_shift_plot.png")

#------------- Apply the overall RV shift -------------------
#  
#  Uses the median bulk RV to shuft all orders.
#  We also shift the continuum, uncertainty, and tellurics array
#

shift_spec = np.ones(target_hdu[1].data['spectrum'].shape)
shift_cont = np.ones(target_hdu[1].data['continuum'].shape)
shift_error = np.ones(target_hdu[1].data['uncertainty'].shape)
shift_tell = np.ones(target_hdu[1].data['tellurics'].shape)

new_hdr.set('VRAD', vrad_bulk, 'Median radial velocity from all orders.')
new_hdr.set('VRAD_PLUS', vrad_errs[0], 'Upper uncertainty on VRAD.')
new_hdr.set('VRAD_MINUS', vrad_errs[1], 'Lower uncertainty on VRAD.')

for order in range(0,hdu[1].data['bary_wavelength'].shape[0]):
    
    targ_wave = target_hdu[1].data['bary_wavelength'][order,:]
    
    wave = hdu[1].data['bary_wavelength'][order,:]
    spec = hdu[1].data['spectrum'][order,:]
    error = hdu[1].data['uncertainty'][order,:]
    cont = hdu[1].data['continuum'][order,:]
    tell = hdu[1].data['tellurics'][order,:]
    
    new_wave = wave - (wave * (vrad_bulk/const.c.to(u.km/u.s).value))
    new_spec = interp1d(new_wave, spec, bounds_error=False,fill_value=(np.nan,np.nan))(targ_wave)
    new_cont = interp1d(new_wave, cont, bounds_error=False,fill_value=(np.nan,np.nan))(targ_wave)
    new_error = interp1d(new_wave, error, bounds_error=False,fill_value=(np.nan,np.nan))(targ_wave)
    new_tell = interp1d(new_wave, tell, bounds_error=False,fill_value=(np.nan,np.nan))(targ_wave)
    
    shift_spec[order,:] = new_spec
    shift_cont[order,:] = new_cont
    shift_error[order,:] = new_error
    shift_tell[order,:] = new_tell
    
# make data dictionary for HDU conversion
data_dict = {'reg_spectrum':shift_spec,
             'reg_uncertainty':shift_error,
             'reg_continuum':shift_cont,
             'reg_wavelength':target_hdu[1].data['bary_wavelength'],'reg_tellurics':shift_tell}

new_ext = fits.BinTableHDU.from_columns(dict_to_hdu(data_dict),name='registered',header=new_hdr)
hdu.append(new_ext)
hdu.writeto(f"{out_dir}/{new_file_name}")











