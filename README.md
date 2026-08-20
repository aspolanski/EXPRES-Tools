# EXPRES-Tools
Tools for observing and reduction with EXPRES @ LDT
====================

## Target List Generator:

To generate a list with requests you will need a CSV file that defines the Star name, RA/Dec, Vmag, expected exposure times, number of exposures, number of exposure meter counts, optional comments, and whether it is an EPRV target. An example can be found [here](https://docs.google.com/spreadsheets/d/e/2PACX-1vRBx1Q26pa51QyiDdIfm-f0kFNy1WNvIBXu73HvWoPqu7Q9luI9av4UNeQ1Id_nfTWrEQv6VHy0KltG/pubhtml?gid=1181914448&single=true "EXPRES Request Example")

```python
from expres_utils import TargetList
import pandas as pd
targets = pd.read_csv('expres_requests.csv',comment='#')
  
tl = TargetList(targets=targets)
tl.make_targetlist(target_list_name='example_requests')
  ```

## Spectral Registration:

`register.py` is a self-contained routine to take reduced EXPRES spectra and find the bulk radial velocity of the star by cross-correlating with a target spectrum, in a typical case would be a LOST Solar spectrum. A new fits extension containing the shift spectra is made and a FITS file is saved. Command useage is given below with the optional diagnostic plot. 

```python
python register.py -h
options:
  -h, --help         show this help message and exit
  --in_file IN_FILE  The input FITS file.
  --out_dir OUT_DIR  Directory of the registered spectra.
  --plot IF_PLOT     Whether to make shift plots.
  --target TARGET    FITS file to register to (defaults to Solar spectrum).
  ```

<img src="https://github.com/aspolanski/EXPRES-Tools/blob/main/HD19994_shift_plot.png" width="650" height="400" />


## Color Correction:

The `color_correction.py` routine corrects EXPRES spectra for chromatic variations due to changes in airmass. Especially useful for Rossiter-McLaughlin datasets where the CCF contrast is highly
sensitive to such effects. In summary, each spectrum is compared to the spectrum of a rapidly rotating B star, which estimates the throughput of the optical system. The differences are corrected for using a high order polynomial. 

```
python color_correction.py -h
options:
  -h, --help            show this help message and exit
  --in_dir IN_DIR       The input FITS file.
  --out_dir OUT_DIR     Directory of the registered spectra.
  --plot IF_PLOT        Make plots (not implemented)
  --target-temp TARGET_TEMP
                        Effective temperature of the target star (K)
```

## CCF Fitting Window:

`ccf_width.py` takes the coadded EXPRES CCF of one or more exposures, fits the CCF function and returns a recommended fitting window to derive the RVs (default is 1.5 times the FWHM).

```
python ccf_width.py -h
options:
  -h, --help            show this help message and exit
  --ccf_dir CCF_DIR     Directory containing CCF fits files. (format = /*/*/)
  --fwhm_factor FWHM_FACTOR
                        Fraction of the FWHM to use as the fitting window.
  --plot                Generate plots.
```
