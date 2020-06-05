import numpy as np
import sip_tpv as st 
from astropy.io import fits

np.warnings.filterwarnings('ignore') # ignore version warnings

fits_image_filename = './test-fits/test-pv.fits'

hdul = fits.open(fits_image_filename)  # open a FITS file

hdul.info() # headers info

hdr = hdul[1].header   # ImageHDU

siphdr = st.pv_to_sip(hdr)  #tpv -> sip

print(repr(hdr))