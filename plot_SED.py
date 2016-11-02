#traslation of Lucas' routine to plot the SED using X-Shooter data
#usage: python -W ignore plot_SED.py -o 'objname' -r  $Rv
 
from astropy.io import fits
from pylab import *
import numpy as np
import itertools
import argparse
from astropy.convolution import convolve, Box1DKernel
from scipy import interpolate

 
path='/Users/macla/Dropbox/M17/SEDs/'
#objname='M17-B331'
# Distance to the star in parsec
d=1980.
 
parser = argparse.ArgumentParser()
 
parser.add_argument('-o', type=str)
parser.add_argument('-r', type=float)
 
objname = parser.parse_args().o
#Rv = parser.parse_args().r
 
if not objname:
        objname = raw_input('Enter name of the object: ')
#if not Rv:
#        Rv = raw_input('Enter the Rv value: ')
try:
    print objname
#    print 'Rv =', Rv 
except IOError, msg:
    parser.error(str(msg))
 
 
 
############################################################
#
################### DEFINE FUNCTIONS #######################
#
############################################################
 
def get_data1D(f):
    flux = f[0].data
    err = f['ERRS'].data
    qual = f['QUAL'].data
 
    # construct wavelength axis from header info
    wave = []
    for i in range(len(f[0].data)):
      currentWlc = (f[0].header["CRVAL1"] + i * f[0].header["CDELT1"])
      wave.append(currentWlc)
    wave = np.array(wave)
       
    return flux, err, qual, wave
 
def normSpec(flux, wl, err, degreePoly):
    left_x = len(flux)/10
    right_x = len(flux) - left_x
    fitData = flux[left_x:right_x]
    fitWlc = wl[left_x:right_x]
    fit = np.polyfit(fitWlc, fitData, degreePoly)
    p = np.poly1d(fit)
    norm = flux / p(wl)
    error = err/p(wl)
    return norm, wl, error
 
def smoothSpec(bin_size,wave,flux):
    new_w = wave
    new_flux = convolve(flux, Box1DKernel(bin_size))
    return new_w, new_flux
     
       
def get_dered(Rv, npix):
 
  wave_start=300. #start wavelength value in nm
  wave_end=4500. #end wavelength value in nm
  wave=np.linspace(wave_start, wave_end, npix)
 
  #Values from Cardelli et al. 1989 for different Rv (Table 3):
  x_cat_C = np.array([2.78,2.27,1.82,1.43,1.11,.8,.63,.46,.29])
  a_cat_C = np.array([.9530,.9982,1.,.8686,.6800,.4008,.2693,.1615,.0800])
  b_cat_C = np.array([1.9090,1.0495,0.,-.3660,-.6239,-.3679,-.2473,-.1483,-.0734])
  AlAv_cat_C = a_cat_C+b_cat_C/Rv
  wave_cat_C = 1000/x_cat_C
 
  #fit a polynomial of degree 5 through the catalog values
  coeffs_C = np.polyfit(wave_cat_C,AlAv_cat_C,5)
 
  a=coeffs_C
  AlAv=np.poly1d(coeffs_C)(wave)
  AlAvarray=AlAv_cat_C
   
  return wave, a, AlAv, AlAvarray
 
 
############################################################
#
######### INSTRUMENT SETTING PARAMETERS PER STAR ###########
#
# DIT       = DETECTOR INTEGRATION TIME
# UVB/VIS LINES = FACTOR TO CORRECT FOR DIFFERENCES IN SLITLOSS BETWEEN ARMS
# slitloss_uvb/vis = EMPIRICAL FACTOR TO MAKE THE RESPONSE CURVES TO FIT EACH OTHER AND THE PHOTOMETRY
#
######### OBSERVATIONAL PARAMETERS PER STAR ################
#
# Teff = FROM HOFFMEISTER ET AL, OR FORM SPECTRAL CLASSIFICATION
# FIXED Rv=3.3
# LOGg = ASSUME 3.5 (GIANTS) OR 4.0 (DWARFS)
 
############################################################
 
 
if objname =='M17-B111':
    conversion = 1 #2.5
    uvbdit=270.
    visdit=300.
    nirdit=30.
 
    uvblines=1
    vislines=1
 
    uvbdit_std=120.  #Feige 110
    visdit_std=190. #Feige 110
    #visdit_std=20. #HD 175644
    nirdit_std= 90.#HD 175644
     
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
    slitloss_uvb=.7
    slitloss_vis=.6
 
    uvbdit=slitloss_uvb
    visdit=slitloss_vis
     
    temp=45000
    logg=4.5
    tempK=44852
    loggK=3.92
    SpT = ' O4.5 V'
#    lum=0 #from Hoffmeister
    BC = -4.4 #from Allen's astrophysical quantities
 
#    Rv=3.2
elif objname == 'M17-B163':
    conversion = 5
    uvbdit=300.
    visdit=300.
    nirdit=300.
 
    uvblines=1
    vislines=1
 
    uvbdit_std=90.  #EG131
    visdit_std=120 #EG131
    #visdit_std=6. #HD 205130
    nirdit_std= 200.#HD 205130
     
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
#    slitloss_uvb=1
#    slitloss_vis=1
 
#    uvbdit=slitloss_uvb
#    visdit=slitloss_vis
     
    temp=8250
    logg=3.5
    tempK=8200
    loggK=3.5
    SpT = 'kA5 III'
#    lum = 00 #
    BC = -1.1 #from Allen's astrophysical quantities
     
#    Rv=3.3
 
elif objname =='M17-B164':
    conversion = 1 #2.5
    uvbdit=270.
    visdit=300.
    nirdit=30.
 
    uvblines=1
    vislines=1
 
    uvbdit_std=120.  #Feige 110
    visdit_std=190. #Feige 110
    #visdit_std=20. #HD 175644
    nirdit_std= 90.#HD 175644
     
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
    slitloss_uvb=.6
    slitloss_vis=.6
 
    uvbdit=slitloss_uvb
    visdit=slitloss_vis
     
    temp=37000
    logg=4.0
    tempK=36872
    loggK=3.92
    SpT = ' O6 V'
#    lum=0 #from Hoffmeister
    BC = -3.9 #from Allen's astrophysical quantities
 
#    Rv=3.2

elif objname == 'M17-IRS15': #B215
    conversion = 3
    uvbdit=900.
    visdit=400.
    nirdit=25
 
    uvblines=1
    vislines=1
 
    uvbdit_std=300.  #EG153
    visdit_std=30 #EG153
    #visdit_std=5. #HD183430 
    nirdit_std=30. #HD183430 
     
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
    slitloss_uvb=1.
    slitloss_vis=1.1
 
    uvbdit=slitloss_uvb
    visdit=slitloss_vis
     
    temp=25000
    logg=4.0
    tempK=25400
    loggK=3.90
    SpT = ' B0 V'
#    lum=3.917 
    BC = -3.2 #from Allen's astrophysical quantities
#    Rv=3.3

elif objname == 'M17-B243':
    conversion = 1.2
    uvbdit=880.
    visdit=450.
    nirdit=10.0
 
    uvblines=0.6
    vislines=1
     
    uvbdit_std=600.  #EG153
    visdit_std=600. #EG153
    #visdit_std=5. #HD175644
    nirdit_std=5. #HD175644
 
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
    slitloss_uvb=1.6
    slitloss_vis=1
 
    uvbdit=slitloss_uvb
    visdit=slitloss_vis
     
    temp=13000
    logg=4.0
    tempK=13000
    loggK=4.04
    SpT = ' B8 V'

#    lum=2.437 
    BC = -0.8 #from Allen's astrophysical quantities
#    Rv=3.3

elif objname =='M17-B253':
    conversion = 1.6 #2.5
    uvbdit=470.
    visdit=500.
    nirdit=20.
 
    uvblines=1
    vislines=1
 
    uvbdit_std=600.  #Feige 110
    visdit_std=14. #
    #visdit_std=20. #HD 175644
    nirdit_std= 40.#HD 175644
     
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
    slitloss_uvb=1
    slitloss_vis=1
 
    uvbdit=slitloss_uvb
    visdit=slitloss_vis
     
    temp=15000
    logg=3.5
    tempK=15000
    loggK=3.49
    SpT = ' B5 III'
#    lum=0 #from Hoffmeister
    BC = -1.9 #from Allen's astrophysical quantities
 
#    Rv=3.2

elif objname == 'M17-B268':
    conversion = 1.2
    uvbdit=870.
    visdit=450.
    nirdit=50.
 
    uvblines=1
    vislines=1
 
    uvbdit_std=600.  #EG153
    visdit_std=600. #EG153
    #visdit_std=5. #HD175644
    nirdit_std=5. #HD175644
 
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
    slitloss_uvb=1
    slitloss_vis=0.93
 
    uvbdit=slitloss_uvb
    visdit=slitloss_vis
 
    temp=12000
    logg=4.0
    tempK=11900
    loggK=4.04
    SpT = ' A0 V'
#    lum=5.593 
    BC = -1.1 #from Allen's astrophysical quantities
#    Rv = 5
     
elif objname == 'M17-B275':
    conversion = 0.5  # factor to fit the spectrum to the photometry
    uvbdit=685.
    visdit=285.
    nirdit=0.9
 
    uvblines=1 # Factor to fit the arms
    vislines=1
     
     
    uvbdit_std=80.  #EG274
    visdit_std=240. #EG274
    #visdit_std=5. #HD180699
    nirdit_std=5. #HD180699
     
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
    slitloss_uvb=0.667
    slitloss_vis=0.25
 
    uvbdit=slitloss_uvb
    visdit=slitloss_vis
     
# Data from Hoffmeister et al. This two parameters are used to select the Kurucz models.
     
    temp=13000
    logg=3.5
    tempK=13000
    loggK=3.5
    SpT = ' B7 III'
    BC = -1.0 #from Allen's astrophysical quantities
#    lum = 3.16933344665 #from Hoffmeister
 
#    Rv=3.3
 
elif objname == 'M17-B289':
     conversion = 7
     uvbdit=900.
     visdit=400.
     nirdit=40.
 
     uvblines=1
     vislines=1
 
     uvbdit_std=10.  #GD153
     visdit_std=0.5 #GD153
     #visdit_std=5. #HD183430 
     nirdit_std=3. #HD183430 
     
     uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
     visdit=visdit/visdit_std  # science DIT / standard DIT
     nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
     slitloss_uvb=1.7
     slitloss_vis=1.9
 
     uvbdit=slitloss_uvb
     visdit=slitloss_vis
     
     temp=34000
     logg=4.0
     tempK=33879
     loggK=3.92
     SpT = ' O9.5 V'
 #    lum=4.720 #from Hoffmeister
     BC = -3.25 #from Allen's astrophysical quantities
 #    Rv=5.5


elif objname == 'M17-B311':
    conversion = 1 #2.5
    uvbdit=270.
    visdit=300.
    nirdit=30.
 
    uvblines=1
    vislines=1
 
    uvbdit_std=120.  #Feige 110
    visdit_std=190. #Feige 110
    #visdit_std=20. #HD 175644
    nirdit_std= 5.#HD 175644
     
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
    slitloss_uvb=0.6
    slitloss_vis=0.45
 
    uvbdit=slitloss_uvb
    visdit=slitloss_vis
     
    temp=36000
    logg=4.0
    tempK=35874
    loggK=3.92
    SpT = ' O8.5 Vz'
#    lum=0 #from Hoffmeister
    BC = -3.35 #from Allen's astrophysical quantities
 
#    Rv=3.5
 
elif objname == 'M17-B331':
    conversion = 0.5
    uvbdit=870.
    visdit=900.
    nirdit=50.
 
    uvblines=0.15
    vislines=0.15
 
    uvbdit_std=600.  #GD153
    visdit_std=600 #GD153
    #visdit_std=5. #HD095884
    nirdit_std= 30.#HD095884
     
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
#    slitloss_uvb=1
#    slitloss_vis=1
 
#    uvbdit=slitloss_uvb
#    visdit=slitloss_vis
     
    temp=22000 # From Hoffmeister+ 2008, ApJ, 686, 310
    logg=3.5
    tempK=22000
    loggK=3.50
    SpT = ' B2 III'
#    lum=3.657 #from Hoffmeister
    BC = -2.4 #from Allen's astrophysical quantities
#    Rv=3.3

 
elif objname == 'M17-B337':
    conversion = 1 #2.5
    uvbdit=870.
    visdit=450.
    nirdit=50.
 
    uvblines=1
    vislines=1
 
    uvbdit_std=600.  #GD153
    visdit_std=600. #GD153
    #visdit_std=6. #HIP 084690 / HD 156293 
    nirdit_std= 140.#HIP 084690 / HD 156293 
     
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
#    slitloss_uvb=1
#    slitloss_vis=1
 
#    uvbdit=slitloss_uvb
#    visdit=slitloss_vis
     
    temp=15000
    logg=3.5
    tempK=15400
    loggK=4.04
    SpT = ' B5 V'
#    lum=4.103 #from Hoffmeister
    BC = -2.2 #from Allen's astrophysical quantities
#    d=1300.
 
#    Rv=3.3
 
     
else:
    conversion = 1 #2.5
    uvbdit=270.
    visdit=300.
    nirdit=30.
 
    uvblines=1
    vislines=1
 
    uvbdit_std=120.  #Feige 110
    visdit_std=190. #Feige 110
    #visdit_std=20. #HD 175644
    nirdit_std= 5.#HD 175644
     
    uvbdit=uvbdit/uvbdit_std  # science DIT / standard DIT
    visdit=visdit/visdit_std  # science DIT / standard DIT
    nirdit=nirdit/nirdit_std # science DIT / standard DIT
 
    slitloss_uvb=0.6
    slitloss_vis=0.45
 
    uvbdit=slitloss_uvb
    visdit=slitloss_vis
     
    temp=15000
    logg=4.0
    tempK=00000
    loggK=0000
    SpT = '??'
#    lum=4.103 #from Hoffmeister
    BC = -1.1 #from Allen's astrophysical quantities
 
#    Rv=3.3
 
########################    Spectral energy distributions      #############################
#
# This program plots the flux calibrated X-shooter spectrum in three arms, with the following input files / parameters:
# * object name
# * tabulated magnitudes of star (in_f)
# * reduced fluxcalibrated 1d data (nir: flux calibrated using xtellcor)
# * slit losses (3 element array)
# 
#
####################################################################
#
#  Create structures to store data
#
####################################################################
 
if objname == 'M17-B243':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B243/xsh_scired_slit_nod/reflex_end_products/2015-05-08T12:43:59/XSHOO.2013-07-17T03:21:47.437_tpl/M17-B243_nod_SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B243/xsh_scired_slit_nod/reflex_end_products/2015-05-08T12:43:59/XSHOO.2013-07-17T03:21:52.597_tpl/M17-B243_nod_SCI_SLIT_FLUX_MERGE1D_VIS.fits')
 
elif objname == 'M17-B268':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B268/UVB/reflex_end_products/2015-05-07T15:43:02/XSHOO.2013-07-17T04:04:07.577_tpl/M17-B268_nod_SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B268/VIS/reflex_end_products/2015-05-07T16:34:12/XSHOO.2013-07-17T04:04:12.718_tpl/M17-B268_nod_SCI_SLIT_FLUX_MERGE1D_VIS.fits')
 
elif objname == 'M17-B337':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B337/xsh_scired_slit_nod/M17-B337_nod_SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B337/xsh_scired_slit_nod/M17-B337_nod_SCI_SLIT_FLUX_MERGE1D_VIS.fits')
    file_nir = fits.open('/Volumes/Macla_data/X-shooter-data/B337/M17-B337_nod_SCI_SLIT_FLUX_MERGE1D_NIR.fits')
 
 
elif objname == 'M17-B275':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B275/xsh_scired_slit_nod/SCI_SLIT_FLUXCAL1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B275/xsh_scired_slit_nod/SCI_SLIT_FLUXCAL1D_VIS.fits')
 
elif objname == 'M17-IRS15':
#    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/M17-IRS15/xsh_scired_slit_nod/SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/IRS15/M17-IRS15_SCI_SLIT_FLUX_MERGE1D_stareA_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/IRS15/M17-IRS15_SCI_SLIT_FLUX_MERGE1D_VIS.fits')
    file_nir = fits.open('/Volumes/Macla_data/X-shooter-data/IRS15/M17-IRS15_SCI_SLIT_FLUX_MERGE1D_NIR.fits')

 
elif objname == 'M17-B289':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B289/xsh_scired_slit_nod/SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B289/xsh_scired_slit_nod/SCI_SLIT_FLUX_MERGE1D_VIS.fits')
#    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B289/xsh_scired_slit_nod/SCI_SLIT_FLUX_MERGE1D_UVB.fits')
#    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B289/xsh_scired_slit_nod/SCI_SLIT_FLUX_MERGE1D_VIS.fits')
 
elif objname == 'M17-B331':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B331/xsh_scired_slit_nod/SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B331/xsh_scired_slit_nod/SCI_SLIT_FLUX_MERGE1D_VIS.fits')
 
elif objname == 'M17-B163':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/B163/BACK-M17-B163_SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/B163/BACK-M17-B163_SCI_SLIT_FLUX_MERGE1D_VIS.fits')
 
elif objname == 'M17-B311':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B311/xsh_scired_slit_nod/SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B311/xsh_scired_slit_nod/SCI_SLIT_FLUX_MERGE1D_VIS.fits')
 
elif objname == 'M17-B111':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B111/xsh_scired_slit_nod/M17-B111_SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B111/xsh_scired_slit_nod/M17-B111_SCI_SLIT_FLUX_MERGE1D_VIS.fits')
 
elif objname == 'M17-B164':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/B164/M17-B164_SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/B164/M17-B164_SCI_SLIT_FLUX_MERGE1D_VIS.fits')
 
elif objname == 'M17-B253':
    file_uvb = fits.open('/Volumes/Macla_data/X-shooter-data/B253/M17-B253_SCI_SLIT_FLUX_MERGE1D_UVB.fits')
    file_vis = fits.open('/Volumes/Macla_data/X-shooter-data/B253/M17-B253_SCI_SLIT_FLUX_MERGE1D_VIS.fits')
 
 
if objname == 'M17-B164':
    file_nir = fits.open('/Volumes/Macla_data/X-shooter-data/B164/M17-B164_SCI_SLIT_FLUX_MERGE1D_NIR.fits')
elif objname == 'M17-B253':
    file_nir = fits.open('/Volumes/Macla_data/X-shooter-data/B253/M17-B253_SCI_SLIT_FLUX_MERGE1D_NIR.fits')
elif objname == 'M17-B111':
    file_nir = fits.open('/Volumes/Macla_data/X-shooter-data/M17-B111/xsh_scired_slit_nod/M17-B111_SCI_SLIT_FLUX_MERGE1D_NIR.fits')
elif objname == 'M17-B163':
    file_nir = fits.open('/Volumes/Macla_data/X-shooter-data/B163/BACK-M17-B163_SCI_SLIT_FLUX_MERGE1D_NIR.fits')
     
else: 
    file_nir_j = np.genfromtxt('/Volumes/Macla_data/X-shooter-data/'+objname+'/xsh_scired_slit_nod/' + objname + '_tell_corr_NIR_j.dat')
    file_nir_h = np.genfromtxt('/Volumes/Macla_data/X-shooter-data/'+objname+'/xsh_scired_slit_nod/' + objname + '_tell_corr_NIR_h.dat')
    file_nir_k = np.genfromtxt('/Volumes/Macla_data/X-shooter-data/'+objname+'/xsh_scired_slit_nod/' + objname + '_tell_corr_NIR_k.dat')
 
 
flux_uvb, err_uvb, q_uvb, wave_uvb = get_data1D(file_uvb)
flux_vis, err_vis, w_vis, wave_vis = get_data1D(file_vis)

if objname == 'M17-IRS15':
        wave_uvb, flux_uvb = np.loadtxt('/Volumes/Macla_data/X-shooter-data/IRS15/M17-IRS15_SCI_SLIT_FLUX_MERGE1D_stare_UVB.dat').T
 
 
if objname == 'M17-B164' or objname == 'M17-B111' or objname == 'M17-IRS15' or objname == 'M17-B253' or objname == 'M17-B337' or objname == 'M17-B163':
    flux_nir, err_nir, w_nir, wave_nir = get_data1D(file_nir)
    wave_nir_j=wave_nir[np.logical_and(wave_nir >= 1010, wave_nir < 1390)]
    wave_nir_h=wave_nir[np.logical_and(wave_nir >= 1390, wave_nir < 1895)]
    wave_nir_k=wave_nir[wave_nir >= 895]
 
    flux_nir_j=flux_nir[np.logical_and(wave_nir >= 1010, wave_nir < 1390)]
    flux_nir_h=flux_nir[np.logical_and(wave_nir >= 1390, wave_nir < 1895)]
    flux_nir_k=flux_nir[wave_nir >= 895]
else:   
    flux_nir_j = file_nir_j[:,1]
    wave_nir_j = file_nir_j[:,0]
 
    flux_nir_h = file_nir_h[:,1]
    wave_nir_h = file_nir_h[:,0]
 
    flux_nir_k = file_nir_k[:,1]
    wave_nir_k = file_nir_k[:,0]
 
flux_uvb=flux_uvb/uvbdit/uvblines
flux_vis=flux_vis/visdit/vislines
 
 
if objname == 'M17-B311':
    flux_nir_j=flux_nir_j/6.5
    flux_nir_h=flux_nir_h/9.
    flux_nir_k=flux_nir_k/3.
else:
    flux_nir_j=flux_nir_j/nirdit
    flux_nir_h=flux_nir_h/nirdit
    flux_nir_k=flux_nir_k/nirdit
         
w_uvb = np.logical_and(wave_uvb > 300, wave_uvb < 560)
w_vis = np.logical_and(wave_vis >= 560, wave_vis < 1018)
w_nir_j = np.logical_and( wave_nir_j >= 1010, wave_nir_j < 1390)
w_nir_h = np.logical_and(wave_nir_h >= 1390, wave_nir_h < 1895)
w_nir_k = wave_nir_k >= 1895
 
wave_uvb = wave_uvb[w_uvb].tolist()
wave_vis = wave_vis[w_vis].tolist()
wave_nir_j = wave_nir_j[w_nir_j].tolist()
wave_nir_h = wave_nir_h[w_nir_h].tolist()
wave_nir_k = wave_nir_k[w_nir_k].tolist()
 
flux_uvb = flux_uvb[w_uvb].tolist()
flux_vis = flux_vis[w_vis].tolist()
flux_nir_j = flux_nir_j[w_nir_j].tolist()
flux_nir_h = flux_nir_h[w_nir_h].tolist()
flux_nir_k = flux_nir_k[w_nir_k].tolist()
 
wave_all = np.array(wave_uvb + wave_vis + wave_nir_j + wave_nir_h + wave_nir_k)
flux_all = np.array(flux_uvb + flux_vis + flux_nir_j + flux_nir_h + flux_nir_k)*conversion
 
#plt.plot(wave_uvb,flux_uvb, wave_vis,flux_vis, wave_nir_j,flux_nir_j, wave_nir_h,flux_nir_h, wave_nir_k,flux_nir_k)
 
 
filt1 = np.logical_and(wave_all >= 1014, wave_all <=1024)
filt2 = np.logical_and(wave_all >= 1330, wave_all <=1450)
filt3 = np.logical_and(wave_all >= 1825, wave_all <=1940)
filt4 = np.logical_and(wave_all > 300, wave_all <= 340)
 
flux_all[filt1] = 'nan'
flux_all[filt2] = 'nan'
flux_all[filt3] = 'nan'
flux_all[filt4] = 'nan'
 
 
####################################################################
#
#  READ MAGNITUDES FROM FILES (OBJNAME+.DAT)
#
# Flux zero points from Spitzer Telescope Handbook
# http://irsa.ipac.caltech.edu/data/SPITZER/docs/spitzermission/missionoverview/spitzertelescopehandbook/19/
#http://casa.colorado.edu/~ginsbura/filtersets.htm
 
# UBVRI magnitudes from Neckel, T. & Chini, R. 1980
# JHK magnitudes from 2MASS Skrutskie et al. 2006
# 3.6 and 5.8 um from Spitzer GLIMSPE Benjamin et al. 2003
# N and Q mag from Nielbock et al. 2001
#
####################################################################
 
#readcol,objname+'.dat',bands,mags,format=('A,F')
#readcol,objname+'.dat',std_type,format='A',numline=1
#readcol,objname+'.dat',std_band,std_mag,format='A,D',skipline=1
 
f_phot = np.genfromtxt(path+objname+'.dat', dtype=None)
std_band = [f_phot[i][0] for i in range(len(f_phot))]
std_mag = [f_phot[i][1] for i in range(len(f_phot))]
nbands=len(std_band)
f0=[]
std_wave=[]
 
# Write the appropriate (wave,flux) values in (nm, erg cm^-2 A^-1)
for i in range(nbands):
  if std_band[i] == 'U': 
    f0.append(4.266e-9)
    std_wave.append(365.)
  elif std_band[i] == 'B':
    f0.append(6.822e-9)
    std_wave.append(433.)
  elif std_band[i] == 'V':
    f0.append(3.802e-9)
    mag_v = std_mag[i]
    std_wave.append(550.)
  elif std_band[i] == 'R': 
    mag_r = std_mag[i]
    f0.append(1.738e-9)
    std_wave.append(700.)
  elif std_band[i] == 'I':
    f0.append(8.318e-10)
    std_wave.append(900.)
  elif std_band[i] == 'J':
    f0.append(3.129e-10)
    std_wave.append(1235.)
  elif std_band[i] == 'H':
    f0.append(1.133e-10)
    std_wave.append(1662.)
  elif std_band[i] == 'K':
    mag_k = std_mag[i]
    f0.append(4.283e-11)
    std_wave.append(2159.)
  elif std_band[i] == '36um':
    f0.append(6.601e-12)
    std_wave.append(3550.)
  elif std_band[i] == '45um':
    f0.append(2.545968e-12)
    std_wave.append(4500.) 
  elif std_band[i] == '58um':
    f0.append(1.064e-12)
    std_wave.append(5731.)
  elif std_band[i] == '80um':
    f0.append(3.004e-13)
    std_wave.append(8000.)
   
# include Spitzer Glimpse:
#nbands=nbands+2
#f0.extend([6.601e-12,1.064e-12])
#std_mag.extend([6.662,5.346])
#std_wave.extend([3550.,5731.])
 
 
std_flux = []
for i in range(len(std_mag)):
  std_flux.append(f0[i]*(10**(-0.4*std_mag[i])))
   
#Append N and Q bands Nielbock et al. 2001 and SPITZER
# conversion to erg s-1 cm-2 A-1 in http://localhost:8888/notebooks/mJy_to_erg_s-1_cm-2_Hz-1.ipynb
 
if objname == 'M17-B275':
  std_wave.extend([10600,20000])
  std_flux.extend([5.06946920612e-15,5.9208999e-15])
elif objname == 'M17-B289':
  std_flux[8] = 1.900E-14 #3.6 micron IRAC (Band 1) magnitude from visier
  std_flux[9] = 6.858E-15 #5.8 micron  IRAC (Band 1) magnitude from visier  
  std_wave.extend([10600,20000])
  std_flux.extend([6.67035421858e-14,6.04081686e-14])
elif objname == 'M17-B331':
  std_wave.extend([10600,20000])
  std_flux.extend([2.29460185119e-14,2.274674835e-14])
elif objname == 'M17-B337':
  std_flux[5] = 3.317E-15 #3.6 micron IRAC (Band 1) magnitude from visier
  std_flux[6] = 5.181E-15 #5.8 micron IRAC (Band 1) magnitude from visier 
  std_wave.extend([10600,20000])
  std_flux.extend([1.73429209683e-15,2.43581325e-15])
elif objname == 'M17-IRS15':
  std_flux[4] = 8.293E-15 #3.6 micron IRAC (Band 1) magnitude from visier
  std_flux[5] = 2.960E-15 #5.8 micron IRAC (Band 1) magnitude from visier 
  std_wave.extend([10600,20000])
  std_flux.extend([8.72482331791e-14,6.88023558e-14])
elif objname == 'M17-B311':
  std_wave.extend([10600,20000])
  std_flux.extend([1.73429209683e-14,1.61138415e-14])
   
 
std_flux = np.array(std_flux)
std_wave = np.array(std_wave) 
 
####################################################################
#
#  Overplot Kurucz model
#
####################################################################
 
 
kurucz=fits.open(path+'Castelli-Kurucz/ckp00_'+str(temp)+'.fits')
 
wave_k = kurucz[1].data['WAVELENGTH']/10.
npix_k=len(wave_k)
 
if logg == 4.: #Fluxes tabulated in units of erg/s/cm^2/A
  flux_k = kurucz[1].data['g40'] #main sequence 
elif logg == 3.5:
  flux_k = kurucz[1].data['g35'] #giant
elif logg == 3.0:
  flux_k = kurucz[1].data['g30'] #
elif logg == 2.5:
  flux_k = kurucz[1].data['g25'] #
elif logg == 2.0:
  flux_k = kurucz[1].data['g20'] #
elif logg == 4.5:
  flux_k = kurucz[1].data['g45'] #
elif logg == 5:
  flux_k = kurucz[1].data['g50'] #giant
   
 
####################################################################
####################################################################
##############         Compute Dereddened flux        ##############
####################################################################
####################################################################
 
 
print ' ------------- calculating Av ------------- \n'

dist=d*3.086e18 #cm
Lsun = 3.84e33 #erg/s

Av_array = []
slopecoeff_array=[]
scalecoeff_array=[]
std_flux_der_array=[]
chi_sq_array=[]
scalefactor_array=[]
rstar_array=[]
Lum_array=[]
chi_sq_notred_array = []
Mo_array = []

# Limit the magnitudes to the range covered by cardelli or to avoid the NIR excess in B243, B268, and B275
if objname == 'M17-B275' or objname == 'M17-B268'  or objname == 'M17-B243':
    max_std = 5
else:
    max_std = len(std_wave[std_wave < 5000])


Rv_array = np.arange(1.9,5.5,0.2) 
for Rv in Rv_array:
    
    print 'Rv = ',Rv
 
    wave, a, AlAv, AlAvarray = get_dered(Rv, len(flux_all))
 
    AlAv=np.poly1d(a)(wave_all)
    AlAv_phot=np.poly1d(a)(std_wave[0:max_std])
 
    ###################
    #obtain best-fit Av
    ###################
 
    Avmin=3
    Avmax=20
    ngrid=1000
    Avgrid=np.linspace(Avmin, Avmax, ngrid)
 
    # define the photospheric domain
    if objname == 'M17-B163':
        ph_k = np.logical_and(wave_k > 900, wave_k < 1000)
        ph = np.logical_and(wave_all > 900, wave_all < 1000)
    elif objname == 'M17-B337':
        ph_k = np.logical_and(wave_k > 600, wave_k < 900)
        ph = np.logical_and(wave_all > 600, wave_all < 900)
    else:
        ph_k = np.logical_and(wave_k > 400, wave_k < 820)
        ph = np.logical_and(wave_all > 400, wave_all < 820)
 
    #plt.loglog(wave_all[ph], flux_der[ph]*10*wave_all[ph])
    #plt.show()
 
 
    coeff_k=np.polyfit(np.log10(wave_k[ph_k]),np.log10(wave_k[ph_k]*10*flux_k[ph_k]),1)
 
    slope_k=coeff_k[0]
    scale_k=coeff_k[1]  # They coincide with Lucas' results
 
 
    # find best-fit Av
    slopegrid=[]
    scalegrid=[]
 
    for i in range(ngrid):
      flux_der = (flux_all*10**(.4*AlAv*Avgrid[i]))
      idx = np.isfinite(np.log10(wave_all[ph])) & np.isfinite(np.log10(wave_all[ph]*10*flux_der[ph]))
      coeff=np.polyfit(np.log10(wave_all[ph])[idx],np.log10(wave_all[ph]*10*flux_der[ph])[idx],1)
 
      slopegrid.append(coeff[0])
      scalegrid.append(coeff[1])
 
    diff = abs(slopegrid-slope_k).tolist()
    ibest = diff.index(min(diff))
    bestslopecoeff_temp=slopegrid[ibest]
    bestscalecoeff_temp=scalegrid[ibest]
    Avbest_temp=Avgrid[ibest]
    
    print 'Temporary best Av = ', Avbest_temp
    
    slopecoeff_array.append(bestslopecoeff_temp)
    scalecoeff_array.append(bestscalecoeff_temp)
    Av_array.append(Avbest_temp)
            
    std_flux_der_temp=std_flux[0:max_std]*10**(.4*AlAv_phot*Avbest_temp)
    
    std_flux_der_array.append(std_flux_der_temp)

####################################################################
#
# find corresponding radius by measuring the distance between tne 
# kurucz model and the data in a region around 500 nm
#
####################################################################
 
# scale flux to distance.
#To convert to observed flux at Earth, multiply by a factor of (R/D)^2
#where R is the stellar radius, and D is the distance to Earth. 
#http://www.stsci.edu/hst/observatory/crds/k93models.html

 
 
    vflux=np.median(np.poly1d([bestslopecoeff_temp,bestscalecoeff_temp])(np.log10([500.,510.,520.,530.,540.,550.,560.])))
    vflux_k=np.median(np.poly1d(coeff_k)(np.log10([500.,510.,520.,530.,540.,550.,560.])))
    
    scalefactor_temp = 10**(vflux-vflux_k)
    scalefactor_array.append(scalefactor_temp) # at lambda=500 nm 
    
    rbest_temp=sqrt(scalefactor_temp)*dist/6.955e10 # in R sun (Rsun = 6.955e10cm)

    rstar_array.append(rbest_temp)
    
    if objname == 'M17-B337':
        Ak_t = 0.056*Avbest_temp #from Chini et al 1998, A&A, 329, 161
        Mk = mag_k + 5 - 5*np.log10(d)
        Ko_t = Mk - Ak_t
        Mbol_t = Ko_t + BC
        L_t =  10**(0.4*(88.72 - Mbol_t))
        lum_temp = np.log10(L_t/Lsun)
        Mo_array.append(Ko_t)
    elif objname == 'M17-B163':
        Ar_t = 0.751*Avbest_temp #from Cardelli 1998
        Mr = mag_r + 5 - 5*np.log10(d)
        Ro_t = Mr - Ar_t
        Mbol_t = Ro_t + BC
        L_t =  10**(0.4*(88.72 - Mbol_t))
        lum_temp = np.log10(L_t/Lsun)
        Mo_array.append(Ro_t)
    else:
        Mv = mag_v + 5 - 5*np.log10(d)
        Mo_t = Mv - Avbest_temp
        Mbol_t = Mo_t + BC
        L_t = 10**(0.4*(88.72 - Mbol_t))
        lum_temp = np.log10(L_t/Lsun)
        Mo_array.append(Mo_t)
        
    Lum_array.append(lum_temp)
    
    
    # Combine lists into list of tuples
    points = zip(wave_k, flux_k)

    # Sort list of tuples by x-value
    points = sorted(points, key=lambda point: point[0])

    # Split list of tuples into two list of x values any y values
    wave_k_sort,flux_k_sort = [],[]
    for i in points:
        wave_k_sort.append(i[0])
        flux_k_sort.append(i[1])
    wave_k_sort = np.array(wave_k_sort)
    flux_k_sort = np.array(flux_k_sort)
    
    # interpolate model in order to have it at the same resolution of the spectra
    funct = interpolate.interp1d(wave_k,wave_k*10*scalefactor_temp*flux_k, bounds_error=False)
    fnew_k = funct(std_wave[0:max_std])
    
    # plt.plot(np.log10(wave_k), np.log10(wave_k*10*scalefactor_temp*flux_k))
    # plt.plot(np.log10(w_std), fnew_k, 'o')
    # plt.show()

    chi_sq_array.append(np.sum((std_wave[0:max_std]*10*std_flux_der_temp - fnew_k)**2/(std_wave[0:max_std]*10*std_flux_der_temp*0.01)**2)/(len(std_wave[0:max_std])-2))
    chi_sq_notred_array.append(np.sum((std_wave[0:max_std]*10*std_flux_der_temp - fnew_k)**2/(std_wave[0:max_std]*10*std_flux_der_temp*0.01)**2))
    
    print 'chi^2 = ', np.sum((std_wave[0:max_std]*10*std_flux_der_temp - fnew_k)**2/(std_wave[0:max_std]*10*std_flux_der_temp*0.01)**2)/(len(std_wave[0:max_std]-2))
    
    # # the kurucz model
    # plt.loglog(wave_k,wave_k*10*scalefactor_temp*flux_k, 'k--', alpha=0.5, linewidth=1)
    # # the magnitudes
    # plt.plot(std_wave,std_wave*10*std_flux, 's', color='0.35', markersize=5, linewidth=3)# fillstyle='none',
    # #The kurucz magnitudes
    # plt.plot(std_wave[0:max_std],fnew_k, 'go', markersize=9, linewidth=3)# fillstyle='none',
    # #The dereddened magnitudes
    # plt.plot(std_wave[0:max_std],std_wave[0:max_std]*10*std_flux_der_temp, 'd',color='royalblue', markersize=5, linewidth=2)#,fillstyle='none'
    # yrange=[10**-15,10**-6]
    # plt.xlim(250,23000)
    # plt.ylim(yrange)
    # plt.xlabel(r'$\lambda$ $\mathrm{(\mu m)}$', fontsize=14)
    # plt.ylabel(r'$\lambda F_{\lambda}$ ($\mathrm{erg}$ $\mathrm{s^{-1} cm^{-2}}$)', fontsize=14)
    # plt.xticks([500,1000,2000, 5000,10000, 20000], ['0.5','1.0','2.0', '5.0', '10.0', '20.0'], fontsize = 12)
    # plt.yticks(fontsize = 12)
    # plt.show()

    
#diff_best = abs(slopecoeff_array - slope_k).tolist()
#index_best = diff_best.index(min(diff_best))
index_best = chi_sq_array.index(min(chi_sq_array))
bestslopecoeff = slopecoeff_array[index_best]
bestscalecoeff = scalecoeff_array[index_best]
Avbest = Av_array[index_best]
Rv = Rv_array[index_best]
std_flux_der = std_flux_der_array[index_best]
scalefactor = scalefactor_array[index_best]
rstar = rstar_array[index_best]
lum = Lum_array[index_best]
Mo = Mo_array[index_best]

out = file('/Users/macla/Dropbox/M17/SEDs/chisq_fit/chi_sq_'+objname+'.dat', 'w')
out.write('# Best parameters for ' + objname +' (Number of points fitted: '+str(len(std_wave[0:max_std]))+')'+ '\n')
out.write('# Rv = '+str(Rv)+'\n')
out.write('# Av = '+str(Avbest)+'\n')
out.write('# Rstar = '+str(rstar)+'\n')
out.write('# Lum = '+str(lum)+'\n')
out.write('# M_0 = '+str(Mo)+'\n')
out.write('# Rv \t \t  chi_red^2 \t \t chi^2 \t\t Av \t \t  Rstar \t \t Lum \t\t M_0 \n')
for i in range(len(Rv_array)):
    out.write(str(Rv_array[i])+'\t'+str(chi_sq_array[i])+'\t'+str(chi_sq_notred_array[i])+'\t'+str(Av_array[i])+'\t'+str(rstar_array[i])+'\t'+str(Lum_array[i])+ '\t' +str(Mo_array[i]) + '\n')

out.close()

#plt.plot(Rv_array,chi_sq_array,'o')
#plt.savefig('/Users/macla/Dropbox/M17/SEDs/chisq_fit/Rv_vs_chi_sq_'+objname+'.pdf', bbox_inches='tight')
# plt.show()



wave, a, AlAv, AlAvarray = get_dered(Rv, len(flux_all))

AlAv=np.poly1d(a)(wave_all)

# rbest=sqrt(scalefactor)*dist/6.955e10 # in R sun (Rsun = 6.955e10cm)
#
# rstar = rbest


print '------------- \n best fit Av=', Avbest
 

####################################################################
#
# Calculate the luminosity based on V mag, BC (from astrophysical 
# quantities), and Av 
# The 88.72 is a cnt = 2.5*log(Lsun) + Mbol_sun
# Lsun = 3.84e33 erg/s
# Mbol_sun = 4.76
####################################################################
# print '\n ------------- calculating luminosity  ------------- \n'
#
# if objname == 'M17-B337':
#     print 'k =', mag_k
#     Ak = 0.056*Avbest #from Chini et al 1998, A&A, 329, 161
#     Mk = mag_k + 5 - 5*np.log10(d)
#     print 'Mk:', mag_k, '+ 5 - 5log10(d) =', Mk
#     print 'Ak =', Ak
#     Ko = Mk - Ak
#     print 'M_0:', Mk, '-', Ak, '=',  Ko
#     Mbol = Ko + BC
#     print 'Mbol:', Ko, '+', BC, '=', Mbol
#     L =  10**(0.4*(88.72 - Mbol))
#     lum = np.log10(L/Lsun)
#     print 'log(L/Lsun) =', lum
# elif objname == 'M17-B163':
#     print 'r =', mag_r
#     Ar = 0.751*Avbest #from Cardelli 1998
#     Mr = mag_r + 5 - 5*np.log10(d)
#     print 'Mr:', mag_r, '+ 5 - 5log10(d) =', Mr
#     print 'Ar =', Ar
#     Ro = Mr - Ar
#     print 'M_0:', Mr, '-', Ar, '=',  Ro
#     Mbol = Ro + BC
#     print 'Mbol:', Ro, '+', BC, '=', Mbol
#     L =  10**(0.4*(88.72 - Mbol))
#     lum = np.log10(L/Lsun)
#     print 'log(L/Lsun) =', lum
# else:
#     print 'v =', mag_v
#     Mv = mag_v + 5 - 5*np.log10(d)
#     print 'Mv:', mag_v, '+ 5 - 5log10(d) =', Mv
#     print 'Av =', Avbest
#     Mo = Mv - Avbest
#     print 'M_0:', Mv, '-', Avbest, '=',  Mo
#     Mbol = Mo + BC
#     print 'Mbol:', Mo, '+', BC, '=', Mbol
#     L = 10**(0.4*(88.72 - Mbol))
#     print 'L: 10^0.4(88.72 -',Mbol, ') = ', L
#     lum = np.log10(L/Lsun)
#     print 'log(L/Lsun) =', lum
#
#  
 
 

#####################
print '\n ------------- calculating radius ------------- \n'
print 'best fit Rstar=',rstar

flux_der=flux_all*10**(.4*AlAv*Avbest)
 
 
# Smooth spectra
wave_all_smooth, flux_all_smooth = smoothSpec(20, wave_all, flux_all) 
 
wave_all_smooth, flux_der_smooth = smoothSpec(20, wave_all, flux_der) 
 
print 'scalefactor ', Rv, '=', scalefactor
 
####################################################################
#
#  Write files
#
####################################################################
 
#out = open(path+objname+'_SED_Rv'+str(Rv)+'_spectra.dat', 'w')
#print '\n ------------- Writing file ------------- \n'
#print ('File written in: '+path+objname+'_SED_Rv'+str(Rv)+'_spectra.dat')
#out.write('##'+ objname + ': Luminosity= '+str(lum)+' radius= '+str(rstar)+' Av= '+str(Avbest)+' Rv= '+str(Rv)+'\n')
#out.write('## Wavelength \t Flux (erg s-1 cm-2 A-1) \t dereddened flux (erg s-1 cm-2 A-1) \n')
#for i in range(len(wave_all_smooth)):
#   out.write(str(wave_all_smooth[i])+'\t'+str(flux_all_smooth[i])+'\t'+str(flux_der_smooth[i])+'\n')
 
#out.close()
####################################################################
#
#  Plot everything
#
####################################################################
 
 
print '\n ------------- Plotting spectrum ------------- \n'
#outname = path+'00_newSEDs/'+objname+'_T_'+str(tempK)+'_Rv_' + str(Rv) + '_logg_' + str(loggK) + '_SED.pdf'
outname = '/Users/macla/Dropbox/M17/SEDs/chisq_fit/'+objname+'_T_'+str(temp)+'_Rv_' + str(Rv) + '_logg_' + str(logg) + '_SED.pdf'
#print 'spectrum plotted in: ', outname+'\n'
 
#yrange=[min(wave_all_smooth*10*flux_all_smooth),10**-6]
#if objname == 'M17-B111':
#    yrange=[10**-11,10**-6]
#elif objname == 'M17-B243' or objname == 'M17-B268' or objname == 'M17-B163':
#    yrange=[10**-14,10**-8]
#elif objname == 'M17-B289' or objname == 'M17-B311':
#    yrange=[10**-13,10**-6]
#elif objname == 'M17-B331':
#    yrange=[10**-14,10**-7]
#else:
#    yrange=[10**-13,10**-7]
yrange=[10**-15,10**-6]
     
#if objname == 'M17-IRS15':
#  for i in range(len(wave_all_smooth)):
#    if wave_all_smooth[i] < 550 or wave_all_smooth[i] > 1020:
#      wave_all_smooth[i] = 'nan'

print objname, ' Photometric Flux:'
for w,f in zip(std_wave,std_flux):
    print w, f
    
# Plot the spectrum
fig = figure(figsize=(5,4))
ax = fig.add_subplot(111)
if objname == 'M17-B163':
    ax.loglog(wave_all_smooth[wave_all_smooth>700],wave_all_smooth[wave_all_smooth>700]*10*flux_all_smooth[wave_all_smooth>700], '0.35')
elif objname == 'M17-B337':
    ax.loglog(wave_all_smooth[wave_all_smooth>600],wave_all_smooth[wave_all_smooth>600]*10*flux_all_smooth[wave_all_smooth>600], '0.35')
elif objname == 'M17-B331':
    ax.loglog(wave_all_smooth[wave_all_smooth>530],wave_all_smooth[wave_all_smooth>530]*10*flux_all_smooth[wave_all_smooth>530], '0.35')
else:
    ax.loglog(wave_all_smooth,wave_all_smooth*10*flux_all_smooth, '0.35')

ax.set_xlim(250,23000)
ax.set_ylim(yrange)
ax.set_xlabel(r'$\lambda$ $\mathrm{(\mu m)}$', fontsize=14)
ax.set_ylabel(r'$\lambda F_{\lambda}$ ($\mathrm{erg}$ $\mathrm{s^{-1} cm^{-2}}$)', fontsize=14)
# the dereddened spectrum
if objname == 'M17-B163':
    ax.loglog(wave_all_smooth[wave_all_smooth>700],wave_all_smooth[wave_all_smooth>700]*10*flux_der_smooth[wave_all_smooth>700], 'royalblue')
elif objname == 'M17-B337':
    ax.loglog(wave_all_smooth[wave_all_smooth>600],wave_all_smooth[wave_all_smooth>600]*10*flux_der_smooth[wave_all_smooth>600], 'royalblue')
elif objname == 'M17-B331':
    ax.loglog(wave_all_smooth[wave_all_smooth>530],wave_all_smooth[wave_all_smooth>530]*10*flux_der_smooth[wave_all_smooth>530], 'royalblue')
else:
    ax.loglog(wave_all_smooth,wave_all_smooth*10*flux_der_smooth, 'royalblue')

# the kurucz model
ax.loglog(wave_k,wave_k*10*scalefactor*flux_k, 'k--', alpha=0.5, linewidth=1)
 
# the magnitudes
ax.plot(std_wave,std_wave*10*std_flux, 's', color='0.35', markersize=5, linewidth=3)# fillstyle='none',
    
plt.xticks([500,1000,2000, 5000,10000, 20000], ['0.5','1.0','2.0', '5.0', '10.0', '20.0'], fontsize = 12)
plt.yticks(fontsize = 12)
 

ax.plot(std_wave[0:max_std],std_wave[0:max_std]*10*std_flux_der, 'd',color='royalblue', markersize=5, linewidth=2)#, fillstyle='none'
 
#The box with radius and Av
textstr = '$\mathrm{%3s}$ \n $T_{\mathrm{eff}}=%5d \mathrm{K}$ \n $\log{g}=%.2f$ \n  \n\n\n \n $\mathrm{R_V=%.1f}$ \n$\mathrm{A_V=%.2f}$ \n $\mathrm{R_{\star}=%.2f R_{\odot}}$ \n $\mathrm{log(L/L_{\odot})=%.2f}$'%(SpT, tempK, loggK, Rv, Avbest, rstar, lum)

ax.text(0.5, 0.5,textstr, ha='left', va='center', fontsize=15, transform=ax.transAxes)
#ax.set_title(objname, fontsize=20)
plt.gcf().subplots_adjust(left=0.2, bottom=0.1)
#fig.savefig(outname, bbox_inches='tight')
#fig.savefig('/Users/macla/Desktop/try.pdf', bbox_inches='tight')
#plt.show()
print '\n ------------- done ------------- \n'
