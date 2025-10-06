#!/usr/local/bin/python3

# Deep2 method:
# see https://ui.adsabs.harvard.edu/abs/2013ApJS..208....5N/abstract
# toward the end

from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.time
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
import os.path
import pickle
import sys

# this should be the part a new user needs to edit
prefs = {'PYZE_DIR':r'C:\PHY199',
         'version':'1.1.6'}
##### end user-editable part ##########################
prefs['TELLURIC_FILE'] = prefs['PYZE_DIR']+'/Telluric/MKtransmission.fits'
prefs['THRUPUT_FILE'] = prefs['PYZE_DIR']+'/Templates/thruput.pickle'
mycheby = None # place to store Chebyshev model of thruput

class Spec:
    def __init__(self, wavelength, flux, ivar):
        self.W = wavelength
        self.I = flux
        self.IVAR = ivar
        return
    
class PypeItSpec:
    def __init__(self, filename):
        hdulist = fits.open(filename)  # open the FITS file
        #print(hdulist)
        data = hdulist[1].data
        hdr = hdulist[0].header # yes, the keywords I need are in the PHDU
        # sanity check: print('Wvlngth range: %f-%f\n' % (data['wave'][0],data['wave'][-1]))
        # store the raw data
        self.raw = Spec(data['wave'],data['flux'],data['ivar'])
        # hack for one weird spectrum
        #hackmask = (self.raw.W > 7150) * (self.raw.W < 7610)
        #self.raw.IVAR[hackmask] = 1e-6
        # hack for another weird spectrum
        #self.raw.IVAR = (self.raw.IVAR**-1+25)**-1
        self.airmass = hdr['AIRMASS']
        self.ra = hdr['RA_OBJ']
        self.dec = hdr['DEC_OBJ']
        self.longitude = hdr['LON-OBS']
        self.latitude = hdr['LAT-OBS']
        self.altitude = hdr['ALT-OBS']
        self.mjd = hdr['MJD']
        hdulist.close()
        self.clean()
        self.undo_telluric_abs()
        self.remove_largescale()
        self.barycentric_correction()
        # barycorr is not usually needed b/c PypeIt already does this at
        # the per-exposure level. However, it's computed and stored here
        # in case some utility function wants to access it (eg to back out
        # duplicative corrections)
        return
    def undo_telluric_abs(self):
        tell_abs = get_telluric_absorption(prefs['TELLURIC_FILE'],self.airmass)
        # telluric absorption model only goes to 9128 A. If the observations
        # go redder, let transmission=1 there
        tell_abs_resamp = griddata(tell_abs.W,tell_abs.I,self.W,fill_value=1.0)
        self.I /= tell_abs_resamp
        self.IVAR *= tell_abs_resamp**2
        return
    def barycentric_correction(self): # just computes this factor
        loc = EarthLocation.from_geodetic(lat=self.latitude*u.deg,lon=self.longitude*u.deg, height=self.altitude*u.m) # could also read the TELESCOP keyword then use EarthLocation.of_site(telescop) 
        sc = SkyCoord(ra=self.ra*u.deg,dec=self.dec*u.deg)
        obstime = astropy.time.Time(self.mjd,format='mjd')
        barycorr = sc.radial_velocity_correction(obstime=obstime,location=loc)
        # https://docs.astropy.org/en/stable/coordinates/velocities.html
        # says to ADD this correction to an observed velocity to get barycentric velocity
        # so I will convert it to a stretch factor that augments the observed s.f.
        self.barycorr_factor = 1+barycorr.to(u.km/u.s).value/299792
        #print('Barycorr stretch factor:',self.barycorr_factor)
        return
    def clean(self):
        # all masks are true if "good"

        # (1) bad skyline at 5580 AA leaves severe artifacts not captured by
        # the IVAR array, so just delete it.
        # 4/12/24: maybe we don't need this anymore? try leaving it out
        #skymask = np.logical_or(self.raw.W<5575,self.raw.W>5583)
        # (this then anded w/goodmask below)
        
        # (2) very negative points are spurious
        goodmask = self.raw.I>-50
        # (3) mask off the last ~15 AA on red edge: usually bad
        goodmask[-40:] = False
        # (4) very positive points can be bad too, but be conservative to
        # not delete strong emission lines!
        goodmask = np.logical_and(goodmask,self.raw.I<20000)
        # in rare cases we have points with wave=flux=ivar=0. Remove those
        goodmask = np.logical_and(goodmask,self.raw.W>0.0)
        self.W = self.raw.W[goodmask] # wvlngth
        self.I = self.raw.I[goodmask] # intensity
        self.IVAR = self.raw.IVAR[goodmask] # ivar
        return
    
    # original divide continuum version
    #def remove_largescale(self,kernel_size=451,hole_radius=15):
     #   kernel = np.ones(kernel_size) / (kernel_size-2*hole_radius-1)
      #  ctr = int(kernel_size/2)
       # kernel[ctr-hole_radius:ctr+hole_radius+1] = 0.0
        #unsharp = np.correlate(self.I,kernel,mode='same')
        # correct for boundary effect
        #boundary_effect = np.correlate(np.ones(len(self.I)),kernel,mode='same')
        #unsharp /= boundary_effect
        #self.U = Spec(self.W,self.I/unsharp,self.IVAR*unsharp**2)
        #return
    
    def remove_largescale(self,kernel_size=451,hole_radius=15): 
    # Removes large-scale variations in the spectrum. 
        kernel = np.ones(kernel_size) / (kernel_size-2*hole_radius-1) 
        ctr = int(kernel_size/2) 
        kernel[ctr-hole_radius:ctr+hole_radius+1] = 0.0 
        unsharp = np.correlate(self.I,kernel,mode='same') 
        self.U = Spec(self.W,self.I-unsharp,self.IVAR) 
        return
    
    
class TemplateSpec:
    def __init__(self, filename):
        self.name = os.path.basename(filename)
        hdulist = fits.open(filename)  # open the FITS file
        data = hdulist[0].data
        I = data[0,:]
        valid = I>0
        self.I = I[valid]
        # wvlngth is stored as an equation, not array
        log_start_wvlngth = hdulist[0].header['COEFF0']
        log_delta_wvlngth = hdulist[0].header['COEFF1']
        log_wvlngth = np.arange(len(I))*log_delta_wvlngth+log_start_wvlngth
        self.W = 10**log_wvlngth[valid]
        self.dlambda = 10**log_delta_wvlngth
        self.remove_largescale()
        return
    def remove_largescale(self,kernel_size=151,hole_radius=5):
        kernel = np.ones(kernel_size) / (kernel_size-2*hole_radius-1)
        ctr = int(kernel_size/2)
        kernel[ctr-hole_radius:ctr+hole_radius+1] = 0.0
        self.U = self.I/np.correlate(self.I,kernel,mode='same')
        return
    
def show_smoothed(wvlngth,I,whichplot=None):
    kernel_size = 8
#    kernel = np.ones(kernel_size) / kernel_size
#    smoothed = np.correlate(I,kernel,mode='same')
    smoothed = gaussian_filter(I,kernel_size,mode='reflect')
    if whichplot:
        whichplot.plot(wvlngth,smoothed,linewidth=0.5)
    else:
        plt.plot(wvlngth,smoothed,linewidth=0.5)
    return np.nanmax(smoothed)

# this fn does not seem to be used!
def broaden_template(I,kms):
    npix = kms/300000/1e-4 # b/c template pix are spaced by 1 part in 1e4
    return scipy.ndimage.gaussian_filter1d(I,npix)

def throughput(W,model=None):
    # W: wavelength grid on which to realize the model
    if not model:
        # use default model: load it if necessary
        if 'thruput' not in prefs:
            with open(prefs['THRUPUT_FILE'],'rb') as f:
                 prefs['thruput'] = pickle.load(f)
        return prefs['thruput'](W)/prefs['thruput'](W).mean()
        # DEPRECATED: super-simple model
        #thruput = 0.25+0.0005*(W-5100)
        #thruput[W>6300] = 0.85+0.00025*(W[W>6300]-6300)
        #thruput[thruput>1.0] = 1.0
        #return thruput/thruput.mean()

    try:
        thruput = model(W)
        return thruput/thruput.mean()
    except:
        sys.stderr.write('Cannot evaluate thruput model!\n')
        return 1


def getlinelist():
    # most of these lines are from the SDSS line list.  I added the
    # FeI lines based on Scott Adler's juxtaposing the NIST Atomic
    # Spectra Database with spectra & templates.  NB: I haven't been
    # consistent in taking all the rest wavelengths from a single
    # consistent source. Eg the 8 UV lines are listed in air here, but
    # I'm not sure all the rest are in air. That's ok in the sense
    # these are not used in z calculations, only to guide the eye, and
    # the eye will notice being off by 1.0003
    return {'2326.00' : 'CII]',
            '2344.21' : 'FeII',
            '2374.46' : 'FeII',
            '2382.76' : 'FeII',
            '2586.65' : 'FeII',
            '2600.17' : 'FeII',
            '2796.35' : 'MgII',
            '2803.53' : 'MgII',
            '3727.0921' : 'OII',
            '3771' : 'Hι',
            '3798' : 'Hθ',
            '3835.40' : 'Hη',
            '3869' : 'NeIII',
            '3889.06' : 'Hζ',
            '3970.08' : 'Hε',
            '4102.89' : 'Hδ',
            '4341.68' : 'Hγ',
            '4862.68' : 'Hβ',
            '4932.603' : 'OIII',
            '4960.295' : 'OIII',
            '5008.240' : 'OIII',
            '6549.86' : 'NII',
            '6564.61' : 'Hα',
            '6585.27' : 'NII',
            '6718.29' : 'SII',
            '6732.67' : 'SII',
            '3934.777' : 'K',
            '3969.588' : 'H',
            '4305.61' : 'G',
            '5176.7' : 'Mg',
            '5269.5' : 'FeI',
            '5328.0' : 'FeI',
            '5371.5' : 'FeI',
            '5895.6' : 'Na',
            '8500.36' : 'CaII',
            '8544.44' : 'CaII',
            '8664.52' : 'CaII',
            '5578.5' : 'Sky',
            '5894.6' : 'Sky',
            '6301.7' : 'Sky',
            '6880' : 'Sabs',
            '7246.0' : 'Sky'}


def mark_lines(z=-1,whichplot=None,ylabel=-10,lambdamin=None,lambdamax=None):
    if whichplot:
        parent = whichplot
    else:
        parent = plt
    linelist = getlinelist()
    if z > -1:
        parent.set_title('z=%.6f' % (z))
    for wvlngthstr in linelist.keys():
        wvlngth = float(wvlngthstr)
        mylabel = linelist[wvlngthstr]
        if mylabel == 'Sky' or mylabel == 'Sabs':
            if lambdamin and wvlngth<lambdamin:
                continue
            if lambdamax and wvlngth>lambdamax:
                continue
            parent.axvline(wvlngth,linewidth=0.5,linestyle=':',color='r')
            parent.text(wvlngth,ylabel+3,mylabel,ha='center',fontdict={'fontsize':8})
        elif z>=0:
            wvlngth *= (1+z)
            parent.axvline(wvlngth,linewidth=0.5,linestyle=':')
            parent.text(wvlngth,ylabel,mylabel,ha='center',fontdict={'fontsize':8})
    return

def open_notes(path,filename='pyzenotes.pickle'):
    mydir = os.path.dirname(path)
    notesfile = mydir+'/'+filename
    if not os.path.exists(notesfile):
        return {}
    with open(notesfile,'rb') as f:
        return pickle.load(f)


def save_notes(notesdict,path,filename):
    mydir = os.path.dirname(path)
    notesfile = mydir+'/' + filename
    try:
        f = open(notesfile,'wb')
        pickle.dump(notesdict,f)
        f.close()
    except:
        # file not writable! Try a backup location
        sys.stderr.write('ERROR: cannot write to %s\n' % (notesfile))
        mydir = os.path.expanduser('~')
        notesfile = mydir+'/' + filename
        sys.stderr.write('will try writing to %s\n' % (notesfile))
        with open(notesfile,'wb') as f:
            pickle.dump(notesdict,f)
    return

# my file is reformatted from
#https://drive.google.com/drive/folders/1FFRWjUZ58HiDuDD33MYqBzMWDQanBRRy
def get_telluric_absorption(filename,airmass):
    with fits.open(filename) as hdulist:
        hdr = hdulist[0].header
        lambdamin = hdr['LAMBDA0']
        stepfac = hdr['STEPFAC']
        data = hdulist[0].data
        (nairmass,nwvlngths) = data.shape
        # next lines are based on Joe's hdr, which is crypytic; I should dbl-chk
        airmass_step = (3.25-1)/(nairmass-1) 
        wvlngthlist = lambdamin*stepfac**np.arange(nwvlngths)
        airmassindex = round((airmass-1)/airmass_step)
        if airmassindex >= nairmass:
            airmassindex = nairmass-1
        # smooth to ~1A resolution
        abs_spec = gaussian_filter(data[airmassindex,:],3,mode='reflect')
        return Spec(wvlngthlist,abs_spec,0)
