#!/usr/bin/env python

# Gaussian fit + GUI stuff stolen from Scott Ransom's pygaussfit.py and subsequently HACKED
# To be used with PSRCHIVE Archive files

# This software lays on the bed of Procrustes all too comfortably.

#Next two lines needed for dispatching on nodes
#import matplotlib
#matplotlib.use("Agg")

from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from psr_utils import rotate, fft_rotate
import numpy as np
import psrchive as pr
import numpy.fft as fft
import numpy.ma as ma
import scipy.optimize as opt
import lmfit as lm
import sys,time

#plt.copper()
#plt.gray()
#plt.bone()
plt.pink()
plt.close('all')
cols = ['b','g','r','c','m','y','b','g','r','c','m','y','b','g','r','c','m','y','b','g','r','c','m','y','b','g','r','c','m','y']

#Dconst = 4.148808e3     #Dispersion constant (e**2/(2*pi*m_e*c)) [MHz**2 pc**-1 cm**3 s], used by PRESTO
Dconst = 0.000241**-1   #"Traditional" dispersion constant, used by psrchive

class GaussianSelector:
    def __init__(self, ax, profile, errs, minspanx=None,
                 minspany=None, useblit=True):
        print "============================================="
        print "Left mouse click to draw a Gaussian component"
        print "Middle mouse click to fit components to data"
        print "Right mouse click to remove a component"
        print "============================================="
        print "Press 'q' or close window when done fitting"
        print "============================================="
        self.ax = ax.axes
        self.profile = profile
        self.proflen = len(profile)
        self.phases = np.arange(self.proflen, dtype='d')/self.proflen
        self.errs = errs
        self.visible = True
        self.DCguess = sorted(profile)[len(profile)/10+1]
        self.init_params = [self.DCguess]
        self.ngauss = 0
        self.canvas = ax.figure.canvas
        self.canvas.mpl_connect('motion_notify_event', self.onmove)
        self.canvas.mpl_connect('button_press_event', self.press)
        self.canvas.mpl_connect('button_release_event', self.release)
        self.canvas.mpl_connect('draw_event', self.update_background)
        self.canvas.mpl_connect('key_press_event',self.keypress)
        self.background = None
        self.rectprops = dict(facecolor='white', edgecolor = 'black',
                              alpha=0.5, fill=False)
        self.to_draw = Rectangle((0,0), 0, 1, visible=False, **self.rectprops)
        self.ax.add_patch(self.to_draw)
        self.useblit = useblit
        self.minspanx = minspanx
        self.minspany = minspany
        # will save the data (position at mouseclick)
        self.eventpress = None
        # will save the data (pos. at mouserelease)
        self.eventrelease = None
        self.plot_gaussians(self.init_params)

    def update_background(self, event):
        'force an update of the background'
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)

    def ignore(self, event):
        'return True if event should be ignored'
        # If no button was pressed yet ignore the event if it was out
        # of the axes
        if self.eventpress == None:
            return event.inaxes!= self.ax
        # If a button was pressed, check if the release-button is the
        # same.
        return (event.inaxes!=self.ax or
                event.button != self.eventpress.button)

    def press(self, event):
        'on button press event'
        # Is the correct button pressed within the correct axes?
        if self.ignore(event): return
        # make the drawed box/line visible get the click-coordinates,
        # button, ...
        self.eventpress = event
        if event.button==1:
            self.to_draw.set_visible(self.visible)
            self.eventpress.ydata = self.DCguess

    def release(self, event):
        'on button release event'
        if self.eventpress is None or self.ignore(event): return
        # release coordinates, button, ...
        self.eventrelease = event
        if event.button==1:
            # make the box/line invisible again
            self.to_draw.set_visible(False)
            self.canvas.draw()
            xmin, ymin = self.eventpress.xdata, self.eventpress.ydata
            xmax, ymax = self.eventrelease.xdata, self.eventrelease.ydata
            # calculate dimensions of box 
            if xmin>xmax: xmin, xmax = xmax, xmin
            if ymin>ymax: ymin, ymax = ymax, ymin
            spanx = xmax - xmin
            spany = ymax - ymin
            xproblems = self.minspanx is not None and spanx<self.minspanx
            yproblems = self.minspany is not None and spany<self.minspany
        # call desired function
        self.onselect()
        self.eventpress = None                # reset the variables to their
        self.eventrelease = None              #   inital values

    def update(self):
        'draw using newfangled blit or oldfangled draw depending on useblit'
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            self.ax.draw_artist(self.to_draw)
            self.canvas.blit(self.ax.bbox)
        else:
            self.canvas.draw_idle()

    def onmove(self, event):
        if self.eventpress is None or self.ignore(event): return
        x,y = event.xdata, event.ydata # actual position with button still pressed
        minx, maxx = self.eventpress.xdata, x # click-x and actual mouse-x
        miny, maxy = self.eventpress.ydata, y # click-y and actual mouse-y
        if minx>maxx: minx, maxx = maxx, minx # get them in the right order
        if miny>maxy: miny, maxy = maxy, miny
        self.to_draw.set_x(minx)             # set lower left of box
        self.to_draw.set_y(miny)
        self.to_draw.set_width(maxx-minx)     # set width and height of box
        self.to_draw.set_height(maxy-miny)
        self.update()

    def keypress(self, event):
        'on key press event'
        if self.ignore(event): return
        self.eventpress = event
        if event.key=='q':
            self.eventpress = None
            self.close()

    def plot_gaussians(self, params):
        plt.subplot(211)
        plt.cla()
        # Re-plot the original profile
        plt.plot(self.phases, self.profile, c='black', lw=3, alpha=0.3)
        plt.xlabel('Pulse Phase')
        plt.ylabel('Pulse Amplitude')
        DC = params[0]
        # Plot the individual gaussians
        for ii in xrange(self.ngauss):
            loc, wid, amp = params[1+ii*3:4+ii*3]
            plt.plot(self.phases, DC + amp*gaussian_profile(self.proflen, loc, wid),'%s'%cols[ii])

    def onselect(self):
        event1 = self.eventpress
        event2 = self.eventrelease
        # Left mouse button = add a gaussian
        if event1.button == event2.button == 1:
            x1, y1 = event1.xdata, event1.ydata
            x2, y2 = event2.xdata, event2.ydata
            loc = 0.5*(x1+x2)
            wid = np.fabs(x2-x1)
            #amp = np.fabs(1.05*(y2-self.init_params[0])*(x2-x1))
            amp = np.fabs(1.05*(y2-self.init_params[0]))
            self.init_params += [loc, wid, amp]
            self.ngauss += 1
            self.plot_gaussians(self.init_params)
            plt.draw()
        # Middle mouse button = fit the gaussians
        elif event1.button == event2.button == 2:
            fit_params, chi_sq, dof, residuals = \
                        fit_gaussian_profile(self.profile, self.init_params,
                                      np.zeros(self.proflen)+self.errs,quiet=True)
            self.fit_params = fit_params
            # scaled uncertainties
            #scaled_fit_errs = fit_errs * np.sqrt(chi_sq / dof)

            # Plot the best-fit profile
            self.plot_gaussians(fit_params)
            fitprof = gen_gaussian_profile(fit_params, self.proflen)
            plt.plot(self.phases, fitprof, c='black', lw=1)
            plt.draw()

            # Plot the residuals
            plt.subplot(212)
            plt.cla()
            residuals = self.profile - fitprof
            plt.plot(self.phases, residuals,'k')
            plt.xlabel('Pulse Phase')
            plt.ylabel('Data-Fit Residuals')
            plt.draw()
        # Right mouse button = remove last gaussian
        elif event1.button == event2.button == 3:
            if self.ngauss:
                self.init_params = self.init_params[:-3]
                self.ngauss -= 1
                self.plot_gaussians(self.init_params)
                plt.draw()
                plt.subplot(212)
                plt.cla()
                plt.xlabel('Pulse Phase')
                plt.ylabel('Data-Fit Residuals')
                plt.draw()

    def close(self):
        plt.close(1)
        plt.close(2)

def gaussian_profile(N, loc, wid, norm=False, abs_wid=False, zeroout=True):
    """
    gaussian_profile(N, loc, wid):
        Return a gaussian pulse profile with 'N' bins and
        an integrated 'flux' of 1 unit (if norm=True; default norm=False and peak ampltiude = 1).
            'N' = the number of points in the profile
            'loc' = the pulse phase (0-1)
            'wid' = the gaussian pulses full width at half-max
        Note:  The FWHM of a gaussian is approx 2.35482 sigma (exactly 2*sqrt(2*ln(2)))
    """
    #Maybe should move these checks to gen_gaussian_portrait?
    if abs_wid: wid = abs(wid)
    if wid > 0.0: pass
    elif wid == 0.0: return np.zeros(N, 'd')
    elif wid < 0.0 and zeroout: return np.zeros(N, 'd')
    elif wid < 0.0 and not zeroout: pass
    else: return 0
    sigma = wid / (2*np.sqrt(2*np.log(2)))
    mean = loc % 1.0
    locval = np.arange(N, dtype='d') / float(N)
    if (mean < 0.5):
        locval = np.where(np.greater(locval, mean+0.5),
                           locval-1.0, locval)
    else:
        locval = np.where(np.less(locval, mean-0.5),
                           locval+1.0, locval)
    try:
        zs = (locval-mean)/sigma
        okzinds = np.compress(np.fabs(zs)<20.0, np.arange(N))   #Why 20?
        okzs = np.take(zs, okzinds)
        retval = np.zeros(N, 'd')
        np.put(retval, okzinds, np.exp(-0.5*(okzs)**2.0)/(sigma*np.sqrt(2*np.pi)))
        if norm: return retval
        #else: return retval/np.max(retval)
        else:
            if np.max(abs(retval)) == 0.0: return retval   #TP hack
            else: return retval/np.max(abs(retval))  #TP hack
    except OverflowError:
        print "Problem in gaussian prof:  mean = %f  sigma = %f" % \
              (mean, sigma)
        return np.zeros(N, 'd')

def gen_gaussian_profile(params, N):
    """
    gen_gaussian_profile(params, N):
        Return a model of a DC-component + M gaussians
            params is a sequence of 1+M*3 values
                the first value is the DC component.  Each remaining
                group of three represents the gaussians loc (0-1),
                wid (FWHM) (0-1), and amplitude (>0.0).
            N is the number of points in the model.
    """
    ngauss = (len(params)-1)/3
    model = np.zeros(N, dtype='d') + params[0]
    for ii in xrange(ngauss):
        loc, wid, amp = params[1+ii*3:4+ii*3]
        model += amp * gaussian_profile(N, loc, wid)
    return model

def gen_gaussian_portrait(params, phases, freqs, nu_ref):
    """
    """
    refparams = np.array([params[0]]+list(params[1::2]))
    locparams = params[2::6]
    widparams = params[4::6]
    ampparams = params[6::6]
    ngauss = len(refparams[1::3])
    nbin = len(phases)
    nchan = len(freqs)
    gport = np.empty([nchan,nbin])
    gparams = np.empty([nchan,len(refparams)])
    gparams[:,0] = refparams[0]    #DC term
    gparams[:,1::3] = np.outer(freqs-nu_ref,locparams)+np.outer(np.ones(nchan),refparams[1::3])    #Locs
    gparams[:,2::3] = np.outer(freqs-nu_ref,widparams)+np.outer(np.ones(nchan),refparams[2::3])    #Wids
    #gparams[:,0::3][:,1:] = np.exp(np.outer(np.log(freqs)-np.log(nu_ref),ampparams)+np.outer(np.ones(nchan),np.log(refparams[0::3][1:])))    #Amps; I am unsure why I needed this fix at some point
    gparams[:,3::3] = np.exp(np.outer(np.log(freqs)-np.log(nu_ref),ampparams)+np.outer(np.ones(nchan),np.log(refparams[3::3])))    #Amps
    for nn in range(nchan):
        gport[nn] = gen_gaussian_profile(gparams[nn],nbin) #Maybe need to contrain so values don't go negative?
    return gport

def powlaw(nu,nu0,A,alpha):
    """
    Power-law spectrum given by:
    F(nu) = A*(nu/nu0)**(alpha)
    """
    return A*((nu/nu0)**(alpha))

def powlawint(nu2,nu1,nu0,A,alpha):
    """
    Returns the integral over a powerlaw of form A*(nu/nu0)**(alpha)
    from nu1 to nu2
    """
    alpha = np.float(alpha)
    if alpha == -1.0:
        return A*nu0*np.log(nu2/nu1)
    else:
        C = A*(nu0**-alpha)/(1+alpha)
        diff = ((nu2**(1+alpha))-(nu1**(1+alpha)))
        return C*diff

def powlawnus(lo,hi,N,alpha,mid=False):
    """
    Returns frequencies such that a bandwidth from lo to hi frequencies
    split into N chunks contains the same amount of power in each chunk
    given a power-law across the band with spectral index alpha.  Default
    behavior returns N+1 frequencies (includes both lo and hi freqs); if
    mid=True, will return N frequencies, corresponding to the middle frequency
    in each chunk.
    """
    alpha = np.float(alpha)
    nus = np.zeros(N+1)
    if alpha == -1.0:
        nus = np.exp(np.linspace(np.log(lo),np.log(hi),N+1))
    else:
        nus = np.power(np.linspace(lo**(1+alpha),hi**(1+alpha),N+1),(1+alpha)**-1)
        #Equivalently:
        #for nn in xrange(N+1):
        #    nus[nn] = ((nn/np.float(N))*(hi**(1+alpha)) + (1-(nn/np.float(N)))*(lo**(1+alpha)))**(1/(1+alpha))
    if mid:
        midnus = np.zeros(N)
        for nn in xrange(N):
            midnus[nn] = 0.5*(nus[nn]+nus[nn+1])
        nus = midnus
    return nus

def fit_gauss_function(params, data=None, errs=None):
    """
    """
    prms = np.array([param.value for param in params.itervalues()])
    return (data - gen_gaussian_profile(prms, len(data))) / errs

def fit_gaussian_portrait_function(params, phases, freqs, nu_ref, fjac=None, data=None, errs=None):
    """
    """
    prms = np.array([param.value for param in params.itervalues()])
    deviates = np.ravel((data - gen_gaussian_portrait(prms, phases, freqs, nu_ref)) / errs)
    return deviates

def fit_powlaws_function(params, freqs, nu0, weights=None, fjac=None, data=None, errs=None):
    """
    """
    prms = np.array([param.value for param in params.itervalues()])
    A = prms[0]
    alpha = prms[1]
    d = []
    f = []
    for ii in xrange(len(weights)):
        if weights[ii]:
            d.append(data[ii])
            f.append(freqs[ii])
        else: pass
    d=np.array(d)
    f=np.array(f)
    return (d - powlaw(f,nu0,A,alpha)) / errs

def fit_portrait_function(params, model=None, p=None, data=None, d=None, errs=None, P=None, freqs=None, nu_ref=np.inf):
    """
    """
    phase = params[0]
    m = 0.0
    if P == None or freqs == None:
        Cdm = 0.0
        freqs = np.inf*np.ones(len(model))
    else: Cdm = Dconst*params[1]/P
    for nn in xrange(len(freqs)):
        err = errs[nn]
        freq = freqs[nn]
        harmind = np.arange(len(model[nn]))
        phasor = np.exp(harmind * 2.0j*np.pi*(phase+(Cdm*(freq**-2.0 - nu_ref**-2.0))))
        mm = np.real(data[nn,:] * np.conj(model[nn,:]) * phasor).sum()
        m += (mm**2.0)*err/p[nn]
    return d-m

def fit_portrait_function_deriv(params, model=None, p=None, data=None, d=None, errs=None, P=None, freqs=None, nu_ref=np.inf):
    """
    """
    phase = params[0]
    Cdm = Dconst*params[1]/P
    d_phi,d_DM = 0.0,0.0
    for nn in xrange(len(freqs)):
        err = errs[nn]
        freq = freqs[nn]
        harmind = np.arange(len(model[nn]))
        phasor = np.exp(harmind * 2.0j*np.pi*(phase+(Cdm*(freq**-2.0 - nu_ref**-2.0))))
        g1 = np.real(data[nn,:]*np.conj(model[nn,:]) * phasor).sum()
        gp2 = np.real(2j*np.pi*harmind * data[nn,:] * np.conj(model[nn,:]) * phasor).sum()
        gd2 = np.real(2j*np.pi*harmind * (freq**-2.0 - nu_ref**-2.0)*(Dconst/P) *data[nn,:] * np.conj(model[nn,:]) * phasor).sum()
        d_phi += -2*g1*gp2*err/p[nn]
        d_DM += -2*g1*gd2*err/p[nn]
    return np.array([d_phi,d_DM])

def fit_portrait_function_2deriv(params, model=None, p=None, data=None, d=None, errs=None, P=None, freqs=None, nu_ref=np.inf):      #Covariance matrix...??
    """
    """
    phase = params[0]
    Cdm = Dconst*params[1]/P
    d2_phi,d2_DM = 0.0,0.0
    for nn in xrange(len(freqs)):
        err = errs[nn]
        freq = freqs[nn]
        harmind = np.arange(len(model[nn]))
        phasor = np.exp(harmind * 2.0j*np.pi*(phase+(Cdm*(freq**-2.0 - nu_ref**-2.0))))
        g1 = np.real(data[nn,:]*np.conj(model[nn,:]) * phasor).sum()
        gp2 = np.real(2.0j*np.pi*harmind * data[nn,:] * np.conj(model[nn,:]) * phasor).sum()
        gd2 = np.real(2.0j*np.pi*harmind * (freq**-2.0 - nu_ref**-2.0)*(Dconst/P) *data[nn,:] * np.conj(model[nn,:]) * phasor).sum()
        gp3 = np.real(pow(2.0j*np.pi*harmind,2.0)*data[nn,:] * np.conj(model[nn,:])* phasor).sum()
        gd3 = np.real(pow(2.0j*np.pi*harmind*(freq**-2.0 - nu_ref**-2.0)*(Dconst/P),2) * data[nn,:] * np.conj(model[nn,:]) * phasor).sum()
        d2_phi += -2.0*err*(pow(gp2,2.0)+(g1*gp3))/p[nn]
        d2_DM += -2.0*err*(pow(gd2,2.0)+(g1*gd3))/p[nn]
    return np.array([d2_phi,d2_DM])

def estimate_portrait(phase, DM, data, scales, P, freqs, nu_ref=np.inf): #here, all vars have additional epoch-index except nu_ref, i.e. all have to be arrays of at least len 1; errs are precision
    """
    """
    #Next lines should be done just as in fit_portrait
    dFFT = fft.rfft(data,axis=2)
    unnorm_errs = np.real(dFFT[:,:,-len(dFFT[0,0])/4:]).std(axis=2)**-2.    #Precision FIX
    #norm_dFFT = np.transpose((unnorm_errs**0.5)*np.transpose(dFFT))
    #norm_errs = np.real(norm_dFFT[:,:,-len(norm_dFFT[0,0])/4:]).std(axis=2)**-2.
    errs = unnorm_errs
    D = Dconst*DM/P
    freqs2 = freqs**-2.0 - nu_ref**-2.0
    phiD = np.outer(D,freqs2)
    phiprime = np.outer(phase,np.ones(len(freqs))) + phiD
    weight = np.sum(pow(scales,2.0)*errs,axis=0)**-1
    phasor = np.array([np.exp(2.0j*np.pi*kk*phiprime) for kk in xrange(len(dFFT[0,0]))]).transpose(1,2,0)
    p = np.sum(np.transpose(np.transpose(scales*errs)*np.transpose(phasor*dFFT)),axis=0)
    wp = np.transpose(weight*np.transpose(p))
    return wp

def wiener_filter(prof,noise):      #FIX does not work
    """
    prof is noisy template
    noise is standard deviation of the gaussian noise in the data
    """
    FFT = fft.rfft(prof)
    pows = np.real(FFT*np.conj(FFT)) / len(prof)        #Check Normalization
    return pows/(pows+(noise**2))
    #return (pows - (noise**2)) / pows

def brickwall_filter(n,kc):
    """
    """
    fk = np.zeros(n)
    fk[:kc] = 1.0
    return fk

def find_kc(prof,noise):
    """
    """
    wf = wiener_filter(prof,noise)
    N = len(wf)
    X2 = np.zeros(N)
    for ii in xrange(N):
        X2[ii] = np.sum((wf-brickwall_filter(N,ii))**2)
    return X2.argmin()

def fit_powlaws(data, freqs, nu0, weights, init_params, errs):
    """
    """
    nparam = len(init_params)
    # Generate the parameter structure
    params = lm.Parameters()
    params.add('amp',init_params[0],vary=True,min=None,max=None)
    params.add('alpha',init_params[1],vary=True,min=None,max=None)
    other_args = {'freqs':freqs, 'nu0':nu0, 'weights':weights, 'data':data, 'errs':errs}
    # Now fit it
    results = lm.minimize(fit_powlaws_function, params, kws=other_args)
    fit_params = np.array([param.value for param in results.params.itervalues()])
    dof = results.nfree
    chi_sq = results.chisqr
    redchi_sq = results.redchi
    residuals = results.residual
    fit_errs = np.array([param.stderr for param in results.params.itervalues()])
    return fit_params, fit_errs, chi_sq, dof, residuals

def fit_gaussian_profile(data, init_params, errs, quiet=True):
    """
    """
    nparam = len(init_params)
    ngauss = (len(init_params)-1)/3
    # Generate the parameter structure
    params = lm.Parameters()
    for ii in xrange(nparam):
        if ii == 0:
            params.add('dc', init_params[ii], vary=True, min=None, max=None, expr=None)
        elif ii in range(nparam)[1::3]:
            params.add('loc%s'%str((ii-1)/3+1), init_params[ii], vary=True, min=None, max=None, expr=None)
        elif ii in range(nparam)[2::3]:
            params.add('wid%s'%str((ii-1)/3+1), init_params[ii], vary=True, min=0.0, max=None, expr=None)
        elif ii in range(nparam)[3::3]:
            params.add('amp%s'%str((ii-1)/3+1), init_params[ii], vary=True, min=0.0, max=None, expr=None)
        else:
            print "Undefined index %d."%ii
            sys.exit()
    other_args = {'data':data, 'errs':errs}
    # Now fit it
    results = lm.minimize(fit_gauss_function, params, kws=other_args)
    #fit_params = results.vars
    fit_params = np.array([param.value for param in results.params.itervalues()])
    dof = results.nfree
    redchi_sq = results.redchi
    residuals = results.residual
    if not quiet:
        print "------------------------------------------------------------------"
        print "Multi-Gaussian Profile Fit Results"
        print "------------------------------------------------------------------"
        print "lmfit status:", results.message
        print "gaussians:", ngauss
        print "DOF:", dof
        print "reduced chi-sq: %.2f" % redchi_sq
        print "residuals mean: %.3g" % np.mean(residuals)
        print "residuals stdev: %.3g" % np.std(residuals)
        print "------------------------------------------------------------------"
    return fit_params, redchi_sq, dof, residuals

def fit_gaussian_portrait(data, errs, init_params, fix_params, phases, freqs, nu_ref, quiet=True):
    """
    """
    nparam = len(init_params)
    ngauss = (len(init_params)-1)/6
    fixloc,fixwid,fixamp = fix_params
    # Generate the parameter structure
    params = lm.Parameters()
    for ii in xrange(nparam):
        if ii == 0:     #DC, limited by 0
            params.add('dc', init_params[ii], vary=True, min=None, max=None, expr=None)
        elif ii%6 == 1:     #loc limits
            params.add('loc%s'%str((ii-1)/6+1), init_params[ii], vary=True, min=None, max=None, expr=None)
        elif ii%6 == 2:     #loc slope limits
            params.add('m_loc%s'%str((ii-1)/6+1), init_params[ii], vary=not(fixloc), min=None, max=None, expr=None)
        elif ii%6 == 3:     #wid limits, limited by 0
            params.add('wid%s'%str((ii-1)/6+1), init_params[ii], vary=True, min=0.0, max=None, expr=None)
        elif ii%6 == 4:     #wid slope limits
            params.add('m_wid%s'%str((ii-1)/6+1), init_params[ii], vary=not(fixwid), min=None, max=None, expr=None)
        elif ii%6 == 5:     #amp limits, limited by 0
            params.add('amp%s'%str((ii-1)/6+1), init_params[ii], vary=True, min=0.0, max=None, expr=None)
        elif ii%6 == 0:     #amp index limits
            params.add('alpha%s'%str((ii-1)/6+1), init_params[ii], vary=not(fixamp), min=None, max=None, expr=None)
        else:
            print "Unfortunate index."
            sys.exit()
    other_args = {'data':data, 'errs':errs, 'phases':phases, 'freqs':freqs, 'nu_ref':nu_ref}
    # Now fit it
    results = lm.minimize(fit_gaussian_portrait_function, params, kws=other_args)
    fit_params = np.array([param.value for param in results.params.itervalues()])
    dof = results.nfree
    chi_sq = results.chisqr
    redchi_sq = results.redchi
    residuals = results.residual
    model = gen_gaussian_portrait(fit_params,phases,freqs,nu_ref)
    if not quiet:
        print "------------------------------------------------------------------"
        print "Gaussian Portrait Fit"
        print "------------------------------------------------------------------"
        print "lmfit status:", results.message
        print "gaussians:", ngauss
        print "DOF:", dof
        print "reduced chi-sq: %.2f" % redchi_sq
        print "residuals mean: %.3g" % np.mean(residuals)
        print "residuals stdev: %.3g" % np.std(residuals)
        print "------------------------------------------------------------------"
    return fit_params, chi_sq, dof

def fit_portrait(data,model,init_params,P=None,freqs=None,nu_ref=np.inf,scales=True):
    """
    """
    #errs = get_noise(data,tau=True,chans=True,fd=True,frac=4) #tau = precision = 1/variance.  FIX Need to use better filtering instead of frac     #FIX get_noise is not right
    dFFT = fft.rfft(data,axis=1)
    mFFT = fft.rfft(model,axis=1)
    unnorm_errs = np.real(dFFT[:,-len(dFFT[0])/4:]).std(axis=1)**-2.    #Precision FIX
    norm_dFFT = np.transpose((unnorm_errs**0.5)*np.transpose(dFFT))
    norm_errs = np.real(norm_dFFT[:,-len(norm_dFFT[0])/4:]).std(axis=1)**-2.
    errs = unnorm_errs
    d = np.real(np.sum(np.transpose(errs*np.transpose(dFFT*np.conj(dFFT)))))
    p = np.real(np.sum(mFFT*np.conj(mFFT),axis=1))
    #other_args = {'model':mFFT,'p':p,'data':dFFT,'d':d,'errs':errs,'P':P,'freqs':freqs,'nu_ref':nu_ref}
    other_args = (mFFT,p,dFFT,d,errs,P,freqs,nu_ref)    #Order matters
    minimize = opt.minimize
    method = 'TNC'  #Seems to work best, fastest
    bounds = [(None,None),(None,None)]    #Bounds on phase and DM
    start = time.time()
    results = minimize(fit_portrait_function,init_params,args=other_args,method=method,jac=fit_portrait_function_deriv,bounds=bounds,options={'disp':False})
    duration = time.time()-start
    phi = results.x[0]
    DM = results.x[1]
    nfeval = results.nfev
    return_code = results.status
    #If the fit fails...????
    if results.success is not True: print "Fit failed.  Return code is %d"%results.status
    param_errs = list(pow(fit_portrait_function_2deriv(np.array([phi,DM]),mFFT,p,dFFT,d,errs,P,freqs,nu_ref),-0.5))
    DoF = len(data.ravel()) - (len(freqs)+2)
    red_chi2 = results.fun / DoF
    if scales:
        scales = get_scales(data,model,phi,DM,P,freqs,nu_ref)
        param_errs += list(pow(2*p*errs,-0.5))  #Errors on scales, if ever needed
        return phi, DM, nfeval, return_code, scales, np.array(param_errs), red_chi2, duration
    else: return phi, DM, nfeval, return_code, np.array(param_errs), red_chi2, duration

def first_guess(data,model,nguess=1000):
    """
    """
    #Get initial guesses for phase, and amplitudes...
    #guessparams = []
    #guesschi2s = []
    crosscorr = np.empty(nguess)
    #phaseguess = np.linspace(0,1.0,nguess)
    phaseguess = np.linspace(-0.5,0.5,nguess)
    for ii in range(nguess):
        phase = phaseguess[ii]
        crosscorr[ii] = np.correlate(fft_rotate(np.sum(data,axis=0),phase*len(np.sum(data,axis=0))),np.sum(model,axis=0))
    phaseguess = phaseguess[crosscorr.argmax()]
    return phaseguess

def make_model(modelfile,phases,freqs,quiet=False):
    """
    """
    nbin = len(phases)
    nchan = len(freqs)
    modeldata = open(modelfile,"r").readlines()
    ngauss = len(modeldata)-3
    params = np.zeros(ngauss*6+1)
    source = modeldata.pop(0)[:-1]
    nu_ref = float(modeldata.pop(0))
    params[0] = float(modeldata.pop(0))
    for gg in xrange(ngauss):
        comp = map(float,modeldata[gg].split())
        params[1+gg*6:7+(gg*6)] = comp
    model = gen_gaussian_portrait(params,phases,freqs,nu_ref)
    if not quiet: print "Made %d component model for %s with %d frequency channels, %d profile bins, %.0f MHz bandwidth centered near %.2f MHz"%(ngauss,source,nchan,nbin,(freqs[-1]-freqs[0])+((freqs[-1]-freqs[-2])),freqs.mean())
    return source,ngauss,model

def get_noise(data,frac=4,tau=False,chans=False,fd=False):     #FIX: Make sure to use on portraits w/o zapped freq. channels, i.e. portxs     FIX: MAKE SIMPLER!!!    FIX: Implement k_max from wiener/brick-wall filter fit        #FIX This is not right 
    """
    """
    shape = data.shape
    if len(shape) == 1:
        prof = data
    elif shape[0] == 1:
        prof = data[0]
    elif shape[1] == 1:
        prof = data[:,0]
    else: pass
    try:
        FFT = fft.rfft(prof)
        if fd:
            if tau: return np.std(np.real(FFT)[-len(FFT)/frac:])**-2
            #if tau: return (np.std(np.real(FFT)[-len(FFT)/frac:])**-2,np.std(np.imag(FFT)[-len(FFT)/frac:])**-2)
            else: return np.std(np.real(FFT)[-len(FFT)/frac:])
            #else: return (np.std(np.real(FFT)[-len(FFT)/frac:]),np.std(np.imag(FFT)[-len(FFT)/frac:]))
        else:
            pows = np.real(FFT*np.conj(FFT))/len(prof)    #!!!CHECK NORMALIZATION
            if tau: return (np.mean(pows[-len(pows)/frac:]))**-1
            else: return np.sqrt(np.mean(pows[-len(pows)/frac:]))
    except(NameError):
        noise = np.zeros(len(data))
        if fd:
            for nn in range(len(noise)):
                    prof = data[nn]
                    FFT = fft.rfft(prof)
                    noise[nn] = np.std(np.real(FFT)[-len(FFT)/frac:])
            if chans:
                if tau: return noise**-2
                else: return noise
            else:
                if tau: return np.median(noise)**-2     #not statistically rigorous
                else: return np.median(noise)
        else:
            for nn in range(len(noise)):
                prof = data[nn]
                FFT = fft.rfft(prof)
                pows = np.real(FFT*np.conj(FFT))/len(prof)    #!!!CHECK NORMALIZATION
                noise[nn] = np.sqrt(np.mean(pows[-len(pows)/frac:]))
            if chans:
                if tau: return noise**-2
                else: return noise
            else:
                if tau: return np.median(noise)**-2     #not statistically rigorous
                else: return np.median(noise)

def get_scales(data,model,phase,DM,P,freqs,nu_ref):
    """
    """
    scales = np.zeros(len(model))
    dFFT = fft.rfft(data,axis=1)
    mFFT = fft.rfft(model,axis=1)
    p = np.real(np.sum(mFFT*np.conj(mFFT),axis=1))
    Cdm = DM*Dconst/P
    for kk in range(len(mFFT[0])):  #FIX vectorize
        scales += np.real(dFFT[:,kk]*np.conj(mFFT[:,kk])*np.exp(2j*np.pi*kk*(phase+(Cdm*(pow(freqs,-2)-pow(nu_ref,-2))))))/p
    return scales

def rotate_portrait(port,phase,DM=None,P=None,freqs=None,nu_ref=np.inf):
    """
    Positive values of phase and DM rotate to earlier phase.
    """
    pFFT = fft.rfft(port,axis=1)
    for nn in xrange(len(pFFT)):
        if DM is None and freqs is None: pFFT[nn,:] *= np.exp(np.arange(len(pFFT[nn])) * 2.0j*np.pi*phase)
        else:
            Cdm = DM*Dconst/P
            freq = freqs[nn]
            phasor = np.exp(np.arange(len(pFFT[nn])) * 2.0j*np.pi*(phase+(Cdm*(freq**-2.0 - nu_ref**-2.0))))
            pFFT[nn,:] *= phasor
    return fft.irfft(pFFT)

def plot_PL_results(M,withprof=1,witherrors=1, negative=0):
    """
    M is a pymc MCMC object (a chain).
    """
    nbin = M.nbin
    alphamu = M.alphas.stats()['mean']
    alpha95 = M.alphas.stats()['95% HPD interval']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    if witherrors:
        if negative: ax1.errorbar((np.arange(nbin)+0.5)/nbin,-alphamu,yerr=[alpha95.transpose()[1]-alphamu,alphamu-alpha95.transpose()[0]],fmt='r+')
        else: ax1.errorbar((np.arange(nbin)+0.5)/nbin,alphamu,yerr=[alphamu-alpha95.transpose()[0],alpha95.transpose()[1]-alphamu],fmt='r+')
    else:
        if negative: ax1.plot((np.arange(nbin)+0.5)/nbin, -alphamu,'r+')
        else: ax1.plot((np.arange(nbin)+0.5)/nbin, alphamu,'r+')
    if withprof:
        plt.axes()
        plt.plot(np.arange(nbin, dtype='d')/nbin,prof,'k-')
    plt.ylim(alphamu.min()-5,plt.ylim()[1])
    plt.show()

def write_model(filenm,source,model_params,nu_ref):
    """
    """
    outfile = open(filenm,"a")
    outfile.write("%s\n"%source)
    outfile.write("%.8f\n"%nu_ref)
    outfile.write("%.8f\n"%model_params[0])
    ngauss = (len(model_params)-1)/6
    for nn in xrange(ngauss):
        comp = model_params[1+nn*6:7+nn*6]
        outfile.write("%.8f\t %.8f\t %.8f\t %.8f\t %.8f\t %.8f\n"%(comp[0],comp[1],comp[2],comp[3],comp[4],comp[5]))
    outfile.close()
    print "%s written."%filenm

def load_data(filenm,dedisperse=False,tscrunch=False,pscrunch=False,quiet=False,rm_baseline=(0,0),Gfudge=1.0):
    """
    Will read data using PSRCHIVE.
    """
    #Load archive
    arch = pr.Archive_load(filenm)
    source = arch.get_source()
    if not quiet: print "Reading data from %s on source %s..."%(filenm,source)
    #Get some metadata
    nu0 = arch.get_centre_frequency()   #Center of the band
    #bw = abs(arch.get_bandwidth())      #For the -200 MHz cases.  Good fix?
    bw = arch.get_bandwidth()      #For the -200 MHz cases.  Good fix?
    nchan = arch.get_nchan()
    #By-hand frequency calculation, equivalent to below from psrchive
    #chanwidth = bw/nchan
    #lofreq = nu0-(bw/2)
    #freqs = np.linspace(lofreq+(chanwidth/2.0),lofreq+bw-(chanwidth/2.0),nchan)     #Centers of frequency channels
    freqs = np.array([arch.get_Integration(0).get_centre_frequency(ii) for ii in range(nchan)])
    #freqs.sort()                #Again, for the negative BW cases.  Good fix?
    nbin = arch.get_nbin()
    phases = np.arange(nbin, dtype='d')/nbin
    #Dedisperse?
    if dedisperse: arch.dedisperse()
    else: arch.dededisperse()
    #This is where I think the bandpass is being removed (needs to be robust...)
    baseline_removed = 0
    if rm_baseline[0] + rm_baseline[1] == 0:
        arch.remove_baseline()
        baseline_removed = 1
    #pscrunchd?
    if pscrunch: arch.pscrunch()
    else:
        print "Full Stokes not ready.  Try again later.  Sorry"
        return 0
    #if arch.get_state != 'Intensity': arch.pscrunch()  #Stokes
    #tscrunch?
#    if tscrunch:
#        arch.tscrunch()
#        nsub = arch.get_nsubint()
#        P = arch.get_Integration(0).get_folding_period()
#        #Get data
#        port = arch.get_data()[0]*Gfudge    #Stokes
#        for ss in range(4):                 #Stokes
#            for ff in range(nchan):         #Stokes
#                port[ss,ff] -= port[ss,ff,260:298].mean()   #Stokes     ##UGLY HARDCODE!!
#                #port[ss,ff] -= port[ss,ff,260/6:298/6].mean()   #Stokes        ##UGLY HARDCODE!!
#        #Get weights !!!Careful about this!!!
#        #weights = arch.get_weights()[0]
#        weights = port[0].sum(axis=1)  #STOKES FAKE WEIGHTS
#        normweights = np.divide(map(int,weights),map(int,weights))
#        maskweights = (normweights+1)%2
#        portweights = np.array([np.ones(nbin)*maskweights[ii] for ii in xrange(nchan)])
#        portx = []    #Stokes
#        for ss in range(4):     #Stokes
#            sp = screen_portrait(port[ss],portweights)  #Does not change port in this case (already "masked")   #Stokes
#            portx.append(sp[1])
#        portx = np.array(portx)
#        #Estimate noise
#        noise_stdev = np.zeros(4)     #Stokes
#        for ss in range(4):
#            noise_stdev[ss] = get_noise(portx[ss])   #Stokes #FIX This is probably not right
#        #Make flux profile
#        #fluxprof = port.sum(1)/nbin  #This is about equal to bscrunch to ~6 places     #Stokes
#        #fluxprofx = ma.masked_array(fluxprof,mask=maskweights).compressed()            #Stokes
#        freqsx = ma.masked_array(freqs,mask=maskweights).compressed()
#        #Get pulse profile
#        arch.fscrunch()
#        prof = arch.get_data()[0][0][0]*Gfudge
#        if not quiet:
#            print "\tcenter freq. (MHz) = %.5f\n\
#            bandwidth (MHz)    = %.1f\n\
#            # bins in prof.    = %d\n\
#            # channels         = %d\n\
#            unzapped chan.     = %d"%(nu0,bw,nbin,nchan,len(portx[0]))
#        arch.refresh()
#        return source,arch,port,portx,noise_stdev,prof,nbin,phases,nu0,bw,nchan,freqs,freqsx,nsub,P,weights,normweights,maskweights,portweights    #Stokes
#    else:
    #tscrunch?
    if tscrunch: arch.tscrunch()
    nsub = arch.get_nsubint()
    #Get data
    ports = arch.get_data()[:,0,:,:]*Gfudge     #FIX Here assumes pscrunched in second index
    Ps = np.array([arch.get_Integration(ii).get_folding_period() for ii in xrange(nsub)],dtype=np.double)
    MJDs = [arch.get_Integration(ii).get_epoch() for ii in xrange(nsub)]
    #Get weights !!!Careful about this!!!
    weights = arch.get_weights()
    normweights = np.array([np.divide(map(int,weights[ii]),map(int,weights[ii])) for ii in xrange(len(weights))])
    maskweights = (normweights+1)%2
    portweights = np.array([np.array([np.ones(nbin)*maskweights[ii,jj] for jj in xrange(nchan)]) for ii in xrange(nsub)])
    ports,portxs = np.array([screen_portrait(ports[ii],portweights[ii])[0] for ii in xrange(nsub)]),np.array([screen_portrait(ports[ii],portweights[ii])[1] for ii in xrange(nsub)])    #FIX latter part may not work if portxs have different sizes
    arch.tscrunch()
    port = ports.mean(axis=0)
    portx = portxs.mean(axis=0)
    #Estimate noise
    noise_stdev = np.zeros(nsub)
    for nn in range(nsub):
        noise_stdev[nn] = get_noise(portxs[nn])     #FIX This is probably not right
    fluxprof = port.sum(1)/nbin  #This is about equal to bscrunch to ~6 places
    fluxprofx = ma.masked_array(fluxprof,mask=np.array(map(round,maskweights.mean(axis=0)))).compressed()
    freqsx = ma.masked_array(freqs,mask=np.array(map(round,maskweights.mean(axis=0)))).compressed()     #FIX will not work if portxs have different sizes/different things zapped
    #Get pulse profile
    arch.fscrunch()
    prof = arch.get_data()[0][0][0]*Gfudge
    if not quiet:
        print "\tcenter freq. (MHz) = %.5f\n\
        bandwidth (MHz)    = %.1f\n\
        # bins in prof.    = %d\n\
        # channels         = %d\n\
        # unzapped chan.   ~ %d\n\
        # sub ints         = %d"%(nu0,bw,nbin,nchan,int(np.array(map(len,portxs)).mean()),nsub) #FIX might not work if subints masked differently
    arch.refresh()      #FIX return as is or as requested scrunched?
    if tscrunch: return source,arch,ports[0],portxs[0],noise_stdev[0],fluxprof,fluxprofx,prof,nbin,phases,nu0,bw,nchan,freqs,freqsx,nsub,Ps[0],MJDs[0],weights[0],normweights[0],maskweights[0],portweights[0]
    else: return source,arch,ports,portxs,noise_stdev,fluxprof,fluxprofx,prof,nbin,phases,nu0,bw,nchan,freqs,freqsx,nsub,Ps,MJDs,weights,normweights,maskweights,portweights

def screen_portrait(port,portweights):
    """
    """
    #nbin = portweights.sum(axis=1).max()
    try: nbin = len(port[0])
    except(TypeError): nbin = len(port)
    normweights = (portweights.sum(axis=1)+1)%(nbin+1)
    nchan = nbin-portweights.sum(axis=0).max()
    maskedport = np.transpose(normweights*np.transpose(port))
    portx = ma.masked_array(port,mask=portweights).compressed()
    portx = portx.reshape(len(portx)/nbin,nbin)
    return maskedport,portx

def plot_lognorm(mu,tau,lo=0.0,hi=5.0,npts=500,plot=1,show=0):
    """
    """
    pts = np.empty(npts)
    xs = np.linspace(lo,hi,npts)
    for ii in xrange(npts):
        pts[ii] = np.exp(pm.lognormal_like(xs[ii],mu,tau))
    if plot:
        plt.plot(xs,pts)
    if show:
        plt.show()
    return xs,pts

def plot_gamma(alpha,beta,lo=0.0,hi=5.0,npts=500,plot=1,show=0):
    """
    """
    pts = np.empty(npts)
    xs = np.linspace(lo,hi,npts)
    for ii in xrange(npts):
        pts[ii] = np.exp(pm.gamma_like(xs[ii],alpha,beta))
    if plot:
        plt.plot(xs,pts)
    if show:
        plt.show()
    return xs,pts

def DM_delay(DM,freq,freq2=None,P=None):
    """
    Calculates the delay [s] of emitted frequency freq [MHz] from dispersion measure DM [cm**-3 pc] relative to infinite frequency.  If freq2 is provided, the relative delay is caluclated.  If a period P [s] is provided, the delay is returned in phase.
    """
    if freq2:
        delay = Dconst*DM*((freq**-2)-(freq2**-2))
    else:
        delay = Dconst*DM*(freq**-2)
    if P:
        return delay/P
    else: return delay

def make_fake():
    noise = np.random.normal(size=64**2)
    noise = noise.reshape(64,64)
    fake = np.zeros([64,64])
    fake[:,25] += 5
    fake[:,26] += 3
    fake[:,27] += 1
    fake[:,24] += 3
    fake[:,23] += 1
    model = 12.1*rotate_portrait(fake,-0.17)
    fake += noise
    return fake,model

def write_princeton_toa(toa_MJDi, toa_MJDf, toaerr, freq, DM, obs='@', name=' '*13):
    """
    RIPPED FROM PRESTO

    Princeton Format

    columns     item
    1-1     Observatory (one-character code) '@' is barycenter
    2-2     must be blank
    16-24   Observing frequency (MHz)
    25-44   TOA (decimal point must be in column 30 or column 31)
    45-53   TOA uncertainty (microseconds)
    69-78   DM correction (pc cm^-3)
    """
    # Splice together the fractional and integer MJDs
    toa = "%5d"%int(toa_MJDi) + ("%.13f"%toa_MJDf)[1:]
    if DM!=0.0:
        print obs+" %13s %8.3f %s %8.3f              %9.5f" % \
              (name, freq, toa, toaerr, DM)
    else:
        print obs+" %13s %8.3f %s %8.3f" % \
              (name, freq, toa, toaerr)

def correct_freqs_doppler(freqs,doppler_factor):
    """
    Input topocentric frequencies, output barycentric frequencies.
    doppler_factor = nu_source / nu_observed = sqrt( (1+beta) / (1-beta))
        for beta = v/c ; v positive for increasing source distance
    NB: PSRCHIVE is defining doppler_factor as the inverse of the above.
    """
    return doppler_factor*freqs