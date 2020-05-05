import numpy as np
from numpy import pi,exp,log,log10,sqrt
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from astropy import units as u
from astropy.units.function import LogQuantity as ulog10
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import Planck15
from astropy.convolution import convolve, Gaussian1DKernel, Gaussian2DKernel, Box1DKernel
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.integrate import simps

cosmoJ = FlatLambdaCDM(H0=70, Om0=0.3, name='Jesper')
lam0   = 1215.67 * u.AA

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax, orientation='horizontal')
#------------------------------------------------------------------------------

def readMoCaLaTA(filename,
                 skipIFU    = False,
                 newz       = None,
                 cosmo      = cosmoJ,
                 endianness = 'big'
                 ):
    """
    Read binary file generated by MoCaLaTA.

    Probably there's a more straightforward/efficient/Pythonic way to do this:

    YES! YES THERE IS! You can just:
        >>> from scipy.io import FortranFile
        >>> f     = FortranFile(filename, 'r' )
        >>> ii,ff = np.int32,np.float64
        >>> rec1  = f.read_record(ii) # n_sofar
        >>> rec2  = f.read_record(ii,ii,ii,ii,ff,ff,ff,ff,ff,ff) # pix,SpecRes1D,...
        etc.

    But now that this function works, let's stick with that. It's pretty fast.

    Anyway...

    Each Fortran-write is a record, which is surrounded by two integers giving
    the size of the record in bytes. These are read using the variables EoR and
    BoR (I think they must be two separate variables, since the fields must have
    different names).

    When a record is read (say in rec1), it is an np.ndarray of length 1.
    Its only element, rec1[0], is an np.void of length n+2, where n is the total
    number of elements in the record (i.e. a 2x3 array counts as 6 elements),
    and the +2 account for the record length integers that wrap the record
    Variables can be referred to by their field name rec1['FieldName'], but this
    is an array, even for scalar values, so to extract a scalar, use
        var1 = rec1['FieldName'][0]

    Keywords
    --------
    filename:   Yes. It's the name of the file to be read.
    skipIFU:    A call to readMoCaLaTA() takes ~33 ms if all data is read. If
                that's too much, reading only the parameters and the 1D
                spectrum takes ~0.1 ms.

    Returns
    -------
    par:    Dictionary with variables from the simulation, with units added.
    spec1D: Spectrum of all photons escaping in this direction, in the observed
            frame, in units of erg/s/cm2/arcsec2/AA.
    IFU:    Datacube
    """

    binfile = open(filename, 'rb')

    # Formats
    if endianness=='big':
        i4,f4,f8 = '>i4','>f4','>f8'
    elif endianness=='little':
        i4,f4,f8 = '<i4','<f4','<f8'
    else:
        assert False, "Keyword `endianness` should be either 'little' (<) or 'big' (>)."

    # Record wrappers
    BoR = ('BoR', i4) # Beginning of record
    EoR = ('EoR', i4) # End of record

    # Record 1-3
    dt1  = np.dtype([BoR, ('n_sofar',  f8), EoR])
    dt2  = np.dtype([BoR, ('pix',      i4),
                          ('SpecRes1D',i4),
                          ('SpecRes2D',i4),
                          ('n_ph',     i4),
                          ('L_tot',     f8),
                          ('BW',       f8, (2,)),
                          ('z',        f8),
                          ('d_L',       f8),
                          ('R_box',    f8), EoR])
    dt3  = np.dtype([BoR, ('obspar',   f8, (4,)), EoR])

    rec1 = np.fromfile(file=binfile, dtype=dt1, count=1)
    rec2 = np.fromfile(file=binfile, dtype=dt2, count=1)
    rec3 = np.fromfile(file=binfile, dtype=dt3, count=1)

    par = {}
    par['n_eff']     = rec1['n_sofar'][0]
    par['pix']       = rec2['pix'][0]
    par['SpecRes1D'] = rec2['SpecRes1D'][0]
    par['SpecRes2D'] = rec2['SpecRes2D'][0]
    par['n_ph']      = rec2['n_ph'][0]
    par['Ltot']      = rec2['L_tot'][0]      * u.erg / u.s
    par['BW']        = rec2['BW'][0]         * u.angstrom
    par['z']         = rec2['z'][0]
    par['dL']        = rec2['d_L'][0]        * u.kpc
    par['Rbox']      = rec2['R_box'][0]      * u.kpc    # center-to-face distance (actually R_obs)

    par['Lbox']      = par['Rbox'] * 2
    par['obspar']    = rec3['obspar'][0]     * np.array([u.cm**-2,1,1,1])
    # Add additional elements:
    if newz is not None: par['z'] = newz
    par['dL']        = cosmo.luminosity_distance(par['z']).to(u.cm)
    par['dlam1']     = (par['BW'][1] - par['BW'][0]) / par['SpecRes1D']
    par['dlam2']     = (par['BW'][1] - par['BW'][0]) / par['SpecRes2D']
    par['as_kpc']    = cosmo.arcsec_per_kpc_proper(par['z'])      # arcsec per kpc
    par['Dtheta']    = par['Lbox'] * par['as_kpc']                # Angle spanned by box
    par['dx']        = par['Lbox'] / par['pix']                   # Length of one pixel
    par['dtheta']    = par['dx'] * par['as_kpc']                  # Angle spanned by pixel
    par['dOmega']    = par['dtheta']**2                           # Solid angle of a pixel
    par['Rlim']      = u.Quantity([-par['Rbox'], par['Rbox']])# Box limits in kpc
    par['anglim']    = par['Rlim'] * par['as_kpc']                # Box limits in arcsec
    par['Rvals']     = np.linspace(par['Rlim'][0] + par['dx']/2,  # Position of each pixel's center in kpc along one dimension
                                   par['Rlim'][1] - par['dx']/2,
                                   par['pix'])
    par['angvals']   = par['Rvals']  * par['as_kpc']              # Position of each pixel in arcsec


    # Record 4
    try:
        dt4    = np.dtype([BoR, ('spec1D', f8, (par['SpecRes1D'],)), EoR])
    except ValueError:
        otherEnd = ('big','little')[endianness=='big']
        raise ValueError("Perhaps the endianness is wrong. Try setting endianness='"+otherEnd+"'.")
    rec4   = np.fromfile(file=binfile, dtype=dt4, count=1)
    spec1D = rec4['spec1D'][0]
    spec1D = spec1D / par['n_eff'] * par['Ltot'] / par['dL']**2 \
             / par['dlam1']
  # plt.clf()
  # plt.plot(spec1D)
  # assert False

    if not skipIFU:
        # Record 5
        CCDshape = (par['pix'],par['pix'],par['SpecRes2D'])
      # CCDshape = (par['SpecRes2D'],par['pix'],par['pix'])
        voxtot   = np.product(CCDshape)               #Total number of elements in IFU
        dt5      = np.dtype([BoR, ('IFU', f8, (voxtot,)), EoR])
        rec5     = np.fromfile(file=binfile, dtype=dt5, count=1)
        IFU      = rec5['IFU'][0].reshape(CCDshape, order='F')
        IFU      = IFU / par['n_eff']                #Now sum(IFU) = fesc
        print('IFU  =', IFU.sum())
        IFU      = IFU * par['Ltot']                 #Scale to luminosity
        IFU      = IFU / (4*pi*par['dL']**2)         #Convert to flux
        IFU      = IFU / par['dlam2'] / par['dOmega']#-"- flux density per area
      # print('IFU  =', IFU.sum())
        print('Transposing IFU from (pix,pix,specres) to (specres,pix,pix).')
        return par,spec1D,IFU.T
    else:
        return par,spec1D

    binfile.close()
#------------------------------------------------------------------------------

class MoCaLaTAclass:
    def __init__(self,MoCabin):
        self.MoCabin = MoCabin

        #Read output from MoCaLaTA
        self.par,self.spec1D,self.IFU = readMoCaLaTA(MoCabin) # or self.MoCabin?

        #set default values for par:
        self.keywords = {}
        self.keywords['convplot']    = False
        self.keywords['IGMfile']     = None
        self.keywords['NBaperture']  = None
        self.keywords['SFRmult']     = 1.
        self.keywords['SpecAp']      = None
        self.keywords['brightest']   = False
        self.keywords['lab']         = ''
#       self.keywords['ncolors']     =
        self.keywords['r0']          = 0.
        self.keywords['rmax']        = None
        self.keywords['seeing']      = .1
        self.keywords['sky']         = 0.
        self.keywords['spmax']       = None
        self.keywords['writeSB']     = False
        self.keywords['writespec']   = False
        self.keywords['xrspec']      = None
        self.keywords['ylog']        = False
        self.keywords['yrSB']        = False

    def RunAnalysis(self,CommandLineKeywords):
        for key in CommandLineKeywords:#.keys()?????
            self.keywords[key] = CommandLineKeywords[key]

        if IGMfile != None:
            self.ApplyIGM()

    def ApplyIGM(self):
        self.IFU_IGMApplied_Mean  = 0.5 * self.IFU
        self.IFU_IGMApplied_Upper = 0.5 * self.IFU
        self.IFU_IGMApplied_Lower = 0.5 * self.IFU
#------------------------------------------------------------------------------

def SB_profile(SBmap,par,ap):
    Rbox  = par['Rbox'].value
    Rvals = par['Rvals'].value
    if ap is None: ap = np.array([0,0,Rbox])

    rmax = Rbox
    n    = int(0.5 * par['pix']) #/ 2
  # d    = [[sqrt((ap[0]-par['Rvals'][i])**2 + (ap[1]-par['Rvals'][j])**2) for i in range(par['pix'])] for j in range(par['pix'])]
    r    = np.linspace(0,rmax,n)
    s    = np.zeros_like(r)
    p    = np.zeros_like(r)
    for i in range(par['pix']):
        for j in range(par['pix']):
            d = sqrt((ap[0]-Rvals[i])**2 + (ap[1]-Rvals[j])**2)
            k = int(d/rmax * n)
            if k<n:
                s[k] += SBmap[i,j]
                p[k] += 1
    s = s / (p*par['dOmega'].value)

    return r,s
#------------------------------------------------------------------------------

def showobs(galdir, view,
    spatial_aperture    = None,
    spec_ap             = None,
    mother              = './',
    cosmo               = Planck15,
    vmin                = None,
    vmax                = None,
    endianness          = 'big',
    seeing              = 0,
    get_value           = None,
    showit              = True,
    comp_spec           = 'spectrum.dat',
    obs2rest            = True,
    fix_amplitude       = 1,
    smoothspec          = 1*u.AA,#0,
    ymax_spec           = 1e-16,#None,
    saveit              = None,     #file extension for saved figure, e.g. 'png'
    window_position     = '+580+580' # x-y position of plot's upper left corner; set to None to let matplotlib decide
    ):
    """
    Show spectrum and SB map in a given direction of MoCaLaTA output.

    Example
    -------

    >>> showobs('gal03/0068/','xm')
    """

    # Load binfile
    binfile        = mother + '/' + galdir + '/' + view + '.bin'
    par,spec1D,IFU = readMoCaLaTA(binfile,endianness=endianness,cosmo=cosmo)
    Rvals,angvals  = par['Rvals'],par['angvals']

    #Fix amplitude
    spec1D = fix_amplitude * spec1D
    IFU    = fix_amplitude * IFU

    # Cut out desired wavelength range
    wavelength = (1+par['z']) * np.linspace(par['BW'][0],par['BW'][1],par['SpecRes2D'])
    if spec_ap is None:
        spec_ap = [(1+par['z'])*par['BW'][0], (1+par['z'])*par['BW'][1]]
    isp  = (spec_ap[0] <= wavelength) & (wavelength <= spec_ap[1])
    IFU  = IFU[isp,:,:]
    wavelength = wavelength[isp]

    # Collapse along spectral direction
    SBmap = np.sum(IFU, axis=0) * par['dlam2']
    SBmap = SBmap.value

    # Blur image by seeing
    if seeing != 0:
        assert seeing.unit != u.dimensionless_unscaled, "Seeing must have units of angle"
        seeing_kpc = (seeing / par['as_kpc']).to(u.kpc)
        seeing_pix = (seeing_kpc / par['dx']).value
        stddev_pix = seeing_pix / 2.355
        kernel     = Gaussian2DKernel(stddev_pix)
        SBmap      = convolve(SBmap,kernel)

    # Extract spectrum from aperture
    if spatial_aperture is not None:
        mask2d = np.array([[LA.norm([x-spatial_aperture[0],y-spatial_aperture[1]]) < spatial_aperture[2]
            for x in Rvals.value] for y in Rvals.value])            #True if inside aperture
        mask3d = np.broadcast_to(mask2d, IFU.shape)
        spec1D = np.sum(IFU*mask3d,axis=(1,2)) * par['dOmega']#Coll. along spatial directions
    else:
        spec1D = np.sum(IFU,       axis=(1,2)) * par['dOmega']#Coll. along spatial directions

    # Redshift-dilute spectral density
    spec1D = spec1D / (1+par['z'])
    # Smooth spectrum
    if smoothspec != 0:
        assert (u.Quantity(smoothspec)).unit.is_equivalent(u.m), '`smoothspec` must have dimensions of length.'
        smoothres = (smoothspec / (par['dlam2'] * (1+par['z']))).decompose()
      # print('par  =', par['dlam2'])
      # print('len  =', len(spec1D))
        kernel    = Gaussian1DKernel(smoothres)
      # print('smoothres    =', smoothres)
        spec1D    = convolve(spec1D.value,kernel) * spec1D.unit

    SBmap[SBmap==0] = SBmap[np.nonzero(SBmap)].min() # adding the min ???
    logSBmap   = log10(SBmap)

    # Plot SB map, SB profile, and spectrum
    if showit:
        dx    = 1 / 13.
        wtot  = 10.
        htot  = 5 * dx * wtot
        dy    = dx / htot * wtot
        mleft = .8 * dx
        wfig  = 3 * dx
        hfig  = 3 * dy
        mbot  = .8 * dy
        mtop  = 1 * dy
        pad   = .39
        hbar  = .3 * dy

        SBlo       = SBmap.min()
        SBhi       = SBmap.max()
        logSBlo    = log10(SBlo)
        logSBhi    = log10(SBhi)
        if vmin is None: vmin = logSBlo
        if vmax is None: vmax = logSBhi
      # print(logSBlo, logSBhi)

        # Trying to color pixels according to red/blue peak ratio
        red_blue = False
        if red_blue:
            assert False, "This doesn't really seem to work"
          # ncol = 256
          # if spec_ap == None: spec_ap = [1213*u.AA,1219*u.AA]
          # cmap = np.zeros((ncol,4))
          # for i in range(ncol):
            ired  = (spec_ap[0] <= wavelength) & (wavelength <= lam0)
            iblue = (lam0       <= wavelength) & (wavelength <= spec_ap[1])
            redmap  = np.sum(IFU[ired, :,:], axis=0).value
            bluemap = np.sum(IFU[iblue,:,:], axis=0).value
          # Fred  = simps(spec1D[ired], x=wavelength[ired])
          # Fblue = simps(spec1D[iblue],x=wavelength[iblue])
          # img   = Fblue / Fred
            img   = log10(redmap / bluemap)
            vmin  = img.min()
            vmax  = img.max()
            cmap  = 'bwr'
        else:
            img  = logSBmap
            cmap = 'hot'

        plt.close('all')
        fig = plt.figure(figsize=(wtot,htot))
        ax1 = fig.add_axes([mleft,mbot,wfig,hfig+.13])
      # alphamap = (logSBmap-logSBlo) / (logSBhi - logSBlo)
        im = ax1.imshow(img, cmap=cmap, origin='lower',
                vmin=vmin, vmax=vmax,
              # alpha=alphamap,
                extent=[par['Rlim'][0].value,
                        par['Rlim'][1].value,
                        par['Rlim'][0].value,
                        par['Rlim'][1].value],
                aspect='auto')

        backend = mpl.get_backend() 
        if backend == 'TkAgg':
            wm = plt.get_current_fig_manager()
            wm.window.wm_geometry(window_position)
        else:
            if window_position is not None:
                print("WARNING: window positioning not possible with backend '"+backend+"'. Try 'TkAgg', or set 'window_position' to None.")
        ax1.set_xlabel('x / kpc')
        ax1.set_ylabel('y / arcsec')
        Lbox = par['Lbox']
        Rtix = u.Quantity([-Lbox/2, -Lbox/4, 0*u.kpc, Lbox/4, Lbox/2])
        ax1.set_yticks(Rtix.value)
        angtix = ('{:4.1f} '*len(Rtix)).format(*(Rtix*par['as_kpc']).value).split()
        ax1.set_yticklabels(angtix)

        # Color bar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax1)
        cax = divider.new_vertical(size="5%", pad=pad, pack_start=False)
        fig.add_axes(cax)
        cb = fig.colorbar(im, cax=cax, orientation="horizontal")
        cb.ax.xaxis.set_ticks_position('bottom')
        cb.ax.xaxis.set_label_position('top')
        cb.ax.tick_params(labelsize=8)
        cb.set_label(label='log(SB / $\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-2}\,\mathrm{arcsec}^{-1}$)',
                size=8)

        # SB profile
        r,SBprof = SB_profile(SBmap,par,spatial_aperture)
        if spatial_aperture is not None:
            apCircle = plt.Circle((spatial_aperture[0],spatial_aperture[1]),spatial_aperture[2], color='lime', fill=False)
            ax1.add_artist(apCircle)
        ax2 = fig.add_axes([mleft+wfig+1.3*dx,mbot,wfig,hfig])
        ax2.set_xlim([0,r.max()])
        ax2.set_ylim([1e-4*SBprof.max(),2*SBprof.max()])
        ax2.plot(r,SBprof)
        ax2.set_yscale('log')
        ax2.set_ylabel('log(SB / $\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-2}\,\mathrm{arcsec}^{-1}$)')

        # Spectrum
        ax3 = fig.add_axes([mleft+2*(wfig+1.3*dx),mbot,wfig,hfig])
        if ymax_spec is None: ymax_spec = 2*spec1D.value.max()
        ax3.set_ylim([0,ymax_spec])
        if comp_spec is not None:
            wavelength_comp,flux_comp = np.loadtxt(comp_spec,unpack=True)
            comp_spec = np.interp(wavelength,wavelength_comp,flux_comp)
            ax3.plot(wavelength,comp_spec,'-k',alpha=.5)
        ax3.plot(wavelength,spec1D,'-b')
      # ax3.scatter(wavelength,spec1D,color='r',s=5)
        ax3.plot(lam0*(1+par['z'])*np.array([1,1]), [0,spec1D.max().value],'k--',alpha=.25)
        ax3.set_ylabel('Flux / $\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-2}\,\mathrm{\AA}^{-1}$' )

        if saveit is not None:
            figname = mother + '/' + galdir + '/' + view + '.' + saveit
            plt.savefig(figname,box_inches='tight',dpi=200)

    # Output
  # print(spec1D[100])
  # print(wavelength[100])
    Ftot   = simps(spec1D.to(u.erg/u.s/u.AA/u.cm**2), x=wavelength.to(u.AA))
  # ftest1 = (spec1D*(1+par['z']) * par['dlam2']).sum().value
  # ftest2 = (SBmap * par['dOmega']).sum().value
    Ltot   = Ftot * 4*pi*(par['dL'].to(u.cm).value)**2
  # SFR    = Ltot / 1.1012e42
  # np.testing.assert_approx_equal(Ftot,ftest1,significant=2) #only if full aperture
  # np.testing.assert_approx_equal(Ftot,ftest2,significant=2) #only if full aperture
  # print('dlam2  =', par['dlam2'])
  # print('dL     =', par['dL'])
  # print('Ftot   =', Ftot)
  # print('ftest1 =', ftest1)
  # print('ftest2 =', ftest2)
  # print('Ltot   =', Ltot)
  # print('SFR    =', SFR)

    if get_value == 'Ftot':
        return Ftot
    elif get_value == 'Ltot':
        return Ltot
#------------------------------------------------------------------------------

def checkBuildAMR(parfile,cellfile,**kwargs):
    """
    Purpose
    -------
    Check that BuildAMRfromParticles.f90 builds the cells around the particles
    created by mkClouds.f90 in the right places.

    Only cloud cells are plotted. If you want to include the field cells, in
    BuildAMRfromParticles.f90's subroutine CountCells() remove the part
    " .and. CurrentCell%phase.eq.1" from the if statement (and recompile).

    Keywords
    --------
    Robs:       Plot only within Robs kpc of center

    Usage
    -----
        >>> checkBuildAMR('clouds_par.dat','clouds_cell.dat',Robs=.1)
    """

    def _iltR(x,y,z,R):
        return np.where((x<=R) & (y<=R) & (z<=R))

    fig   = plt.figure(figsize=(5,5))
    ax    = fig.add_subplot(111, projection='3d')

    x,y,z = np.loadtxt(parfile,unpack=True)
    print('n_par  =', len(x))
    if len(x) > 100000: print('WARNING: This is going to be slow!')

  # ax.set_alpha(.01)

    if 'Robs' in kwargs:
        Robs = kwargs['Robs']
        ax.set_xlim([-Robs,Robs])
        ax.set_ylim([-Robs,Robs])
        ax.set_zlim([-Robs,Robs])
        ind = _iltR(x,y,z,Robs)
        ax.scatter(x[ind],y[ind],z[ind],s=1,label='Particles')
    else:
        ax.scatter(x,y,z,s=1,label='Particles')

    x,y,z = np.loadtxt(cellfile,unpack=True)
    print('n_cell =', len(x))

    if 'Robs' in kwargs:
        ind = _iltR(x,y,z,Robs)
        ax.scatter(x[ind],y[ind],z[ind],s=1,label='Cells (excluding ICM cells)')
    else:
        ax.scatter(x,y,z,s=1,label='Cells (excluding ICM cells)')

    ax.legend()
#------------------------------------------------------------------------------

def showAnisotropy(galid,
                   label            = 'sample_9kpc',   # or 'sample_Rbox'
                   nside            = 16,       # Number of pixels per side in map
                   nearest_neighbor = False,    # If True, interpolate from four nearest pixels
                 # fesctype         = 'Lya',    # or 'FUV' or 'unity'
                   flux             = False,
                   showmap          = False,
                   showhist         = False,
                   cosmo            = cosmoJ,
                   endianness       = 'big'
                   ):
    """
    Example from https://stackoverflow.com/questions/31573572/healpy-from-data-to-healpix-map

    Example
    -------
    >>> import analyzeMoCaLaTA as am
    >>> am.showAnisotropy('030068',nearest_neighbor=True, showmap=True)
    """
    import healpy as hp
    from mymath import car2sph

    assert nearest_neighbor, "FIX THIS: nearest_neighbor=False doesn't work!"

    mother      = '/Users/pela/Projects/z8p8/CHAclz8p77rz12_512x64_Alex/'
    galdir      = mother + 'gal' + galid[:2] + '/' + galid[2:] + '/'
    samplefile  = galdir + label + '.out'

    def get_nph_emitted():
        """Get number of emitted photons"""
        d       = mother + 'gal' + galid[:2] + '/'
        logfile = d + galid[2:] + '_' + label + '.log'
        with open(logfile) as f:
            lines = f.readlines()
        i = 0
        while True:
            i -= 1
            lastline = lines[i]
            try:
                if lastline.split()[0] == 'Launched':
                    nph = int(lastline.split()[3]) # 4th item of ['Launched', 'photon', '#', '157800', 'of', '1000000']
                    break
            except:
                pass
            assert i > -18, "Logfile doesn't contain the word 'Launched'"
        return nph # + 5 # Since photons are written for every ten. Better: Subtract penultim. number and div. by 2.
    #-------------------------------

    def get_isoflux():
        """Get measured flux if escape were isotropic"""
        binfile  = galdir + 'xm' + label + '.bin'     #Arbitrarily using 'xm'
        par,_ = readMoCaLaTA(binfile,skipIFU=True,endianness=endianness,cosmo=cosmo) # '_' is just a dummy
        dL       = par['dL']
        Ltot     = par['Ltot']
        return Ltot / (4*pi*dL**2)
    #-------------------------------

    nhat = np.loadtxt(samplefile,unpack=True,usecols=(3,4,5))
    nph  = np.shape(nhat)[1]                    #Number of escaping photons

    npix        = hp.nside2npix(nside)          #Total number of pixels in map
    r,theta,phi = car2sph(nhat)                 #Convert to spherical coords

    hpmap = np.zeros(npix, dtype=np.float)      #Initiate map
    for i in range(nph):
        t,p   = theta[i], phi[i]
        near4 = hp.get_interp_weights(nside,t,p)#Get four nearest pixel and their weights
        if nearest_neighbor:
            pix        = near4[0][near4[1].argmax()]#Pick the one with largest weight
            hpmap[pix] = hpmap[pix] + 1             #Add 1 to corresponding pixel
        else:
            hpmap[near4[0]] = near4[1]

    nph_em = get_nph_emitted()
  # print('nph_em   =', nph_em )
    hpmap  = hpmap / float(nph_em) * npix # * fesc
    fesc   = nph / nph_em

    if flux:
        few

    if showmap:
        hp.mollview(hpmap)                          #Show map

    if showhist:
        nmin = 0#hpmap.min()
        nmax = 3#hpmap.max()
        Dn   = nmax - nmin
        dn   = .1
        bins = int(Dn/dn)
        plt.clf()
        P,edges = np.histogram(hpmap, bins=bins, range=(nmin,nmax), density=True)
        x  = edges[0:-1] + dn/2.

        plt.close()
        plt.plot(x,P,'-b',ls='steps-post',lw=2,alpha=.5,label='Full sky distribution')
        plt.xlim([0,3.])
        plt.ylim([0,1.5])
        plt.xlabel('Relative escape fraction')
        plt.ylabel('Normalized probability')

        Ffullmed,Ffulllo,Ffullhi = np.percentile(hpmap, [50,15.87,84.17])
        plt.errorbar([Ffullmed],[max(P)/2], xerr=[[Ffullmed-Ffulllo],[Ffullhi-Ffullmed]],
                fmt='ro',mec='r',capsize=10,capthick=4,lw=4,markersize=15,alpha=.5)

        F6dirmed,F6dirlo,F6dirhi = get_6dirs(galid,label,endianness=endianness)
        plt.errorbar([F6dirmed],[max(P)/2], xerr=[[F6dirmed-F6dirlo],[F6dirhi-F6dirmed]],
                fmt='co',mec='c',capsize=10,capthick=4,lw=4,markersize=15,alpha=.5)

        plt.errorbar([666],[666], xerr=[[666],[666]], fmt='ro',mec='r',capsize=5,capthick=2,lw=2,markersize=4,alpha=.5,
                label='Full sky [50,16,84] percentiles')

        plt.errorbar([666],[666], xerr=[[666],[666]], fmt='co',mec='c',capsize=5,capthick=2,lw=2,markersize=4,alpha=.5,
                label='Six dirs [50,16,84] percentiles')

        plt.legend()
        plt.show()

    print(np.mean(hpmap))
    return np.percentile(hpmap, [50,15.87,84.17])
#------------------------------------------------------------------------------

def anemall(precalc=True,nside=20):
    import sys
    sys.path.append('/Users/pela/Projects/z8p8/CHAclz8p77rz12_512x64_Alex')
    import quickplot as qp

    gals = ['030061', '030068', '030104', '050026', '070048', '090066',
            '090075', '090095', '120055', '130076', '130082', '150088',
            '150091', '330092', '330113', '480062', '520066', '520067',
            '530090', '530113', '620075', '620100', '650107', '650110',
            '650132', '650153', '710058', '760138', '890070', '890079']

    label = 'sample_9kpc'
    n     = len(gals)

    if precalc:
        Ffullmed = np.array([ 0.88343558, 0.85171103, 0.94398093, 0.80267559, 0.76750879, 0.78527607, 0.86375321, 0.91428571, 0.88006112, 0.89719626, 0.88157061, 0.86413502, 0.912     , 0.5780956 , 0.91534884, 0.79283315, 0.87697929, 0.96      , 0.9218543 , 0.90566038, 0.336     , 0.82572924, 0.81428571, 0.81355932, 0.9316277 , 0.92347157, 0.6937371 , 0.83054893, 0.74366197, 0.92676431])
        Ffulllo  = np.array([ 0.47852761, 0.51711027, 0.45768772, 0.32107023, 0.53725616, 0.47116564, 0.37017995, 0.45714286, 0.47669977, 0.48598131, 0.48638379, 0.54008439, 0.432     , 0.41550621, 0.52465116, 0.537514  , 0.52618758, 0.58285714, 0.57218543, 0.36226415, 0.2592    , 0.43081526, 0.51428571, 0.48813559, 0.61402735, 0.51303976, 0.41293875, 0.45823389, 0.33802817, 0.51131824])
        Ffullhi  = np.array([ 1.39877301, 1.33840304, 1.51609058, 1.60535117, 1.04381196, 1.0601227 , 1.72750643, 1.6       , 1.50343774, 1.42056075, 1.33755541, 1.17468354, 1.44      , 0.74971773, 1.52930233, 1.0212766 , 1.52009744, 1.33714286, 1.36688742, 1.63018868, 0.408     , 1.57965595, 1.62857143, 1.10939908, 1.35509484, 1.33390338, 0.95801789, 1.46062053, 1.75774648, 1.34221039])
        F6dirmed = np.array([ 0.81819018, 0.57089531, 1.06261205, 0.80124869, 0.82996426, 0.50042148, 0.66221193, 1.12872692, 0.81600982, 0.83267545, 0.60234834, 0.48175574, 0.6506801 , 0.30997793, 0.79730022, 0.68953002, 1.06653949, 0.64774038, 0.72683011, 1.30270788, 0.27403854, 0.91497625, 0.81724196, 0.77081942, 0.83649245, 0.67873479, 0.35241104, 0.4957306 , 0.68068191, 0.94662543])
        F6dirlo  = np.array([ 0.63665815, 0.36662507, 0.93133318, 0.38119017, 0.54096622, 0.28310982, 0.54143849, 0.92691883, 0.41202662, 0.72735711, 0.38460092, 0.28142063, 0.41070265, 0.17081527, 0.52083775, 0.49598953, 0.62399436, 0.51186395, 0.49477837, 0.44460862, 0.18024051, 0.5010263 , 0.64515015, 0.25917197, 0.64098035, 0.43973012, 0.31944992, 0.3731934 , 0.41421792, 0.60183702])
        F6dirhi  = np.array([ 1.0084844 , 0.99781449, 1.69370721, 1.47584211, 0.92911222, 0.62697899, 1.27170998, 1.24515481, 1.09977817, 1.19060657, 1.08385582, 0.70637284, 1.20298018, 0.41768997, 1.43059404, 0.78012809, 1.35182201, 0.96334103, 0.97111814, 1.48375399, 0.3117092 , 1.44014353, 1.46503143, 0.86842268, 1.06361962, 1.10040028, 0.66428342, 1.15937711, 1.70080653, 1.16001772])
    else:
        from stuff import counter

        Ffullmed = np.empty(n)
        Ffulllo  = np.empty(n)
        Ffullhi  = np.empty(n)
        F6dirmed = np.empty(n)
        F6dirlo  = np.empty(n)
        F6dirhi  = np.empty(n)

        for i,gal in enumerate(gals):
            counter(i,n)
            Ffullmed[i],Ffulllo[i],Ffullhi[i] = \
                showAnisotropy(gal,label=label,nside=nside,nearest_neighbor=True)
            F6dirmed[i],F6dirlo[i],F6dirhi[i] = \
                get_6dirs(gal,label)

    galids,logMvir = qp.readGalData(1.4,(20,10))

    idx     = [i for i,g in enumerate(galids) if g in gals] # Indices of gals in galids
    logMvir = logMvir[idx]                      #Remove the other elements
    Mvir    = 10**logMvir

    Mvir,Ffullmed,Ffulllo,Ffullhi,F6dirmed,F6dirlo,F6dirhi = \
        [np.array(i) for i in zip(*sorted(zip(Mvir,Ffullmed,Ffulllo,Ffullhi,F6dirmed,F6dirlo,F6dirhi)))]

    xax = range(n)

    Ffullmed[-6:-4] = .6 * Ffullmed[-6:-4]
    Ffulllo [-6:-4] = .6 * Ffulllo [-6:-4]
    Ffullhi [-6:-4] = .6 * Ffullhi [-6:-4]

    plt.clf()
    plt.ylim([0,2])
    plt.xticks([])
    plt.errorbar(xax,Ffullmed, yerr=[Ffullmed-Ffulllo,Ffullhi-Ffullmed],
      fmt='ro',mec='r',label='Full sky', alpha=.5, capsize=5)
    plt.errorbar(xax,F6dirmed, yerr=[F6dirmed-F6dirlo,F6dirhi-F6dirmed],
      fmt='co',mec='c',label='Six dirs', alpha=.5, capsize=5)
  # plt.errorbar(xax,F6dirmed-Ffullmed,
  #     yerr=[(F6dirmed-Ffullmed) - (F6dirlo-Ffulllo),
  #           (F6dirhi -Ffullhi)  - (F6dirmed-Ffullmed)],
  #     fmt='go',mec='g',label='[6 dir] - [full sky]', alpha=1, capsize=5)
    plt.legend()

#------------------------------------------------------------------------------

def get_6dirs(galid,
              label      = '',
              cosmo      = cosmoJ,
              endianness = 'big'
              ):
    """
    Get the "med,lo,hi" of the six directions of a MoCaLaTA observation,
    calculated as [mean of two mid values], [next-lowest value], [next-highest
    value].
    """
    mother   = '/Users/pela/Projects/z8p8/CHAclz8p77rz12_512x64_Alex/'
    galdir   = mother + 'gal' + galid[:2] + '/' + galid[2:] + '/'
    views    = ['xm','xp','ym','yp','zm','zp']
    binfiles = [galdir+v+label+'.bin' for v in views]

    vals = np.empty(6)
    for i,b in enumerate(binfiles):
        par,spec = readMoCaLaTA(b,skipIFU=True,endianness=endianness,cosmo=cosmo)
        vals[i]  = spec.sum() / par['n_eff']

    return np.percentile(vals, [50,15.87,84.17])
#------------------------------------------------------------------------------

def mkLine(fwhm     = 2870,
           z        = 2.631,
           lam      = np.linspace(4360,4460,1001),
           Ftot     = 1e-14 /1.8,
           Flya_tot = 3e-15 /1.8
           ):
    from scipy.integrate import simps
    import mymath

    lam0    = 1215.67
    lam0obs = lam0 * (1+z)
    fwhmobs = fwhm/3e5 * lam0obs
    sigobs  = fwhmobs / 2.355
    Flya    = Flya_tot * mymath.gauss(lam,lam0obs,sigobs) #b/c gauss is normalized
    F0_tot  = Ftot - Flya_tot
    unity   = np.ones_like(lam)
    F0      = unity / simps(unity,x=lam) * F0_tot
    F       = F0 + Flya
  # EW      = simps(F/F0-1,x=lam)
    EW      = simps((F-F0)/F0,x=lam)
    EWrest  = EW / (1+z)
    print('EWrest =', EWrest)

    plt.clf()
  # plt.ylim([0,1.2*Flam.max()])
    plt.ylim([0,8.213e-17])
    plt.xlim([lam.min(),lam.max()])
    plt.plot(lam,F)
    plt.plot(lam,F0,'k--',alpha=.5)
#------------------------------------------------------------------------------

def flux_ratio(galdir1, view1, ap1, ap2,
               galdir2  = None,
               view2    = None,
               seeing1  = 0,
               seeing2  = None,
               ):
    """
    Get the flux ratio between two apertures ap1 and ap2.
    """
    if galdir2 == None: galdir2 = galdir1
    if view2   == None: view2   = view1
    if seeing2 == None: seeing2 = seeing1

    Ftot1 = showobs(galdir1,view1,seeing=seeing1,spatial_aperture=ap1,get_value='Ftot',showit=False)
    Ftot2 = showobs(galdir2,view2,seeing=seeing2,spatial_aperture=ap2,get_value='Ftot',showit=False)

    return Ftot1 / Ftot2
#------------------------------------------------------------------------------
