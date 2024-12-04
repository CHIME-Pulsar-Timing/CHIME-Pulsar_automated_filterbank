from builtins import range
from builtins import object
import sys
import numpy as Num
import copy, random, struct
from presto_without_presto import psr_utils, infodata, polycos
from matplotlib import pyplot as plt
from matplotlib import colors as mplc
from matplotlib import ticker as mpltick
import six
import numbers
from presto_without_presto.bestprof import bestprof
from numba import jit

try:
    from presto.presto import chi2_sigma
except ImportError:
    print("Couldn't import chi2_sigma from presto. Implemented workaround may give slight differences")
    from scipy.stats import chi2 as scipychi2 

    def chi2_sigma(chi2, dof):
        """
        Return the approximate significance in Gaussian sigmas
        sigmas of a chi^2 value of chi2 given dof degrees of freedom.
        """
        if chi2 <= 0.0:
            return 0.0
        logp = scipychi2.logsf(chi2, dof)
        return psr_utils.extended_equiv_gaussian_sigma(logp)


def set_labels(ax, labx="", laby=""):
    if labx:
        ax.set_xlabel(labx)
    if laby:
        ax.set_ylabel(laby)

# set some default plotting parameters
params = {
    "xtick.top": True,
    "ytick.right": True,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.labelsize": "small",
    "ytick.labelsize": "small",
    "axes.labelsize": "small",
}
plt.rcParams.update(params)

class PGPalette(object):
    def __init__(self, palette):
        self.setpalette(palette)

    # Set the color palette
    # copied from presto Pgplot.py
    def setpalette(self, palette):
        """
        setpalette(self, palette):
            Set the color palette for imag-style routines
        """
        if (palette == 'rainbow'):
            self.l = Num.array([0.0, 0.015, 0.225, 0.4, 0.59,
                                0.6, 0.775, 0.955, 0.965, 1.0])
            self.r = Num.array([1.0, 1.0, 1.0, 0.0, 0.0,
                                0.0, 0.0, 0.947, 1.0, 1.0])
            self.g = Num.array([0.0, 0.0, 1.0, 1.0, 1.0,
                                0.946, 0.0, 0.8, 0.844, 1.0])
            self.b = Num.array([0.0, 0.0, 0.0, 0.0, 0.95,
                                1.0, 1.0, 1.0, 1.0, 1.0])
        elif (palette == 'antirainbow'):
            self.l = Num.array([0.0, 0.035, 0.045, 0.225, 0.4,
                                0.41, 0.6, 0.775, 0.985, 1.0])
            self.r = Num.array([1.0, 1.0, 0.947, 0.0, 0.0,
                                0.0, 0.0, 1.0, 1.0, 1.0])
            self.g = Num.array([1.0, 0.844, 0.8, 0.0, 0.946,
                                1.0, 1.0, 1.0, 0.0, 0.0])
            self.b = Num.array([1.0, 1.0, 1.0, 1.0, 1.0,
                                0.95, 0.0, 0.0, 0.0, 0.0])
        elif (palette == 'astro'):
            self.l = Num.array([0.0, 0.167, 0.333, 0.5,
                                0.667, 0.833, 1.0])
            self.r = Num.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
            self.g = Num.array([0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0])
            self.b = Num.array([0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0])
        elif (palette == 'hue'):
            self.l = Num.array([0.0, 0.167, 0.333, 0.5,
                                0.667, 0.833, 1.0])
            self.r = Num.array([1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0])
            self.g = Num.array([0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
            self.b = Num.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0])
        elif (palette == 'heat'):
            self.l = Num.array([0.0, 0.48, 0.7, 0.75, 1.0])
            self.r = Num.array([0.0, 1.0, 1.0, 1.0, 1.0])
            self.g = Num.array([0.0, 0.0, 0.423, 0.519, 1.0])
            self.b = Num.array([0.0, 0.0, 0.0, 0.0, 1.0])
        elif (palette == 'gamma'):
            self.l = Num.array([0.0, 0.33, 0.66, 1.0])
            self.r = Num.array([0.3, 1.0, 0.0, 0.0])
            self.g = Num.array([0.0, 0.3, 1.0, 0.0])
            self.b = Num.array([0.0, 0.0, 0.3, 1.0])
        elif (palette == 'antigray' or palette == 'antigrey'):
            self.l = Num.array([0.0, 1.0])
            self.r = Num.array([1.0, 0.0])
            self.g = Num.array([1.0, 0.0])
            self.b = Num.array([1.0, 0.0])
        elif (palette == 'apjgray' or palette == 'apjgrey'):
            self.l = Num.array([0.0, 1.0])
            self.r = Num.array([1.0, 0.25])
            self.g = Num.array([1.0, 0.25])
            self.b = Num.array([1.0, 0.25])
        else:  # altered this so can try a matplotlib cmap instead
            raise AttributeError(f"{palette} is not a defined pgplot palette in PGPalette")

    @property
    def cmap(self):
        """Convert to matplotlib colormap"""
        cdict = { 
            "red": [(l, r, r) for l,r in zip(self.l, self.r)],
            "green": [(l, g, g) for l,g in zip(self.l, self.g)],
            "blue": [(l, b, b) for l,b in zip(self.l, self.b)]
        }

        return mplc.LinearSegmentedColormap("", cdict)
    
    @property
    def colors(self):
        """Colors as a list of (r, g, b) tuples"""
        return [(r,g,b) for r,g,b in zip(self.r, self.g, self.b)]

class pfd(object):

    def __init__(self, filename):
        self.pfd_filename = filename
        infile = open(filename, "rb")
        # See if the .bestprof file is around
        try:
            self.bestprof = bestprof(filename+".bestprof")
        except IOError:
            self.bestprof = 0
        swapchar = '<' # this is little-endian
        data = infile.read(5*4)
        testswap = struct.unpack(swapchar+"i"*5, data)
        # This is a hack to try and test the endianness of the data.
        # None of the 5 values should be a large positive number.
        if (Num.fabs(Num.asarray(testswap))).max() > 100000:
            swapchar = '>' # this is big-endian
        (self.numdms, self.numperiods, self.numpdots, self.nsub, self.npart) = \
                      struct.unpack(swapchar+"i"*5, data)
        (self.proflen, self.numchan, self.pstep, self.pdstep, self.dmstep, \
         self.ndmfact, self.npfact) = struct.unpack(swapchar+"i"*7, infile.read(7*4))

        self.filenm = infile.read(struct.unpack(swapchar+"i", infile.read(4))[0])
        self.candnm = infile.read(struct.unpack(swapchar+"i", infile.read(4))[0]).decode("utf-8")
        self.telescope = infile.read(struct.unpack(swapchar+"i", infile.read(4))[0]).decode("utf-8")
        self.pgdev = infile.read(struct.unpack(swapchar+"i", infile.read(4))[0])
        test = infile.read(16)
        if not test[:8]==b"Unknown" and b':' in test:
            self.rastr = test[:test.find(b'\0')]
            test = infile.read(16)
            self.decstr = test[:test.find(b'\0')]
        else:
            self.rastr = "Unknown"
            self.decstr = "Unknown"
            if ':' not in test:
                infile.seek(-16, 1) # rewind the file before the bad read
        (self.dt, self.startT) = struct.unpack(swapchar+"dd", infile.read(2*8))
        (self.endT, self.tepoch, self.bepoch, self.avgvoverc, self.lofreq, \
         self.chan_wid, self.bestdm) = struct.unpack(swapchar+"d"*7, infile.read(7*8))
        # The following "fixes" (we think) the observing frequency of the Spigot
        # based on tests done by Ingrid on 0737 (comparing it to GASP)
        # The same sorts of corrections should be made to WAPP data as well...
        # The tepoch corrections are empirically determined timing corrections
        # Note that epoch is only double precision and so the floating
        # point accuracy is ~1 us!
        if self.telescope=='GBT':
            if (Num.fabs(Num.fmod(self.dt, 8.192e-05) < 1e-12) and \
                ("spigot" in filename.lower() or "guppi" not in filename.lower()) and \
                (self.tepoch < 54832.0)):
                sys.stderr.write("Assuming SPIGOT data...\n")
                if self.chan_wid==800.0/1024: # Spigot 800 MHz mode 2
                    self.lofreq -= 0.5 * self.chan_wid
                    # original values
                    #if self.tepoch > 0.0: self.tepoch += 0.039334/86400.0
                    #if self.bestprof: self.bestprof.epochf += 0.039334/86400.0
                    # values measured with 1713+0747 wrt BCPM2 on 13 Sept 2007
                    if self.tepoch > 0.0: self.tepoch += 0.039365/86400.0
                    if self.bestprof: self.bestprof.epochf += 0.039365/86400.0
                elif self.chan_wid==800.0/2048:
                    self.lofreq -= 0.5 * self.chan_wid
                    if self.tepoch < 53700.0:  # Spigot 800 MHz mode 16 (downsampled)
                        if self.tepoch > 0.0: self.tepoch += 0.039352/86400.0
                        if self.bestprof: self.bestprof.epochf += 0.039352/86400.0
                    else:  # Spigot 800 MHz mode 14
                        # values measured with 1713+0747 wrt BCPM2 on 13 Sept 2007
                        if self.tepoch > 0.0: self.tepoch += 0.039365/86400.0
                        if self.bestprof: self.bestprof.epochf += 0.039365/86400.0
                elif self.chan_wid==50.0/1024 or self.chan_wid==50.0/2048: # Spigot 50 MHz modes
                    self.lofreq += 0.5 * self.chan_wid
                    # Note: the offset has _not_ been measured for the 2048-lag mode
                    if self.tepoch > 0.0: self.tepoch += 0.039450/86400.0
                    if self.bestprof: self.bestprof.epochf += 0.039450/86400.0
        (self.topo_pow, tmp) = struct.unpack(swapchar+"f"*2, infile.read(2*4))
        (self.topo_p1, self.topo_p2, self.topo_p3) = struct.unpack(swapchar+"d"*3, \
                                                                   infile.read(3*8))
        (self.bary_pow, tmp) = struct.unpack(swapchar+"f"*2, infile.read(2*4))
        (self.bary_p1, self.bary_p2, self.bary_p3) = struct.unpack(swapchar+"d"*3, \
                                                                   infile.read(3*8))
        (self.fold_pow, tmp) = struct.unpack(swapchar+"f"*2, infile.read(2*4))
        (self.fold_p1, self.fold_p2, self.fold_p3) = struct.unpack(swapchar+"d"*3, \
                                                                   infile.read(3*8))
        # Save current p, pd, pdd
        # NOTE: Fold values are actually frequencies!
        self.curr_p1, self.curr_p2, self.curr_p3 = \
                psr_utils.p_to_f(self.fold_p1, self.fold_p2, self.fold_p3)
        self.pdelays_bins = Num.zeros(self.npart, dtype='d')
        (self.orb_p, self.orb_e, self.orb_x, self.orb_w, self.orb_t, self.orb_pd, \
         self.orb_wd) = struct.unpack(swapchar+"d"*7, infile.read(7*8))
        self.dms = Num.asarray(struct.unpack(swapchar+"d"*self.numdms, \
                                             infile.read(self.numdms*8)))
        if self.numdms==1:
            self.dms = self.dms[0]
        self.periods = Num.asarray(struct.unpack(swapchar+"d"*self.numperiods, \
                                                 infile.read(self.numperiods*8)))
        self.pdots = Num.asarray(struct.unpack(swapchar+"d"*self.numpdots, \
                                               infile.read(self.numpdots*8)))
        self.numprofs = self.nsub*self.npart
        if (swapchar=='<'):  # little endian
            self.profs = Num.zeros((self.npart, self.nsub, self.proflen), dtype='d')
            for ii in range(self.npart):
                for jj in range(self.nsub):
                    self.profs[ii,jj,:] = Num.fromfile(infile, Num.float64, self.proflen)
        else:
            self.profs = Num.asarray(struct.unpack(swapchar+"d"*self.numprofs*self.proflen, \
                                                   infile.read(self.numprofs*self.proflen*8)))
            self.profs = Num.reshape(self.profs, (self.npart, self.nsub, self.proflen))
        if (self.numchan==1):
            try:
                idata = infodata.infodata(self.filenm[:self.filenm.rfind(b'.')]+b".inf")
                try:
                    if idata.waveband=="Radio":
                        self.bestdm = idata.DM
                        self.numchan = idata.numchan
                except:
                        self.bestdm = 0.0
                        self.numchan = 1
            except IOError:
                print("Warning!  Can't open the .inf file for "+filename+"!")
        self.binspersec = self.fold_p1*self.proflen
        self.chanpersub = self.numchan // self.nsub
        self.subdeltafreq = self.chan_wid*self.chanpersub
        self.hifreq = self.lofreq + (self.numchan-1)*self.chan_wid
        self.losubfreq = self.lofreq + self.subdeltafreq - self.chan_wid
        self.subfreqs = Num.arange(self.nsub, dtype='d')*self.subdeltafreq + \
                        self.losubfreq
        self.subdelays_bins = Num.zeros(self.nsub, dtype='d')
        # Save current DM
        self.currdm = 0
        self.killed_subbands = []
        self.killed_intervals = []
        self.pts_per_fold = []
        # Note: a foldstats struct is read in as a group of 7 doubles
        # the correspond to, in order:
        #    numdata, data_avg, data_var, numprof, prof_avg, prof_var, redchi
        self.stats = Num.zeros((self.npart, self.nsub, 7), dtype='d')
        for ii in range(self.npart):
            currentstats = self.stats[ii]
            for jj in range(self.nsub):
                if (swapchar=='<'):  # little endian
                    currentstats[jj] = Num.fromfile(infile, Num.float64, 7)
                else:
                    currentstats[jj] = Num.asarray(struct.unpack(swapchar+"d"*7, \
                                                                 infile.read(7*8)))
            self.pts_per_fold.append(self.stats[ii][0][0])  # numdata from foldstats
        self.start_secs = Num.add.accumulate([0]+self.pts_per_fold[:-1])*self.dt
        self.pts_per_fold = Num.asarray(self.pts_per_fold)
        self.mid_secs = self.start_secs + 0.5*self.dt*self.pts_per_fold
        if (not self.tepoch==0.0):
            self.start_topo_MJDs = self.start_secs/86400.0 + self.tepoch
            self.mid_topo_MJDs = self.mid_secs/86400.0 + self.tepoch
        if (not self.bepoch==0.0):
            self.start_bary_MJDs = self.start_secs/86400.0 + self.bepoch
            self.mid_bary_MJDs = self.mid_secs/86400.0 + self.bepoch
        self.Nfolded = Num.add.reduce(self.pts_per_fold)
        self.T = self.Nfolded*self.dt
        self.avgprof = (self.profs/self.proflen).sum()
        self.varprof = self.calc_varprof()
        # nominal number of degrees of freedom for reduced chi^2 calculation
        self.DOFnom = float(self.proflen) - 1.0
        # corrected number of degrees of freedom due to inter-bin correlations
        self.dt_per_bin = self.curr_p1 / self.proflen / self.dt
        self.DOFcor = self.DOFnom * self.DOF_corr()
        infile.close()
        self.barysubfreqs = None
        # KC: changed this so doppler corrects if self.avgvoverc is initially nonzero
        if self.avgvoverc != 0:
            # Make the Doppler correction
            #self.barysubfreqs = None
            self.barysubfreqs = self.subfreqs*(1.0+self.avgvoverc)
        elif self.candnm.startswith("PSR_"):
            try:
                psrname = self.candnm[4:]
                self.polycos = polycos.polycos(psrname,
                                               filenm=self.pfd_filename+".polycos")
                midMJD = self.tepoch + 0.5*self.T/86400.0
                self.avgvoverc = self.polycos.get_voverc(int(midMJD), midMJD-int(midMJD))
                #sys.stderr.write("Approximate Doppler velocity (in c) is:  %.4g\n"%self.avgvoverc)
                # Make the Doppler correction
                self.barysubfreqs = self.subfreqs*(1.0+self.avgvoverc)
            except IOError:
                self.polycos = 0
        if self.barysubfreqs is None:
            self.barysubfreqs = self.subfreqs

    def __str__(self):
        out = ""
        for k, v in list(self.__dict__.items()):
            if k[:2]!="__":
                if isinstance(self.__dict__[k], six.string_types):
                    out += "%10s = '%s'\n" % (k, v)
                elif isinstance(self.__dict__[k], numbers.Integral):
                    out += "%10s = %d\n" % (k, v)
                elif isinstance(self.__dict__[k], numbers.Real):
                    out += "%10s = %-20.15g\n" % (k, v)
        return out

    def dedisperse(self, DM=None, interp=False, doppler=False):
        """
        dedisperse(DM=self.bestdm, interp=False, doppler=False):
            Rotate (internally) the profiles so that they are de-dispersed
                at a dispersion measure of DM.  Use FFT-based interpolation if
                'interp' is non-zero (NOTE: It is off by default!).
                Doppler shift subband frequencies if doppler is non-zero.
                (NOTE: It is off by default. However, if you fold raw data
                for a search candidate (and let it search), prepfold *does*
                doppler correct the frequencies! It does *not* doppler
                correct if you fold with polycos for timing, for instance.)
        """
        if DM is None:
            DM = self.bestdm
        # Note:  Since TEMPO Doppler corrects observing frequencies, for
        #        TOAs, at least, we need to de-disperse using topocentric
        #        observing frequencies.
        if doppler:
            freqs = psr_utils.doppler(self.subfreqs, self.avgvoverc)
        else:
            freqs = self.subfreqs
        self.subdelays = psr_utils.delay_from_DM(DM, freqs)
        self.hifreqdelay = self.subdelays[-1]
        self.subdelays = self.subdelays-self.hifreqdelay
        delaybins = self.subdelays*self.binspersec - self.subdelays_bins
        if interp:
            new_subdelays_bins = delaybins
            for ii in range(self.npart):
                for jj in range(self.nsub):
                    tmp_prof = self.profs[ii,jj,:]
                    self.profs[ii,jj] = psr_utils.fft_rotate(tmp_prof, delaybins[jj])
            # Note: Since the rotation process slightly changes the values of the
            # profs, we need to re-calculate the average profile value
            self.avgprof = (self.profs/self.proflen).sum()
        else:
            new_subdelays_bins = Num.floor(delaybins+0.5)
            for ii in range(self.nsub):
                rotbins = int(new_subdelays_bins[ii])%self.proflen
                if rotbins:  # i.e. if not zero
                    subdata = self.profs[:,ii,:]
                    self.profs[:,ii] = Num.concatenate((subdata[:,rotbins:],
                                                        subdata[:,:rotbins]), 1)
        self.subdelays_bins += new_subdelays_bins
        self.sumprof = self.profs.sum(0).sum(0)
        if Num.fabs((self.sumprof/self.proflen).sum() - self.avgprof) > 1.0:
            print("self.avgprof is not the correct value!")
        self.currdm = DM

    def freq_offsets(self, p=None, pd=None, pdd=None):
        """
        freq_offsets(p=*bestp*, pd=*bestpd*, pdd=*bestpdd*):
            Return the offsets between given frequencies
            and fold frequencies.

            If p, pd or pdd are None use the best values.

            A 3-tuple is returned.
        """
        if self.fold_pow == 1.0:
            bestp = self.bary_p1
            bestpd = self.bary_p2
            bestpdd = self.bary_p3
        else:
            if self.topo_p1 == 0.0:
                bestp = self.fold_p1
                bestpd = self.fold_p2
                bestpdd = self.fold_p3
            else:
                bestp = self.topo_p1
                bestpd = self.topo_p2
                bestpdd = self.topo_p3
        if p is not None:
            bestp = p
        if pd is not None:
            bestpd = pd
        if pdd is not None:
            bestpdd = pdd

        # self.fold_p[123] are actually frequencies, convert to periods
        foldf, foldfd, foldfdd = self.fold_p1, self.fold_p2, self.fold_p3
        foldp, foldpd, foldpdd = psr_utils.p_to_f(self.fold_p1, \
                                        self.fold_p2, self.fold_p3)

        # Get best f, fd, fdd
        # Use folding values to be consistent with prepfold_plot.c
        bestfdd = psr_utils.p_to_f(foldp, foldpd, bestpdd)[2]
        bestfd = psr_utils.p_to_f(foldp, bestpd)[1]
        bestf = 1.0/bestp

        # Get frequency and frequency derivative offsets
        f_diff = bestf - foldf
        fd_diff = bestfd - foldfd

        # bestpdd=0.0 only if there was no searching over pdd
        if bestpdd != 0.0:
            fdd_diff = bestfdd - foldfdd
        else:
            fdd_diff = 0.0

        return (f_diff, fd_diff, fdd_diff)

    def DOF_corr(self):
        """
        DOF_corr():
            Return a multiplicative correction for the effective number of
            degrees of freedom in the chi^2 measurement resulting from a
            pulse profile folded by PRESTO's fold() function
            (i.e. prepfold).  This is required because there are
            correlations between the bins caused by the way that prepfold
            folds data (i.e. treating a sample as finite duration and
            smearing it over potenitally several bins in the profile as
            opposed to instantaneous and going into just one profile bin).
            The correction is semi-analytic (thanks to Paul Demorest and
            Walter Brisken) but the values for 'power' and 'factor' have
            been determined from Monte Carlos.  The correction is good to
            a fractional error of less than a few percent as long as
            dt_per_bin is > 0.5 or so (which it usually is for pulsar
            candidates).  There is a very minimal number-of-bins
            dependence, which is apparent when dt_per_bin < 0.7 or so.
            dt_per_bin is the width of a profile bin in samples (a float),
            and so for prepfold is pulse period / nbins / sample time.  Note
            that the sqrt of this factor can be used to 'inflate' the RMS
            of the profile as well, for radiometer eqn flux density estimates,
            for instance.
        """
        power, factor = 1.806, 0.96  # From Monte Carlo
        return self.dt_per_bin * factor * \
               (1.0 + self.dt_per_bin**(power))**(-1.0/power)

    def use_for_timing(self):
        """
        use_for_timing():
            This method returns True or False depending on whether
            the .pfd file can be used for timing or not.  For this
            to return true, the pulsar had to have been folded with
            a parfile and -no[p/pd]search (this includes -timing), or
            with a p/pdot/pdotdot and a corresponding -no[p/pd]search.
            In other words, if you let prepfold search for the best
            p/pdot/pdotdot, you will get bogus TOAs if you try timing
            with it.
        """
        T = self.T
        bin_dphi = 1.0/self.proflen
        # If any of the offsets causes more than a 0.1-bin rotation over
        # the obs, then prepfold searched and we can't time using it
        # Allow up to a 0.5 bin shift for pdd/fdd since the conversions
        # back and forth can cause float issues.
        offsets = Num.fabs(Num.asarray(self.freq_offsets()))
        dphis = offsets * Num.asarray([T, T**2.0/2.0, T**3.0/6.0])
        if max(dphis[:2]) > 0.1 * bin_dphi or dphis[2] > 0.5 * bin_dphi:
            return False
        else:
            return True

    def time_vs_phase(self, p=None, pd=None, pdd=None, interp=0):
        """
        time_vs_phase(p=*bestp*, pd=*bestpd*, pdd=*bestpdd*):
            Return the 2D time vs. phase profiles shifted so that
                the given period and period derivative are applied.
                Use FFT-based interpolation if 'interp' is non-zero.
                (NOTE: It is off by default as in prepfold!).
        """
        # Cast to single precision and back to double precision to
        # emulate prepfold_plot.c, where parttimes is of type "float"
        # but values are upcast to "double" during computations.
        # (surprisingly, it affects the resulting profile occasionally.)
        parttimes = self.start_secs.astype('float32').astype('float64')

        # Get delays
        f_diff, fd_diff, fdd_diff = self.freq_offsets(p, pd, pdd)
        #print "DEBUG: in myprepfold.py -- parttimes", parttimes
        delays = psr_utils.delay_from_foffsets(f_diff, fd_diff, fdd_diff, parttimes)

        # Convert from delays in phase to delays in bins
        bin_delays = Num.fmod(delays * self.proflen, self.proflen) - self.pdelays_bins

        # Rotate subintegrations
        # subints = self.combine_profs(self.npart, 1)[:,0,:] # Slower than sum by ~9x
        subints = Num.sum(self.profs, axis=1).squeeze()
        if interp:
            new_pdelays_bins = bin_delays
            for ii in range(self.npart):
                tmp_prof = subints[ii,:]
                # Negative sign in num bins to shift because we calculated delays
                # Assuming +ve is shift-to-right, psr_utils.rotate assumes +ve
                # is shift-to-left
                subints[ii,:] = psr_utils.fft_rotate(tmp_prof, -new_pdelays_bins[ii])
        else:
            new_pdelays_bins = Num.floor(bin_delays+0.5)
            indices = Num.outer(Num.arange(self.proflen), Num.ones(self.npart))
            indices = Num.mod(indices-new_pdelays_bins, self.proflen).T
            indices += Num.outer(Num.arange(self.npart)*self.proflen, \
                                    Num.ones(self.proflen))
            subints = subints.flatten('C')[indices.astype('i8')]
        return subints

    def adjust_period(self, p=None, pd=None, pdd=None, interp=0):
        """
        adjust_period(p=*bestp*, pd=*bestpd*, pdd=*bestpdd*):
            Rotate (internally) the profiles so that they are adjusted to
                the given period and period derivatives.  By default,
                use the 'best' values as determined by prepfold's seaqrch.
                This should orient all of the profiles so that they are
                almost identical to what you see in a prepfold plot which
                used searching.  Use FFT-based interpolation if 'interp'
                is non-zero.  (NOTE: It is off by default, as in prepfold!)
        """
        if self.fold_pow == 1.0:
            bestp = self.bary_p1
            bestpd = self.bary_p2
            bestpdd = self.bary_p3
        else:
            bestp = self.topo_p1
            bestpd = self.topo_p2
            bestpdd = self.topo_p3
        if p is None:
            p = bestp
        if pd is None:
            pd = bestpd
        if pdd is None:
            pdd = bestpdd

        # Cast to single precision and back to double precision to
        # emulate prepfold_plot.c, where parttimes is of type "float"
        # but values are upcast to "double" during computations.
        # (surprisingly, it affects the resulting profile occasionally.)
        parttimes = self.start_secs.astype('float32').astype('float64')

        # Get delays
        f_diff, fd_diff, fdd_diff = self.freq_offsets(p, pd, pdd)
        delays = psr_utils.delay_from_foffsets(f_diff, fd_diff, fdd_diff, parttimes)

        # Convert from delays in phase to delays in bins
        bin_delays = Num.fmod(delays * self.proflen, self.proflen) - self.pdelays_bins
        if interp:
            new_pdelays_bins = bin_delays
        else:
            new_pdelays_bins = Num.floor(bin_delays+0.5)

        # Rotate subintegrations
        for ii in range(self.nsub):
            for jj in range(self.npart):
                tmp_prof = self.profs[jj,ii,:]
                # Negative sign in num bins to shift because we calculated delays
                # Assuming +ve is shift-to-right, psr_utils.rotate assumes +ve
                # is shift-to-left
                if interp:
                    self.profs[jj,ii] = psr_utils.fft_rotate(tmp_prof, -new_pdelays_bins[jj])
                else:
                    self.profs[jj,ii] = psr_utils.rotate(tmp_prof, \
                                            -new_pdelays_bins[jj])
        self.pdelays_bins += new_pdelays_bins
        if interp:
            # Note: Since the rotation process slightly changes the values of the
            # profs, we need to re-calculate the average profile value
            self.avgprof = (self.profs/self.proflen).sum()

        self.sumprof = self.profs.sum(0).sum(0)
        if Num.fabs((self.sumprof/self.proflen).sum() - self.avgprof) > 1.0:
            print("self.avgprof is not the correct value!")

        # Save current p, pd, pdd
        self.curr_p1, self.curr_p2, self.curr_p3 = p, pd, pdd

    def adjust_period_more_efficient(self, p=None, pd=None, pdd=None):
        """
        adjust_period(p=*bestp*, pd=*bestpd*, pdd=*bestpdd*):
            Rotate (internally) the profiles so that they are adjusted to
                the given period and period derivatives.  By default,
                use the 'best' values as determined by prepfold's seaqrch.
                This should orient all of the profiles so that they are
                almost identical to what you see in a prepfold plot which
                used searching.  
                
                No FFT-based interpolation option (default option in prepfold is for
                it to be off)
        """
        if self.fold_pow == 1.0:
            bestp = self.bary_p1
            bestpd = self.bary_p2
            bestpdd = self.bary_p3
        else:
            bestp = self.topo_p1
            bestpd = self.topo_p2
            bestpdd = self.topo_p3
        if p is None:
            p = bestp
        if pd is None:
            pd = bestpd
        if pdd is None:
            pdd = bestpdd

        # Cast to single precision and back to double precision to
        # emulate prepfold_plot.c, where parttimes is of type "float"
        # but values are upcast to "double" during computations.
        # (surprisingly, it affects the resulting profile occasionally.)
        parttimes = self.start_secs.astype('float32').astype('float64')

        # Get delays
        f_diff, fd_diff, fdd_diff = self.freq_offsets(p, pd, pdd)
        delays = psr_utils.delay_from_foffsets(f_diff, fd_diff, fdd_diff, parttimes)

        # Convert from delays in phase to delays in bins
        bin_delays = Num.fmod(delays * self.proflen, self.proflen) - self.pdelays_bins

        new_pdelays_bins = Num.floor(bin_delays+0.5)
        new_pdelays_bins = new_pdelays_bins.astype(int) % self.proflen

        for jj in range(self.npart):
            self.profs[jj,:,:] = Num.roll(self.profs[jj,:,:], new_pdelays_bins[jj], axis=1)
        self.pdelays_bins += new_pdelays_bins

        self.sumprof = self.profs.sum(0).sum(0)
        if Num.fabs((self.sumprof/self.proflen).sum() - self.avgprof) > 1.0:
            print("self.avgprof is not the correct value!")

        # Save current p, pd, pdd
        self.curr_p1, self.curr_p2, self.curr_p3 = p, pd, pdd

    def combine_profs(self, new_npart, new_nsub):
        """
        combine_profs(self, new_npart, new_nsub):
            Combine intervals and/or subbands together and return a new
                array of profiles.
        """
        if (self.npart % new_npart):
            print("Warning!  The new number of intervals (%d) is not a" % new_npart)
            print("          divisor of the original number of intervals (%d)!"  % self.npart)
            print("Doing nothing.")
            return None
        if (self.nsub % new_nsub):
            print("Warning!  The new number of subbands (%d) is not a" % new_nsub)
            print("          divisor of the original number of subbands (%d)!"  % self.nsub)
            print("Doing nothing.")
            return None

        dp = self.npart // new_npart
        ds = self.nsub // new_nsub

        newprofs = Num.zeros((new_npart, new_nsub, self.proflen), 'd')
        for ii in range(new_npart):
            # Combine the subbands if required
            if (self.nsub > 1):
                for jj in range(new_nsub):
                    subprofs = Num.add.reduce(self.profs[:,jj*ds:(jj+1)*ds], 1)
                    # Combine the time intervals
                    newprofs[ii][jj] = Num.add.reduce(subprofs[ii*dp:(ii+1)*dp])
            else:
                newprofs[ii][0] = Num.add.reduce(self.profs[ii*dp:(ii+1)*dp,0])
        return newprofs

    def kill_intervals(self, intervals):
        """
        kill_intervals(intervals):
            Set all the subintervals (internally) from the list of
                subintervals to all zeros, effectively 'killing' them.
        """
        for part in intervals:
            self.profs[part,:,:] *= 0.0
            self.killed_intervals.append(part)
        # Update the stats
        self.avgprof = (self.profs/self.proflen).sum()
        self.varprof = self.calc_varprof()

    def kill_subbands(self, subbands):
        """
        kill_subbands(subbands):
            Set all the profiles (internally) from the list of
                subbands to all zeros, effectively 'killing' them.
        """
        for sub in subbands:
            self.profs[:,sub,:] *= 0.0
            self.killed_subbands.append(sub)
        # Update the stats
        self.avgprof = (self.profs/self.proflen).sum()
        self.varprof = self.calc_varprof()

    def plot_sumprof(self, ax=None):
        """
        plot_sumprof(self, device='/xwin'):
            Plot the dedispersed and summed profile.
        """
        show = False
        if ax is None:
            show = True
            fig, ax = plt.subplots()

        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()
        normprof = self.sumprof - min(self.sumprof)
        normprof /= max(normprof)
        self.plotxy(normprof, labx="Phase Bins", laby="Normalized Flux",
                    ax=ax)

        if show:
            plt.show()
            plt.close()

    def plot_sumprof_prepfold(self, ax=None):
        """
        KC: enough things needed to be tweaked that for now easier to make this
        plot_sumprof(self, device='/xwin'):
            Plot the dedispersed and summed profile.
        """
        show = False
        if ax is None:
            show = True
            fig, ax = plt.subplots()

        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()
        normprof = self.sumprof - min(self.sumprof)
        normprof /= max(normprof)

        lenprof = normprof.shape[0]

        normprof_doubled = Num.zeros((lenprof*2))
        normprof_doubled[:lenprof] = normprof
        normprof_doubled[lenprof:] = normprof

        ax.axhline(normprof.mean(), linestyle=(0, (3, 5)), c="black", linewidth=1)
        # cannot figure out how to offset this like in the prepfold plot
        # hopefully being light, a different colour, and in the background works?
        ax.errorbar(0.1, normprof.mean(), yerr=normprof.std(), capsize=4, ecolor="red", alpha=0.5)

        self.plotxy(normprof_doubled, x=Num.linspace(0,2,lenprof*2), ax=ax)

        if show:
            plt.show()
            plt.close()

    def plot_chi2_vs_int(self, ax=None):
        """
        Based off plot_chi2_vs_sub, but using cumulative profiles
        """
        show = False
        if ax is None:
            show = True
            fig, ax = plt.subplots()

        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()

        # sum over subbands (=>(npart, nbin)), then take sum up to each interval
        # this is a tad dodgy since it's not ignoring killed intervals or subbands when making the profiles
        profs = Num.cumsum(self.profs.sum(1), 0)

        # Compute the averages and variances for the intervals
        avgs = profs.sum(1)/self.proflen
        vars = []
        var = 0.0
        for part in range(self.npart):
            if part in self.killed_intervals:
                vars.append(var)
                continue
            for sub in range(self.nsub):
                if sub in self.killed_subbands:
                    continue
                var += self.stats[part][sub][5]
            vars.append(var)
        chis = Num.zeros(self.npart, dtype='f')
        for ii in range(self.npart):
            chis[ii] = self.calc_redchi2(prof=profs[ii], avg=avgs[ii], var=vars[ii])
 
        # Now plot it
        self.plotxy(chis, labx="Time Intervals", laby=r"Reduced $\chi^2$", rangey=[0.0, max(chis)*1.1], ax=ax)

        if show:
            plt.show()
            plt.close()

        return chis

    def plot_chi2_vs_int_prepfold(self, ax=None):
        """
        based off plot_chi2_vs_sub, but using cumulative profiles
        """
        show = False
        if ax is None:
            show = True
            fig, ax = plt.subplots()

        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()

        # sum over subbands (=>(npart, nbin)), then take sum up to each interval
        # this is a tad dodgy since it's not ignoring killed intervals or subbands when making the profiles
        profs = Num.cumsum(self.profs.sum(1), 0)
        # Compute the averages and variances for the intervals
        avgs = profs.sum(1)/self.proflen
        vars = []
        var = 0.0
        for part in range(self.npart):
            if part in self.killed_intervals:
                vars.append(var)
                continue
            for sub in range(self.nsub):
                if sub in self.killed_subbands:
                    continue
                var += self.stats[part][sub][5]
            vars.append(var)
        chis = Num.zeros(self.npart, dtype='f')
        for ii in range(self.npart):
            chis[ii] = self.calc_redchi2(prof=profs[ii], avg=avgs[ii], var=vars[ii])

        # Now plot it
        self.plotxy(Num.linspace(0,1,self.npart), chis, labx=r"Reduced $\chi^2$", laby="Fraction of Observation", rangey=[0.0, 1.0], ax=ax)
        ax.set_xlim(max(chis)*1.1, 0)
        ax.tick_params(labelleft=False, labelright=True)

        if show:
            plt.show()
            plt.close()

        return chis


    # use the same scaling as in prepfold_plot.c
    def scale2d(self, array2d, indiv_scaling=True):
        """
        scale2d(array2d):
            Make a rescaled version of array2d using the same scaling as prepfold_plot.c
            indiv_scaling=True: scale each profile independently
        """
        global_max = Num.maximum.reduce(Num.maximum.reduce(array2d))
        if (global_max==0.0):  global_max = 1.0
        if indiv_scaling:
            min_parts = Num.minimum.reduce(array2d, 1)
            rescaled = (array2d-min_parts[:,Num.newaxis])/Num.fabs(global_max)
        else:
            global_min = Num.minimum.reduce(Num.minimum.reduce(array2d))
            rescaled = (array2d - global_min)/Num.fabs(global_max)
        return rescaled

    def plot2d(self, array2d, rangex=[0,1], rangey=[0,1], ax=None, cmap_name='antigrey', rangex2=None, rangey2=None, labx="", laby="", labx2="", laby2="", imshow=False, indiv_scaling=True):
        """
        plot2d(array2d, rangex=[0,1], rangey=[0,1], ax=None, cmap='gist_yarg', rangex2=None, rangey2=None, labx="", laby="", labx2="", laby2=""):
            Plot a 2d array (on ax if given)
            range<x/y> defines the outer edges of the plot
            lab<x/y> is the axis label
            range<x/y>2 defines the outer edges of the second x/y axis
            lab<x/y>2 is the label of the second x/y axis
            cmap_name is the colormap to use. It must be a named cmap from either from presto.Pgplot.py or matplotlib
            indiv_scaling determines whether profiles are scaled individually 

            returns a list of any duplicate axes which were made in this order:
            [duplicate_x_axis, duplicate_y_axis]
            if a duplicate axis was not needed it is replaced with None
        """
        array2d = self.scale2d(array2d, indiv_scaling=indiv_scaling)

        show = False
        if ax is None:
            show = True
            fig, ax = plt.subplots()

        # I believe that for pgplot these range values are the centers of the cells
        # for pcolormesh they're the edges so some transformation is required
        x_center = Num.linspace(rangex[0], rangex[-1], array2d.shape[1])
        dx = x_center[1] - x_center[0]
        y_center = Num.linspace(rangey[0], rangey[-1], array2d.shape[0])
        dy = y_center[1] - y_center[0]
        x = Num.linspace(x_center[0]-dx/2, x_center[-1]+dx/2, array2d.shape[1]+1)
        y = Num.linspace(y_center[0]-dy/2, y_center[-1]+dy/2, array2d.shape[0]+1)
        xx, yy = Num.meshgrid(x, y)

        try:
            cmap = PGPalette(cmap_name).cmap
        except AttributeError:
            cmap = cmap_name

        if imshow:
            extent = (x[0], x[-1], y[0], y[-1])
            ax.imshow(array2d, aspect="auto", origin="lower", extent=extent, cmap=cmap)
        else:
            ax.pcolormesh(xx, yy, array2d, cmap=cmap)
            # pclormesh automatically sorts your coordinates into ascending order
            if Num.argmax(x) == 0:
                ax.invert_xaxis()
            if Num.argmax(y) == 0:
                ax.invert_yaxis()
        set_labels(ax, labx=labx, laby=laby)

        # rotate y axis tick labels
        ax.tick_params(axis='y', labelrotation = 90)

        dupe_axes = []
        if rangex2 is not None:
            ax.xaxis.set_ticks_position('bottom')
            dupe_ax_x = ax.twiny()
            dupe_ax_x.xaxis.set_ticks_position('top')
            set_labels(dupe_ax_x, labx=labx2)
            dupe_ax_x.set_xlim(rangex2[0], rangex2[-1])
            dupe_axes.append(dupe_ax_x)
        else:
            dupe_axes.append(None)

        if rangey2 is not None:
            ax.yaxis.set_ticks_position('left')
            dupe_ax_y = ax.twinx()
            dupe_ax_y.yaxis.set_ticks_position('right')
            set_labels(dupe_ax_y, laby=laby2)
            dupe_ax_y.set_ylim(rangey2[0], rangey2[-1])
            dupe_ax_y.tick_params(axis='y', labelrotation = 90)
            dupe_axes.append(dupe_ax_y)
        else:
            dupe_axes.append(None)

        if show:
            plt.show()
            plt.close()
            return None, None

        # I have no idea how these'll work in scope
        # want to keep them around and not overwrite them
        # other than keeping a global list this is what ocurred to me
        return dupe_axes


    def plotxy(self, y, x=None, labx="", laby="", rangey=None, rangex=None, ax=None):
        show = False
        if ax is None:
            show = True
            fig, ax = plt.subplots()

        plot_kwargs = {
            "c": "black",
            "linewidth": 0.6
        }

        if x is not None:
            ax.plot(x, y, **plot_kwargs)
        else:
            ax.plot(y, **plot_kwargs)

        set_labels(ax, labx=labx, laby=laby)

        ax.tick_params(axis='y', labelrotation = 90)

        if rangey is not None:
            ax.set_ylim(rangey[0], rangey[-1])
        if rangex is not None:
            ax.set_xlim(rangex[0], rangex[-1])

        if show:
            plt.show()
            plt.close()

#    def greyscale(self, array2d, **kwargs):
#        """
#        greyscale(array2d, **kwargs):
#            Plot a 2D array as a greyscale image using the same scalings
#                as in prepfold.
#        """
#        # Use the same scaling as in prepfold_plot.c
#        global_max = Num.maximum.reduce(Num.maximum.reduce(array2d))
#        if (global_max==0.0):  global_max = 1.0
#        min_parts = Num.minimum.reduce(array2d, 1)
#        array2d = (array2d-min_parts[:,Num.newaxis])/Num.fabs(global_max)
#        Pgplot.plot2d(array2d, image='antigrey', **kwargs)

    def plot_intervals(self, phasebins='All', ax=None):
        """
        plot_intervals(self, phasebins='All', ax=None):
            Plot the subband-summed profiles vs time.  Restrict
                the bins in the plot to the (low:high) slice defined
                by the phasebins option if it is a tuple (low,high)
                instead of the string 'All'.
        """
        show = False
        if ax is None:
            fig, ax = plt.subplots()
            show=True

        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()
        if phasebins != 'All':
            lo, hi = phasebins
            profs = self.profs[:,:,lo:hi].sum(1)
        else:
            lo, hi = 0.0, self.proflen
            profs = self.profs.sum(1)
        dupe_axes = self.plot2d(
            profs, rangex=[lo, hi], rangey=[0.0, self.npart],
            labx="Phase Bins", labx2="Pulse Phase", laby="Time Intervals",
            rangex2=Num.asarray([lo, hi])*1.0/self.proflen,
            laby2="Time (s)", rangey2=[0.0, self.T],
            ax=ax
        )

        if show:
            plt.show()
            plt.close()
            return None, None
        return dupe_axes
    
    def plot_intervals_prepfold(self, ax=None):
        """
        KC: enough things needed to be tweaked that for now easier to make this
        plot_intervals(self, phasebins='All', ax=None):
            Plot the subband-summed profiles vs time.  Restrict
                the bins in the plot to the (low:high) slice defined
                by the phasebins option if it is a tuple (low,high)
                instead of the string 'All'.
        """
        show = False
        if ax is None:
            fig, ax = plt.subplots()
            show = True

        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()

        profs = self.profs.sum(1)
        profs_doubles = Num.zeros((self.npart, 2*self.proflen))
        profs_doubles[:,:self.proflen] = profs
        profs_doubles[:,self.proflen:] = profs
        _, _ = self.plot2d(
            profs_doubles, rangex=[0, 2], rangey=[0.0, self.T],
            labx="Pulse Phase", laby="Time (s)",
            ax=ax
        )

        # remove the tick marks from the right y-axis
        ax.tick_params(right=False, which="both")

        if show:
            plt.show()
            plt.close()

    def plot_subbands(self, phasebins='All', ax=None):
        """
        plot_subbands(self, phasebins='All', device='/xwin'):
            Plot the interval-summed profiles vs subband.  Restrict
                the bins in the plot to the (low:high) slice defined
                by the phasebins option if it is a tuple (low,high)
                instead of the string 'All'.
        """
        show = False
        if ax is None:
            fig, ax = plt.subplots()
            show=True

        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()
        if phasebins != 'All':
            lo, hi = phasebins
            profs = self.profs[:,:,lo:hi].sum(0)
        else:
            lo, hi = 0.0, self.proflen
            profs = self.profs.sum(0)
        lof = self.lofreq - 0.5*self.chan_wid
        hif = lof + self.chan_wid*self.numchan
        dupe_axes = self.plot2d(profs, rangex=[lo, hi], rangey=[0.0, self.nsub],
                      labx="Phase Bins", labx2="Pulse Phase", laby="Subbands",
                      rangex2=Num.asarray([lo, hi])*1.0/self.proflen,
                      laby2="Frequency (MHz)", rangey2=[lof, hif],
                      ax=ax)

        if show:
            plt.show()
            plt.close()
            return None, None
        return dupe_axes
    
    def plot_subbands_prepfold(self, ax=None, new_nsub=None):
        show = False
        if ax is None:
            fig, ax = plt.subplots()
            show = True

        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()

        profs = self.profs.sum(0)  # now (nsub, nbin)

        nsub = self.nsub
        if new_nsub is not None and new_nsub < self.nsub:
            if self.nsub % new_nsub:
                print("Warning!  The new number of subbands (%d) is not a" % new_nsub)
                print("          divisor of the original number of subbands (%d)!"  % self.nsub)
                print("Not fscrunching to %d." % new_nsub)
            else:
                ds = self.nsub // new_nsub
                print(f"fscrunching by a factor of {ds} ({self.nsub} -> {new_nsub})")
                # didn't want to use combine_profs as it loops through all the intervals
                profs = Num.add.reduce(profs.reshape(-1,ds,profs.shape[-1]), 1)
                nsub = new_nsub

        profs_doubled = Num.zeros((profs.shape[0], profs.shape[1]*2))
        profs_doubled[:,:profs.shape[1]] = profs
        profs_doubled[:,profs.shape[1]:] = profs

        lof = self.lofreq - 0.5*self.chan_wid
        hif = lof + self.chan_wid*self.numchan
        _, dupe_axis_y = self.plot2d(
            profs_doubled, rangex=[0, 2], rangey=[0.0, nsub],
            labx="Phase", laby="Sub-band",
            laby2="Frequency (MHz)", rangey2=[lof, hif],
            ax=ax,
        )

        # do some tick gymnastics to make it match presto
        ax.tick_params(top=False, which="both")
        ax.tick_params(axis="x", which="both", direction="out")
        ax.tick_params(axis="y", which="both", direction="out")
        dupe_axis_y.tick_params(axis="y", which="both", direction="out")

        
        # huh this should work, but if direction="in" the ticks are plotted in the background
        # and I can't seem to stop that, even using zorder
        #dupe_axis_x = ax.secondary_xaxis("top")
        #dupe_axis_x.tick_params(axis="x", which="both", direction="in", color="orange", width=10, length=10, zorder=1000)
        #dupe_axis_x.tick_params(labeltop=False)

        # this came from a suggested solution for matplotlib < 3.1, but works better
        dupe_axis_x = ax.twiny()
        dupe_axis_x.tick_params(which="both", axis="x", top=True, bottom=False, direction="in", labeltop=False)
        dupe_axis_x.sharex(ax)

        if show:
            plt.show()
            plt.close()
            return None
        return dupe_axis_x, dupe_axis_y

    def calc_varprof(self):
        """
        calc_varprof(self):
            This function calculates the summed profile variance of the
                current pfd file.  Killed profiles are ignored.
        """
        varprof = 0.0
        for part in range(self.npart):
            if part in self.killed_intervals: continue
            for sub in range(self.nsub):
                if sub in self.killed_subbands: continue
                varprof += self.stats[part][sub][5] # foldstats prof_var
        return varprof
    
    def calc_vardata(self):
        # NB stats[part][sub] corresponds to the following from foldstats
        # numdata, data_avg, data_var, numprof, prof_avg, prof_var, redchi
        """
        WRONG: does not give the same value as the prepfold plot
        calc_vardata(self):
            This function calculates the summed data variance of the
                current pfd file.  Killed profiles are ignored.
        """
        vardata = 0.0
        for part in range(self.npart):
            if part in self.killed_intervals: continue
            for sub in range(self.nsub):
                if sub in self.killed_subbands: continue
                vardata += self.stats[part][sub][2] # foldstats data_var
        return vardata
    
    def calc_avgdata(self):
        # NB stats[part][sub] corresponds to the following from foldstats
        # numdata, data_avg, data_var, numprof, prof_avg, prof_var, redchi
        """
        WRONG: does not give the same value as the prepfold plot
        (removing the /nsum doesn't do it either)
        calc_avgdata(self):
            This function calculates the summed data average of the
                current pfd file.  Killed profiles are ignored.
        """
        avgdata = 0.0
        nsum = 0
        for part in range(self.npart):
            if part in self.killed_intervals: continue
            for sub in range(self.nsub):
                if sub in self.killed_subbands: continue
                avgdata += self.stats[part][sub][1] # foldstats data_var
                nsum += 1
        return avgdata / nsum


    def calc_redchi2(self, prof=None, avg=None, var=None):
        """
        calc_redchi2(self, prof=None, avg=None, var=None):
            Return the calculated reduced-chi^2 of the current summed profile.
        """
        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()
        if prof is None:  prof = self.sumprof
        if avg is None:  avg = self.avgprof
        if var is None:  var = self.varprof
        # Note:  use the _corrected_ DOF for reduced chi^2 calculation
        return ((prof-avg)**2.0/var).sum() / self.DOFcor

    def calc_sigma(self):
        """
        calc_sigma(self):
            Return the calculated sigma (equivalent gaussian sig) of the summed profile.
        """
        return chi2_sigma(self.calc_redchi2() * self.DOFcor, self.DOFcor)

    def plot_chi2_vs_DM(self, loDM, hiDM, N=100, interp=0, ax=None):
        """
        plot_chi2_vs_DM(self, loDM, hiDM, N=100, interp=0, device='/xwin'):
            Plot (and return) an array showing the reduced-chi^2 versus
                DM (N DMs spanning loDM-hiDM).  Use sinc_interpolation
                if 'interp' is non-zero.
        # KC this is using doppler-corrected freqs
        """
        show = False
        if ax is None:
            fig, ax = plt.subplots()
            show = True

        # Sum the profiles in time
        sumprofs = self.profs.sum(0)
        if not interp:
            profs = sumprofs
        else:
            profs = Num.zeros(Num.shape(sumprofs), dtype='d')
        DMs = psr_utils.span(loDM, hiDM, N)
        chis = Num.zeros(N, dtype='f')
        subdelays_bins = self.subdelays_bins.copy()
        for ii, DM in enumerate(DMs):
            subdelays = psr_utils.delay_from_DM(DM, self.barysubfreqs)
            hifreqdelay = subdelays[-1]
            subdelays = subdelays - hifreqdelay
            delaybins = subdelays*self.binspersec - subdelays_bins
            if interp:
                interp_factor = 16
                for jj in range(self.nsub):
                    profs[jj] = psr_utils.interp_rotate(sumprofs[jj], delaybins[jj],
                                                        zoomfact=interp_factor)
                # Note: Since the interpolation process slightly changes the values of the
                # profs, we need to re-calculate the average profile value
                avgprof = (profs/self.proflen).sum()
            else:
                new_subdelays_bins = Num.floor(delaybins+0.5)
                for jj in range(self.nsub):
                    profs[jj] = psr_utils.rotate(profs[jj], int(new_subdelays_bins[jj]))
                subdelays_bins += new_subdelays_bins
                avgprof = self.avgprof
            sumprof = profs.sum(0)
            chis[ii] = self.calc_redchi2(prof=sumprof, avg=avgprof)
        # Now plot it
        self.plotxy(chis, DMs, labx="DM", laby=r"Reduced $\chi^2$", ax=ax, rangey=[0.0, max(chis)*1.1])

        if show:
            plt.show()
            plt.close()
        return (chis, DMs)

    def plot_chi2_vs_sub(self, ax=None):
        """
        plot_chi2_vs_sub(self, device='/xwin'):
            Plot (and return) an array showing the reduced-chi^2 versus
                the subband number.
        """
        show = False
        if ax is None:
            fig, ax = plt.subplots()
            show = True

        # Sum the profiles in each subband
        profs = self.profs.sum(0)
        # Compute the averages and variances for the subbands
        avgs = profs.sum(1)/self.proflen
        vars = []
        for sub in range(self.nsub):
            var = 0.0
            if sub in self.killed_subbands:
                vars.append(var)
                continue
            for part in range(self.npart):
                if part in self.killed_intervals:
                    continue
                var += self.stats[part][sub][5] # foldstats prof_var
            vars.append(var)
        chis = Num.zeros(self.nsub, dtype='f')
        for ii in range(self.nsub):
            chis[ii] = self.calc_redchi2(prof=profs[ii], avg=avgs[ii], var=vars[ii])
        # Now plot it
        self.plotxy(chis, labx="Subband Number", laby=r"Reduced $\chi^2$",
                    rangey=[0.0, max(chis)*1.1], ax=ax)

        if show:
            plt.show()
            plt.close()
        return chis

    def estimate_offsignal_redchi2(self, numtrials=20):
        """
        estimate_offsignal_redchi2():
            Estimate the reduced-chi^2 off of the signal based on randomly shifting
                and summing all of the component profiles.
        """
        redchi2s = []
        for count in range(numtrials):
            prof = Num.zeros(self.proflen, dtype='d')
            for ii in range(self.npart):
                for jj in range(self.nsub):
                    tmpprof = copy.copy(self.profs[ii][jj])
                    prof += psr_utils.rotate(tmpprof, random.randrange(0,self.proflen))
            redchi2s.append(self.calc_redchi2(prof=prof))
        return Num.mean(redchi2s)

    def adjust_fold_frequency(self, phasebins, profs=None, shiftsubs=False):
        """
        adjust_fold_frequency(phasebins, profs=None, shiftsubs=False):
            Linearly shift the intervals by phasebins over the course of
                the observation in order to change the apparent folding
                frequency.  Return a 2D array containing the de-dispersed
                profiles as a function of time (i.e. shape = (npart, proflen)),
				and the reduced chi^2 of the resulting summed profile.
                If profs is not None, then use profs instead of self.profs.
				If shiftsubs is not False, then actually correct the subbands
				instead of a 2D projection of them.
        """
        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()
        if shiftsubs:
            print("Shifting all the subbands...")
            if profs is None:
                profs = self.profs
            for ii in range(self.npart):
                bins_to_shift = int(round(float(ii)/self.npart * phasebins))
                for jj in range(self.nsub):
                    profs[ii,jj] = psr_utils.rotate(profs[ii,jj], bins_to_shift)
            redchi = self.calc_redchi2(prof=profs.sum(0).sum(0))
        else:
            print("Shifting just the projected intervals (not individual subbands)...")
            if profs is None:
                profs = self.profs.sum(1)
            for ii in range(self.npart):
                bins_to_shift = int(round(float(ii)/self.npart * phasebins))
                profs[ii] = psr_utils.rotate(profs[ii], bins_to_shift)
            redchi = self.calc_redchi2(prof=profs.sum(0))
        print("New reduced-chi^2 =", redchi)
        return profs, redchi

    def dynamic_spectra(self, onbins, combineints=1, combinechans=1,
                        calibrate=True, plot=True):
        """
        dynamic_spectra(onbins, combineints=1, combinechans=1,
                        calibrate=True, plot=True, device='/xwin'):
            Return (and plot) the dynamic spectrum (DS) resulting
                from the folds in the .pfd assuming that the pulsar
                is 'on' during the bins specified in 'onbins' and
                off elsewhere (ON-OFF).  If calibrate is True, the
                DS will be (ON-OFF)/OFF.  combineints and combinechans
                describe how many adjacent intervals or frequency
                channels will be combined when making the DS.
        """
        # Determine the indices of the off-pulse region
        indices = Num.arange(self.proflen)
        Num.put(indices, Num.asarray(onbins), -1)
        offbins = Num.compress(indices >= 0, Num.arange(self.proflen))
        numon = len(onbins)
        numoff = len(offbins)
        # De-disperse if required first
        if 'subdelays' not in self.__dict__:
            print("Dedispersing first...")
            self.dedisperse()
        # The following is the average offpulse level
        offpulse = Num.sum(Num.take(self.profs, offbins, 2), 2)/float(numoff)
        # The following is the average onpulse level
        onpulse  = Num.sum(Num.take(self.profs,  onbins, 2), 2)/float(numon)
        # Now make the DS
        self.DS = onpulse - offpulse
        self.DSnpart = self.npart
        self.DSstart_secs = self.start_secs
        self.DSintdt = self.DSstart_secs[1] - self.DSstart_secs[0]
        self.DSnsub = self.nsub
        self.DSsubfreqs = self.subfreqs
        self.DSsubdeltafreq = self.subdeltafreq
        if (calibrate):
            # Protect against division by zero
            offpulse[offpulse==0.0] = 1.0
            self.DS /= offpulse
        # Combine intervals if required
        if (combineints > 1):
            # First chop off any extra intervals
            if (self.npart % combineints):
                self.DSnpart = (self.npart // combineints) * combineints
                self.DS = self.DS[:self.DSnpart,:]
            # Now reshape and add the neighboring intervals
            self.DS = Num.reshape(self.DS, (self.DSnpart // combineints,
                                            combineints, self.DSnsub))
            print(Num.shape(self.DS))
            self.DS = Num.sum(self.DS, 1)
            self.DSstart_secs = self.DSstart_secs[::combineints]
            self.DSintdt *= combineints
            self.DSnpart //= combineints
        # Combine channels if required
        if (combinechans > 1):
            # First chop off any extra channels
            if (self.nsub % combinechans):
                self.DSnsub = (self.nsub // combinechans) * combinechans
                self.DS = self.DS[:,:self.DSnsub]
            # Now reshape and add the neighboring intervals
            self.DS = Num.reshape(self.DS, (self.DSnpart,
                                            self.DSnsub // combinechans, combinechans))
            self.DS = Num.sum(self.DS, 2)
            self.DSsubfreqs = psr_utils.running_avg(self.subfreqs[:self.DSnsub], combinechans)
            self.DSsubdeltafreq *= combinechans
            self.DSnsub //= combinechans
        print("DS shape = ", Num.shape(self.DS))
        # Plot it if required
        if plot:
            lof = self.subfreqs[0]-0.5*self.DSsubdeltafreq
            hif = self.subfreqs[-1]+0.5*self.DSsubdeltafreq
            lot = 0.0
            hit = self.DSstart_secs[-1] + self.DSintdt
            dupe_axes = self.plot2d(self.DS, rangex=[lof, hif], rangey=[lot, hit],
                           labx="Frequency (MHz)", labx2="Subband Number",
                           laby="Time (s)", laby2="Interval Number",
                           rangex2=[0, self.DSnsub], rangey2=[0, self.DSnpart],
                           ax=None)
            # didn't want to return dupe_axes as this function already returns something
            plt.show()
            plt.close()
        return self.DS

    def compute_periods_pdots(self, pstep, pdstep, npfact):
        """Compute periods and pdots to search over for a given pstep, pdstep, npfact"""
        # weirdly this is not the same as self.dt*self.Nfolded
        # this is the T prepfold uses to calculate the periods etc
        tot_time = self.start_secs[-1] + self.start_secs[1] - self.start_secs[0]

        numtrials = 2 * self.proflen * npfact + 1
        foldf, foldfd = self.fold_p1, self.fold_p2
        periods = Num.zeros((numtrials))
        pdots = Num.zeros((numtrials))
        fdots = Num.zeros((numtrials))
        fdotdots = Num.zeros((numtrials))

        for ii in range(numtrials):
            totpdelay = ii - (numtrials - 1) / 2
            dtmp = (totpdelay * pstep) / self.proflen
            periods[ii] = 1.0 / (foldf + dtmp / tot_time)
            dtmp = (totpdelay * pdstep) / self.proflen
            fdots[ii] = phasedelay2fdot(dtmp, tot_time)
            pdots[ii] = psr_utils.p_to_f(foldf, foldfd + fdots[ii])[1]
            fdotdots[ii] = phasedelay2fdotdot(dtmp, tot_time)

        return periods, pdots


    def ppdot_grid(self, dm=None, periods=None, pdots=None, search=False, pdd=None, doppler=False):
        """
        ppdot_grid(self, dm=None, periods=None, pdots=None, search=False, pdd=None):
            Make the p-pdot grid.
            search=True: Find the best p and pdot in the grid and update the pfd with 
                these values.
            search=False: Make the grid but afterwards return the pfd to its original
                state.
            If periods, pdots is None, use values from prepfold search.
            If dm, pdd is None, use best values already stored in pfd.

            NB no search is done over dm or pdd.
            Note this function uses adjust_period_more_efficient.
            The original adjust_period notes in its docstring that the
            results are *almost* identical to those from the prepfold search.
            adjust_period_more_efficient has 2 differences from adjust_period
              - it uses Num.roll rather than psr_utils.rotate
              - there is not option to use fft interpolation

            grid has shape (npdots, nperiods)

            returns grid, periods, pdots, dm, pdd, best_p, best_pd
        """
        if periods is None:
            periods = self.periods
        if pdots is None:
            pdots = self.pdots
        if pdd is None:
            if self.fold_pow == 1.0:
                pdd = self.bary_p3
            else:
                pdd = self.topo_p3
        if dm is None:
            dm = self.bestdm

        # check if already have an identical one stored
        if 'ppdotgrid' in self.__dict__:
            print("Found existing ppdot grid, checking for compatibility")
            if self.ppdotgrid.same_input_pars(periods, pdots, pdd, dm, doppler):
                print("Same input parameters. Using existing ppdotgrid")
                return self.ppdotgrid.redchi2s, self.ppdotgrid.periods, self.ppdotgrid.pdots, self.ppdotgrid.dm, self.ppdotgrid.pdd, self.ppdotgrid.best_p, self.ppdotgrid.best_pd
            else:
                print("Not a match. Proceeding with the search; ppdotgrid will be overwritten")

        # so can return the pfd to the state it started in if search=False
        starting_values = [self.currdm, self.curr_p1, self.curr_p2, self.curr_p3]
        print(f"For p-pdot search, using dm: {dm}, pdd: {pdd}")
        self.dedisperse(dm, doppler=doppler)
        redchi2s = Num.zeros((len(pdots), len(periods)))
        pp, pdpd = Num.meshgrid(periods, pdots)
        for idx, p in Num.ndenumerate(pp):
            self.adjust_period_more_efficient(p=p, pd=pdpd[idx], pdd=pdd)
            redchi2s[idx] = self.calc_redchi2()

        best_index = Num.unravel_index(Num.argmax(redchi2s), redchi2s.shape)
        best_p = periods[best_index[1]]
        best_pd = pdots[best_index[0]]

        if not search:
            self.dedisperse(starting_values[0])
            p, pd, pdd = starting_values[1:]
            #print(f"Resetting dm, p, pd, pdd to their starting values: {starting_values}")
            self.dedisperse(starting_values[0])
            # I think this is necessary so other functions behave as expected
            if starting_values[0] == 0:
                del self.subdelays
            self.adjust_period_more_efficient(p=p, pd=pd, pdd=pdd)
        else:
            self.adjust_period_more_efficient(p=best_p, pd=best_pd, pdd=pdd)

        self.ppdotgrid = RedChi2s(redchi2s, periods, pdots, pdd, dm, doppler, best_p=best_p, best_pd=best_pd)

        return redchi2s, periods, pdots, dm, pdd, best_p, best_pd

    def ppdot_grid2(self, dm=None, pdd=None, periods=None, pdots=None, doppler=False):
        # this should leave the profiles etc unchanged
        # BUT it also assumes that the curr_p1 etc values are those derived from fold_p1 etc
        if periods is None:
            periods = self.periods
        if pdots is None:
            pdots = self.pdots
        if pdd is None:
            if self.fold_pow == 1.0:
                pdd = self.bary_p3
            else:
                pdd = self.topo_p3
        if dm is None:
            dm = self.bestdm

        # check if already have an identical one stored
        if 'ppdotgrid' in self.__dict__:
            print("Found existing ppdot grid, checking for compatibility")
            if self.ppdotgrid.same_input_pars(periods, pdots, pdd, dm, doppler):
                print("Same input parameters. Using existing ppdotgrid")
                return self.ppdotgrid.redchi2s, self.ppdotgrid.periods, self.ppdotgrid.pdots, self.ppdotgrid.dm, self.ppdotgrid.pdd, self.ppdotgrid.best_p, self.ppdotgrid.best_pd
            else:
                print("Not a match. Proceeding with the search; ppdotgrid will be overwritten")

        print(f"For p-pdot search, using dm: {dm}, pdd: {pdd}")
        startingdm = self.currdm
        self.dedisperse(dm, doppler=doppler)

        # pre-calculate some stuff
        parttimes = self.start_secs.astype('float32').astype('float64')
        fold_ps = self.fold_p1, self.fold_p2, self.fold_p3
        # not searching in dm so can sum all subbands
        tmp_profs = self.profs.sum(1)
        pp, pdpd = Num.meshgrid(periods, pdots)

        redchi2s = calc_many_redchis(pp, pdpd, pdd, parttimes, fold_ps, tmp_profs, self.proflen, self.npart, self.pdelays_bins, self.avgprof, self.varprof, self.DOFcor)
        best_index = Num.unravel_index(Num.argmax(redchi2s), redchi2s.shape)
        best_p = periods[best_index[1]]
        best_pd = pdots[best_index[0]]

        self.dedisperse(startingdm, doppler=doppler)
        if startingdm == 0:
                del self.subdelays

        self.ppdotgrid = RedChi2s(redchi2s, periods, pdots, pdd, dm, doppler, best_p=best_p, best_pd=best_pd)

        return redchi2s, periods, pdots, dm, pdd, best_p, best_pd

    def plot_ppdot_grid(self, dm=None, pdd=None, periods=None, pdots=None, axs=None, search=True, use_fold=False, imshow=False, doppler=False):
        """
            if search=True the best p and pdot in the grid will be used for the p/pdot vs redchi2 plots
                BUT pfd will still be in its original state (aka ppdot_grid is run with search=False)
            if search=False 
                and use_fold=False the best p and pdot stored in the pfd (aka from the prepfold search) will be used
                and use_fold=True the original p and pdot from the prepfold command will be used
                (if you leave periods, pdots as None these should both exist in the grid
                if you pass in your own, and there's no match, the closest points in the grid will be used)

            axs can be one axis, in which case just the grid will be plotted
            or a list of 3 axes for the grid, period, pdot in that order.
            If axs is None, all 3 will be plotted by default.

            returns:
            dupe_axes_grid, dm, best_p, best_pd, pdd
            
            where dupe_axes_grid is a list of duplicate axes used for plotting
            [dupe_asis_x, dupe_axis_y]
            (if passed in axs=None, these will both be None)
        """
        show = False
        if axs is None:
            show = True
            fig, [[ax_p, ax_delete],[ax_grid, ax_pdot]] = plt.subplots(
                2,2, 
                sharex="col", sharey="row",
                gridspec_kw={"width_ratios": [2, 1], "height_ratios": [1, 2]},
            )
            ax_delete.set_axis_off()
        elif isinstance(axs, list):
            ax_grid, ax_p, ax_pdot = axs
        else:
            ax_grid = axs
            ax_p, ax_pdot = None, None

        pfold, pdfold, pddfold = psr_utils.p_to_f(self.fold_p1, self.fold_p2, self.fold_p3)
        ffold, fdfold = self.fold_p1, self.fold_p2
        redchi2s, periods, pdots, dm, pdd, best_p, best_pd = self.ppdot_grid2(dm=dm, periods=periods, pdots=pdots, doppler=doppler)

        # get p,pd to use for 1D plots, and to mark on chisq grid
        mark_bestppd = True
        if search:
            best_index = Num.unravel_index(Num.argmax(redchi2s), redchi2s.shape)
            # check best values are what expect
            assert best_p == periods[best_index[1]]
            assert best_pd == pdots[best_index[0]]
            print(f"Found best period:\n{best_p}")
            print(f"Found best pdot:\n{best_pd}")
        else:
            if use_fold:
                best_p = pfold
                best_pd = pdfold
                # it'll be at [0,0] so don't need to mark it
                mark_bestppd = False
            else:
                if self.fold_pow == 1.0:
                    best_p = self.bary_p1
                    best_pd = self.bary_p2
                else:
                    best_p = self.topo_p1
                    best_pd = self.topo_p2

            best_index_p = Num.argmin(abs(periods - best_p))
            best_index_pd = Num.argmin(abs(pdots - best_pd))
            best_index = [best_index_pd, best_index_p]

        p_slice = redchi2s[best_index[0],:]
        pd_slice = redchi2s[:, best_index[1]]

        prange = (Num.array([periods[0], periods[-1]]) - pfold) * 1000
        frange = (1/Num.array([periods[0], periods[-1]]) - ffold)
        pdrange = (Num.array([pdots[0], pdots[-1]]) - pdfold)
        fdrange = [
            psr_utils.p_to_f(pfold, pdots[0])[1] - fdfold,
            psr_utils.p_to_f(pfold, pdots[-1])[1] - fdfold,
        ]
        labx = f"Period - {pfold*1000:.8f} (ms)"
        labx2 = f"Freq - {ffold:.6f} (Hz)"
        if pdfold == 0:
            laby = "P-dot (s/s)"
        elif pdfold < 0:
            laby = f"P-dot + {-pdfold:.5g} (s/s)"
        else:
            laby = f"P-dot - {pdfold:.5g} (s/s)"
        if fdfold == 0:
            laby2 = "F-dot (Hz/s)"
        elif fdfold < 0:
            laby2 = f"F-dot + {-fdfold:.5g} (Hz/s)"
        else:
            laby2 = f"F-dot - {fdfold:.5g} (Hz/s)"

        dupe_axes_grid = self.plot2d(redchi2s, rangex=prange, rangey=pdrange, ax=ax_grid, cmap_name='antirainbow', rangex2=frange, rangey2=fdrange, labx=labx, laby=laby, labx2=labx2, laby2=laby2, imshow=imshow, indiv_scaling=False)
        if mark_bestppd:
            ax_grid.scatter((best_p-pfold)*1000, best_pd-pdfold, marker="+", linewidth=2)

        if ax_p is not None:
            if mark_bestppd:
                ax_p.axvline((best_p-pfold)*1000, c='lightsalmon')
            # test if it's in a corner plot
            if ax_grid in ax_p.get_shared_x_axes().get_siblings(ax_p):
                _labx = None
            else:
                _labx = labx
            rangey=[0.0, max(p_slice)*1.1]
            self.plotxy(p_slice, x=(periods-pfold)*1000, labx=_labx, laby=r"Reduced $\chi^2$", rangex=prange, rangey=rangey, ax=ax_p)

        if ax_pdot is not None:
            if ax_grid in ax_pdot.get_shared_y_axes().get_siblings(ax_pdot):
                if mark_bestppd:
                    ax_pdot.axhline((best_pd-pdfold), c='lightsalmon')
                x = pd_slice
                y = pdots - pdfold
                _laby = ""
                _labx = r"Reduced $\chi^2$"
                rangex = [0.0, max(pd_slice)]
                rangey = pdrange
            else:
                if mark_bestppd:
                    ax_pdot.axvline((best_pd-pdfold), c='lightsalmon')
                x = pdots - pdfold
                y = pd_slice
                _labx = laby
                _laby = r"Reduced $\chi^2$"
                rangex = pdrange
                rangey = [0.0, max(pd_slice)]
            self.plotxy(y, x=x, labx=_labx, laby=_laby, rangex=rangex, rangey=rangey, ax=ax_pdot)

        if show:
            plt.show()
            plt.close()
            return None, None, dm, best_p, best_pd, pdd

        return dupe_axes_grid, dm, best_p, best_pd, pdd


    def plot_intervals_corner(self, axs=None):
        """
            axs, if passed in, must be 2x2 which can be indexed like, e.g., axs[1,1]
        """
        show = False
        if axs is None:
            show = True
            fig, axs = plt.subplots(
                2,2, sharex="col", 
                gridspec_kw=dict(hspace=0, wspace=0, height_ratios=[3,8], width_ratios=[2,1])
            )
        else:
            # share axes if not already set up
            if axs[1,0] not in axs[0,0].get_shared_x_axes().get_siblings(axs[0,0]):
                axs[0,0].sharex(axs[1,0])
#            if axs[1,0] not in axs[1,1].get_shared_y_axes().get_siblings(axs[1,1]):
#                axs[1,1].sharey(axs[1,0])

        # delete some axis spines etc
        axs[0,1].set_axis_off()
        axs[0,0].spines['left'].set_visible(False)
        axs[0,0].spines['top'].set_visible(False)
        axs[0,0].spines['right'].set_visible(False)
        axs[0,0].set_yticks([])

        self.plot_intervals_prepfold(ax=axs[1,0])  # no duplicate axes
        self.plot_chi2_vs_int_prepfold(ax=axs[1,1])
        self.plot_sumprof_prepfold(ax=axs[0,0])
        
        # set some tickmark things
        axs[1,1].yaxis.set_label_position("right")
        axs[0,0].xaxis.set_ticks_position('bottom')
        axs[0,0].set_xticks([1])

        axs[0,0].set_title("2 Pulses of Best Profile")

        if show:
            plt.show()
            plt.close()


    def prepfold_plot(self, dm=None, search=False, use_fold=False, new_nsub=None, lodm=None, hidm=None, ndm=None, reset_afterwards=True, doppler=True):
        starting_values = [self.currdm, self.curr_p1, self.curr_p2, self.curr_p3]
        #print("DEBUG: barysubfreqs==subfreqs?",(self.barysubfreqs==self.subfreqs).all())
        fig = plt.figure(figsize=(10.11,7.56), constrained_layout=False)
        #fig = plt.figure(figsize=(12,8))
        outer_grid = fig.add_gridspec(11, 3, width_ratios=[3,2,2], wspace=0.4)

        # set up axes layout
        grid_text = outer_grid[:3,1:].subgridspec(1,1)
        ax_text = grid_text.subplots()
        grid_int = outer_grid[:,0].subgridspec(2,2, wspace=0, hspace=0, height_ratios=[3,8], width_ratios=[2,1])
        axs_int = grid_int.subplots(sharex="col")
        axs_int[0,1].spines['top'].set_visible(False)
        axs_int[0,1].spines['right'].set_visible(False)
        grid_subband_dmchi = outer_grid[3:,1].subgridspec(2,1,height_ratios=[5,3], hspace=0.3)
        axs_subband_dmchi = grid_subband_dmchi.subplots()
        # OK because the spacing is difficult think I want these to be two separate gridspecs, le sigh
        grid_ppdotchi = outer_grid[3:,2].subgridspec(2,1, hspace=0.5, height_ratios=[1,1.05])  # 2D plot and 2 1D plots
        grid_ppdotchi1d = grid_ppdotchi[0].subgridspec(2,1, hspace=0.8)  # 2 1D plots
        axs_ppdotchi1d = grid_ppdotchi1d.subplots()
        grid_ppdotchi2d = grid_ppdotchi[1].subgridspec(1,1)
        axs_ppdotchi2d = grid_ppdotchi2d.subplots()

        # Need to do the ppdot grid first since it finds the "best" p and pd values
        # and we want those applied for the other plots
        # if search == False:
        #    if use_fold == True:
        #        best_p, best_pd are the fold values (arguments passed into prepfold)
        #    else:
        #        best_p, best_pd are the best values from the original prepfold search
        print("Plotting P-Pdot")
        dupe_axs_ppdot, dm, best_p, best_pd, pdd = self.plot_ppdot_grid(dm=dm, search=search, use_fold=use_fold, axs=[axs_ppdotchi2d, axs_ppdotchi1d[1], axs_ppdotchi1d[0]], doppler=doppler)
        # adjust the pfd to the values to be used for the other plots
        self.adjust_period_more_efficient(p=best_p, pd=best_pd, pdd=pdd)

        # DM-chi2 plot is next because I think it expects the profiles to be at DM0
        assert self.currdm == 0
        axs_subband_dmchi[1].axvline(dm, c='lightsalmon')
        print("Plotting DM")
        if lodm is None:
            lodm = self.dms.min()
        if hidm is None:
            hidm = self.dms.max()
        if ndm is None:
            ndm = len(self.dms)
        _ , _ = self.plot_chi2_vs_DM(lodm, hidm, ndm, ax=axs_subband_dmchi[1])

        # Now we dedisperse before plotting the intervals and subbands
        self.dedisperse(DM=dm, doppler=doppler)
        print("Plotting intervals")
        dupe_axs_intervals = self.plot_intervals_corner(axs=axs_int)
        print("Plotting subbands")
        dupe_axs_subband = self.plot_subbands_prepfold(ax=axs_subband_dmchi[0], new_nsub=new_nsub)

        # will need to set some tick marks and such
        # move P, P-dot vs chi plot labels to the right
        axs_ppdotchi1d[0].tick_params(labelleft=False, labelright=True)
        axs_ppdotchi1d[1].tick_params(labelleft=False, labelright=True)
        axs_ppdotchi1d[0].yaxis.set_label_position("right")
        axs_ppdotchi1d[1].yaxis.set_label_position("right")

        # Sometimes too many yticks and the labels run into each other
        # tried MaxNLocator but it gave weird results with just 1 tick
        for axx in [axs_ppdotchi2d, dupe_axs_ppdot[1]]:
            if (len(axx.yaxis.get_ticklabels()) - 2) > 4:  # added a -2 since it often trims the ends
                resample_ticks(axx, 'y', factor=2, keep=0)

        # The text box bit
        ax_text.set_axis_off()

        tepoch_str = r"$\rm{Epoch_{topo}}$ = "
        if self.tepoch <= 0:
            tepoch_str += "N/A"
        else:
            tepoch_str += f"{self.tepoch:.11f}"

        tbary_str = r"$\rm{Epoch_{bary}}$ = "
        if self.bepoch == 0:
            tbary_str += "N/A"
        else:
            tbary_str += f"{self.bepoch:.11f}"

        cand_text = [
            f"Candidate: {self.candnm}",
            f"Telescope: {self.telescope}",
            tepoch_str,
            tbary_str,
            r"$\rm{T_{sample}}$",  # can't find an attribute which indicates whether an Events were used
            f"Data Folded",  # can't find an attribute which indicates whether an Events were used
            f"Data Avg", # beststats.data_avg
            f"Data StdDev", # sqrt(beststats.data_var)
            f"Profile Bins",
            f"Profile Avg", # beststats.prof_avg
            f"Profile StdDev", # sqrt(beststats.prof_var)
        ]
        cand_text2 = [
            f" = {self.dt:.5g}",
            f" = {int(self.Nfolded)}",
            f" = {self.calc_avgdata():.4g}",
            f" = {Num.sqrt(self.calc_vardata()):.4g}",
            f" = {self.proflen:d}",
            f" = {self.avgprof:.4g}",
            f" = {Num.sqrt(self.varprof):.4g}",
        ]
        offset = 0.2
        cand_text_x = -0.36
        linespacing=1.1
        ax_text.text(cand_text_x, 0.01, "\n".join(cand_text), transform=ax_text.transAxes, fontsize='smaller', linespacing=linespacing)
        ax_text.text(cand_text_x+offset, 0.01, "\n".join(cand_text2), transform=ax_text.transAxes, fontsize='smaller', linespacing=linespacing)

        # set orbital parameter strings
        Porb_str = r"$\rm{P_{orb}}$" + " (s) = "
        asini_str = r"$\rm{a_1 sin(i)/c}$ = "
        Tperi_str = r"$\rm{T_{peri}}$ = "
        ecc_str = r"$\rm{e}$ = "
        omega_str = r"$\omega$" + " (rad) = "
        if self.orb_p == 0:
            Porb_str += "N/A"
            asini_str += "N/A"
            Tperi_str += "N/A"
            ecc_str += "N/A"
            omega_str += "N/A"
        else:
            Porb_str += f"{self.orb_p}"
            asini_str += f"{self.orb_x}"
            Tperi_str += f"{self.orb_t}"
            ecc_str += f"{self.orb_e}"
            omega_str += f"{self.orb_w}"

        search_text_lhs = [
            r"$\rm{RA_{J2000}}$" + f" = {self.rastr.decode('utf-8')}",
            "           Folding Parameters",
            r"$\rm{DOF_{eff}}$" + f" = {self.DOFcor:.2f}" + r"  $\chi^2_{\rm{red}}$" + f" = {self.calc_redchi2():.3f}",
            r"Dispersion Meansure (DM; $\rm{pc/cm^3}$)" + f" = {self.currdm}",
            r"$\rm{P_{x}}$" + f" (ms) = {self.curr_p1*1000}",  # I have no idea if these are bary or topo! It maybe depends. I think topo in my file but not sure if that's always true
            r"$\rm{P'_{x}}$" + f" (s/s) = {self.curr_p2}",
            r"$\rm{P''_{x} (s/s^2)}$" + f" = {self.curr_p3}",
            "           Binary Parameters",
            Porb_str,
            asini_str,
            Tperi_str,
        ]

        # this won't work if did the actual presto thing
        logp = scipychi2.logsf(self.calc_redchi2() * self.DOFcor, self.DOFcor)

        search_text_rhs = [
            r"$\rm{DEC_{J2000}}$" + f" = {self.decstr.decode('utf-8')}",
            "",
            f"P(Noise) < {Num.exp(logp):.2g}  ({self.calc_sigma():.1f})" + r"$\sigma$)",
            "",
            r"$\rm{ }$",
            r"$\rm{P_{bary}}$",
            r"$\rm{P'_{bary}}$",
            r"$\rm{P''_{bary}}$",
            "",
            ecc_str,
            omega_str,
            r"$\rm{ }$"
        ]

        offset = 0.6
        fold_text_x = 0.1
        linespacing=0.9
        ax_text.text(fold_text_x, 0.01, "\n".join(search_text_lhs), transform=ax_text.transAxes, fontsize='smaller', linespacing=linespacing)
        ax_text.text(fold_text_x+offset, 0.01, "\n".join(search_text_rhs), transform=ax_text.transAxes, fontsize='smaller', linespacing=linespacing)

        if reset_afterwards:
            #print(f"Resetting dm, p, pd, pdd to their starting values: {starting_values}")
            self.dedisperse(starting_values[0], doppler=doppler)
            if starting_values[0] == 0:
                del self.subdelays
            self.adjust_period_more_efficient(p=starting_values[1], pd=starting_values[2], pdd=starting_values[3])

        return fig

# added this so I could store it in the pfd and didn't have to recompute it every time
class RedChi2s:
    def __init__(self, redchi2s, periods, pdots, pdd, dm, doppler,  best_p=None, best_pd=None):
        self.redchi2s = redchi2s
        self.periods = periods
        self.pdots = pdots
        self.pdd = pdd
        self.dm = dm
        self.doppler = doppler
        # these are derived quantities but this is simpler
        self.best_p = best_p
        self.best_pd = best_pd
        best_index = Num.unravel_index(Num.argmax(redchi2s), redchi2s.shape)
        if self.best_p is None:
            self.best_p = periods[best_index[1]]
        if self.best_pd is None:
            self.best_pd = pdots[best_index[0]]

    def __eq__(self, other): 
        if not isinstance(other, RedChi2s):
            # don't attempt to compare against unrelated types
            return NotImplemented

        condition = (
            (self.redchi2s == other.redchi2s).all() 
            and (self.periods == other.periods).all()
            and (self.pdots == other.pdots).all()
            and self.pdd == other.pdd
            and self.dm == other.dm
            and self.doppler == other.doppler
        )

        return condition

    def same_input_pars(self, periods, pdots, pdd, dm, doppler):
        condition = (
            (self.periods == periods).all()
            and (self.pdots == pdots).all()
            and self.pdd == pdd
            and self.dm == dm
            and self.doppler == doppler
        )

        return condition

def phasedelay2fdot(phasedelay, time):
    """Translated from prepfold.c"""
    if (time == 0.0):
        return 0.0
    else:
        return 2.0 * phasedelay / (time * time)

def phasedelay2fdotdot(phasedelay, time):
    """Translated from prepfold.c"""
    if (time == 0.0):
        return 0.0
    else:
        return 6.0 * phasedelay / (time * time * time)

@jit(nopython=True)
def p_to_f(p, pd, pdd=None):
    """
    p_to_f(p, pd, pdd=None):
       Convert period, period derivative and period second
       derivative to the equivalent frequency counterparts.
       Will also convert from f to p.
    """
    f = 1.0 / p
    fd = -pd / (p * p)
    if pdd is None:
        return [f, fd]
    else:
        if pdd == 0.0:
            fdd = 0.0
        else:
            fdd = 2.0 * pd * pd / (p**3.0) - pdd / (p * p)
        return [f, fd, fdd]

@jit(nopython=True)
def delay_from_foffsets(df, dfd, dfdd, times):
    """
    Return the delays in phase caused by offsets in
    frequency (df), and two frequency derivatives (dfd, dfdd)
    at the given times in seconds.
    """
    f_delays = df * times
    fd_delays = dfd * times**2 / 2.0
    fdd_delays = dfdd * times**3 / 6.0
    return f_delays + fd_delays + fdd_delays

@jit(nopython=True)
def calc_many_redchis(pp, pdpd, pdd, parttimes, fold_ps, tmp_profs, proflen, npart, pdelays_bins, avgprof, varprof, DOFcor):
    redchi2s = Num.zeros(pp.shape)
    foldf, foldfd, foldfdd = fold_ps
    foldp, foldpd, foldpdd = p_to_f(*fold_ps)

    for idx, p in Num.ndenumerate(pp):
        pd = pdpd[idx]
        new_profs = Num.zeros_like(tmp_profs)
        fdd = p_to_f(foldp, foldpd, pdd)[2]
        fd = p_to_f(foldp, pd)[1]
        f = 1.0/p
        f_diff = f - foldf
        fd_diff = fd - foldfd
        fdd_diff = fdd - foldfdd
        delays = delay_from_foffsets(f_diff, fd_diff, fdd_diff, parttimes)
        # Convert from delays in phase to delays in bins
        bin_delays = Num.fmod(delays * proflen, proflen) - pdelays_bins
        new_pdelays_bins = Num.floor(bin_delays+0.5)
        new_pdelays_bins = Num.array([int(x) % proflen for x in new_pdelays_bins])

        for jj in range(npart):
            new_profs[jj,:] = Num.roll(tmp_profs[jj,:], new_pdelays_bins[jj])
        sumprof = new_profs.sum(0)
        redchi2s[idx] = ((sumprof - avgprof)**2.0/varprof).sum() / DOFcor

    return redchi2s

def resample_ticks(ax, xory, factor=2, keep=0):
    """
    Resample ticks by a factor of <factor>, making sure to keep the value <keep> if it appears in the ticks
    xory = 'x' or 'y'
    """
    if xory == 'x':
        axxy = ax.xaxis
        set_fn = ax.set_xticks
    elif xory == 'y':
        axxy = ax.yaxis
        set_fn = ax.set_yticks
    else:
        raise RuntimeError(f"xory ({xory}) must be either 'x' or 'y'")

    current_ticks = axxy.get_majorticklocs()
    offset=0
    if keep in current_ticks:
        preserve = Num.where(current_ticks==keep)[0][0]
        offset = preserve % factor

    current_ticks = list(current_ticks)
    custom_ticks = [t for t in current_ticks if not (current_ticks.index(t)-offset)%factor] 
    # set as new labels
    set_fn(custom_ticks)


if __name__ == "__main__":
    #testpfd = "/home/ransom/tmp_pfd/M5_52725_W234_PSR_1518+0204A.pfd"
    #testpfd = "/home/ransom/tmp_pfd/M13_52724_W234_PSR_1641+3627C.pfd"
    testpfd = "M13_53135_W34_rficlean_DM30.10_PSR_1641+3627C.pfd"

    tp = pfd(testpfd)

    if (0):
        print(tp.start_secs)
        print(tp.mid_secs)
        print(tp.start_topo_MJDs)
        print(tp.mid_topo_MJDs)
        print(tp.T)

    #tp.kill_subbands([6,7,8,9,30,31,32,33])
    #tp.kill_intervals([2,3,4,5,6])

    #tp.plot_chi2_vs_sub()
    #(chis, DMs) = tp.plot_chi2_vs_DM(0.0, 50.0, 501, interp=1)
    #best_index = Num.argmax(chis)
    #print "Best DM = ", DMs[best_index]

    (chis, DMs) = tp.plot_chi2_vs_DM(0.0, 50.0, 501)
    best_index = Num.argmax(chis)
    print("Best DM = ", DMs[best_index])

    tp.dedisperse()
    tp.plot_subbands()
    tp.plot_sumprof()
    print("DM =", tp.bestdm, "gives reduced chi^2 =", tp.calc_redchi2())

    tp.dedisperse(27.0)
    tp.plot_subbands()
    tp.plot_sumprof()
    print("DM = 27.0 gives reduced chi^2 =", tp.calc_redchi2())

    tp.dedisperse(33.0)
    tp.plot_subbands()
    tp.plot_sumprof()
    print("DM = 33.0 gives reduced chi^2 =", tp.calc_redchi2())

    tp.plot_intervals()
