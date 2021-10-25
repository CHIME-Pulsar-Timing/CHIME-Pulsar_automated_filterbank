from __future__ import print_function
from __future__ import absolute_import

from matplotlib import pyplot as plt
from matplotlib import colorbar
import attr
from attr import attrs, attrib
from typing import List
import os

from builtins import zip
import numpy as np

import argparse


######################################################
### Normal DDplan.py things
# Major changes
# changed classes to use attrs
# added the option to include scattering
# added supporting functions to allow ddplan's to be written to file
# took out all dependence on presto
# changes pgplot to matplotlib

@attrs
class observation(object):
    """
    An observation

    Attributes:
    dt (float): sampling time/resolution in sec
    f_ctr (float): central frequency in MHz
    BW (float): bandwidth in MHz
    numchan (int): number of frquency channels
    cDM (float): coherent dedispersion, default is 0

    chanwidth (float): Calculated on initialization BW/numchan
    """
    dt: float = attrib()
    f_ctr: float = attrib()
    BW: float = attrib()
    numchan: int = attrib()
    cDM: float = attrib(default=0)
    chanwidth: float = attrib(init=False)

    @chanwidth.default
    def _chanwidth(self):
        return self.BW / self.numchan

    def guess_dDM(self, DM):
        """
        guess_dDM(self, DM):
            Choose a reasonable dDM by setting the maximum smearing across the
                'BW' to equal the sampling time 'dt'.
        """
        return self.dt*0.0001205*self.f_ctr**3.0/(0.5*self.BW)

    def as_dict(self):
        """
        Return a dictionary representation for easy storage
        """
        ret_dict = dict(
            dt=self.dt,
            f_ctr=self.f_ctr,
            BW=self.BW,
            numchan=self.numchan,
            cDM=self.cDM
        )
        return ret_dict


@attrs
class dedisp_method(object):
    obs: observation = attrib()
    downsamp: int = attrib()
    numsub: int = attrib()
    numprocs: int = attrib()
    smearfact: float = attrib()
    scatteringfact: float = attrib()
    include_scattering: bool = attrib()
    DMs: np.ndarray = attrib(eq=attr.cmp_using(eq=np.array_equal))
    DMs_per_prepsub: int = attrib()
    dsubDM: int = attrib()
    numprepsub: int = attrib()
    # did make these properties but lead to annoying long floats
    dDM: float = attrib()
    loDM: float = attrib()
    hiDM: float = attrib()

    @property
    def numDMs(self):
        return len(self.DMs)


    @property
    def BW_smearing(self):
        return BW_smear(self.dDM, self.obs.BW, self.obs.f_ctr)

    @property
    def sub_smearing(self):
        return subband_smear(self.dsubDM, self.numsub, self.obs.BW, self.obs.f_ctr)

    def chan_smear(self, DM):
        """
        Return the smearing (in ms) in each channel at the specified DM
        """
        try:
            DM = np.where(DM-self.obs.cDM==0.0, self.obs.cDM+self.dDM/2.0, DM)
        except TypeError:
            if (DM-self.obs.cDM==0.0): DM = self.obs.cDM+self.dDM/2.0
        return dm_smear(DM, self.obs.chanwidth, self.obs.f_ctr, self.obs.cDM)

    def scattering_smear(self, DM):
        """
        Return the smearing (in ms) due to Bhat(2004) scattering at the top
        of the band.
        """
        return scattering_smearing(self.obs.f_ctr+self.obs.BW/2, DM)

    def total_smear_without_scattering(self, DM):
        """
        Return the totals smearing in ms without including scattering.
        Includes smearing due to:
            the sampling rate,
            the smearing over each channel,
            the smearing over each subband
            (if numsub > 0), the smearing over the full BW assuming the
        worst-case DM error
        This is useful because assume less scattering to determine a
        conservative DDplan, but for a conservative sensitivity estimate
        you want to assume more scattering.
        """
        return np.sqrt((1000.0*self.obs.dt)**2.0 +
                       (1000.0*self.obs.dt*self.downsamp)**2.0 +
                       self.BW_smearing**2.0 +
                       self.sub_smearing**2.0 +
                       self.chan_smear(DM)**2.0)

    def total_smear(self, DM):
        """
        Return the total smearing in ms due to the sampling rate,
        the smearing over each channel, the smearing over each subband
        (if numsub > 0), the smearing over the full BW assuming the
        worst-case DM error, and, if included, the smearing due to
        scattering at the top of the band.
        """
        if self.include_scattering:
            scatter = self.scattering_smear(DM)
        else:
            scatter = 0
        return np.sqrt((1000.0*self.obs.dt)**2.0 +
                       (1000.0*self.obs.dt*self.downsamp)**2.0 +
                       scatter**2.0 +
                       self.BW_smearing**2.0 +
                       self.sub_smearing**2.0 +
                       self.chan_smear(DM)**2.0)

    def DM_for_smearfact(self):
        """
        Return the DM where the smearing in a single channel is a factor smearfact
        larger than all the other smearing causes combined.
        (except scattering)
        """
        return DM_for_smearfact(self.obs, self.dDM, self.dsubDM, self.numsub, self.downsamp, self.smearfact)

    def DM_for_scatteringfact(self):
        """
        Return the DM where the scattering at the top of the band is a factor scatteringfact
        larger than all the other smearing causes combined.
        (except in-channel DM smearing)
        """
        return DM_for_scatteringfact(self.obs, self.dDM, self.dsubDM, self.numsub, self.downsamp, self.scatteringfact)

    def plot(self, work_fract, axs):
        DMspan = self.DMs[-1]-self.DMs[0]
        loDM  = self.DMs[0]  + DMspan*0.02
        hiDM  = self.DMs[-1] - DMspan*0.02
        midDM = self.DMs[0]  + DMspan*0.5
        dt_ms = 1000.0*self.obs.dt*self.downsamp
        axs.plot(self.DMs, self.total_smear(self.DMs), linewidth=1, color='black')
        axs.text(midDM, 1.1*self.total_smear(midDM), f"{self.numDMs} ({100*work_fract:.1f}%)", horizontalalignment='left', verticalalignment='bottom', rotation=90)

        # Sample time
        axs.plot(self.DMs, np.zeros(self.numDMs)+dt_ms, linewidth=1, color='green')
        axs.text(loDM, 0.85*dt_ms, f"{dt_ms}", color='green')

        # DM stepsize smearing
        axs.plot(self.DMs, np.zeros(self.numDMs)+self.BW_smearing, linewidth=1, color='red')
        axs.text(hiDM, 0.85*self.BW_smearing, f"{self.dDM}", color='red')

        # channel smearing
        axs.plot(self.DMs, self.chan_smear(self.DMs), linewidth=1, color='blue')

        # subband smearing
        if (self.numsub):
            axs.plot(self.DMs, np.zeros(self.numDMs)+self.sub_smearing, linewidth=1, color='purple')
            axs.text(midDM, 0.85*self.sub_smearing, f"{self.dsubDM} ({self.numprepsub})")

        # scattering smearing
        if self.include_scattering:
            axs.plot(self.DMs, self.scattering_smear(self.DMs), linewidth=1, color='yellow')


    def __str__(self):
        if (self.numsub):
            return "%9.3f  %9.3f  %6.2f    %4d  %6.2f  %6d  %6d  %6d " % \
                   (self.loDM, self.hiDM, self.dDM, self.downsamp, self.dsubDM,
                    self.numDMs, self.DMs_per_prepsub, self.numprepsub)
        else:
            return "%9.3f  %9.3f  %6.2f    %4d  %6d" % \
                   (self.loDM, self.hiDM, self.dDM, self.downsamp, self.numDMs)

    def as_dict(self):
        """dictionary containing all information for the ddplan - for easy saving"""
        ret_dict = dict(
            obs=self.obs.as_dict(),
            downsamp=self.downsamp,
            numsub=self.numsub,
            numprocs=self.numprocs,
            smearfact=self.smearfact,
            scatteringfact=self.scatteringfact,
            include_scattering=self.include_scattering,
            DMs=self.DMs,
            DMs_per_prepsub=self.DMs_per_prepsub,
            dsubDM=self.dsubDM,
            numprepsub=self.numprepsub,
            dDM=self.dDM,
            hiDM=self.hiDM,
            loDM=self.loDM,
        )
        return ret_dict

    @classmethod
    def from_dict(cls, ddplan_dict):
        """make a instance from a dictionary in the same format as returned by as_dict
        (popping the obs dict out and making it separately was getting annoying)"""
        obs_dict = ddplan_dict.pop('obs')
        return cls(obs=observation(**obs_dict), **ddplan_dict)

    @classmethod
    def make(cls, obs, downsamp, loDM, hiDM, dDM, numDMs=0, numsub=0, numprocs=1, smearfact=2.0,
            scatteringfact=2.0, include_scattering=False):
        if hiDM and numDMs:
            print("Warning: setting numDMs will override hiDM like so:\nhiDM = loDM + numDMs*dDM")

        if (numsub):  # Calculate the maximum subband smearing we can handle
            DMs_per_prepsub = 2
            BW_smearing = BW_smear(dDM, obs.BW, obs.f_ctr)
            while(1):
                next_dsubDM = (DMs_per_prepsub+2) * dDM
                next_ss = subband_smear(next_dsubDM, numsub, obs.BW, obs.f_ctr)
                # The 0.8 is a small fudge factor to make sure that the subband
                # smearing is always the smallest contribution
                if (next_ss > 0.8*min(BW_smearing, 1000.0*obs.dt*downsamp)):
                    dsubDM = DMs_per_prepsub*dDM
                    break
                DMs_per_prepsub += 2
        else:
            dsubDM = dDM
            # these two aren't really accurate
            # (DMs_per_prepsub would be just he number of DMs and numprepsub would be 1)
            # but in the previous version they just aren't assigned and I wanted to keep in line with that
            DMs_per_prepsub = 0
            numprepsub = 0

        # Calculate the nominal DM to move to the next method
        DM_for_smearf = DM_for_smearfact(obs, dDM, dsubDM, numsub, downsamp, smearfact)
        if include_scattering:
            DM_for_scatteringf = DM_for_scatteringfact(obs, dDM, dsubDM, numsub, downsamp, scatteringfact)
            cross_DM = min(DM_for_smearf, DM_for_scatteringf)
        else:
            cross_DM = DM_for_smearf
        if (cross_DM > hiDM):
            cross_DM = hiDM

        if (numDMs==0):
            numDMs = int(np.ceil((cross_DM-loDM)/dDM))
            if (numsub):
                numprepsub = int(np.ceil(numDMs*dDM / dsubDM))
                if (numprocs > 1 and numprepsub % numprocs):
                    # Make sure the number of "calls" is a multiple of numprocs
                    numprepsub = (numprepsub // numprocs + 1) * numprocs
                    # Now adjust DMs_per_prepsub in case numprepsub increased a lot
                    while (DMs_per_prepsub > 1 and
                           numprepsub * DMs_per_prepsub > numDMs):
                        DMs_per_prepsub -= 1
                numDMs = numprepsub * DMs_per_prepsub

        # Make sure the number of DMs is divisible by the number of processors
        if (numprocs > 1 and numDMs % numprocs):
            numDMs = (numDMs // numprocs + 1) * numprocs
        hiDM = loDM + numDMs*dDM
        DMs = np.arange(numDMs, dtype='d')*dDM + loDM

        init_dict = dict(
            obs=obs,
            downsamp=downsamp,
            numsub=numsub,
            numprocs=numprocs,
            smearfact=smearfact,
            scatteringfact=scatteringfact,
            include_scattering=include_scattering,
            DMs=DMs,
            DMs_per_prepsub=DMs_per_prepsub,
            dsubDM=dsubDM,
            numprepsub=numprepsub,
            dDM=dDM,
            hiDM=hiDM,
            loDM=loDM,
        )
        #print(f"{init_dict}")
        return cls(**init_dict)


def DM_for_smearfact(obs, dDM, dsubDM, numsub, downsamp, smearfact):
    """
    Return the DM where the smearing in a single channel is a factor smearfact
    larger than all the other smearing causes combined.
    (except scattering)
    """
    BW_smearing = BW_smear(dDM, obs.BW, obs.f_ctr)
    sub_smearing = subband_smear(dsubDM, numsub, obs.BW, obs.f_ctr)
    other_smear = np.sqrt((1000.0*obs.dt)**2.0 +
                          (1000.0*obs.dt*downsamp)**2.0 +
                          BW_smearing**2.0 +
                          sub_smearing**2.0)
    return smearfact*0.001*other_smear/obs.chanwidth*0.0001205*obs.f_ctr**3.0 + obs.cDM

def DM_for_scatteringfact(obs, dDM, dsubDM, numsub, downsamp, scatteringfact):
    """
    Return the DM where the scattering at the top of the band is a factor scatteringfact
    larger than all the other smearing causes combined.
    (except in-channel DM smearing)
    """
    BW_smearing = BW_smear(dDM, obs.BW, obs.f_ctr)
    sub_smearing = subband_smear(dsubDM, numsub, obs.BW, obs.f_ctr)
    other_smear = np.sqrt((1000.0*obs.dt)**2.0 +
                          (1000.0*obs.dt*downsamp)**2.0 +
                          BW_smearing**2.0 +
                          sub_smearing**2.0)
    # from Bhat 2004
    A, B, C, D = -6.46, 0.154, 1.07, -3.86
    f_MHz =  obs.f_ctr + obs.BW/2

    det = B**2 - 4*C*(A+D*np.log10(f_MHz/1000) - np.log10(scatteringfact*other_smear))
    return 10**((-B + np.sqrt(det))/(2*C))


def choose_downsamps(blocklen):
    """
    choose_downsamps(blocklen):
        Return a good list of possible downsample sizes given a
        block of data of length blocklen spectra.
    """
    # This is first cut.  We will then remove redundant ones.
    x = np.asarray([n for n in np.arange(1, 260) if blocklen%n==0])
    if len(x)==1: return x
    # Now only choose those where the ratio is between 1.5 and 2, if possible
    if (x[1:]/x[:-1]).min() < 1.5:
        newx = [1]
        if 2 in x: newx.append(2)
        if 3 in x: newx.append(3)
        maxnewx = newx[-1]
        while maxnewx < x[-1]:
            if round(1.5*maxnewx+1e-7) in x:
                newx.append(round(1.5*maxnewx+1e-7))
            elif 2*maxnewx in x:
                newx.append(2*maxnewx)
            else:
                if x[-1] > 1.5*maxnewx:
                    newx.append(int(x[x>1.5*maxnewx].min()))
                else:
                    return newx
            maxnewx = newx[-1]
        return newx
    else:
        return x

def dm_smear(DM, BW, f_ctr, cDM=0.0):
    """
    dm_smear(DM, BW, f_ctr, cDM=0.0):
        Return the smearing in ms caused by a 'DM' over a bandwidth
        of 'BW' MHz centered at 'f_ctr' MHz.
    """
    return 1000.0*np.fabs(DM-cDM)*BW/(0.0001205*f_ctr**3.0)

def BW_smear(DMstep, BW, f_ctr):
    """
    BW_smear(DMstep, BW, f_ctr):
        Return the smearing in ms caused by a search using a DM stepsize of
        'DMstep' over a bandwidth of 'BW' MHz centered at 'f_ctr' MHz.
    """
    maxDMerror = 0.5*DMstep
    return dm_smear(maxDMerror, BW, f_ctr)

def guess_DMstep(dt, BW, f_ctr):
    """
    guess_DMstep(dt, BW, f_ctr):
        Choose a reasonable DMstep by setting the maximum smearing across the
        'BW' to equal the sampling time 'dt'.
    """
    return dt*0.0001205*f_ctr**3.0/(0.5*BW)

def subband_smear(subDMstep, numsub, BW, f_ctr):
    """
    subband_smear(subDMstep, numsub, BW, f_ctr):
        Return the smearing in ms caused by a search using a subband
        DM stepsize of 'subDMstep' over a total bandwidth of 'BW' MHz
        centered at 'f_ctr' MHz, and having numsub subbands.
    """
    if (numsub==0): return 0.0
    subBW = BW/numsub
    maxsubDMerror = 0.5*subDMstep
    return dm_smear(maxsubDMerror, subBW, f_ctr)

def scattering_smearing(f_MHz, DM):
    """returns scattering broadening at f_MHz in ms
    based on Bhat et al. (2004, ApJ,605, 759),"""
    if 'ndarray' in DM.__class__.__name__:
        if (DM < 1E-2).any():  # was just DM=0 but discovered it doesn't play nice with DMs < 1E-2
            scatter_smear = np.zeros_like(DM)
            scatter_smear[DM >= 1E-2] = calc_smear(f_MHz, DM[DM >= 1E-2])
            scatter_smear[DM < 1E-2] = 0
        else:
            scatter_smear = calc_smear(f_MHz, DM)
    elif DM < 1E-2:
        scatter_smear = 0
    else:
        scatter_smear = calc_smear(f_MHz, DM)
    return scatter_smear

def calc_smear(f_MHz, DM):
    log10DM = np.log10(DM)
    log10tau = -6.46 + 0.154*log10DM + 1.07*log10DM**2 -3.86*np.log10(f_MHz/1000)
    return 10**(log10tau)

def total_smear(DM, DMstep, dt, f_ctr, BW, numchan, subDMstep, cohdm=0.0, numsub=0, include_scattering=False):
    """
    total_smear(DM, DMstep, dt, f_ctr, BW, numchan, subDMstep, cohdm=0.0, numsub=0):
        Return the total smearing in ms due to the sampling rate,
        the smearing over each channel, the smearing over each subband
        (if numsub > 0), the smearing over the full BW assuming the
        worst-case DM error, and (if included) the smearing due to scattering at the
        top of the band.
    """
    if include_scattering:
        scatter = scattering_smearing(f_ctr+BW/2, DM)
    else:
        scatter = 0
    return np.sqrt(2 * (1000.0*dt)**2.0 +
                   scatter**2.0 +
                   dm_smear(DM, BW/numchan, f_ctr, cohdm)**2.0 +
                   subband_smear(subDMstep, numsub, BW, f_ctr)**2.0 +
                   BW_smear(DMstep, BW, f_ctr)**2.0)

def dm_steps(loDM, hiDM, obs, cohdm=0.0, numsub=0, numprocs=1,
             ok_smearing=0.0, blocklen=None, ax=None, include_scattering=False):
    """
    dm_steps(loDM, hiDM, obs, cohdm=0.0, numsub=0, numprocs=1,
             ok_smearing=0.0, blocklen=None, device="/XWIN"):
        Return the optimal DM stepsizes (and subband DM stepsizes if
        numsub>0) to keep the total smearing below 'ok_smearing' (in ms),
        for the DMs between loDM and hiDM.  If 'ok_smearing'=0.0, then
        use the best values based only on the data.  If the blocklen is
        not None, use it to determine possible downsampling values.
        And if device is not None, use it as the PGPLOT device for plotting.

        KC changed device to a matplotlib axis
    """
    # Allowable DM stepsizes
    allow_dDMs = [0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0,
    #allow_dDMs = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0,
                  2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0]

    # Allowable number of downsampling factors
    allow_downsamps = choose_downsamps(blocklen)

    # Initial values
    index_downsamps = index_dDMs = 0
    downsamp = allow_downsamps[index_downsamps]
    dDM = allow_dDMs[index_dDMs]
    dtms = 1000.0*obs.dt

    # Fudge factor that "softens" the boundary defining
    # if 2 time scales are equal or not
    ff = 1.2

    # This is the array that will hold the de-dispersion plans
    methods = []

    # Minimum possible smearing
    min_tot_smearing = total_smear(loDM+0.5*dDM, dDM, obs.dt, obs.f_ctr,
                                   obs.BW, obs.numchan, allow_dDMs[0], cohdm, 0)
    # Minimum channel smearing
    min_chan_smearing = dm_smear(np.linspace(loDM, hiDM, 10000),
                                 obs.chanwidth, obs.f_ctr, cohdm).min()
    # Minimum smearing across the obs.BW
    min_BW_smearing = BW_smear(dDM, obs.BW, obs.f_ctr)

    print()
    print("Minimum total smearing     : %.3g ms" % min_tot_smearing)
    print("--------------------------------------------")
    print("Minimum channel smearing   : %.3g ms" % min_chan_smearing)
    print("Minimum smearing across BW : %.3g ms" % min_BW_smearing)
    print("Minimum sample time        : %.3g ms" % dtms)
    print()

    ok_smearing = max([ok_smearing, min_chan_smearing, min_BW_smearing, dtms])
    print("Setting the new 'best' resolution to : %.3g ms" % ok_smearing)

    # See if the data is too high time resolution for our needs
    if (ff*min_chan_smearing > dtms or
        ok_smearing > dtms):
        if (ok_smearing > ff*min_chan_smearing):
            print("   Note: ok_smearing > dt (i.e. data is higher resolution than needed)")
            okval = ok_smearing
        else:
            print("   Note: min_chan_smearing > dt (i.e. data is higher resolution than needed)")
            okval = ff*min_chan_smearing

        while (dtms*allow_downsamps[index_downsamps+1] < okval):
            index_downsamps += 1
        downsamp = allow_downsamps[index_downsamps]
        print("         New dt is %d x %.12g ms = %.12g ms" % \
              (downsamp, dtms, dtms*downsamp))

    # Calculate the appropriate initial dDM
    dDM = guess_DMstep(obs.dt*downsamp, obs.BW, obs.f_ctr)
    print("Best guess for optimal initial dDM is %.3f" % dDM)
    while (allow_dDMs[index_dDMs+1] < ff*dDM):
        index_dDMs += 1

    # If numprocs > 1, we are using mpiprepsubband, so let the
    # user know
    if (numprocs > 1):
        print("\nAssuming we are using mpiprepsubband with %d dedispersing CPUs:"
              % numprocs)
        print("Each line of the dedispersion plan is one or more distinct runs of")
        print("mpiprepsubband, and each 'call' is the work that a single CPU is doing.")

    # Create the first method
    methods = [dedisp_method.make(obs, downsamp, loDM, hiDM,
                             allow_dDMs[index_dDMs], numsub=numsub,
                             numprocs=numprocs, include_scattering=include_scattering)]
    numDMs = [methods[-1].numDMs]

    # Calculate the next methods
    while(methods[-1].hiDM < hiDM):

        # Determine the new downsample factor
        index_downsamps += 1
        downsamp = allow_downsamps[index_downsamps]
        eff_dt = dtms*downsamp

        # Determine the new DM step
        while (BW_smear(allow_dDMs[index_dDMs+1], obs.BW, obs.f_ctr) < ff*eff_dt):
            index_dDMs += 1
        dDM = allow_dDMs[index_dDMs]

        # Get the next method
        methods.append(dedisp_method.make(obs, downsamp, methods[-1].hiDM,
                                     hiDM, dDM, numsub=numsub,
                                     numprocs=numprocs, include_scattering=include_scattering))
        numDMs.append(methods[-1].numDMs)

    # Calculate the DMs to search and the smearing at each
    total_numDMs = sum(numDMs)
    DMs = np.zeros(total_numDMs, dtype='d')
    total_smears = np.zeros(total_numDMs, dtype='d')

    # Calculate the DMs and optimal smearing for all the DMs
    for ii, offset in enumerate(np.add.accumulate([0]+numDMs[:-1])):
        DMs[offset:offset+numDMs[ii]] = methods[ii].DMs
        total_smears[offset:offset+numDMs[ii]] = methods[ii].total_smear(methods[ii].DMs)

    # Calculate the predicted amount of time that will be spent in searching
    # this batch of DMs as a fraction of the total
    work_fracts = [meth.numDMs/float(meth.downsamp) for meth in methods]
    work_fracts = np.asarray(work_fracts)/sum(work_fracts)

    # The optimal smearing
    tot_smear = total_smear(DMs, allow_dDMs[0], obs.dt, obs.f_ctr,
                            obs.BW, obs.numchan, allow_dDMs[0], cohdm, 0)

    if ax is not None:
        # Plot them
        ax.plot(DMs, tot_smear, c="orange")
        ax.set_yscale('log')
        ax.set_xlim([loDM, hiDM])
        ax.set_ylim([0.3*min(tot_smear), 2.5*max(tot_smear)])
        ax.set_xlabel("Dispersion Measure (pc/cm$^3$)")
        ax.set_ylabel("Smearing (ms)")

        # add some info as the title
        title = f"$f_{{{'ctr'}}}$ = {obs.f_ctr:d} MHz \t "
        if dtms < 0.1:
            title += f"dt = {dtms*1000} ms \t "
        else:
            title += f"dt = {dtms} ms \t "
        title += f"BW = {obs.BW:d} MHz \t $N_{{{'chan'}}}$ = {obs.numchan:d}"
        if numsub:
            title += f" \t $N_{{{'sub'}}}$ = {numsub:d}"

        ax.set_title(title.expandtabs())

        x_txt = 0.98
        y_txt = 0.0
        dy_txt = +0.05
        ax.text(x_txt, y_txt + 7*dy_txt, "Total Smearing", color='black', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes,)
        ax.text(x_txt, y_txt + 6*dy_txt, "Optimal Smearing", color='orange', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes,)
        if cohdm:
            ax.text(x_txt, y_txt + 5*dy_txt, "Chan Smearing (w/ coherent dedisp)", color='blue', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes,)
        else:
            ax.text(x_txt, y_txt + 5*dy_txt, "Channel Smearing", color='blue', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes,)
        ax.text(x_txt, y_txt + 4*dy_txt, "Sample Time (ms)", color='green', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes,)
        ax.text(x_txt, y_txt + 3*dy_txt, "DM Stepsize Smearing", color='red', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes,)
        if include_scattering:
            ax.text(x_txt, y_txt + 2*dy_txt, f"Scattering Smearing ({int(obs.f_ctr+obs.BW/2)} MHz)", color='yellow', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes,)
        if numsub:
            ax.text(x_txt, y_txt + 1*dy_txt, "Subband Stepsize Smearing (# passes)", color='purple', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes,)

        if (numsub):
            print("\n  Low DM    High DM     dDM  DownSamp  dsubDM   #DMs  DMs/call  calls  WorkFract")
        else:
            print("\n  Low DM    High DM     dDM  DownSamp   #DMs  WorkFract")
        for method, fract in zip(methods, work_fracts):
            print(method, "  %.4g" % fract)
            method.plot(fract, ax)
        print("\n\n")

    return methods


############### Writing dedisp.py ####################

def truncate(n, decimals=0):
    """truncation version of rounding
    truncate(7.398493, 2) would give 7.39"""
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

# ugh this is gross, but I can't think of a nicer way to do it right now
dedisp_template1 = """
from __future__ import print_function
from builtins import zip
from builtins import range
import os
def myexecute(cmd):
    print("'%s'"%cmd)
    os.system(cmd)
# By default, do not output subbands
outsubs = False
"""

dedisp_template2 = """
# Loop over the DDplan plans
for dDM, dsubDM, dmspercall, downsamp, subcall, startDM in zip(dDMs, dsubDMs, dmspercalls, downsamps, subcalls, startDMs):
    # Loop over the number of calls
    for ii in range(subcall):
        subDM = startDM + (ii+0.5)*dsubDM
        loDM = startDM + ii*dsubDM
        if outsubs:
            # Get our downsampling right
            subdownsamp = downsamp // 2
            datdownsamp = 2
            if downsamp < 2: subdownsamp = datdownsamp = 1
            # First create the subbands
            myexecute("prepsubband -sub -subdm %.2f -nsub %d -downsamp %d -o %s %s" %
                      (subDM, nsub, subdownsamp, basename, rawfiles))
            # And now create the time series
            subnames = basename+"_DM%.2f.sub[0-9]*"%subDM
            myexecute("prepsubband -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -o %s %s" %
                      (loDM, dDM, dmspercall, datdownsamp, basename, subnames))
        else:
            myexecute("prepsubband -nsub %d -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -o %s %s" %
                      (nsub, loDM, dDM, dmspercall, downsamp, basename, rawfiles))
"""

dedisp_template2_with_mask = """
# Loop over the DDplan plans
for dDM, dsubDM, dmspercall, downsamp, subcall, startDM in zip(dDMs, dsubDMs, dmspercalls, downsamps, subcalls, startDMs):
    # Loop over the number of calls
    for ii in range(subcall):
        subDM = startDM + (ii+0.5)*dsubDM
        loDM = startDM + ii*dsubDM
        if outsubs:
            # Get our downsampling right
            subdownsamp = downsamp // 2
            datdownsamp = 2
            if downsamp < 2: subdownsamp = datdownsamp = 1
            # First create the subbands
            myexecute("prepsubband -sub -subdm %.2f -nsub %d -downsamp %d -mask %s_rfifind.mask -o %s %s" %
                      (subDM, nsub, subdownsamp, basename, basename, rawfiles))
            # And now create the time series
            subnames = basename+"_DM%.2f.sub[0-9]*"%subDM
            myexecute("prepsubband -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -o %s %s" %
                      (loDM, dDM, dmspercall, datdownsamp, basename, subnames))
        else:
            myexecute("prepsubband -nsub %d -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -mask %s_rfifind.mask -o %s %s" %
                      (nsub, loDM, dDM, dmspercall, downsamp, basename, basename, rawfiles))
"""

dedisp_template2_nosubbanding_withmask = """
# No subbanding
# Loop over the DDplan plans
for dDM, dmpermethod, downsamp, startDM in zip(dDMs, dmspermethod, downsamps, startDMs):
    myexecute("prepsubband -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -mask %s_rfifind.mask -o %s %s" %
                      (loDM, dDM, dmpermethod, downsamp, basename, basename, rawfiles))
"""


dedisp_template2_nosubbanding = """
# No subbanding
# Loop over the DDplan plans
for dDM, dmpermethod, downsamp, startDM in zip(dDMs, dmspermethod, downsamps, startDMs):
    myexecute("prepsubband -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -o %s %s" %
                      (loDM, dDM, dmpermethod, downsamp, basename, rawfiles))
"""

def print_dedisp_methods(methods, rawfiles, mask=True):
    """
    Print a dedisp.py for methods

    if mask=True, include the --mask flag in prepsubband,
        assuming mask of form filename_rfifind.mask
    """
    dDMs = [m.dDM for m in methods]
    startDMs = [m.loDM for m in methods]
    downsamps = [m.downsamp for m in methods]
    try:
        dsubDMs = [m.dsubDM for m in methods]
        dmspercalls = [m.DMs_per_prepsub for m in methods]
        subcalls = [m.numprepsub for m in methods]
        subbanding = True
    except AttributeError:
        dmspermethod = [len(m.DMs) for m in methods]
        subbanding = False

    basename, ext = os.path.splitext(rawfiles)
    print(dedisp_template1)
    print(f"nsub = {numsubbands}\n\n")
    print(f"basename = {basename}\n")
    print(f"rawfiles = {rawfiles}\n\n")

    print(f"""# dDM steps from DDplan.py
dDMs        = {dDMs}\n""")
    if subbanding:
        print(f"""# dsubDM steps
dsubDMs     = {dsubDMs}\n""")
    print(f"""# downsample factors
downsamps   = {downsamps}\n""")
    if subbanding:
        print(f"""# number of calls per set of subbands
subcalls    = {subcalls}\n""")
    # truncated this to avoid 0.435000000000002-esque numbers
    # may come back to bite me in the arse
    print(f"""# The low DM for each set of DMs
startDMs    = {[truncate(x, 3) for x in startDMs]}\n""")
    if subbanding:
        print(f"""# DMs/call
dmspercalls = {dmspercalls}\n""")
    else:
        print(f"""# DMs/plan
dmspermethod = {dmspermethod}\n""")

    if subbanding:
        if mask:
            print(dedisp_template2_with_mask)  # untested
        else:
            print(dedisp_template2)
    else:
        if mask:
            print(dedisp_template2_nosubbanding_withmask)
        else:
            print(dedisp_template2_nosubbanding)

def write_dedisp_methods_to_file(methods, rawfiles, mask=True):
    """
    Write methods to dedisp_filename.py

    if mask=True, include the --mask flag in prepsubband,
        assuming mask of form filename_rfifind.mask
    """
    dDMs = [m.dDM for m in methods]
    startDMs = [m.loDM for m in methods]
    downsamps = [m.downsamp for m in methods]
    try:
        dsubDMs = [m.dsubDM for m in methods]
        dmspercalls = [m.DMs_per_prepsub for m in methods]
        subcalls = [m.numprepsub for m in methods]
        subbanding = True
    except AttributeError:
        dmspermethod = [len(m.DMs) for m in methods]
        subbanding = False

    basename, ext = os.path.splitext(rawfiles)

    with open('dedisp_%s.py'%basename, 'w') as f:
        f.write(dedisp_template1)
        f.write(f"nsub = {numsubbands}\n\n")
        f.write(f"basename = {basename}\n")
        f.write(f"rawfiles = {rawfiles}\n\n")

        f.write(f"""# dDM steps from DDplan.py
dDMs        = {dDMs}\n""")
        if subbanding:
            print(f"""# dsubDM steps
dsubDMs     = {dsubDMs}\n""")
        f.write(f"""# downsample factors
downsamps   = {downsamps}\n""")
        if subbanding:
            f.write(f"""# number of calls per set of subbands
subcalls    = {subcalls}\n""")
        # truncated this to avoid 0.435000000000002-esque numbers
        # may come back to bite me in the arse
        f.write(f"""# The low DM for each set of DMs
startDMs    = {[truncate(x, 3) for x in startDMs]}\n""")
        if subbanding:
            f.write(f"""# DMs/call
dmspercalls = {dmspercalls}\n""")
        else:
            f.write(f"""# DMs/plan
dmspermethod = {dmspermethod}\n""")

        if subbanding:
            if mask:
                f.write(dedisp_template2_with_mask)
            else:
                f.write(dedisp_template2)
        else:
            if mask:
                f.write(dedisp_template2_nosubbanding_withmask)
            else:
                f.write(dedisp_template2_nosubbanding)

def write_dedisp_methods(methods, filfilename, mask=True, to_file=True):
    if to_file:
        write_dedisp_methods_to_file(methods, filfilename, mask=mask)
    else:
        print_dedisp_methods(methods, filfilename, mask=mask)

######################################################

################## ddplan classes ####################

# made a ddplan class
# stores a set of ddplan_methods, can write and read them to file

@attrs
class ddplan:
    """A dedispersion plan - a collection of dedisp_method's"""
    methods: list = attrib(factory=list)

    @classmethod
    def read_from_npz(filename):


        init_dict = dict(
            obs=observation(**data['obs']),
            downsamp=data['downsamp'],
            numsub=data['numsub'],
            numprocs=data['numprocs'],
            smearfact=data['smearfact'],
            scatteringfact=data['scatteringfact'],
            include_scattering=data['include_scattering'],
            DMs=data['DMs'],
            DMs_per_prepsub=data['DMs_per_prepsub'],
            dsubDM=data['dsubDM'],
            numprepsub=data['numprepsub'],
            dDM=data['dDM'],
            hiDM=data['hiDM'],
            loDM=data['loDM'],
        )

        return cls(**init_dict)

    def write_dedisp_py(filfilename, to_file=True, mask=True):
        """
        Write a dedisp.py

        If to_file=True it will be writted to a file of the form
        dedisp_filfile.py, presuming filfilename was filfile.fil
        If to_file=False, output is printed.

        mask=True includes the --mask option in the prepsubband call,
        assuming it's of the form filfile_rfifind.mask"""
        write_dedisp_methods(self.methods, filfilename, mask=mask, to_file=to_file)

    def write_to_npz(self, filename):
        """Save ddplan to npz file"""
        method_dicts = []
        for method in self.methods:
            method_dicts.append(method.as_dict())
        np.savez(filename, method_dicts=method_dicts)
        print(f"Wrote {len(method_dicts)} dedisp_methods to {filename}")

    @classmethod
    def read_from_npz(cls, filename):
        data = np.load(filename, allow_pickle=True)
        methods = []
        for method_dict in data['method_dicts']:
            methods.append(dedisp_method.from_dict(method_dict))
        print(f"Read {len(methods)} dedisp_methods from {filename}")
        return cls(methods=methods)

@attrs
class frankenstein_ddplan(ddplan):
    """A ddplan but with evenly spaced coherently dedispersed observations up to some limit"""
    methods_by_cdm: dict = attrib(factory=dict)

    @property
    def cdms(self):
        return list(self.methods_by_cdm.keys())

    @property
    def hiDM(self):
        # NB if you make it yourself or add methods manually this may not be true
        return self.methods[-1].hiDM

    @property
    def loDM(self):
        # NB if you make it yourself or add methods manually this may not be true
        return self.methods[0].loDM


    def write_dedisp_for(self, fname, to_file=True, mask=True):
        """
        When have a coherently dedispersed filterbank (with a filename like J2119+49_DM273_59488_pow.fil)
        you only want the part of the frankenstein plan at that coherent DM.
        This selects the relevent ddplan arc, and writes out the appropriate dedisp.py
        (NB the DM in the filename must match (to 0dp) a cdm in the methods_by_cdm keys)

        if to_file=False, result will be printed
        otherwise it's written to dedisp_filename.py
        """
        splt_fname = fname.split("_")
        trunc_DM = int(splt_fname[1][2:])
        cDM = None
        for key in list(self.methods_by_cdm.keys()):
            if key // 1 == trunc_DM:
                cDM = key
                break
        if cDM is None:
            print(f"Match for {trunc_DM} not found among keys {list(self.methods_by_cdm.keys())}")
        else:
            write_dedisp_methods(self.methods_by_cdm[cDM], fname, to_file=to_file, mask=mask)

    @classmethod
    def make(cls, ddplan_methods, max_meth, hi_DM, troubleshoot=False, extend_up_to=0):
        """
        takes ddplan_methods[:max_meth] and repeats them across the DM range up to hi_DM, producing
        ddplan_methods spaced out at a series of coherent DMs

        At the moment this will give multiple coherently dedsipersed ddplan arcs, until cDM > hi_DM
        extend_up_to let's you extend the DDplan above this
        e.g. if you only want to coherently dedisperse up to 300, but want a plan up to 500,
        pass in 300 as hi_DM and set extend_up_to to 500
        """
        # NB assumes lo_DM is 0
        base_methods = ddplan_methods[0:max_meth]
        half_span_per_cDM = base_methods[-1].hiDM

        if extend_up_to == 0:
            extend_up_to = hi_DM

        def print_message(message, troubleshoot=troubleshoot):
            if troubleshoot:
                print(message)

        coh_DMs = []
        amalgamated_ddplan_methods = []
        plans_split_by_cdm = {}

        cDM = base_methods[-1].hiDM

        while cDM <= hi_DM:
            cDM_plans = []
            coh_DMs.append(cDM)
            curr_meth = base_methods[0]
            tmp_obs = observation(curr_meth.obs.dt, curr_meth.obs.f_ctr, curr_meth.obs.BW, curr_meth.obs.numchan, cDM)

            print_message("## LEFT ##")
            for i in range(max_meth)[::-1]:
                curr_meth = ddplan_methods[i]
                new_lo = cDM - curr_meth.hiDM
                new_hi = min(cDM - curr_meth.loDM, hi_DM)
                new_meth = dedisp_method.make(
                    tmp_obs,
                    curr_meth.downsamp,
                    new_lo,
                    new_hi,
                    curr_meth.dDM,
                    numsub=numsubbands,
                    numprocs=numprocs,
                    include_scattering=False
                )
                print_message(f"{i}: {new_meth.loDM} -> {new_meth.hiDM} in steps of {new_meth.dDM}")
                amalgamated_ddplan_methods.append(new_meth)
                cDM_plans.append(new_meth)

            print_message("## RIGHT ##")
            # should always be in the right part of the arc when you hit the end of the range
            for i in range(max_meth):
                curr_meth = ddplan_methods[i]
                new_lo = cDM + curr_meth.loDM
                new_hi = min(cDM + curr_meth.hiDM, extend_up_to)
                new_meth = dedisp_method.make(
                    tmp_obs,
                    curr_meth.downsamp,
                    new_lo,
                    new_hi,
                    curr_meth.dDM,
                    numsub=numsubbands,
                    numprocs=numprocs,
                    include_scattering=False
                )
                print_message(f"{i}: {new_meth.loDM} -> {new_meth.hiDM} in steps of {new_meth.dDM}")
                amalgamated_ddplan_methods.append(new_meth)
                cDM_plans.append(new_meth)
                if new_hi == extend_up_to:
                    plans_split_by_cdm[cDM] = cDM_plans
                    return cls(methods=amalgamated_ddplan_methods, methods_by_cdm=plans_split_by_cdm)

            plans_split_by_cdm[cDM] = cDM_plans
            cDM += 2*base_methods[-1].hiDM
            print_message(f"cDM is now {cDM}")

        # extend the last ddplan arc, up to the top of the desired range
        #reset cDM to last value
        cDM = amalgamated_ddplan_methods[-1].obs.cDM
        i+=1
        while new_hi < extend_up_to:
            try:
                curr_meth = ddplan_methods[i]
                new_lo = cDM + curr_meth.loDM
                new_hi = min(cDM + curr_meth.hiDM, extend_up_to)
                new_meth = dedisp_method.make(
                    tmp_obs,
                    curr_meth.downsamp,
                    new_lo,
                    new_hi,
                    curr_meth.dDM,
                    numsub=numsubbands,
                    numprocs=numprocs,
                    include_scattering=False
                )
                print_message(f"{i}: {new_meth.loDM} -> {new_meth.hiDM} in steps of {new_meth.dDM}")
                amalgamated_ddplan_methods.append(new_meth)
                plans_split_by_cdm[cDM].append(new_meth)
            except IndexError:
                print(
                    f"Ran out of ddplan methods at {i}\n"
                    f"Got up to a DM of {amalgamated_ddplan_methods[-1].hiDM}\n"
                    f"Make ddplan_methods with a larger DM range to get up to wanted DM of {extend_up_to}"
                )
                break

            i += 1

        return cls(methods=amalgamated_ddplan_methods, methods_by_cdm=plans_split_by_cdm)

    def write_to_npz(self, filename):
        """Save ddplan to npz file"""
        method_dicts = []
        for method in self.methods:
            method_dicts.append(method.as_dict())


        cdms=[]

        save_dict = dict(
            method_dicts=method_dicts,
        )

        for key in list(self.methods_by_cdm.keys()):
            tmp_meth_list = [meth.as_dict() for meth in self.methods_by_cdm[key]]
            save_dict[str(key)] = tmp_meth_list
            cdms.append(key)


        save_dict['cdms']=cdms

        np.savez(filename, **save_dict)
        print(f"Wrote {len(method_dicts)} frankensteined dedisp_methods to {filename}")

    @classmethod
    def read_from_npz(cls, filename):
        """Read ddplan from npz"""
        data = np.load(filename, allow_pickle=True)
        methods = []
        for method_dict in data['method_dicts']:
            obs_dict = method_dict.pop('obs')
            methods.append(dedisp_method(obs=observation(**obs_dict), **method_dict))

        methods_split_by_cdm = {}
        for cdm in data['cdms']:
            tmp_meth_list = [dedisp_method.from_dict(method_dict) for method_dict in data[str(cdm)]]
            methods_split_by_cdm[cdm] = tmp_meth_list

        print(f"Read {len(methods)} frankensteined dedisp_methods from {filename}")
        return cls(methods=methods, methods_by_cdm=methods_split_by_cdm)

###############################################################
