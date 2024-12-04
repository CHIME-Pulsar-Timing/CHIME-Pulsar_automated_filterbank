# straight copy 14/01/2022
from builtins import object

## Automatically adapted for numpy Apr 14, 2006 by convertcode.py


class infodata(object):
    def __init__(self, filenm):
        self.breaks = 0
        for line in open(filenm, encoding="latin-1"):
            if line.startswith(" Data file name"):
                self.basenm = line.split("=")[-1].strip()
                continue
            if line.startswith(" Telescope"):
                self.telescope = line.split("=")[-1].strip()
                continue
            if line.startswith(" Instrument"):
                self.instrument = line.split("=")[-1].strip()
                continue
            if line.startswith(" Object being observed"):
                self.object = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Right Ascension"):
                self.RA = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Declination"):
                self.DEC = line.split("=")[-1].strip()
                continue
            if line.startswith(" Data observed by"):
                self.observer = line.split("=")[-1].strip()
                continue
            if line.startswith(" Epoch"):
                self.epoch = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Barycentered?"):
                self.bary = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of bins"):
                self.N = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Width of each time series bin"):
                self.dt = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Any breaks in the data?"):
                self.breaks = int(line.split("=")[-1].strip())
                if self.breaks:
                    self.onoff = []
                continue
            if line.startswith(" On/Off bin pair"):
                vals = line.split("=")[-1].strip().split(",")
                self.onoff.append((int(vals[0]), int(vals[1])))
                continue
            if line.startswith(" Type of observation"):
                self.waveband = line.split("=")[-1].strip()
                continue
            if line.startswith(" Beam diameter"):
                self.beam_diam = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Dispersion measure"):
                self.DM = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Central freq of low channel"):
                self.lofreq = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Total bandwidth"):
                self.BW = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of channels"):
                self.numchan = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Channel bandwidth"):
                self.chan_width = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Data analyzed by"):
                self.analyzer = line.split("=")[-1].strip()
                continue
            if line.startswith(" Orbit removed?"):
                self.deorbited = int(line.split("=")[-1].strip())
                continue

    def to_file(self, inffn, notes=None):
        if not inffn.endswith(".inf"):
            raise ValueError(
                "PRESTO info files must end with '.inf'. " "Got: %s" % inffn
            )
        with open(inffn, "w") as ff:
            if hasattr(self, "basenm"):
                ff.write(
                    " Data file name without suffix          =  %s\n" % self.basenm
                )
            if hasattr(self, "telescope"):
                ff.write(
                    " Telescope used                         =  %s\n" % self.telescope
                )
            if hasattr(self, "instrument"):
                ff.write(
                    " Instrument used                        =  %s\n" % self.instrument
                )
            if hasattr(self, "object"):
                ff.write(
                    " Object being observed                  =  %s\n" % self.object
                )
            if hasattr(self, "RA"):
                ff.write(" J2000 Right Ascension (hh:mm:ss.ssss)  =  %s\n" % self.RA)
            if hasattr(self, "DEC"):
                ff.write(" J2000 Declination     (dd:mm:ss.ssss)  =  %s\n" % self.DEC)
            if hasattr(self, "observer"):
                ff.write(
                    " Data observed by                       =  %s\n" % self.observer
                )
            if hasattr(self, "epoch"):
                ff.write(
                    " Epoch of observation (MJD)             =  %05.15f\n" % self.epoch
                )
            if hasattr(self, "bary"):
                ff.write(" Barycentered?           (1=yes, 0=no)  =  %d\n" % self.bary)
            if hasattr(self, "N"):
                ff.write(
                    " Number of bins in the time series      =  %-11.0f\n" % self.N
                )
            if hasattr(self, "dt"):
                ff.write(" Width of each time series bin (sec)    =  %.15g\n" % self.dt)
            if hasattr(self, "breaks") and self.breaks:
                ff.write(" Any breaks in the data? (1 yes, 0 no)  =  1\n")
                if hasattr(self, "onoff"):
                    for ii, (on, off) in enumerate(self.onoff, 1):
                        ff.write(
                            " On/Off bin pair #%3d                   =  %-11.0f, %-11.0f\n"
                            % (ii, on, off)
                        )
            else:
                ff.write(" Any breaks in the data? (1 yes, 0 no)  =  0\n")
            # These two were left out, I assume because there's some variation? e.g. beam_diam might not always be in arcsec
            if hasattr(self, "waveband"):
                ff.write(
                    " Type of observation (EM band)          =  %s\n" % self.waveband
                )
            if hasattr(self, "beam_diam"):
                ff.write(
                    " Beam diameter (arcsec)                 =  %d\n" % self.beam_diam
                )
            if hasattr(self, "DM"):
                ff.write(" Dispersion measure (cm-3 pc)           =  %.12g\n" % self.DM)
            if hasattr(self, "lofreq"):
                ff.write(
                    " Central freq of low channel (Mhz)      =  %.12g\n" % self.lofreq
                )
            if hasattr(self, "BW"):
                ff.write(" Total bandwidth (Mhz)                  =  %.12g\n" % self.BW)
            if hasattr(self, "numchan"):
                ff.write(
                    " Number of channels                     =  %d\n" % self.numchan
                )
            if hasattr(self, "chan_width"):
                ff.write(
                    " Channel bandwidth (Mhz)                =  %.12g\n"
                    % self.chan_width
                )
            if hasattr(self, "analyzer"):
                ff.write(
                    " Data analyzed by                       =  %s\n" % self.analyzer
                )
            if hasattr(self, "deorbited"):
                ff.write(
                    " Orbit removed?          (1=yes, 0=no)  =  %d\n" % self.deorbited
                )
            ff.write(" Any additional notes:\n")
            if notes is not None:
                ff.write("    %s\n" % notes.strip())


# added tweaked version so can make it from a dict
class infodata2(object):
    def __init__(self, init_dict):
        if init_dict.get("basenm", None) is not None:
            self.basenm = init_dict["basenm"]
        if init_dict.get("telescope", None) is not None:
            self.telescope = init_dict["telescope"]
        if init_dict.get("instrument", None) is not None:
            self.instrument = init_dict["instrument"]
        if init_dict.get("object", None) is not None:
            self.object = init_dict["object"]
        if init_dict.get("RA", None) is not None:
            self.RA = init_dict["RA"]
        if init_dict.get("DEC", None) is not None:
            self.DEC = init_dict["DEC"]
        if init_dict.get("observer", None) is not None:
            self.observer = init_dict["observer"]
        if init_dict.get("epoch", None) is not None:
            self.epoch = init_dict["epoch"]
        if init_dict.get("bary", None) is not None:
            self.bary = init_dict["bary"]
        if init_dict.get("N", None) is not None:
            self.N = init_dict["N"]
        if init_dict.get("dt", None) is not None:
            self.dt = init_dict["dt"]
        if init_dict.get("breaks", None) is not None:
            self.breaks = init_dict["breaks"]
        if init_dict.get("onoff", None) is not None:
            self.onoff = init_dict["onoff"]
        if init_dict.get("waveband", None) is not None:
            self.waveband = init_dict["waveband"]
        if init_dict.get("beam_diam", None) is not None:
            self.beam_diam = init_dict["beam_diam"]
        if init_dict.get("DM", None) is not None:
            self.DM = init_dict["DM"]
        if init_dict.get("lofreq", None) is not None:
            self.lofreq = init_dict["lofreq"]
        if init_dict.get("BW", None) is not None:
            self.BW = init_dict["BW"]
        if init_dict.get("numchan", None) is not None:
            self.numchan = init_dict["numchan"]
        if init_dict.get("chan_width", None) is not None:
            self.chan_width = init_dict["chan_width"]
        if init_dict.get("analyzer", None) is not None:
            self.analyzer = init_dict["analyzer"]
        if init_dict.get("deorbited", None) is not None:
            self.deorbited = init_dict["deorbited"]

    @classmethod
    def from_file(cls, filenm):
        init_dict = {}
        init_dict["breaks"] = 0
        for line in open(filenm, encoding="latin-1"):
            if line.startswith(" Data file name"):
                init_dict["basenm"] = line.split("=")[-1].strip()
                continue
            if line.startswith(" Telescope"):
                init_dict["telescope"] = line.split("=")[-1].strip()
                continue
            if line.startswith(" Instrument"):
                init_dict["instrument"] = line.split("=")[-1].strip()
                continue
            if line.startswith(" Object being observed"):
                init_dict["object"] = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Right Ascension"):
                init_dict["RA"] = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Declination"):
                init_dict["DEC"] = line.split("=")[-1].strip()
                continue
            if line.startswith(" Data observed by"):
                init_dict["observer"] = line.split("=")[-1].strip()
                continue
            if line.startswith(" Epoch"):
                init_dict["epoch"] = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Barycentered?"):
                init_dict["bary"] = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of bins"):
                init_dict["N"] = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Width of each time series bin"):
                init_dict["dt"] = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Any breaks in the data?"):
                init_dict["breaks"] = int(line.split("=")[-1].strip())
                if init_dict["breaks"]:
                    init_dict["onoff"] = []
                continue
            if line.startswith(" On/Off bin pair"):
                vals = line.split("=")[-1].strip().split(",")
                init_dict["onoff"].append((int(vals[0]), int(vals[1])))
                continue
            if line.startswith(" Type of observation"):
                init_dict["waveband"] = line.split("=")[-1].strip()
                continue
            if line.startswith(" Beam diameter"):
                init_dict["beam_diam"] = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Dispersion measure"):
                init_dict["DM"] = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Central freq of low channel"):
                init_dict["lofreq"] = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Total bandwidth"):
                init_dict["BW"] = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of channels"):
                init_dict["numchan"] = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Channel bandwidth"):
                init_dict["chan_width"] = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Data analyzed by"):
                init_dict["analyzer"] = line.split("=")[-1].strip()
                continue
            if line.startswith(" Orbit removed?"):
                init_dict["deorbited"] = int(line.split("=")[-1].strip())
                continue
        return cls(init_dict)

    def to_file(self, inffn, notes=None):
        if not inffn.endswith(".inf"):
            raise ValueError(
                "PRESTO info files must end with '.inf'. " "Got: %s" % inffn
            )
        with open(inffn, "w") as ff:
            if hasattr(self, "basenm"):
                ff.write(
                    " Data file name without suffix          =  %s\n" % self.basenm
                )
            if hasattr(self, "telescope"):
                ff.write(
                    " Telescope used                         =  %s\n" % self.telescope
                )
            if hasattr(self, "instrument"):
                ff.write(
                    " Instrument used                        =  %s\n" % self.instrument
                )
            if hasattr(self, "object"):
                ff.write(
                    " Object being observed                  =  %s\n" % self.object
                )
            if hasattr(self, "RA"):
                ff.write(" J2000 Right Ascension (hh:mm:ss.ssss)  =  %s\n" % self.RA)
            if hasattr(self, "DEC"):
                ff.write(" J2000 Declination     (dd:mm:ss.ssss)  =  %s\n" % self.DEC)
            if hasattr(self, "observer"):
                ff.write(
                    " Data observed by                       =  %s\n" % self.observer
                )
            if hasattr(self, "epoch"):
                ff.write(
                    " Epoch of observation (MJD)             =  %05.15f\n" % self.epoch
                )
            if hasattr(self, "bary"):
                ff.write(" Barycentered?           (1=yes, 0=no)  =  %d\n" % self.bary)
            if hasattr(self, "N"):
                ff.write(
                    " Number of bins in the time series      =  %-11.0f\n" % self.N
                )
            if hasattr(self, "dt"):
                ff.write(" Width of each time series bin (sec)    =  %.15g\n" % self.dt)
            if hasattr(self, "breaks") and self.breaks:
                ff.write(" Any breaks in the data? (1 yes, 0 no)  =  1\n")
                if hasattr(self, "onoff"):
                    for ii, (on, off) in enumerate(self.onoff, 1):
                        ff.write(
                            " On/Off bin pair #%3d                   =  %-11.0f, %-11.0f\n"
                            % (ii, on, off)
                        )
            else:
                ff.write(" Any breaks in the data? (1 yes, 0 no)  =  0\n")
                # NB exploredat won't work with no On/Off pair lines

            # These two were left out, I assume because there's some variation? e.g. beam_diam might not always be in arcsec
            if hasattr(self, "waveband"):
                ff.write(
                    " Type of observation (EM band)          =  %s\n" % self.waveband
                )
            if hasattr(self, "beam_diam"):
                ff.write(
                    " Beam diameter (arcsec)                 =  %d\n" % self.beam_diam
                )
            if hasattr(self, "DM"):
                ff.write(" Dispersion measure (cm-3 pc)           =  %.12g\n" % self.DM)
            if hasattr(self, "lofreq"):
                ff.write(
                    " Central freq of low channel (Mhz)      =  %.12g\n" % self.lofreq
                )
            if hasattr(self, "BW"):
                ff.write(" Total bandwidth (Mhz)                  =  %.12g\n" % self.BW)
            if hasattr(self, "numchan"):
                ff.write(
                    " Number of channels                     =  %d\n" % self.numchan
                )
            if hasattr(self, "chan_width"):
                ff.write(
                    " Channel bandwidth (Mhz)                =  %.12g\n"
                    % self.chan_width
                )
            if hasattr(self, "analyzer"):
                ff.write(
                    " Data analyzed by                       =  %s\n" % self.analyzer
                )
            if hasattr(self, "deorbited"):
                ff.write(
                    " Orbit removed?          (1=yes, 0=no)  =  %d\n" % self.deorbited
                )
            ff.write(" Any additional notes:\n")
            if notes is not None:
                ff.write("    %s\n" % notes.strip())
