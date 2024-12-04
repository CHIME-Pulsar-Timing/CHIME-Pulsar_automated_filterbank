from builtins import map
import re
import glob

# import presto.sifting as sifting
import sifting_DMprec as sifting
from operator import itemgetter, attrgetter
import argparse
import yaml

parser = argparse.ArgumentParser()
parser.add_argument(
    "--zmax",
    help="zmax used in search, used in the glob to find candidate files",
    type=str,
    default="0",
)
parser.add_argument(
    "--basename",
    help="Base filename to look for in glob, e.g. J2119+49_59487_pow",
    type=str,
    default="*",
)
parser.add_argument(
    "--min_num_DMs",
    help="In how many DMs must a candidate be detected to be considered 'good'",
    type=int,
    default=2,
)
parser.add_argument(
    "--low_DM_cutoff",
    help="Lowest DM to consider as a 'real' pulsar",
    type=float,
    default=2.0,
)
parser.add_argument(
    "--r_err",
    help="How close a candidate has to be to another candidate to consider it the same candidate (in Fourier bins)",
    type=float,
    default=1.1,
)
parser.add_argument(
    "--sigma_threshold",
    help="Ignore candidates with a sigma (from incoherent power summation) less than this",
    type=float,
    default=4.0,
)
parser.add_argument(
    "--c_pow_threshold",
    help="Ignore candidates with a coherent power less than this",
    type=float,
    default=100.0,
)
parser.add_argument(
    "--short_period",
    help="Shortest period candidates to consider (s)",
    type=float,
    default=0.0005,
)
parser.add_argument(
    "--long_period",
    help="Longest period candidates to consider (s)",
    type=float,
    default=15.0,
)
parser.add_argument(
    "--harm_pow_cutoff",
    help="Ignore any candidates where at least one harmonic does exceed this power",
    type=float,
    default=8.0,
)
parser.add_argument(
    "--dm_iqr_limit",
    help="Ignore any candidate where the DMs of its hits have an interquartile range above this value",
    type=float,
    default=20,
)
parser.add_argument(
    "--dm_precision",
    help="Precision in DM (aka number of decimal places)",
    default=2,
    type=int,
)
parser.add_argument(
    "--infyaml",
    help="yaml file containing a list of all of the inf files used in the search",
    type=str,
)
parser.add_argument(
    "--outfile", help="Where to write sifted candidates", default=None, type=str
)
args = parser.parse_args()

# glob for ACCEL files
globaccel = args.basename + "_DM*ACCEL_" + args.zmax
# glob for .inf files
globinf = args.basename + "_DM*.inf"


# set DM precision
prec = args.dm_precision
sifting.dmprec = args.dm_precision
sifting.DM_re = sifting.refresh_DM_re()
# set other sifting parameters
min_num_DMs = args.min_num_DMs
low_DM_cutoff = args.low_DM_cutoff
sifting.sigma_threshold = args.sigma_threshold
sifting.c_pow_threshold = args.c_pow_threshold
sifting.r_err = args.r_err
sifting.short_period = args.short_period
sifting.long_period = args.long_period
sifting.harm_pow_cutoff = args.harm_pow_cutoff
sifting.dm_iqr_limit = args.dm_iqr_limit

# If the birds file works well, the following shouldn't
# be needed at all...  If they are, add tuples with the bad
# values and their errors.
#                (ms, err)
sifting.known_birds_p = []
#                (Hz, err)
sifting.known_birds_f = []

# --------------------------------------------------------------

# Try to read the .inf files first, as _if_ they are present, all of
# them should be there.  (if no candidates are found by accelsearch
# we get no ACCEL files...
# I have no inffiles . . .

# inffiles = glob.glob(globinf)
if args.infyaml:
    with open(args.infyaml, "r") as f:
        inffiles = yaml.safe_load(f)
else:
    inffiles = []
candfiles = glob.glob(globaccel)
# Check to see if this is from a short search
# if len(re.findall("_[0-9][0-9][0-9]M_" , inffiles[0])):
if inffiles:
    dmstrs = [x.split("DM")[-1].split("_")[0] for x in inffiles]
else:
    dmstrs = [x.split("DM")[-1].split("_")[0] for x in candfiles]

dms = list(map(float, dmstrs))
dms.sort()
dmstrs = [f"{x:.{prec}f}" for x in dms]

# Read in all the candidates
cands = sifting.read_candidates(candfiles)

# Remove candidates that are duplicated in other ACCEL files
if len(cands):
    cands = sifting.remove_duplicate_candidates(cands)

# Remove candidates with DM problems
if len(cands):
    cands = sifting.remove_DM_problems(cands, min_num_DMs, dmstrs, low_DM_cutoff)

# Remove candidates that are harmonically related to each other
# Note:  this includes only a small set of harmonics
if len(cands):
    cands = sifting.remove_harmonics(cands)

# Write candidates to STDOUT
if len(cands):
    cands.sort(key=attrgetter("sigma"), reverse=True)
    sifting.write_candlist(cands, candfilenm=args.outfile)
