from modeller import *
import os
import mdt
import mdt.features

env = environ()
mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
restyp = mdt.features.ResidueType(mlib)
chi2 = mdt.features.Chi2Dihedral(mlib, bins=mdt.uniform_bins(144, -180, 2.5))

m = mdt.Table(mlib, file='mdt.mdt')

# remove the bins corresponding to undefined values for each of the 3 variables:
m = m.reshape(features=(xray, restyp, chi2), offset=(0,0,0), shape=(1,-2,-1))

# Let's get rid of the resolution variable from the output MDT table:
m = m.integrate(features=(restyp, chi2))

# Process the raw histograms to get appropriate pdf 1D splines for restraints:

# Start by smoothing with a uniform prior (equal weight when 10 points per bin),
# producing a normalized distribution that sums to 1 (not a pdf when dx != 1):
m = m.smooth(dimensions=1, weight=10)

# Normalize it to get the true pdf (Integral p(x) dx = 1):
# (the scaling actually does not matter, because I am eventually taking the
#  log and subtracting the smallest element of the final pdf, so this command
#  could be omitted without impact):
m = m.normalize(to_pdf=True, dimensions=1, dx_dy=2.5, to_zero=True)

# Take the logarithm of the smoothed frequencies
# (this is safe: none of bins is 0 because of mdt.smooth()):
m = m.log_transform(offset=0., multiplier=1.)

# Reverse the sign:
m = m.linear_transform(offset=0., multiplier=-1.)

# Offset the final distribution so that the lowest value is at 0:
m = m.offset_min(dimensions=1)

mdt.write_splinelib(file("chi2.py", "w"), m, "chi2", density_cutoff=0.1)

text = """
SET TICK_FONT = 5, CAPTION_FONT = 5
SET Y_TICK = -999 -999 -999
SET WORLD_WINDOW = -999 -999 -999 -999
"""
m.write_asgl(asglroot='modlib-a', plot_type='PLOT2D', every_x_numbered=20,
             text=text, dimensions=1, plot_position=1, plots_per_page=8)
os.system('asgl modlib-a')
