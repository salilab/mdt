from modeller import *
import mdt

env = environ()
mlib = mdt.mdt_library(env, '../lib/mdt2.bin')

m = mdt.mdt(mlib, file='mdt.mdt')

# remove the bins corresponding to undefined values for each of the 3 variables:
m = m.reshape(features=(35,1,9), offset=(0,0,0), shape=(-1,-2,-1))

# Let's get rid of the resolution variable (feature 35) from the output
# MDT table:
m = m.integrate(features=(1,9))

# Process the raw histograms to get appropriate pdf 1D splines for restraints:

# Start by smoothing with a uniform prior (equal weight when 10 points per bin),
# producing a normalized distribution that sums to 1 (not a pdf when dx <> 1):
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

mdt.write_splinelib(file("psi.py", "w"), m, "psi", density_cutoff=0.1)

text = """
SET TICK_FONT = 5, CAPTION_FONT = 5
SET Y_TICK = -999 -999 -999
SET WORLD_WINDOW = -999 -999 -999 -999
"""
m.write_asgl(asglroot='modlib-a', plot_type='PLOT2D', every_x_numbered=20,
             text=text)
env.system('asgl modlib-a')
