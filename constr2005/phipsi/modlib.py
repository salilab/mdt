from modeller import *
import mdt
import mdt.features

env = environ()
mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
restyp = mdt.features.ResidueType(mlib)
psi = mdt.features.PsiDihedral(mlib, bins=mdt.uniform_bins(72, -180, 5.0))
phi = mdt.features.PhiDihedral(mlib, bins=mdt.uniform_bins(72, -180, 5.0))

m = mdt.Table(mlib, file='mdt.mdt')

# Eliminate the bins corresponding to undefined values:
m = m.reshape(features=(xray, restyp, psi, phi), offset=(0,0,0,0),
              shape=(1,-2,-1,-1))

# Let's get rid of the resolution variable from the output MDT table:
m = m.integrate(features=(restyp, psi, phi))

# Process the raw histograms to get appropriate pdf 1D splines for restraints:

# Start by smoothing with a uniform prior (equal weight when 10 points per bin),
# producing a normalized distribution that sums to 1 (not a pdf when dx != 1):
m = m.smooth(dimensions=2, weight=10)

# Normalize it to get the true pdf (Integral p(x) dx = 1):
# (the scaling actually does not matter, because I am eventually taking the
#  log and subtracting the smallest element of the final pdf, so this command
#  could be omitted without impact):
m = m.normalize(to_pdf=True, dimensions=2, dx_dy=(5., 5.), to_zero=True)

# Take the logarithm of the smoothed frequencies
# (this is safe: none of bins is 0 because of mdt.smooth()):
m = m.log_transform(offset=0., multiplier=1.)

# Reverse the sign:
m = m.linear_transform(offset=0., multiplier=-1.)

# Offset the final distribution so that the lowest value is at 0:
m = m.offset_min(dimensions=2)

mdt.write_2dsplinelib(file("phipsi.py", "w"), m, density_cutoff=0.1)
