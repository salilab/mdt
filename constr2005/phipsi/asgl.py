from modeller import *
import os
import mdt
import mdt.features

env = Environ()
mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
restyp = mdt.features.ResidueType(mlib)
psi = mdt.features.PsiDihedral(mlib, bins=mdt.uniform_bins(72, -180, 5.0))
phi = mdt.features.PhiDihedral(mlib, bins=mdt.uniform_bins(72, -180, 5.0))

m = mdt.Table(mlib, file='mdt.mdt')
m = m.reshape(features=(xray, restyp, psi, phi),
              offset=(0,0,0,0), shape=(1,-2,-1,-1))

text = """
SET TICK_FONT = 5, CAPTION_FONT = 5
SET WORLD_WINDOW = -999 -999 -999 -999
SET NO_XY_SCOLUMNS = 0 0, DPLOT_BOUNDS 0.0 -999
TRANSFORM TRF_TYPE=LOGARITHMIC4, TRF_PARAMETERS=10 1
"""
m.write_asgl(asglroot='asgl1-a', plots_per_page=3, dimensions=2,
             plot_position=9, every_x_numbered=12, every_y_numbered=12,
             text=text, x_decimal=0, y_decimal=0)

os.system("asgl asgl1-a")
os.system("ps2pdf asgl1-a.ps")
