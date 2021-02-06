from modeller import *
import os
import mdt.features

env = Environ()
mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
restyp = mdt.features.ResidueType(mlib)
chi3 = mdt.features.Chi3Dihedral(mlib, bins=mdt.uniform_bins(144, -180, 2.5))

m = mdt.Table(mlib, file='mdt.mdt')
m = m.reshape(features=(xray, restyp, chi3), offset=(0,0,0), shape=(1,-2,-1))

text = """
SET X_LABEL_STYLE = 2, X_TICK_LABEL = -999 -999
SET X_TICK = -999 -999 -999
SET TICK_FONT = 5, CAPTION_FONT = 5
SET Y_TICK = -999 -999 -999
SET WORLD_WINDOW = -999 -999 -999 -999
SET NO_XY_SCOLUMNS = 1 1, XY_SCOLUMNS = 2 1
FILL_COLUMN COLUMN = 2, COLUMN_PARAMETERS = -180. 2.5
SET BAR_XSHIFT = 1.25
"""
m.write_asgl(asglroot='asgl1-a', plots_per_page=8, dimensions=1,
             plot_position=1, every_x_numbered=20, text=text, x_decimal=0)

os.system("asgl asgl1-a")
os.system("ps2pdf asgl1-a.ps")
