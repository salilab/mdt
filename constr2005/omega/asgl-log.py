from modeller import *
import os
import mdt
import mdt.features

env = Environ()
mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
restyp_1 = mdt.features.ResidueType(mlib, delta=1)
omega = mdt.features.OmegaDihedral(mlib, bins=mdt.uniform_bins(720, -180, 0.5))

m = mdt.Table(mlib, file='mdt.mdt')
m = m.reshape(features=(xray, restyp_1, omega), offset=(0,0,0), shape=(1,-2,-1))

text = """
SET X_LABEL_STYLE = 2, X_TICK_LABEL = -999 -999
SET X_TICK = -999 -999 -999
SET TICK_FONT = 5, CAPTION_FONT = 5
SET Y_TICK = -999 -999 -999
SET WORLD_WINDOW = -999 -999 -999 -999
SET NO_XY_SCOLUMNS = 1 1, XY_SCOLUMNS = 2 1
FILL_COLUMN COLUMN = 2, COLUMN_PARAMETERS = -180. 0.5
SET BAR_XSHIFT = 0.25
TRANSFORM TRF_TYPE = LOGARITHMIC4, ;
          TRF_PARAMETERS = 1 1, NO_XY_SCOLUMNS = 0 1, XY_SCOLUMNS = 1
"""
m.write_asgl(asglroot='asgl2-a', plots_per_page=8, dimensions=1,
             plot_position=1, every_x_numbered=20, text=text, x_decimal=0)

os.system("asgl asgl2-a")
os.system("ps2pdf asgl2-a.ps")
