from modeller import *
import os
import mdt
import mdt.features

env = Environ()
mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 1.4, "under 1.4"),
                                               (1.4,  1.6, "1.4-1.6"),
                                               (1.6,  1.8, "1.6-1.8"),
                                               (1.8,  2.001, "1.8-2.0")])
restyp = mdt.features.ResidueType(mlib)
chi1 = mdt.features.Chi1Dihedral(mlib, bins=mdt.uniform_bins(144, -180, 2.5))


m = mdt.Table(mlib, file='mdt.mdt')
# Remove undefined bins (and gap residue type)
m = m.reshape(features=(xray, restyp, chi1), offset=m.offset, shape=(0,-2,-1))

text = """
SET TICK_FONT = 5, CAPTION_FONT = 5
SET Y_TICK = -999 -999 -999, WORLD_WINDOW = -999 -999 -999 -999
"""
m.write_asgl(asglroot='asgl2-a', plots_per_page=8, dimensions=1,
             plot_position=1, every_x_numbered=20, text=text, x_decimal=1)

os.system("asgl asgl2-a")
os.system("ps2pdf asgl2-a.ps")
