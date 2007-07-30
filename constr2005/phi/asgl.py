from modeller import *
import mdt

env = environ()
mlib = mdt.mdt_library(env, '../lib/mdt2.bin')

m = mdt.mdt(mlib, file='mdt.mdt')
m = m.reshape(features=(35,1,7), offset=(0,0,0), shape=(-1,-2,-1))

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

env.system("asgl asgl1-a")
env.system("ps2pdf asgl1-a.ps")
