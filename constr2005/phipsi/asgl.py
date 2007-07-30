from modeller import *
import mdt

env = environ()
mlib = mdt.mdt_library(env, '../lib/mdt2.bin')

m = mdt.mdt(mlib, file='mdt.mdt')
m = m.reshape(features=(35,1,9,7), offset=(0,0,0,0), shape=(1,-2,-1,-1))

text = """
SET TICK_FONT = 5, CAPTION_FONT = 5
SET WORLD_WINDOW = -999 -999 -999 -999
SET NO_XY_SCOLUMNS = 0 0, DPLOT_BOUNDS 0.0 -999
TRANSFORM TRF_TYPE=LOGARITHMIC4, TRF_PARAMETERS=10 1
"""
m.write_asgl(asglroot='asgl1-a', plots_per_page=3, dimensions=2,
             plot_position=9, every_x_numbered=12, every_y_numbered=12,
             text=text, x_decimal=0, y_decimal=0)

env.system("asgl asgl1-a")
env.system("ps2pdf asgl1-a.ps")
