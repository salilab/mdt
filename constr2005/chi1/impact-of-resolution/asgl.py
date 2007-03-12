from modeller import *
import mdt

env = environ()
mlib = mdt.mdt_library(env, '${LIB}/mdt.ini', '../../lib/mdt1.bin')

m = mdt.mdt(mlib, file='mdt.mdt')

text = """
SET TICK_FONT = 5, CAPTION_FONT = 5
SET Y_TICK = -999 -999 -999, WORLD_WINDOW = -999 -999 -999 -999
"""
m.write_asgl(asglroot='asgl2-a', plots_per_page=8, dimensions=1,
             plot_position=1, every_x_numbered=20, text=text, x_decimal=1)

env.system("asgl asgl2-a")
env.system("ps2pdf asgl2-a.ps")
