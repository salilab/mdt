from modeller import *
import os
import mdt

env = environ()
mlib = mdt.Library(env, 'mdt_test_xray.bin')

m = mdt.Table(mlib, file='mdt2.mdt')

text = """
SET TICK_FONT = 5, CAPTION_FONT = 5
SET Y_TICK = -999 -999 -999, WORLD_WINDOW = -999 -999 -999 -999
"""
m.write_asgl(asglroot='asgl2-a', plots_per_page=8, dimensions=1,
             plot_position=1, every_x_numbered=2, text=text, x_decimal=1)

os.system("asgl asgl2-a")
os.system("ps2pdf asgl2-a.ps")
