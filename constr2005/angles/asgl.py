from modeller import *
import os
import mdt

env = environ()
mlib = mdt.Library(env, '../lib/mdt2.bin')
mlib.angle_classes.read('${LIB}/anggrp.lib')

m = mdt.Table(mlib, file='mdt.mdt')
m = m.reshape(features=(35,111,112), offset=(0,0,0), shape=(1,-1,-1))

text = """
SET X_LABEL_STYLE = 2, X_TICK_LABEL = -999 -999
SET X_TICK = -999 -999 -999
SET TICK_FONT = 5, CAPTION_FONT = 5
SET Y_TICK = -999 -999 -999
SET WORLD_WINDOW = -999 -999 -999 -999
SET NO_XY_SCOLUMNS = 1 1, XY_SCOLUMNS = 2 1
FILL_COLUMN COLUMN = 2, COLUMN_PARAMETERS = 0. 0.25
SET BAR_XSHIFT = 0.125
ZOOM SCALE_WORLDX = 0.08
"""
m.write_asgl(asglroot='asgl1-a', plots_per_page=8, dimensions=1,
             plot_position=1, every_x_numbered=999, text=text, x_decimal=0)

os.system("asgl asgl1-a")
os.system("ps2pdf asgl1-a.ps")
