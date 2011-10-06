from __future__ import division
import sys, os, math, copy, gc
import vtk
import numpy
import pygtk
pygtk.require('2.0')
import gtk
import gobject
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import axes
from matplotlib.widgets import Cursor
import mpl_toolkits.mplot3d.axes3d as p3
from events import Observer
from pbrainlib.gtkutils import error_msg, simple_msg, make_option_menu,\
     get_num_value, get_num_range, get_two_nums, str2int_or_err,\
     OpenSaveSaveAsHBox, ButtonAltLabel, str2num_or_err, exception_to_str

class WaveletRunner(gtk.Window, Observer):
    def __init__(self,eoi,freq,t,data):
        gtk.Window.__init__(self)
        Observer.__init__(self)
        self.eoi = eoi
        self.eegfreq = freq
        self.t = t
        self.data = data

        self.fig = Figure(figsize=(15,15), dpi=72)
        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.canvas.show()
        
        self.resize(700,512)
        self.set_title('Wavelet Runner')
        
