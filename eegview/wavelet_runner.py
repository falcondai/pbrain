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
import wavelet_creator
from scipy.stats.stats import pearsonr

class WaveletRunner(gtk.Window, Observer):
    def __init__(self,eoi,freq,t,data,trial_length,offset,window_length):
        gtk.Window.__init__(self)
        Observer.__init__(self)
        self.choose_file()
        self.channels = eoi
        self.selected_channels = []
        self.wavelets = {} # wavelets is a dict from name to np array
        self.eegfreq = freq
        self.t = t # an array of ms exactly indexed into self.data[0]
        self.trial_length = trial_length
        self.offset = offset
        print t
        print data.shape
        print self.channels
        self.data = data

        self.fig = Figure(figsize=(15,15), dpi=72)
        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.canvas.show()
        
        self.resize(700,512)
        self.set_title('Wavelet Runner')

        vbox = gtk.VBox()
        vbox.show()
        self.add(vbox)
        
        buttonChans = gtk.Button("Chans")
        buttonChans.show()
        buttonChans.connect('clicked', self.load_chans)
        
        buttonPlot = gtk.Button(stock=gtk.STOCK_EXECUTE)
        buttonPlot.show()
        buttonPlot.connect('clicked', self.execute)

        lwindow = gtk.Label()
        lwindow.set_text("window in ms")
        lwindow.show()
        self.window_length_entry = gtk.Entry()
        self.window_length_entry.set_text(str(window_length))
        self.window_length_entry.show()

        lmod = gtk.Label()
        lmod.set_text("modval")
        lmod.show()
        self.modval_entry = gtk.Entry()
        self.modval_entry.set_text('7.0')
        self.modval_entry.show()

        hbox = gtk.HBox()
        hbox.show()
        hbox.set_spacing(3)
        vbox.pack_start(hbox, False, False)
        hbox.pack_start(buttonChans, False, False)
        hbox.pack_start(buttonPlot, False, False)
        hbox.pack_start(lwindow, False, False)
        hbox.pack_start(self.window_length_entry, False, False)
        hbox.pack_start(lmod, False, False)
        hbox.pack_start(self.modval_entry, False, False)


        self.statBar = gtk.Label()
        self.statBar.set_alignment(0,0)
        self.statBar.show()
        self.progBar = gtk.ProgressBar()
        self.progBar.set_orientation(0)  # bottom-to-top
        self.progBar.set_fraction(0)
        self.progBar.show()
        
        vbox.pack_start(self.canvas, True, True)
        vbox.pack_start(self.statBar, False, False)
        vbox.pack_start(self.progBar, False, False)
    

    def choose_file(self):
        chooser = gtk.FileChooserDialog(title="please create dump file", action=gtk.FILE_CHOOSER_ACTION_SAVE, buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_SAVE,gtk.RESPONSE_OK))
        response = chooser.run()
        if response == gtk.RESPONSE_OK:
            filename = chooser.get_filename()
        else:
            chooser.destroy()
            return
        # try and write a dummy file to fname to make sure the dir
        # is writable
        tmpfile = filename + 'tmp'
        try: file(tmpfile, 'wb').write('123')
        except IOError:
            error_msg('Basepath %s does not appear to be writable' % filename,
                      parent=self)
            return
        else:
            os.remove(tmpfile)
        chooser.destroy()
        self.save_file = filename
        return


    def write_line(self,f, channel_name, wavelet_name, time_start, window_length, result):
        print >> f, "%s,%s,%f,%d,%f" %(channel_name,wavelet_name,time_start,window_length,result)

    def execute(self,button):
        file(self.save_file,'wb').write('')
        f = open(self.save_file,'ab')
        self.window_length = float(self.window_length_entry.get_text())
        self.modval = float(self.modval_entry.get_text())
        self.wavelets = wavelet_creator.create_all(self.window_length,self.eegfreq,self.modval)
        if self.selected_channels == []:
            channels = range(len(self.channels))
        else:
            channels = self.selected_channels
        print "CHANNELS: ", channels
        numRows,numCols = self.data.shape
        ind = range(int(self.offset), int(numRows - self.window_length +1), int(self.window_length))
        num_slices = len(ind)
        for entry in ind:
            for channel in channels:
                thisSlice = self.data[entry:(entry+self.window_length),channel]
                # print "CHANNEL: ", channel, "THISSLICE: \n", thisSlice,
                for wavelet in self.wavelets:
                    assert(len(self.wavelets[wavelet]) == len(thisSlice))
                    result,p = pearsonr(thisSlice,self.wavelets[wavelet])
                    self.write_line(f,self.channels[channel], wavelet, self.t[entry],self.window_length,result)
                    print result
        f.close

    def load_chans(self, button):
        dlg = gtk.Dialog("Channel Manipulation")
        dlg.connect("destroy", dlg.destroy)
        dlg.set_size_request(400,400)
        scrolled_window = gtk.ScrolledWindow(None, None)
        scrolled_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        dlg.vbox.pack_start(scrolled_window, True, True, 0)
        scrolled_window.show()

        table = gtk.Table(2,(1+len(self.channels)))
        table.set_row_spacings(8)
        table.set_col_spacings(8)
        scrolled_window.add_with_viewport(table)
        table.show()
        #attach format: obj, beg end x, beg end y
        l1 = gtk.Label("show            channel")
        l1.show()

        table.attach(l1,0,1,0,1)
        #an array to control the check boxes
        chanbuts = []
        for i in range(0, len(self.channels)):
            s1 = "                %s" % (self.channels[i],)
            chanbuts.append(gtk.CheckButton(s1))
            chanbuts[i].show()
            if i in self.selected_channels:
                chanbuts[i].set_active(True) #reactivate previously active channels
            chanbuts[i].connect("toggled", self.chanswitch, self.channels[i])
            table.attach(chanbuts[i], 0,1,i+1,i+2)
        butOK = gtk.Button("OK")    
        butOK.connect('clicked', (lambda b, x: x.destroy()), dlg)
        butOK.show()
        dlg.vbox.pack_start(butOK, False, False)
        dlg.show()

    def chanswitch(self, widget, channelnum):
        print channelnum
        if widget.get_active():
            self.selected_channels.append(self.channels.index(channelnum))
        else:
            self.selected_channels.remove(self.channels.index(channelnum))
