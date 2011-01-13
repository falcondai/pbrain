from __future__ import division
import sys, os, math, copy, gc
import pylab as p
import vtk

from scipy import array, zeros, ones, sort, absolute, sqrt, divide,\
     argsort, take, arange
import numpy
import pygtk
pygtk.require('2.0')
import gtk
import gobject
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.widgets import Cursor
from events import Observer
from shared import fmanager
from pbrainlib.gtkutils import error_msg, simple_msg, make_option_menu,\
     get_num_value, get_num_range, get_two_nums, str2int_or_err,\
     OpenSaveSaveAsHBox, ButtonAltLabel, str2num_or_err

class CohExplorer(gtk.Window, Observer):
    """this class will implement an entirely new view, launching from view3, allowing exploration and data discovery based on coherence dumps from the dump coherences button in the view3 window. """
    def __init__(self, eoi):
        gtk.Window.__init__(self)
        Observer.__init__(self)
        self.bandlist = ['delta', 'theta', 'alpha', 'beta', 'gamma', 'high']
        self.t_data = {} #this is the dictionary used in read_data, from time points to channel cohs unaveraged
        self.x_data = {} #this is the dictionary used in make_fig, from channel labels to xdata
        self.channels = eoi
        self.fig = Figure(figsize=(15,15), dpi=72)
        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.canvas.show()
        self.canvas.mpl_connect('motion_notify_event', self.motion_notify_event)
        self.canvas.mpl_connect('button_press_event', self.button_press_event)
        self.canvas.mpl_connect('button_release_event', self.button_release_event)   
        print "explorer channels are: ", self.channels
        
        self.resize(700,512)
        self.set_title('Coherence Explorer')
        
        self.band = self.bandlist[0] #start with delta
        self.length = 10 #start with a tiny amount of lines
        vbox = gtk.VBox()
        vbox.show()
        self.add(vbox)
        saveEntry = gtk.Entry()
        self.dumpfile = None
        self.oldlength = 10
        self.opt = 'multi'
        self.optchanged = 0
        
        def plot(button):
            self.length = int(lengthEntry.get_text())
            if (self.length > self.oldlength or self.optchanged == 1):
                del self.t_data
                self.t_data = {}
                self.read_data(self.dumpfile, self.length, self.opt)
            del self.x_data
            self.x_data = {}
            self.make_fig(self.band)
            self.canvas.draw()
            self.oldlength = copy.deepcopy(self.length)
            self.optchanged = 0
        
        def load_file(button):
            dumpfile = fmanager.get_filename(title="Select dump file:")
            if dumpfile is None: return
            if not os.path.exists(dumpfile):
                print 'File %s does not exist' % dumpfile
                return

            try: fh = file(dumpfile)
            except IOError, msg:
                msg = exception_to_str('Could not open %s' % dumpfile)
                error_msg(msg)
                return
            fh.close() 
            self.dumpfile = dumpfile
            saveEntry.set_text(str(self.dumpfile))
            self.read_data(self.dumpfile, self.length, self.opt) #the muscle
            return
            
    
        buttonSave = gtk.Button(stock=gtk.STOCK_OPEN)
        buttonSave.show()
        buttonSave.connect('clicked', load_file)
        saveEntry.show()
        
        def set_active_band(combobox):
            model = combobox.get_model()
            index = combobox.get_active()
            label = model[index][0]
            self.band = label
            return
        def set_opts(combobox):
            self.optchanged = 1
            model = combobox.get_model()
            index = combobox.get_active()
            label = model[index][0]
            self.opt = label
            
        bandMenu = make_option_menu(self.bandlist, func=set_active_band)
        optMenu = make_option_menu(['multi', 'forward', 'backward', 'multiphase', 'forwardphase'], func=set_opts)
        
        lengthEntry = gtk.Entry()
        lengthEntry.set_text(str(self.length))
        lengthEntry.show()
        buttonPlot = gtk.Button(stock=gtk.STOCK_EXECUTE) #the execute button replots with changes
        buttonPlot.show()
        buttonPlot.connect('clicked', plot)
        
        hbox = gtk.HBox()
        hbox.show()
        hbox.set_spacing(3)
        vbox.pack_start(hbox, False, False)
        hbox.pack_start(saveEntry, False, False)
        hbox.pack_start(buttonSave, False, False)
        hbox.pack_start(bandMenu, False, False)
        hbox.pack_start(lengthEntry, False, False)
        hbox.pack_start(optMenu, False, False)
        hbox.pack_start(buttonPlot, False, False)
        
        self.progBar = gtk.ProgressBar()
        #self.progBar.set_size_request(10, 100)
        self.progBar.set_orientation(0)  # bottom-to-top
        self.progBar.set_fraction(0)
        self.progBar.show()
        load_file(buttonSave) #load a file on startup
        self.make_fig(self.band) #and plot it
        
        vbox.pack_start(self.canvas, True, True)
        vbox.pack_start(self.progBar, False, False)
   
   
    
    def motion_notify_event(self, event):
        return False
    def button_press_event(self, event):
        return False
    def button_release_event(self, event):
        return False
            
    def read_data(self, df, length, opt):
    
        def progress_callback(frac,  msg):
                if frac<0 or frac>1: return
                self.progBar.set_fraction(frac)
                while gtk.events_pending(): gtk.main_iteration()
        #file looks like:
        #'E1,E2,delta 1-4,theta 4-8,alpha 8-12,beta 12-30,gamma 30-50,high gamma 70-100,delta phase,theta phase,alpha phase,beta phase,gamma phase,high gamma phase'    
        #the spec is a dict from name of trode to list containing numsamples (for each band!) and sumsamples
        #and a dict from time to above dict
        f = open(df, 'rb')
        nl = length
        i = 0
        trodelist = []
        tr = {}
        skip = 0
        tnext = 0
        cc = 0
        while 1:
            if (i != cc and i%10 == 0):
                progress_callback(i/nl, "almost done..")
            line = f.readline()
            if (i == nl):
                break
            if (skip == 1):
                skip = 0
                continue
            if (tnext == 1):
                tnext = 0
                tcur = int(line)
                #print "current time is ", tcur
                tr = copy.deepcopy(tr) #trodes is a dict from name of trode to tuple containing numsamples and sumsamples. copies are important.
                tr = {} #get ready
                trodelist.append(tr) #add  tr to the list of data
                
                #print "appended new list at index ", i
                
                self.t_data[tcur] = trodelist[i] #add the time to the dict, get ready to add the trodes data at that time
                i += 1
                continue
            cc = copy.deepcopy(i)    
            if (line[0] == '['):
                if (line[1:5] == 'swee'):
                    skip = 1
                    continue
                elif (line[1:5] == 'leng'):
                    skip = 1
                    continue
                elif (line[1:5] == 'offs'):
                    tnext = 1
                    continue
            line = line.split(',')
            
            if (opt == 'multi'):
                zstart = 0
                zend = 1
                cstart = 2
                cend = 8
            if (opt == 'forward'):
                zstart = 0
                zend = 0
                cstart = 2
                cend = 8
            if (opt == 'backward'):  
                zstart = 1
                zend = 1
                cstart = 2
                cend = 8
            if (opt == 'multiphase'):  
                zstart = 0
                zend = 1
                cstart = 9
                cend = 14
            if (opt == 'forwardphase'):
                zstart = 0
                zend = 0
                cstart = 9
                cend = 14
            
            for col in arange(cstart,cend): #change nan to 0
                if (math.isnan(float(line[col]))):
                    line[col] = 0 
            for z in (zstart,zend):
                if (line[z] in tr):
                    item = tr[line[z]]
                    item[0] += 1
                    item[1] = [ab+ac for ab,ac in zip(item[1], map(float,line[cstart:cend]))] #add the new data, for all bands
                else:
                    tr[line[z]] = [1, map(float,line[cstart:cend])]
            
        f.close()
        
    def make_fig(self, band):
        #this function modifies self.fig
        self.fig.clear()
        self.lines = []
        col = self.bandlist.index(band)#find the band's index for the t_data
        N = len(self.channels)
        self.ax = self.fig.add_subplot(1,1,1) #new singleplot configuration
        self.cursor = Cursor(self.ax, useblit=True, linewidth=1, color='white')
        self.cursor.horizOn = True
        self.cursor.vertOn = True
        colordict = ['#A6D1FF','green','red','cyan','magenta','yellow','white']
        keys = self.t_data.keys()
        keys = sorted(keys)
        del keys[-1] #for some reason the last time key is coming out null sometimes in the data?
        dumpchans = self.t_data[0].keys()
        #print "dumpchans!!", dumpchans
        #print "EOI: ", self.channels
        
        counter = 0
        xdata = []
        for i in self.t_data[0]: #go through the ordered channel list one by one
            xdata.append([]) #keep all the channel data seperate
            #print "this is the channelkey: ", i
            isplit = i.split(" ")
            isplit = isplit[0],int(isplit[1])
            if not isplit in self.channels: #make sure the channel we plot is in the view3 display
                print "channel error! on channel ", i
            color = colordict[((counter)%7)]
            for t in keys[0:self.length-2]: #at each time point. we don't want to use more of t_data than is asked for.
                (ns,ss) = self.t_data[t][i] #get the data out
                ssband = ss[col] #choose the band
                x1 = ssband/ns #find the average
                xdata[counter].append(x1)
            self.x_data[isplit] = xdata[counter] #save the channel's data into an organized dict of channels
            newp = self.ax.plot(keys[0:self.length-2], xdata[counter], color, label=(str(i))) #plot a channel """arange(len(xdata[counter]))"""
            #print "at time ", counter, "xdata has ", len(xdata[counter]), " channels."
            counter += 1
        self.ax.patch.set_facecolor('black') #black bg
        self.progBar.set_fraction(0)
        return
