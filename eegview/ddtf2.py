import numpy as np
import scikits.statsmodels.tsa.stattools as gtest
import math
import matplotlib.cbook as cbook
import matplotlib
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_agg import FigureCanvasAgg

from matplotlib.backends.backend_gtkagg import NavigationToolbar
from matplotlib.figure import Figure


import signal_gen
import pygtk
pygtk.require('2.0')
import gtk, gobject



class DDTF():

    def __init__(self):

        ttle1='Spectrum el1'
        ttle2=' el2 - . el1'
        ttle3=' el3 - . el1'
        ttle4=' el1 - . el2'
        ttle5=' Spectrum el2'
        ttle6=' el3 - . el2'
        ttle7=' el1 - . el3'
        ttle8='el2 - . el3'
        ttle9='Spectrum el3'


        self.win = gtk.Window()
        self.win.set_border_width(5)
        self.win.resize(800,400)
        vbox = gtk.VBox(spacing=3)
        self.win.add(vbox)
        vbox.show()
        self.fig = Figure(figsize=(7,5), dpi=72)

        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.canvas.show()
        vbox.pack_start(self.canvas, True, True)
        (self.e1,self.e2,self.e3) = signal_gen.signal_gen(3,.1,.001)

        self.ax1 = self.fig.add_subplot(431, title=ttle1)
        self.ax1.set_xlim(0,60)
        self.ax1.set_ylim(0,1)

        self.ax2 = self.fig.add_subplot(432, title=ttle2)
        self.ax2.set_xlim(0,60)
        self.ax2.set_ylim(0,1)

        self.ax3 = self.fig.add_subplot(433, title=ttle3)
        self.ax3.set_xlim(0,60)
        self.ax3.set_ylim(0,1)

        self.ax4 = self.fig.add_subplot(434, title=ttle4)
        self.ax4.set_xlim(0,60)
        self.ax4.set_ylim(0,1)


        self.ax5 = self.fig.add_subplot(435, title=ttle5)
        self.ax5.set_xlim(0,60)
        self.ax5.set_ylim(0,1)


        self.ax6 = self.fig.add_subplot(436, title=ttle6)
        self.ax6.set_xlim(0,60)
        self.ax6.set_ylim(0,1)

        self.ax7 = self.fig.add_subplot(437, title=ttle7)
        self.ax7.set_xlim(0,60)
        self.ax7.set_ylim(0,1)

        self.ax8 = self.fig.add_subplot(438, title=ttle8)
        self.ax8.set_xlim(0,60)
        self.ax8.set_ylim(0,1)

        self.ax9 = self.fig.add_subplot(439, title=ttle9)
        self.ax9.set_xlim(0,60)
        self.ax9.set_ylim(0,1)

        self.ax10 = self.fig.add_subplot(4,3,10, title="el1")
        self.ax11 = self.fig.add_subplot(4,3,11, title="el2")
        

    def ddtf(self,el1,el2,sample_rate=500,duration=20,step=128,increment=5):

        # self.ax10.plot(el1)
        # self.ax11.plot(el2)
        # self.ax12.plot(el3)

        # notes: duration is the length of a window in seconds
        # increment is the length of a step in seconds
        # step is the num points in an fft-analysis epoch
        N = len(el1)
        dt = 1/float(sample_rate)
        fNyq = sample_rate/2
        df = 1/(step*dt)
        f = np.arange(0,fNyq,df) #Frequency axis for the FFT

        count = 0
        end_step = N - duration*sample_rate
        print "end_step ", end_step
        print "stepping by ", increment * sample_rate
        for w in np.arange(0,end_step, increment * sample_rate):
            x=el1[w:w+duration*sample_rate] # should this be - 1 or 2?
            y=el2[w:w+duration*sample_rate]
            # z=el3[w:w+duration*sample_rate]
            # Initialize the Cross-Spectral arrays for averaging
            print "step first is : ", step
            Sxx=np.zeros((1,step - 1)); # - 1 here?
            print "Sxx: " , Sxx.shape
            Syy=Sxx
            # Szz=Sxx
            Sxy=Sxx
            # Sxz=Sxx
            # Syz=Sxx
            # Szy=Sxx
            # print "xshape : ", x.shape
            # print "Sxx shape : ", Sxx.shape
            xtemp=np.arange(0,step-1)
            xtemp_ones = np.ones(len(xtemp))
            # print "xtempshape: ", xtemp.shape
            A = np.vstack([xtemp,xtemp_ones]).T
            # print "A shape: ", A.shape
            inner_end_step = sample_rate*duration - step
            # print "inner_end_step ", inner_end_step
            # print "step ", step
            for i in np.arange(0,inner_end_step - 1,step):
                m,b = np.linalg.lstsq(A,x[i:i+step-1])[0] # the minus 1?
                print "TESTING LINALG SHAPE: ", A.shape, x[i:i+step-1].shape 
                # print "m, b: ", m, b
                trend = m*xtemp + b
                # print "istep : ", (i+step-1)
                x[i:i+step-1] = x[i:i+step-1] - trend # detrend
                x[i:i+step-1] = x[i:i+step-1] - np.mean(x[i:i+step-1]) # demean
                fx = np.fft.fft(x[i:i+step-1] * np.hanning(step-1).T) # windowed fft

                m,b = np.linalg.lstsq(A,y[i:i+step-1])[0] # the minus 1?
                trend = m*xtemp + b
                y[i:i+step-1] = y[i:i+step-1] - trend # detrend
                y[i:i+step-1] = y[i:i+step-1] - np.mean(y[i:i+step-1]) # demean
                fy = np.fft.fft(y[i:i+step-1] * np.hanning(step-1).T) # windowed fft

                # m,b = np.linalg.lstsq(A,z[i:i+step-1])[0] # the minus 1?
                # trend = m*xtemp + b
                # z[i:i+step-1] = z[i:i+step-1] - trend # detrend
                # z[i:i+step-1] = z[i:i+step-1] - np.mean(z[i:i+step-1]) # demean
                # fz = np.fft.fft(z[i:i+step-1] * np.hanning(step-1).T) # windowed fft

                # print "fs are ", fx, fy, fz
                # print "fxconf ", fx.conj()
                # print "Sxx ", Sxx.shape, Sxx.shape
                # print "fxstuff ", ((fx * fx.conj())).shape

                Sxx=Sxx+(fx * fx.conj())
                # print "Sxx2 ", Sxx.shape
                Syy=Syy+(fy * fy.conj())
                # Szz=Szz+(fz * fz.conj())
                Sxy=Sxy+(fx * fy.conj())
                # Sxz=Sxz+(fx * fz.conj())
                # Syz=Syz+(fy * fz.conj())

                # print "Sxx shape: ", Sxx.shape
                # print "Sxy shape: ", Sxy.shape
                # print "Szy shape: ", Sxx.shape
                # print "Syz shape: ", Syz.shape

                Syx = Sxy.conj()
                # Szx = Sxz.conj()
                # Szy = Syz.conj()

            S11=abs(Sxx)**2
            S12=abs(Sxy)**2
            # S13=abs(Sxz)**2
            S21=abs(Syx)**2
            S22=abs(Syy)**2
            # S23=abs(Syz)**2
            # S31=abs(Szx)**2
            # S32=abs(Szy)**2
            # S33=abs(Szz)**2

            sumS = S11 + S12 #  + S13
            sumS2 = S21 + S22 #  + S23
            # sumS3 = S31 + S32 + S33
            # NS11 = S11 / S11.max()
            NS12 = S12 / sumS
            # NS13 = S13 / sumS
            NS21 = S21 / sumS2
            # NS22 = S22 / S22.max()
            # NS23 = S23 / sumS2
            # NS31 = S31 / sumS3
            # NS32 = S32 / sumS3
            # NS33 = S33 / S33.max()

            count += 1
            
            # print "NS13: ", NS13
            # self.ax1.plot(f[0:step/4],NS11[0][0:step/4])
            # self.ax2.plot(f[0:step/4],NS12[0][0:step/4])
            # self.ax3.plot(f[0:step/4],NS13[0][0:step/4])
            # self.ax4.plot(f[0:step/4],NS21[0][0:step/4])
            # self.ax5.plot(f[0:step/4],NS22[0][0:step/4])
            # self.ax6.plot(f[0:step/4],NS23[0][0:step/4])
            # self.ax7.plot(f[0:step/4],NS31[0][0:step/4])
            # self.ax8.plot(f[0:step/4],NS32[0][0:step/4])
            # self.ax9.plot(f[0:step/4],NS33[0][0:step/4])






        return (f,NS12[0], NS21[0])


def start_ddtf():
    d = DDTF()
    f,a,b = d.ddtf(d.e1,d.e2)
    print a, b
    # d.win.show()


def test_do():
    e1,e2,e3 = signal_gen.signal_gen(3,.1,.001)
    do_ddtf_loop(e1,e2)

def do_ddtf_single_loop(el1,el2,sample_rate=500,duration=20):
    
    # notes: duration is the length of a window in seconds - note that this is a moving window
    # increment is the length of a step in seconds
    step = 32
    # should step just default to 64 or should it be increment in points??
    # step is the num points in an fft-analysis epoch
    N = len(el1)
    dt = 1/float(sample_rate)
    fNyq = sample_rate/2
    df = 1/(step*dt)
    f = np.arange(0,fNyq,df) #Frequency axis for the FFT
    # print "duration, sample_rate, step, increment ", duration, sample_rate,step,increment
    x=el1 #[w:w+new_duration*sample_rate] # should this be - 1 or 2?
    y=el2 #[w:w+new_duration*sample_rate]
    # Initialize the Cross-Spectral arrays for averaging
    Sxx=np.zeros((1,step)); # - 1 here?
    # print "Sxx: " , Sxx.shape
    Syy=Sxx
    # Szz=Sxx
    Sxy=Sxx
    # print "xshape : ", x.shape
    # print "Sxx shape : ", Sxx.shape
    xtemp=np.arange(0,step)
    xtemp_ones = np.ones(len(xtemp))
    # print "xtempshape: ", xtemp.shape
    A = np.vstack([xtemp,xtemp_ones]).T
    # print "A shape: ", A.shape, x[0:0+step-1].shape
    inner_end_step = sample_rate*duration - step
    # print "INDS: ", np.append(np.arange(0,inner_end_step,4),inner_end_step)
    for i in np.append(np.arange(0,inner_end_step,4),inner_end_step):
        # print "X SHAPE ", x[i:i+step-1].shape, i
        m,b = np.linalg.lstsq(A,x[i:i+step])[0] # the minus 1?
        # print "m, b: ", m, b
        trend = m*xtemp + b
        # print "istep : ", (i+step-1)
        x[i:i+step] = x[i:i+step] - trend # detrend
        x[i:i+step] = x[i:i+step] - np.mean(x[i:i+step]) # demean
        fx = np.fft.fft(x[i:i+step] * np.hanning(step).T) # windowed fft

        m,b = np.linalg.lstsq(A,y[i:i+step])[0] # the minus 1?
        trend = m*xtemp + b
        y[i:i+step] = y[i:i+step] - trend # detrend
        y[i:i+step] = y[i:i+step] - np.mean(y[i:i+step]) # demean
        fy = np.fft.fft(y[i:i+step] * np.hanning(step).T) # windowed fft

        Sxx=Sxx+(fx * fx.conj())
        Syy=Syy+(fy * fy.conj())
        Sxy=Sxy+(fx * fy.conj())
        Syx = Sxy.conj()

    S11=abs(Sxx)**2
    S12=abs(Sxy)**2
    S21=abs(Syx)**2
    S22=abs(Syy)**2

    sumS = S11 + S12 #  + S13
    sumS2 = S21 + S22 #  + S23

    NS12 = S12 / sumS
    if sumS2[0][0] == 0.0:
        NS12 = NS12 * 0
    NS21 = S21 / sumS2 # possibly do something with this?

    return f,NS12[0]



def do_ddtf_loop(el1,el2,sample_rate=500,duration=20):
    
    # notes: duration is the length of a window in seconds - note that this is a moving window
    # increment is the length of a step in seconds
    new_duration = duration/4
    increment = new_duration/4

    step = 32
    # should step just default to 64 or should it be increment in points??
    # step is the num points in an fft-analysis epoch
    N = len(el1)
    dt = 1/float(sample_rate)
    fNyq = sample_rate/2
    df = 1/(step*dt)
    f = np.arange(0,fNyq,df) #Frequency axis for the FFT
    # print "duration, sample_rate, step, increment ", duration, sample_rate,step,increment
    count = 0
    end_step = N - new_duration*sample_rate

    steps = np.arange(0,end_step,increment*sample_rate)
    total_steps = len(steps)
    # print "STEPS: ", steps, total_steps
    final = [] # np.zeros((total_steps, step - 1))
    # print "end_step ", end_step
    # print "stepping by ", increment * sample_rate
    for w in np.arange(0,end_step, increment * sample_rate):
        x=el1[w:w+new_duration*sample_rate] # should this be - 1 or 2?
        y=el2[w:w+new_duration*sample_rate]

        # z=el3[w:w+duration*sample_rate]
        # Initialize the Cross-Spectral arrays for averaging
        # print "step first is : ", step
        Sxx=np.zeros((1,step - 1)); # - 1 here?
        # print "Sxx: " , Sxx.shape
        Syy=Sxx
        # Szz=Sxx
        Sxy=Sxx
        # print "xshape : ", x.shape
        # print "Sxx shape : ", Sxx.shape
        xtemp=np.arange(0,step-1)
        xtemp_ones = np.ones(len(xtemp))
        # print "xtempshape: ", xtemp.shape
        A = np.vstack([xtemp,xtemp_ones]).T
        # print "A shape: ", A.shape, x[0:0+step-1].shape
        inner_end_step = sample_rate*new_duration - step
        # print "inner_end_step ", inner_end_step
        # print "step ", step
        for i in np.arange(0,inner_end_step - 1,step):
            # print "X SHAPE ", x[i:i+step-1].shape, i
            m,b = np.linalg.lstsq(A,x[i:i+step-1])[0] # the minus 1?
            # print "m, b: ", m, b
            trend = m*xtemp + b
            # print "istep : ", (i+step-1)
            x[i:i+step-1] = x[i:i+step-1] - trend # detrend
            x[i:i+step-1] = x[i:i+step-1] - np.mean(x[i:i+step-1]) # demean
            fx = np.fft.fft(x[i:i+step-1] * np.hanning(step-1).T) # windowed fft

            m,b = np.linalg.lstsq(A,y[i:i+step-1])[0] # the minus 1?
            trend = m*xtemp + b
            y[i:i+step-1] = y[i:i+step-1] - trend # detrend
            y[i:i+step-1] = y[i:i+step-1] - np.mean(y[i:i+step-1]) # demean
            fy = np.fft.fft(y[i:i+step-1] * np.hanning(step-1).T) # windowed fft

            # m,b = np.linalg.lstsq(A,z[i:i+step-1])[0] # the minus 1?
            # trend = m*xtemp + b
            # z[i:i+step-1] = z[i:i+step-1] - trend # detrend
            # z[i:i+step-1] = z[i:i+step-1] - np.mean(z[i:i+step-1]) # demean
            # fz = np.fft.fft(z[i:i+step-1] * np.hanning(step-1).T) # windowed fft

            # print "fs are ", fx, fy, fz
            # print "fxconf ", fx.conj()
            # print "Sxx ", Sxx.shape, Sxx.shape
            # print "fxstuff ", ((fx * fx.conj())).shape

            Sxx=Sxx+(fx * fx.conj())
            # print "Sxx2 ", Sxx.shape
            Syy=Syy+(fy * fy.conj())
            # Szz=Szz+(fz * fz.conj())
            Sxy=Sxy+(fx * fy.conj())
            # Sxz=Sxz+(fx * fz.conj())
            # Syz=Syz+(fy * fz.conj())

            # print "Sxx shape: ", Sxx.shape
            # print "Sxy shape: ", Sxy.shape
            # print "Szy shape: ", Sxx.shape
            # print "Syz shape: ", Syz.shape

            Syx = Sxy.conj()
            # Szx = Sxz.conj()
            # Szy = Syz.conj()

        S11=abs(Sxx)**2
        S12=abs(Sxy)**2
        # S13=abs(Sxz)**2
        S21=abs(Syx)**2
        S22=abs(Syy)**2
        # S23=abs(Syz)**2
        # S31=abs(Szx)**2
        # S32=abs(Szy)**2
        # S33=abs(Szz)**2

        sumS = S11 + S12 #  + S13
        sumS2 = S21 + S22 #  + S23
        # sumS3 = S31 + S32 + S33
        # NS11 = S11 / S11.max()
        NS12 = S12 / sumS
        # NS13 = S13 / sumS
        NS21 = S21 / sumS2
        # NS22 = S22 / S22.max()
        # NS23 = S23 / sumS2
        # NS31 = S31 / sumS3
        # NS32 = S32 / sumS3
        # NS33 = S33 / S33.max()
        # print count
        # print "finalshape: ", final.shape, NS12[0].shape, final[count].shape, count
        final.append(NS12[0]) # [count] = NS12[0]
        # count += 1
        
    # print "finalshape: ", final.shape
    # final = np.mean(final, axis=0)
    # print "finalshape: ", final.shape

    return f,final


