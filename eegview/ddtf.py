import numpy as np
import signal_gen

import pygtk
pygtk.require('2.0')
import gtk
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas

class DDTF(gtk.Window):

    def ddtf(self,el1,el2,el3,sample_rate=400,duration=20,step=128,increment=5):



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
            z=el3[w:w+duration*sample_rate]
            # Initialize the Cross-Spectral arrays for averaging
            print "step first is : ", step
            Sxx=np.zeros((1,step - 1)); # - 1 here?
            print "Sxx: " , Sxx.shape
            Syy=Sxx
            Szz=Sxx
            Sxy=Sxx
            Sxz=Sxx
            Syz=Sxx
            Szy=Sxx
            print "xshape : ", x.shape
            print "Sxx shape : ", Sxx.shape
            xtemp=np.arange(0,step-1)
            xtemp_ones = np.ones(len(xtemp))
            print "xtempshape: ", xtemp.shape
            A = np.vstack([xtemp,xtemp_ones]).T
            print "A shape: ", A.shape
            inner_end_step = sample_rate*duration - step
            print "inner_end_step ", inner_end_step
            print "step ", step
            for i in np.arange(0,inner_end_step - 1,step):
                m,b = np.linalg.lstsq(A,x[i:i+step-1])[0] # the minus 1?
                print "m, b: ", m, b
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

                m,b = np.linalg.lstsq(A,z[i:i+step-1])[0] # the minus 1?
                trend = m*xtemp + b
                z[i:i+step-1] = z[i:i+step-1] - trend # detrend
                z[i:i+step-1] = z[i:i+step-1] - np.mean(z[i:i+step-1]) # demean
                fz = np.fft.fft(z[i:i+step-1] * np.hanning(step-1).T) # windowed fft

                # print "fs are ", fx, fy, fz
                # print "fxconf ", fx.conj()
                # print "Sxx ", Sxx.shape, Sxx.shape
                # print "fxstuff ", ((fx * fx.conj())).shape

                Sxx=Sxx+(fx * fx.conj())
                # print "Sxx2 ", Sxx.shape
                Syy=Syy+(fy * fy.conj())
                Szz=Szz+(fx * fz.conj())
                Sxy=Sxy+(fx * fx.conj())
                Sxz=Sxz+(fx * fy.conj())
                Syz=Syz+(fy * fy.conj())

                # print "Sxx shape: ", Sxx.shape
                # print "Sxy shape: ", Sxy.shape
                # print "Szy shape: ", Sxx.shape
                # print "Syz shape: ", Syz.shape

                Syx = Sxy.conj().T
                Szx = Sxz.conj().T
                Szy = Syz.conj().T

            S11=abs(Sxx)**2
            S12=abs(Sxy)**2
            S13=abs(Sxz)**2
            S21=abs(Syx)**2
            S22=abs(Syy)**2
            S23=abs(Syz)**2
            S31=abs(Szx)**2
            S32=abs(Szy)**2
            S33=abs(Szz)**2

            sumS = S11 + S12 + S13
            sumS2 = S21 + S22 + S23
            sumS3 = S31 + S32 + S33
            NS11 = S11 / S11.max()
            NS12 = S12 / sumS
            NS13 = S13 / sumS
            NS21 = S21 / sumS2
            NS22 = S22 / S22.max()
            NS23 = S23 / sumS2
            NS31 = S31 / sumS3
            NS32 = S32 / sumS3
            NS33 = S33 / S33.max()

            count += 1

            ttle1='Spectrum el1'
            ttle2=' el2 - . el1'
            ttle3=' el3 - . el1'
            ttle4=' el1 - . el2'
            ttle5=' Spectrum el2'
            ttle6=' el3 - . el2'
            ttle7=' el1 - . el3'
            ttle8='el2 - . el3'
            ttle9='Spectrum el3'

            # print "ns11 shape ", NS11.shape
            # print "f shape ", f.shape
            # print "f is: " , f
            # print "step is: ", step

            # print "shape x, y ", f[1:step/4].shape, NS11[1:step/4].shape
            # plot.subplot(211)

            # plot.axis([0, 60, 0, 1]) 



        # print (NS12, NS13, NS21, NS22, NS23, NS31, NS32, NS33)
        return (f ,step,NS11, NS12, NS13, NS21, NS22, NS23, NS31, NS32, NS33)

    def delete_event(self, widget, event, data=None):
        return False
    def destroy(self,widget, data=None):
        gtk.main_quit()
    
    def __init__(self):
        super(DDTF,self).__init__()
        
        self.connect("delete_event", self.delete_event)
        self.connect("destroy", self.destroy)
        e1,e2,e3 = signal_gen.signal_gen(.2,.01,.001)

        (f ,step,NS11, NS12, NS13, NS21, NS22, NS23, NS31, NS32, NS33) = self.ddtf(e1,e2,e3)
        
        # gtk.Window.__init__(self)
        self.fig = Figure(figsize = (20,15), dpi=72)
        
        self.canvas = FigureCanvas(self.fig)
        self.canvas.set_size_request(1800, 640)
        
        t = np.arange(0.0,50.0, 0.01)
        xlim = np.array([0,10])

        
        self.axes = self.fig.add_axes([0.075, 0.25, 0.9, 0.725], axisbg='#FFFFCC')

        self.axes.plot(t, np.sin(2*0.32*np.pi*t) * np.sin(2*2.44*np.pi*t) )
        self.axes.set_xlim([0.0,10.0])
        self.axes.set_xticklabels([])

        self.axesSpec = self.fig.add_axes([0.075, 0.05, 0.9, 0.2])
        t = self.axesSpec.text(
            0.5, 0.5,
            'Click on EEG channel for spectrogram (scroll mouse to expand)',
            verticalalignment='center',
            horizontalalignment='center',
            )
        t.set_transform(self.axes.transAxes)
        self.axesSpec.set_xlim([0.0,10.0])
        self.axesSpec.set_xticklabels([])
        self.axesSpec.set_yticklabels([])
        self.canvas.show()
        self.show()
        # self.axes.plot(f[step/4],NS11[:,0:step/4],'k')
        
        # plot.plot([1,2,3,4])
        # plot.show()
    def main(self):
        gtk.main()


if __name__=='__main__':
    ddtfmain = DDTF()
    ddtfmain.main()
    
