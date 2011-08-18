import numpy as np
import scikits.statsmodels.tsa.stattools as gtest
import math
import matplotlib.cbook as cbook
import matplotlib
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_agg import FigureCanvasAgg

from matplotlib.backends.backend_gtkagg import NavigationToolbar
from matplotlib.figure import Figure
import ddtf2


import signal_gen
import pygtk
pygtk.require('2.0')
import gtk, gobject



def window_hanning(x):
    "return x times the hanning window of len(x)"
    return np.hanning(len(x))*x


def start_test(maxlag, col1,col2, start, stop):

    Xnew = []
    Ynew = []

    example = open('/home/lelz/towlelab/examples/test.rest.neuroscanascii', 'rb')
    counter = 0
    for l in example:
        if counter < start:
            counter += 1
            continue
        counter += 1
        if counter > stop:
            break
        l = l.strip().split()
    #print "t2", len(l[0])
        if len(l) > 0:
    
            if len(l[0]) > 0:
            
                if l[0][0] == '[':
                    continue
        
            # print "appending to x: ", l[0]
            Xnew.append(float(l[col1]))
            # print "appending to y: ", l[1]
            Ynew.append(float(l[col2]))

    Xnew = np.array(Xnew)
    Ynew = np.array(Ynew)
    
    d = np.vstack((Xnew,Ynew)).T
    res = gtest.grangercausalitytests(d, maxlag, verbose=False)
    return res

def start_test_signal_gen(lag,SNR, K1,K2,col1,col2):
    data = signal_gen.signal_gen(SNR, K1,K2)
    Xnew = data[col1]
    Ynew = data[col2]
    d = np.vstack((Xnew,Ynew)).T
    res = gtest.grangercausalitytests(d, lag, verbose=True)
    return res

def donothing_callback(*args):
    pass

def detrend_none(x):
    "Return x: no detrending"
    return x


def old_granger_test(X, ij, newLength=256, NFFT=256, offset=0, Fs=2, maxlag=3, progressCallback=donothing_callback, window=window_hanning, noverlap=0, detrend = detrend_none, gv1=0,gv2=0):
    threshold = .05/(len(ij)*2)
    oldNFFT = NFFT
    NFFT = newLength
    numRows, numCols = X.shape

    if numRows < NFFT:
        tmp = X
        X = np.zeros( (NFFT, numCols), X.dtype)
        X[:numRows,:] = tmp
        del tmp

    numRows, numCols = X.shape
    # get all the columns of X that we are interested in by checking
    # the ij tuples
    allColumns = set()
    for i,j in ij:
        allColumns.add(i); allColumns.add(j)
    Ncols = len(allColumns)

    # for real X, ignore the negative frequencies
    if np.iscomplexobj(X): numFreqs = NFFT
    else: numFreqs = NFFT//2+1

    if cbook.iterable(window):
        assert(len(window) == NFFT)
        windowVals = window
    else:
        windowVals = window_hanning(np.ones(NFFT, X.dtype)) #I changed this from window to window_hanning. window was not doing anything!! -eli

    ind = range(offset, int(numRows-newLength+1), int(oldNFFT-noverlap)) #coherence calcs on each sweep start at offset
    numSlices = len(ind)
    FFTSlices = {}
    Pxx = {}
    if newLength > oldNFFT: #make sure that the newlength is shorter than the segment (epoch) length (see below)
        newLength = oldNFFT
    slices = range(numSlices)
    # normVal = np.linalg.norm(windowVals)**2

    for iCol in allColumns:
        progressCallback(i/Ncols, 'Cacheing FFTs')
        Slices = np.zeros( (numSlices, newLength))
        for iSlice in slices:
            #thisSlice = X[ind[iSlice]:ind[iSlice]+NFFT, iCol] #this is the line that reads the data normally
            thisSlice = X[ind[iSlice]:ind[iSlice]+newLength, iCol] #this is the line that reads sections of epochs
            print "GRANGER TESTING: ", ind[iSlice], " to ", ind[iSlice] + newLength
            print "shape of all: ", Slices.shape
            print "shape of slice: ", thisSlice.shape
            #thisSlice = windowVals*detrend(thisSlice)
            Slices[iSlice] = thisSlice # = np.fft.fft(thisSlice)[:numFreqs]
        #FFTSlices[iCol] = Slices
        Pxx[iCol] = np.mean(Slices,axis=0) # / normVal
        print "shape of pxx one col: ", Pxx[iCol].shape
        # print Pxx[iCol]
    del Slices, ind, windowVals

    Cxy = {}
    Phase = {}
    count = 0

    N = len(ij)

    typedict = ['params_ftest', 'ssr_chi2test']
    for i,j in ij:
        count += 1
        if count%10==0:
            progressCallback(count/N, 'Computing coherences')
        
        d = np.vstack((Pxx[i],Pxx[j])).T
        drev = np.vstack((Pxx[j],Pxx[i])).T
        res = gtest.grangercausalitytests(d,maxlag,verbose=False)
        resrev = gtest.grangercausalitytests(drev,maxlag,verbose=False)
        # print "looking at: ", res[maxlag][0]

        # this gets the p value
        result = res[maxlag][0][typedict[gv1]][not gv2]
        rev_result = resrev[maxlag][0][typedict[gv1]][not gv2]
        if result <= threshold and rev_result > threshold:
            Cxy[i,j] = 1 - (result/threshold)
            Phase[i,j] = 90
        elif result > threshold and rev_result <= threshold:
            Cxy[i,j] = 1 - (rev_result/threshold)    
            Phase[i,j] = -90
        elif result <= threshold and rev_result <= threshold:
            Cxy[i,j] = min((result, rev_result))
            direction = (result,rev_result).index(Cxy[i,j])
            if direction == 0:
                Phase[i,j] = 10
            else:
                Phase[i,j] = -10
        else: 
            Cxy[i,j] = 0
            Phase[i,j] = 0

    
        # print typedict[gv1],gv2
        print "RESIS: ", Cxy[i,j], "RESULT ", result, "REVRESULT ", rev_result

        
        
    print "NFFT, NUMFREQS: ", NFFT, numFreqs
    freqs = Fs/NFFT*np.arange(numFreqs)
    print "FREQS ARE: ", freqs
    return Cxy, Phase, freqs


def ddtf_test(X, ij, newLength=256, NFFT=256, offset=0, Fs=2, progressCallback=donothing_callback, window=window_hanning, noverlap=0, detrend = detrend_none):
    # note that Fs is the frequency of the eeg spectrum, and should never actually be 2
    threshold = .05/(len(ij)*2)
    oldNFFT = NFFT
    NFFT = newLength
    numRows, numCols = X.shape
    print "DDTF: oldNFFT, NFFT, newLength, xshape ", oldNFFT,NFFT,newLength,X.shape
    duration = (NFFT/Fs) # duration in seconds of each slice
    print "DDTF: trial length is ", duration, " Fs is ", Fs

    if numRows < NFFT:
        tmp = X
        X = np.zeros( (NFFT, numCols), X.dtype)
        X[:numRows,:] = tmp
        del tmp

    numRows, numCols = X.shape
    # get all the columns of X that we are interested in by checking
    # the ij tuples
    allColumns = set()
    for i,j in ij:
        allColumns.add(i); allColumns.add(j)
    Ncols = len(allColumns)

    # for real X, ignore the negative frequencies
    if np.iscomplexobj(X): numFreqs = NFFT
    else: numFreqs = NFFT//2+1

    if cbook.iterable(window):
        assert(len(window) == NFFT)
        windowVals = window
    else:
        windowVals = window_hanning(np.ones(NFFT, X.dtype)) #I changed this from window to window_hanning. window was not doing anything!! -eli

    ind = range(offset, int(numRows-newLength+1), int(oldNFFT-noverlap)) #coherence calcs on each sweep start at offset
    numSlices = len(ind)
    FFTSlices = {}
    Pxx = {}
    slices = range(numSlices)
    # normVal = np.linalg.norm(windowVals)**2

    for iCol in allColumns:
        progressCallback(i/Ncols, 'Cacheing FFTs')
        Slices = np.zeros( (numSlices, newLength))
        for iSlice in slices:
            thisSlice = X[ind[iSlice]:ind[iSlice]+newLength, iCol] #this is the line that reads sections of epochs
            Slices[iSlice] = thisSlice # = np.fft.fft(thisSlice)[:numFreqs]
        Pxx[iCol] = np.mean(Slices,axis=0) # / normVal
        # print "shape of pxx one col: ", Pxx[iCol].shape
        
    del Slices, ind, windowVals

    Cxy = {}
    Phase = {}
    count = 0

    N = len(ij)
    All_final = {}
    All = {}
    iters = 4
    for i,j in ij:
        count += 1
        if count%10==0:
            progressCallback(count/N, 'Computing coherences')
        
        f, final = ddtf2.do_ddtf_loop(Pxx[i],Pxx[j],sample_rate=Fs,duration=duration)
        # Cxy[i,j] = result # max((result, rev_result))
        # Phase[i,j] = result
        counter = 0
        for entry in final:
            try:
                All_final[counter][i,j] = entry
            except:
                All_final[counter] = {}
                All_final[counter][i,j] = entry
            counter += 1

        print "FREQUENCIES!!! ", f
    
                
    print "NFFT, NUMFREQS: ", NFFT, numFreqs, oldNFFT
    freqs = f # Fs/NFFT*np.arange(numFreqs)
    print "FREQS ARE: ", freqs
    return All_final, freqs




def granger_test2(X, ij, newLength=256, NFFT=256, offset=0, Fs=2, maxlag=3, progressCallback=donothing_callback, window=window_hanning, noverlap=0, detrend = detrend_none, gv1=0,gv2=0):

    oldNFFT = NFFT
    NFFT = newLength
    numRows, numCols = X.shape

    if numRows < NFFT:
        tmp = X
        X = np.zeros( (NFFT, numCols), X.dtype)
        X[:numRows,:] = tmp
        del tmp

    numRows, numCols = X.shape
    # get all the columns of X that we are interested in by checking
    # the ij tuples
    allColumns = set()
    for i,j in ij:
        allColumns.add(i); allColumns.add(j)
    Ncols = len(allColumns)

    # for real X, ignore the negative frequencies
    if np.iscomplexobj(X): numFreqs = NFFT
    else: numFreqs = NFFT//2+1

    if cbook.iterable(window):
        assert(len(window) == NFFT)
        windowVals = window
    else:
        windowVals = window_hanning(np.ones(NFFT, X.dtype)) #I changed this from window to window_hanning. window was not doing anything!! -eli

    ind = range(offset, int(numRows-newLength+1), int(oldNFFT-noverlap)) #coherence calcs on each sweep start at offset
    numSlices = len(ind)
    threshold = .05/(numSlices*2)
    FFTSlices = {}
    Pxx = {}
    Cxy = {}
    Phase = {}
    
    if newLength > oldNFFT: #make sure that the newlength is shorter than the segment (epoch) length (see below)
        newLength = oldNFFT
    slices = range(numSlices)
    typedict = ['params_ftest', 'ssr_ftest']
    
    # normVal = np.linalg.norm(windowVals)**2
    counter = 0
    for i,j in ij: # FOR EACH ELECTRODE PAIR
        progressCallback(counter/len(ij), 'Cacheing FFTs')
        Slices = np.zeros((numSlices, 1))
        RevSlices = np.zeros((numSlices, 1))
        counter += 1
        counter2 = 0
        print i,j
        for iSlice in slices: # FOR EACH TRIAL
            thisSlice = X[ind[iSlice]:ind[iSlice]+newLength, i] #this is the line that reads sections of epochs
            thisSlice2 = X[ind[iSlice]:ind[iSlice]+newLength, j] 
            #print "GRANGER TESTING: ", ind[iSlice], " to ", ind[iSlice] + newLength
            #print "shape of all: ", Slices.shape
            #print "shape of slice: ", thisSlice.shape
            thisSlice = windowVals*detrend(thisSlice)
            thisSlice2 = windowVals*detrend(thisSlice2)
            d = np.vstack((thisSlice,thisSlice2)).T
            drev = np.vstack((thisSlice2,thisSlice)).T
            res = gtest.grangercausalitytests(d,maxlag,verbose=False)
            resrev = gtest.grangercausalitytests(drev,maxlag,verbose=False)
            print "RES: ", res 
            presult = res[maxlag][0][typedict[gv1]][not gv2]
            prev_result = resrev[maxlag][0][typedict[gv1]][not gv2]
            fresult = np.log(res[maxlag][0][typedict[gv1]][gv2])
            frev_result = np.log(resrev[maxlag][0][typedict[gv1]][gv2])
            print "RESULTS: ", fresult, frev_result
            if presult <= threshold and prev_result > threshold:
                Slices[counter2] = fresult
                RevSlices[counter2] = 0.
            elif presult > threshold and prev_result <= threshold:
                RevSlices[counter2] = frev_result
                Slices[counter2] = 0.
            elif presult <= threshold and prev_result <= threshold:
                Slices[counter2] = fresult
                RevSlices[counter2] = frev_result
            else: 
                Slices[counter2] = 0.
                RevSlices[counter2] = 0.
            
            counter2 += 1
            
        avg_forward = np.mean(Slices)
        avg_backwards = np.mean(RevSlices)
        Cxy[i,j] = max((avg_forward, avg_backwards))
        
        direction = (avg_forward,avg_backwards).index(Cxy[i,j])
        if direction == 0:
            Phase[i,j] = 90
        else:
            Phase[i,j] = -90
        print "TRODE RESULTS: ", Cxy[i,j], Phase[i,j]
        del Slices, RevSlices
        
        
    freqs = Fs/NFFT*np.arange(numFreqs)
    print "FREQS ARE: ", freqs
    return Cxy, Phase, freqs





def run_all(res_dict, lag, start, stop):
    for i in range(40):
        for j in range(40):
            if i != j:
                print "\nRUNNING WITH " + str(i) + " AND " + str(j)
                res_dict[(i,j)] = start_test(lag, i,j, start, stop)
    return res_dict

def start_run(lag, start, stop):
    res_dict = {}
    res_dict = run_all(res_dict, lag, start, stop)
    simple_results = {}
    for i in range(lag):
        print "now calculating for lag " + str(lag)
        lag_dict = {}
        for key in res_dict.keys():
            lag_dict[key] = (res_dict[key][i+1][0]['lrtest'], res_dict[key][i+1][0]['params_ftest'])
        # lag_dict is loaded for the given i lag
        simple_results[i+1] = lag_dict
    return simple_results

def start_all_runs(lag):
    time_simple_results = {}
    start = 0
    stop = 10000
    step = 20
    for ind in range(start, stop, step):
        print "index is now: " + str(ind)
        time_simple_results[ind] = start_run(lag, ind, ind+step)
    return time_simple_results    

def plot_results(res,curlag, ind1, ind2):
    x = []
    y = {}
    for pair in sorted(set(res[0][1].keys())):
        y[pair] = []
    for time in sorted(set(res.keys())):
        x.append(time)
        for lag in res[time]:
            if lag == curlag:
                for pair in sorted(set(res[time][lag].keys())):
                    y[pair].append(res[time][lag][pair][ind1][ind2])
    
    xdim = len(sorted(set(y.keys())))
    ydim = len(y[1,2])
    x = np.array(x)
    z = np.zeros((xdim, ydim))
    counter = 0
    for key in sorted(set(y.keys())):
        z[counter] = y[key]
        counter += 1
    print x.shape, z.shape
    return x,z


def go_granger(lag):
    res = start_all_runs(lag)
    return res
def plot_granger(res, lag, ind1):
    x,z = plot_results(res, lag, ind1, 0)
    plt.plot(x,z.T)
    plt.show()
    return x,z


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

    
        c = 0
        for i in self.e1:
            print i, self.e2[c]
            c += 1
    
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
        self.ax12 = self.fig.add_subplot(4,3,12, title="el3")

    def ddtf(self,el1,el2,el3,sample_rate=400,duration=20,step=128,increment=5):

        self.ax10.plot(el1)
        self.ax11.plot(el2)
        self.ax12.plot(el3)

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
                Szz=Szz+(fz * fz.conj())
                Sxy=Sxy+(fx * fy.conj())
                Sxz=Sxz+(fx * fz.conj())
                Syz=Syz+(fy * fz.conj())

                # print "Sxx shape: ", Sxx.shape
                # print "Sxy shape: ", Sxy.shape
                # print "Szy shape: ", Sxx.shape
                # print "Syz shape: ", Syz.shape

                Syx = Sxy.conj()
                Szx = Sxz.conj()
                Szy = Syz.conj()

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

            self.ax1.plot(f[0:step/4],NS11[0][0:step/4])
            self.ax2.plot(f[0:step/4],NS12[0][0:step/4])
            self.ax3.plot(f[0:step/4],NS13[0][0:step/4])
            self.ax4.plot(f[0:step/4],NS21[0][0:step/4])
            self.ax5.plot(f[0:step/4],NS22[0][0:step/4])
            self.ax6.plot(f[0:step/4],NS23[0][0:step/4])
            self.ax7.plot(f[0:step/4],NS31[0][0:step/4])
            self.ax8.plot(f[0:step/4],NS32[0][0:step/4])
            self.ax9.plot(f[0:step/4],NS33[0][0:step/4])






        return (f ,step,NS11, NS12, NS13, NS21, NS22, NS23, NS31, NS32, NS33)


def start_ddtf():
    d = DDTF()
    d.ddtf(d.e1,d.e2,d.e3)
    d.win.show()
