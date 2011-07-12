import numpy as np
import scikits.statsmodels.tsa.stattools as gtest
import math
import matplotlib.cbook as cbook
import signal_gen


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


def granger_test(X, ij, newLength=256, NFFT=256, offset=0, Fs=2, maxlag=3, progressCallback=donothing_callback, window=window_hanning, noverlap=0, detrend = detrend_none, gv1=0,gv2=0):
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

