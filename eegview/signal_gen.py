import numpy as np
import math
import scipy
import sys

def signal_gen(SNR, K2, K3):

    sample_rate = 400
    freq = 30.
    tim = 40
    dly2 = 5
    dly3 = 10
    # SNR = 3.
    # K2 = .1
    # K3 = .001

    NL = 0.
    f_Nyq = sample_rate / 2.
    dt = 1./sample_rate
    N = tim * sample_rate
    #print "dt is: ", dt
    t = np.arange(0,tim,dt)
    noise = np.random.rand(len(t))
    #print "noise is len: ", len(noise)
    null = np.zeros(len(t) - 1)
    t = np.sin(t)
    S1 = 2*math.pi*freq*t * SNR * scipy.std(noise,0)
    
    S2 = np.zeros(len(S1))
    S3 = np.zeros(len(S1))
    #print "n is : ", N
    np.random.seed(2)
    for k in range(dly2+1, N):
        if (NL == 0):
            np.random.seed(612)
            
            S2[k] = S1[k-dly2]*K2 + np.random.rand(1)
        else:
            np.random.seed(k)
            S2[k] = math.pow(S1[k-dly2], 2) * K2 + np.random.rand(1)

    for k in range(dly3+1, N):
        if (NL == 0):
            np.random.seed(k)
            S3[k] = S1[k-dly3]*K3 + np.random.rand(1)
        else:
            np.random.seed(k+1)
            S3[k] = math.pow(S1[k-dly3], 2) * K3 + np.random.rand(1)

    S1 = S1 + noise
    S1 = (S1 - np.mean(S1))/scipy.std(S1, 0)
    S2 = (S2 - np.mean(S2))/scipy.std(S2, 0)
    S3 = (S3 - np.mean(S3))/scipy.std(S3, 0)

    #print "signal 1: ", S1
    #print "signal 2: ", S2
    #print "signal 3: ", S3

    #print "signal 1: ", S1.shape
    #print "signal 2: ", S2.shape
    #print "signal 3: ", S3.shape

    return (S1, S2, S3)
#(s1, s2, s3) = signal_gen(3,.1,.001)
#for val in s3:
#    sys.stdout.write(str(val) + ",")
