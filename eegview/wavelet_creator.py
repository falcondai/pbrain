import numpy as np

def decreasing_sin(newLength,eegfreq,modval):
    srate = eegfreq
    newLength_t = newLength * 1.0/eegfreq
    # want the freq of the wave to be 4 per length in points
    t = np.arange(0,newLength_t,1.0/srate) # time axis in seconds
    modifier = modval/newLength_t
    freq_array = np.arange(modifier,modifier/2,-((modifier)/(newLength*2)))
    wave = np.sin(2.0*np.pi*freq_array*t)
    return wave

def simple_sin(newLength,eegfreq,modval):
    srate = eegfreq
    newLength_t = newLength * 1.0/eegfreq
    # want the freq of the wave to be 4 per length in points
    t = np.arange(0,newLength_t,1.0/srate) # time axis in seconds
    modifier = modval/newLength_t
    freq_array = np.arange(modifier,modifier/2,-((modifier)/(newLength*2)))
    wave = np.sin(2.0*np.pi*(modifier - modifier/2)*t)
    return wave

def increasing_sin(newLength,eegfreq,modval):
    srate = eegfreq
    newLength_t = newLength * 1.0/eegfreq
    # want the freq of the wave to be 4 per length in points
    t = np.arange(0,newLength_t,1.0/srate) # time axis in seconds
    modifier = modval/newLength_t
    freq_array = np.arange(modifier/2,modifier,((modifier)/(newLength*2)))
    wave = np.sin(2.0*np.pi*freq_array*t)
    return wave

def create_all(wl,ef,mv):
    a = {}
    a['increasing_sin'] = increasing_sin(wl,ef,mv)
    a['simple_sin'] = simple_sin(wl,ef,mv)
    a['decreasing_sin'] = decreasing_sin(wl,ef,mv)
    return a
             
