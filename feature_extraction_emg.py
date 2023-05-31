#!/usr/bin/env python3

import sys

# window of EMG data for processing 
emg_window_size = 128 
emg_window = [0] * emg_window_size

previous_sample = 0
emg = 0

alpha = -0.5

# mean value calcution for filter 
sum = 0
mean = 0

# time-domain features
# integrated EMG and mean absolute value 
iEMG = 0 
MAV = 0

# simple-square integral and variance 
SSI = 0
var = 0 

# myopulse percentage 
myop = 0  

def update_window(new_sample):
    """
    :param new_sample: current EMG sample 
    Adds EMG sample to window, removes oldest sample 
    """
    global previous_sample, emg
    emg = new_sample

    # Glitch filter
    if emg > 10000:
        print("GLITCH:", emg, file=sys.stderr)
        emg = previous_sample
        
    # First Order IIR Filter
    # emg = (1 - alpha)*emg + alpha*previous_sample

    # append new sample into window 
    emg_window.append(emg)
    # remove old sample from window 
    emg_window.pop(0)

    previous_sample = emg

def calc_time():
    """
    :param None:
    Calculate features in time-domain, store in global variables
    """
    global sum, mean, iEMG, MAV, SSI, var 

    # reset variables before assignment 
    sum = 0
    iEMG = 0
    SSI = 0 

    # for-loop for summing over all samples 
    for i in range(128):
        # sum for high-pass filter 
        sum += emg_window[i]
        # iEMG: sum of all samples 
        iEMG += abs(emg_window[i])
        # SSI:  sum of all square samples
        SSI += (emg_window[i]**2)

    # mean  = sum / num_samples
    mean = sum / emg_window_size    
    # MAV   = iEMG / num_samples 
    MAV = iEMG / emg_window_size 
    # var   = SSI / (num_samples-1)
    var = SSI / (emg_window_size - 1)
