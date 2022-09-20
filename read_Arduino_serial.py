# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:47:50 2022

@author: Ryan.Larson
"""

import serial
import serial.tools.list_ports
import time
import pandas as pd
import matplotlib.pyplot as plt
from os.path import exists
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import butter, lfilter, freqz
from scipy import signal


def find_Arduino_and_connect():
    # Determine which port is connected to Arduino
    print("Searching for Arduino on COM ports...")
    ports = list(serial.tools.list_ports.comports())
    Arduino_port = None
    for p in ports:
        # print(p)
        if "Arduino" in p.description:
            Arduino_port = str(p)
            # find the substring with "COM" and the number
            COMstart = Arduino_port.find('COM')
            COMnum = Arduino_port[COMstart+3]
            # print(COMstart)
            Arduino_port = Arduino_port[COMstart:(COMstart+5)]
            # print(Arduino_port)
            # input("Press Enter to continue...")
    try:
        ser = serial.Serial(Arduino_port, 115200, timeout=1)
        time.sleep(2)
        print("STATUS: Connected to Arduino on COM port {}\n".format(COMnum))
    except:
        # if ser:
        #     disconnect_Arduino(ser)
        # else:
        ser = None
        raise Exception("STATUS: Failed to connect to Arduino")
    
    return ser


def disconnect_Arduino(ser):
    if (ser.isOpen() == True):
        ser.close()
        print("STATUS: Serial communication with Arduino CLOSED")


# def read_data_nlines(ser, nlines):
#     data = []
#     for i in range(nlines):
#         line = ser.readline()   # read a byte string
#         if line:
#             string = line.decode()  # convert the byte string to a unicode string
#             print(string)
#             data.append(string)
#     disconnect_Arduino(ser)
    
#     return data

# def read_data_nseconds(ser, nseconds):
#     data = []
#     times = []
#     tstart = time.time()
#     print("START:\t{}".format(tstart))
#     while time.time() < (tstart + nseconds):
#         line = ser.readline()
#         if line:
#             string = line.decode()  # convert the byte string to a unicode string
#             # print(string)
#             data.append(string)
#             times.append(time.time())
#     tend = time.time()
#     print("END:\t{}".format(tend))
#     print("ELAPSED:\t{}".format(tend - tstart))
#     disconnect_Arduino(ser)
    
#     return data, times


def watch_serial(ser):
    if (ser.isOpen() == False):
        ser.open()
    
    data = []
    times = []
    data_captured = False
    while True:
        # print("Looking for switch to turn on...")
        line = ser.readline()
        string = line.decode()
        if "HIGH" in string:
            print("Start of data capture")
            while True:
                current_line = ser.readline().decode()
                data.append(current_line)
                times.append(time.time())
                # print(current_line)
                if "LOW" in current_line:
                    # data.append(string)
                    # times.append(time.time())
                    data_captured = True
                    print("End of data capture")
                    break
        if data_captured == True:
            break
    
    # Cleanup for next run
    disconnect_Arduino(ser)
    
    return data, times


def clean_data(data, times):
    if len(data) != len(times):
        raise ValueError("Data and times are not the same length")
    
    # Subtract t0 from all times
    t0 = times[0]
    for i in range(len(times)):
        times[i] -= t0
    
    # Join the two lists as columns
    data_times = []
    for i in range(len(data)):
        data_times.append([data[i], times[i]])
    
    # Clean data_times
    data_times_cleaned = [x for x in data_times if "HIGH" not in x[0]]
    data_times_cleaned = [x for x in data_times_cleaned if "LOW" not in x[0]]
    
    # Separate data and times again
    data_cleaned = []
    times_cleaned = []
    for i in range(len(data_times_cleaned)):
        data_cleaned.append(data_times_cleaned[i][0])
        times_cleaned.append(data_times_cleaned[i][1])
        
    # Separate data into x, y, and z lists
    x = []
    y = []
    z = []
    for i,line in enumerate(data_cleaned):
        split_line = line.split(',')
        x.append(float(split_line[0]))
        y.append(float(split_line[1]))
        z.append(float(split_line[2]))
        
    return x, y, z, times_cleaned
    

def send_to_csv(x, y, z, times_cleaned, vel, xnew, ynew, znew, times_regular, vel_filtered, filepath):
    d_measured = {"time": times_cleaned, "x": x, "y": y, "z": z, "velocity": vel}
    d_processed = {"x_interp": xnew, "y_interp": ynew, "z_interp": znew,
                   "time_regular": times_regular, "velocity_drift_filtered": vel_filtered}
    df_measured = pd.DataFrame(d_measured)
    df_processed = pd.DataFrame(d_processed)
    df = pd.concat([df_measured, df_processed], axis=1)
    df.to_csv(filepath)
    
    return df


def prepare_accel_data(x,y,z,times,vel):
    # Make sure time data is evenly sampled
    time_diffs = np.diff(times)
    # time_diffs = []
    # for i in range(1,len(times)):
    #     time_diffs.append(times[i] - times[i-1])
        
    tstep = np.mean(time_diffs)
    n_tsteps = int(np.max(times) // tstep)
    times_regular = list(np.arange(0, n_tsteps*tstep, tstep))
    while np.min(times_regular) < np.min(times):
        del times_regular[0]
    
    # Interpolate x, y, and z acceleration at times_regular
    fx = interp1d(times, x)
    fy = interp1d(times, y)
    fz = interp1d(times, z)
    fvel = interp1d(times, vel)
    
    xnew = fx(times_regular)
    ynew = fy(times_regular)
    znew = fz(times_regular)
    velnew = fvel(times_regular)
    
    tdelta = times_regular[1] - times_regular[0] 
    
    # X = np.fft.fft(xnew)
    # Y = np.fft.fft(ynew)
    # Z = np.fft.fft(znew)
    
    # freqx = np.fft.fftfreq(len(xnew), tdelta)
    # freqy = np.fft.fftfreq(len(ynew), tdelta)
    # freqz = np.fft.fftfreq(len(znew), tdelta)
    
    # # Filter requirements.
    # order = 6
    # fs = 1.0/tdelta       # sample rate, Hz
    # cutoff = 7.0  # desired cutoff frequency of the filter, Hz
    
    # # Get the filter coefficients so we can check its frequency response.
    # b, a = butter_lowpass(cutoff, fs, order)
    
    # xfiltered = butter_lowpass_filter(xnew, cutoff, fs, order)
    # yfiltered = butter_lowpass_filter(ynew, cutoff, fs, order)
    # zfiltered = butter_lowpass_filter(znew, cutoff, fs, order)
    
    return xnew, ynew, znew, times_regular, velnew, tdelta
    
    
    
    
def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

    
    

def integrate_acceleration(time, accel):
    vel = []
    for i in range(len(time)):
        if i == 0:
            vel.append(0)
        else:
            current_vel = (time[i] - time[i-1])*((accel[i-1] + accel[i]) / 2) + vel[i-1]
            vel.append(current_vel)
            
    return vel


if __name__ == "__main__":
    ser = find_Arduino_and_connect()
    directory = "C:/Users/Ryan.Larson.ROCKWELLINC/OneDrive - Rockwell Inc/Desktop/"
    # directory = input("Enter the path to the data directory, ending in /")
    
    while True:
        try:
            filename = input("\nEnter the name of the output csv file:\t")
            filepath = directory + filename + ".csv"
            while exists(filepath):
                filename = input("Chosen file name already exists. Enter a new file name:\t")
                filepath = directory + filename + ".csv"
                if exists(filepath) == False:
                    break
            print("Ready")
            data, times = watch_serial(ser)     # The main loop should hold here while the box switch is on
            print("Cleaning data...")
            x,y,z,times = clean_data(data, times)
            print("Integrating acceleration...")
            vel = integrate_acceleration(times, z)  # Assumes z axis is where velocity will be measured
            
            print("Interpolating data...")
            xnew, ynew, znew, times_regular, velnew, tdelta = prepare_accel_data(x,y,z,times,vel)
            fs= 1.0/tdelta
            print("Filtering drift out of velocity...")
            vel_filtered = butter_highpass_filter(velnew, 0.5, fs, order=6)
            print("Sending data to csv file...")
            df = send_to_csv(x, y, z, times, vel, xnew, ynew, znew, times_regular, vel_filtered, filepath)
            
            # # Determine the value of the velocity just before impact
            # for i,accel in enumerate(df["z"]):
                
            
            plt.figure()
            plt.plot(times,vel, label="velocity")
            plt.plot(times_regular,vel_filtered, label="velocity drift filtered")
            plt.legend()
            plt.xlabel("Time (s)")
            plt.ylabel("Velocity (m/s)")
            plt.title(filename)
            
        except KeyboardInterrupt:
            disconnect_Arduino(ser)
            print("\nExiting program")
            sys.exit(0)