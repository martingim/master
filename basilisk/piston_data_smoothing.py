import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import leastsq
from scipy import fft
# smoothes the voltage data from the reading of the piston position
# takes the derivative and converts to m/s
# saves the result to piston_speed.dat in the same folder as fil1.dat
# run with: find . -name "fil1.dat" | xargs python3 piston_data_smoothing.py
# to run in all subdirs where fil1.dat is present
function_fit = 0
fft_lowpass = 1
savgol_fit = 0
fil3 = True # use fil3 instead of the volt reading from the piston
window_length = 20
polyorder = 4
f = 1.425
guess_amp = 0.1
guess_tanh_factor = 1
guess_phase_offset = 1
guess_tanh_offset = 1
wave_time = 70
ramptime = 3
verbose = 1
if len(sys.argv)>1:
    filenames = sys.argv[1:]
else:
    print("no filename supplied in command line")
    filenames = ['2d_piston/boundary-piston/piston_files/1/fil1.dat']
    if fil3:
        filenames = ['2d_piston/boundary-piston/piston_files/1/fil3.dat']
file_samplerate = 100
if fil3:
    convert_to_meters = 0.01
else:
    convert_to_meters = -0.044


for filename in filenames:
    data = np.loadtxt(filename)
    if len(data)%2:
        data = data[:-1]

    if wave_time*file_samplerate>len(data):
        wave_time = int(len(data)/100)
    position = (data-np.mean(data))*convert_to_meters
    mean_guess = np.mean(position)
    print(f"loaded: {filename}")
    t = np.linspace(0,len(position)/file_samplerate, len(position))
    
    # calculate the least squares amplitude
    optimize_func = lambda x: x[0]*np.sin(f*(t[file_samplerate*ramptime:-file_samplerate*(wave_time-ramptime)]+x[1])*2*np.pi) + x[2] - position[file_samplerate*ramptime:-file_samplerate*(wave_time-ramptime)] 
    res = leastsq(optimize_func, [guess_amp, guess_phase_offset, mean_guess])
    amp = res[0][0]
    phase_offset = res[0][1]
    offset = res[0][2]
    
    
    # calulate the starting rampup tanh function
    optimize_func = lambda x: amp*np.sin(f*(t[:file_samplerate*ramptime]+phase_offset)*2*np.pi)*(np.tanh(x[0]*(t[:file_samplerate*ramptime]+x[1]))).clip(min=0)+offset - position[:file_samplerate*ramptime]
    res = leastsq(optimize_func, [guess_tanh_factor, guess_tanh_offset])    
    tanh_factor = res[0][0]
    tanh_offset = res[0][1]
    
    # calculate the rampdown function
    optimize_func = lambda x: amp*np.sin(f*(t[(wave_time-ramptime)*file_samplerate:file_samplerate*wave_time]+phase_offset)*2*np.pi)*(np.tanh(tanh_factor*(t[(wave_time-ramptime)*file_samplerate:file_samplerate*wave_time]+tanh_offset))).clip(min=0)*(np.tanh((x[0]-t[(wave_time-ramptime)*file_samplerate:file_samplerate*wave_time])*x[1])).clip(min=0)+offset - position[(wave_time-ramptime)*file_samplerate:file_samplerate*wave_time]
    res = leastsq(optimize_func, [wave_time, 1])
    rampdown_time = res[0][0]
    rampdown_tanh = res[0][1]
    
    if verbose:
        print('Function fit parameters')
        print('amplitude', amp)
        print('sin phase offset', phase_offset)
        print('tanh_factor', tanh_factor)
        print('tanh_offset', tanh_offset)
        print('rampdown_time:',rampdown_time)
        print('rampdown_tanh', rampdown_tanh)    
    
    fitted_position = amp*np.sin(f*(t+phase_offset)*2*np.pi)*np.tanh((t+tanh_offset)*tanh_factor).clip(min=0)*np.tanh((rampdown_time-t)*rampdown_tanh).clip(min=0)+offset

    ## SAVGOL filter
    smoothed_derivative = savgol_filter(data, window_length, polyorder, deriv=1, delta=1)
    piston_speed = -smoothed_derivative*convert_to_meters*file_samplerate

    smoothed_voltage = savgol_filter(data, window_length, polyorder, delta=1)
    piston_position = smoothed_voltage*convert_to_meters


    ##fft
    real_fft = np.fft.rfft(position)
    real_fft[500:] = 0
    fft_smoothed = np.fft.irfft(real_fft)

    print(f"saving results in: {filename[:-8]}")    

    speed_file = filename[:-8]+'piston_speed.dat'
    position_file = filename[:-8]+'piston_position.dat'

    # save results to file
    if function_fit:
        np.savetxt(speed_file, np.gradient(fitted_position)*file_samplerate, fmt='%.8f')
        np.savetxt(position_file, fitted_position, fmt='%.8f')
    elif fft_lowpass:
        np.savetxt(speed_file, np.gradient(fft_smoothed)*file_samplerate, fmt='%.8f')
        np.savetxt(position_file, fft_smoothed, fmt='%.8f')
    elif savgol_fit:
        np.savetxt(speed_file, np.gradient(fitted_position)*file_samplerate, fmt='%.8f')
        np.savetxt(position_file, fitted_position, fmt='%.8f')

    plt.figure()
    plt.plot(t, position)
    plt.plot(t, fitted_position)
    plt.plot(t, piston_position-piston_position.mean())
    plt.plot(t, fft_smoothed)
    plt.legend(['data', 'func fit', 'savgol smoothing', 'fft_smoothing'])
    plt.title(filename[-10:-9])
    # plt.plot(t[-100:], position[-100:])
    # plt.show()
    #plt.plot(t, savgol_data*np.tanh(t*100))
    #plt.plot(t, np.gradient(data))
    #plt.show()
#plt.show()