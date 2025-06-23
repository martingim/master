import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import fft
## fft of piston file
# filename = '2d_piston/boundary-piston/piston_files/1/fil3.dat'
# real_fft = np.fft.rfft(position_data)
# real_fft[500:] = 0
# fft_smoothed = np.fft.irfft(real_fft)
# data = np.loadtxt(filename)
# if len(data)%2:
#     data = data[:-1]
# position_data = (data-np.mean(data))*0.01

number_of_pistons = 14
speed_file = 'piston_speed.dat'
position_file = 'piston_position.dat'
file_samplerate = 100
end_time = 20

f = 1.425
wave_time = 1/f*10 #how many amplitudes to generate
piston_amplitude = 0.025
piston_phase_offset = 0.2 # how much the pistons lag behind each other


generated_position = np.zeros((end_time*file_samplerate,number_of_pistons))
t = np.linspace(0,wave_time, int(wave_time*100)+1)

offset = 0

for i in range(number_of_pistons):
    generated_position[int(offset*file_samplerate):len(t)+int(offset*file_samplerate),i] = piston_amplitude*np.sin((t)*f*2*np.pi)*np.tanh(t)*np.tanh(wave_time-t)
    if i<6:
        offset += piston_phase_offset
    if i>6:
        offset -= piston_phase_offset

np.savetxt(speed_file, np.gradient(generated_position, axis=0)*file_samplerate, fmt='%.8f')
np.savetxt(position_file, generated_position, fmt='%.8f')
