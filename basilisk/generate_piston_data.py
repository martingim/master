import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import fft
filename = '2d_piston/boundary-piston/piston_files/1/fil3.dat'
speed_file = 'piston_speed.dat'
position_file = 'piston_position.dat'
file_samplerate = 100
f = 1.425
wave_time = 1/f*5
end_time = 20
data = np.loadtxt(filename)
if len(data)%2:
    data = data[:-1]
position_data = (data-np.mean(data))*0.01
real_fft = np.fft.rfft(position_data)
real_fft[500:] = 0
fft_smoothed = np.fft.irfft(real_fft)


t = np.linspace(0,wave_time, int(wave_time*100+1))
generated_position = 0.0125*np.sin(t*f*2*np.pi)*np.tanh(t)*np.tanh(wave_time-t)
t = np.linspace(0, end_time, end_time*file_samplerate+1)
generated_position = np.concatenate((generated_position, np.zeros(len(t)+len(generated_position))))

np.savetxt(speed_file, np.gradient(generated_position)*file_samplerate, fmt='%.8f')
np.savetxt(position_file, generated_position, fmt='%.8f')



#plt.plot(fft_smoothed)
#plt.plot(generated_position)
#plt.show()