import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
# smoothes the voltage data from the reading of the piston position
# takes the derivative and converts to m/s
# saves the result to piston_speed.dat in the same folder as fil1.dat
# run with: find . -name "fil1.dat" | xargs python3 piston_data_smoothing.py
# to run in all subdirs where fil1.dat is present
if len(sys.argv)>1:
    filenames = sys.argv[1:]
else:
    print("no filename supplied in command line")
    filenames = ['2d_piston/boundary-piston/piston_files/2/fil1.dat']
file_samplerate = 100
volt_to_meters = 0.044


for filename in filenames:
    data = np.loadtxt(filename)
    print(f"loaded: {filename}")
    t = np.linspace(0,len(data)/100, len(data))


    smoothed_derivative = savgol_filter(data, 40, 3, deriv=1, delta=1)
    piston_speed = smoothed_derivative*volt_to_meters*file_samplerate

    smoothed_position = savgol_filter(data, 40, 3, delta=1)
    piston_position = smoothed_position*volt_to_meters

    speed_file = filename[:-8]+'piston_speed.dat'
    position_file = filename[:-8]+'piston_position.dat'

    print(f"saving results in: {filename[:-8]}")
    np.savetxt(speed_file, piston_speed, fmt='%.8f')
    np.savetxt(position_file, piston_position, fmt='%.8f')


    #plt.plot(t, savgol_data*np.tanh(t*100))
    #plt.plot(t, np.gradient(data))
    #plt.legend(['smoothed', 'original'])
    #plt.show()
