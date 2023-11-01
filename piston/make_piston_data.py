import numpy as np

filename  = "waveinput.dat"

ntimesteps = 24000

f = 1 #frequency of the piston in hz
a = 0.1 #amplitude of the piston motion
T = 10
t = np.linspace(0,T,ntimesteps)
x = a*np.cos(t*2*np.pi*f)
v = a*np.sin(t*2*np.pi*f)

with open(filename, 'w') as f:
    f.write("@v213")
    f.write(f"[piston wavemaker properties]\n")
    f.write(f"ntimesteps {ntimesteps}\n")
    for i in range(ntimesteps):
        f.write(f"{t[i]:.5f} {x[i]:.7f} {v[i]:.7f}\n")

