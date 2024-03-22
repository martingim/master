import numpy as np

filename  = "waveinput.dat"

ntimesteps = 24000

f = 0.5 #frequency of the piston in hz
a = 0.01 #amplitude of the piston motion
T = 20
t = np.linspace(0,T,ntimesteps)
x = a*np.sin(t*2*np.pi*f)
v = a*np.pi*2*f*np.sin(t*2*np.pi*f)

with open(filename, 'w') as f:
    f.write("@v213\n")
    f.write("""[wave type]
wavemaker
[general input data]
depth  0.6
swl  0.
[pistonwavemaker wave properties]
#alpha_z, alpha_u
0. 0.\n""")
    f.write(f"#ntimestep\n{ntimesteps}\n")
    for i in range(ntimesteps):
        f.write(f"{t[i]:.4f} 0.0 {v[i]:.7f} {x[i]:.7f}\n")