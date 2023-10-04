import pandas as pd
import matplotlib.pyplot as plt
from time import sleep
vel = pd.read_csv("velocities.txt", names=['t', 'y', 'ux', 'eta'])
print(vel[vel['t']==5])

    

plt.figure()
for t in vel['t'].unique():
    df = vel[vel['t']==t]
    print(df['ux'].values)
    plt.plot(df['ux'].values, df['y'].values)
    plt.pause(1)

plt.show()