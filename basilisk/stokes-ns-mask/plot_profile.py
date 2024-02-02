import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("profiles/res_0.0000.txt", names=['x', 'y'], sep=' ')
x = df['x'].to_numpy()
y = df['y'].to_numpy()

plt.plot(x, y, 'x')
plt.show()

