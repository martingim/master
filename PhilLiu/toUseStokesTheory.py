from StokesDispSolver import StokesDispSolver
from StokesU import stokes_u
from StokesEta import StokesEta

import numpy as np
import matplotlib.pyplot as plt

# Input region
h0 = 10
modeNo = 1
a = 1.6
T = 5
theta0 = np.arange(0, 2 * np.pi, 0.01)

# Call StokesDispSolver
Result = StokesDispSolver('h', h0, 'T', T, 'a', a, 'mode', modeNo)

if modeNo == 1:
    a = Result.a
    H = Result.H
elif modeNo == 2:
    a = Result.aw
    H = Result.Hw

# Call StokesEta function
theta = theta0  # Assuming theta0 is the same as theta in the MATLAB code
eta, eta1, eta2, eta3, eta4, eta5 = StokesEta(Result.k, h0, Result.a, theta)

# Plot the result
plt.plot(theta0, eta)
c = Result.L / T  # phase speed

# Velocity and potential function
phi = np.zeros(11)
u1 = np.zeros(11)
u2 = np.zeros(11)
u3 = np.zeros(11)
u4 = np.zeros(11)
u5 = np.zeros(11)
w1 = np.zeros(11)
w2 = np.zeros(11)
w3 = np.zeros(11)
w4 = np.zeros(11)
w5 = np.zeros(11)

for i in range(len(theta)):
    z = np.concatenate([np.linspace(-h, eta[i], 10), [-2.5]])
    phi[i], u, u1[i], u2[i], u3[i], u4[i], u5[i], w, w1[i], w2[i], w3[i], w4[i], w5[i] = StokesU(Result.k, h0, Result.a, theta[i], z)

Result.u = u
Result.w = w

# To-Do: Implement FentonEta function if needed
# etaFen = FentonEta(Result.k, h0, a, theta0)

# For the second mode, a is combined first harmonic amplitude
# modeNo = 2
# Result2 = StokesDispSolver('h', h0, 'T', T0, 'a', a, 'mode', 2)

# Then you can use StokesEta and StokesU functions for Result2
