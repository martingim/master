import numpy as np
from numpy import sin, cos, pi, sqrt, exp



def eta_stokes(theta):
    """free surface Stokes 3rd order deep water waves"""
    return a*((1-1/16*(k*a)**2)*cos(theta) +  1/2*k*a*cos(2*theta) + 3/8*(k*a)**2*cos(3*theta))


def phi_stokes(x, z, t):
    """velocity potenital for Stokes 3rd order deep water"""
    return sqrt(g/k)*exp(k*z)*sin(theta(x, t))

def U_stokes(x, z, t):
    """horizontal velocity for Stokes 3rd order deep water"""
    return -sqrt(g*k)*exp(k*z)*cos(theta(x, t))

def V_stokes(x, z, t):
    """Vertical velocity for Stokes 3rd order deep water waves"""
    return phi(x, z, t)*k

def alpha_kz_stokes(k, z):
    return exp(k*z)

if __name__ == "__main__":
    import matplotlib.pyplot as plt


    a = 0.1
    k = 1
    tau = 1
    omega = 2*pi/tau 
    g = 9.81

    N = 1000
    N_quiver = 20
    def theta(x, t):
        return x*k-omega*t
    x = np.linspace(0,11,N)
    plt.plot(x, eta(theta(x, 0)))
    plt.ylim(-1, 0.2)


    x0 = 1.7*pi
    z = np.linspace(-1,eta(theta(x0,0)), N_quiver)
    x0 = np.zeros(N_quiver) + x0
    plt.quiver(x0, z, U(x0,z,0), V(x0, z, 0))
    plt.show()

