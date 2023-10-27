# This is a solver for nonlinear dispersion relationship of Stokes theory.
# The full nonlinear dispersion relationship has h (mean water depth), either T
# (wave period) or omega0 (leading order frequency), and 
# either a or H is given. Results include k, H, a, omega, omega0. 
# Mode 1, the first order wave amplitude/height is given as a, H.
# Mode 2, the first harmonic term's coefficient is given as aw, Hw. 
# Example 1, 
# Results =StokesDispSolver(h=30,T=5, a=2.5, mode=1)
# Example 2, 
# Results =StokesDispSolver(h=30,T=5, a=2.5, mode=2)
# Example 3, 
# Results =StokesDispSolver(h=30,T=5, H=5, mode=1)
# Example 4, 
# Results =StokesDispSolver(h=1, omega0=1.2566, a=0.15, mode=1)
# Example 5, 
# Results =StokesDispSolver(h=1, omega0=1.2566, H=0.3, mode=1)
import sympy as sp
import numpy as np

def StokesDispSolver(**kwargs):
    h = kwargs['h']
    a = kwargs.get('a', None)
    H = kwargs.get('H', None)
    T = kwargs.get('T', None)
    omega0 = kwargs.get('omega0', None)
    modeNo = kwargs.get('mode', None)
    
    g = 9.81 
    #Make a guess at k using stokes second order for the numerical solver
    if T is not None:
        sigma = 2*np.pi/T
    elif omega0 is not None:
        sigma = omega0
    k0 = 0.0001
    k1 = sigma**2/(g*np.tanh(k0*h))
    error = 1
    while abs(error)>0.0001:
        k0 = k1
        k1 = sigma**2/(g*np.tanh(k0*h))
        error = k0-k1

    stokes2_k = k1 #guess at k




    Results = {}
    Results['h'] = h

    k = sp.Symbol('k')
    alpha1 = sp.cosh(2*k*h)
    sigma = sp.tanh(k*h)
    if a is not None and T is not None:
        omega0 = sp.sqrt(g*k*sp.tanh(k*h))
        Results['T'] = T
        omega = 2*sp.pi/T
        Results['omega'] = omega
        L0 = g*T**2/2/sp.pi
        if modeNo == 1: 
            omega2 = k**2*a**2*(2*alpha1**2+7)/4/(alpha1 -1)**2
            omega4 = (a**4*k**4*(20*alpha1**5 + 112*alpha1**4 - 100*alpha1**3 - 68*alpha1**2 - 
                        211*alpha1 + 328))/(32*(alpha1 - 1)**5)
            omegaFun = omega0*(1+omega2+omega4) - omega  
            kResults = sp.nsolve(omegaFun, stokes2_k)
            Results['a'] = a
            sigma = sp.tanh(kResults*h)
            B31 = (3+8*sigma**2-9*sigma**4)/16/sigma**4
            B51 = (121*sp.cosh(2*h*kResults)**5 + 263*sp.cosh(2*h*kResults)**4 +376*sp.cosh(2*h*kResults)**3
            - 1999*sp.cosh(2*h*kResults)**2 + 2509*sp.cosh(2*h*kResults) - 1108)/(192*(sp.cosh(2*h*kResults) - 1)**5)
            Results['aw'] = a*(1+kResults**2.*a**2*B31+kResults**4.*a**4*B51)
        elif modeNo == 2: # in case first harmonic wave amplitude aw is given
            a0 = sp.Symbol('a0')

            B31 = (3+8*sigma**2-9*sigma**4)/16/sigma**4
            B51 = (121*sp.cosh(2*h*k)**5 + 263*sp.cosh(2*h*k)**4 +376*sp.cosh(2*h*k)**3 
            - 1991*sp.cosh(2*h*k)**2 + 2509*sp.cosh(2*h*k) - 1108)/(192*(sp.cosh(2*h*k) - 1)**5)
            a0Fun = a - a0*(1+k**2.*a0**2*B31+k**4.*a0**4*B51)
            omega2 = k**2*a0**2*(2*alpha1**2+7)/4/(alpha1 -1)**2
            omega4 = (a0**4*k**4*(20*alpha1**5 + 112*alpha1**4 
            - 100*alpha1**3 - 68*alpha1**2 - 211*alpha1 + 328))/(32*(alpha1 - 1)**5)
            omegaFun = omega0*(1+omega2+omega4) - omega 
            S = sp.nsolve([a0Fun, omegaFun], [a0, k],[stokes2_k, a])
            kResults = abs(S[0])
            aResults = abs(S[1])
            Results['a'] = aResults
            Results['aw'] = a
            a = Results['a']
        
        Results['k'] = kResults
        Results['omega0'] = sp.sqrt(g*kResults*sp.tanh(kResults*h))
        sigma = sp.tanh(kResults*h)
        B31 = (3+8*sigma**2-9*sigma**4)/16/sigma**4
        B33 = (27-9*sigma**2+9*sigma**4-3*sigma**6)/64/sigma**6   
        B51 = (121*sp.cosh(2*h*kResults)**5 + 263*sp.cosh(2*h*kResults)**4 +376*sp.cosh(2*h*kResults)**3 
            - 1999*sp.cosh(2*h*kResults)**2 + 2509*sp.cosh(2*h*kResults) - 1108)/(192*(sp.cosh(2*h*kResults) - 1)**5)
        B53 = (57*sp.cosh(2*h*kResults)**7 + 204*sp.cosh(2*h*kResults)**6 - 53*sp.cosh(2*h*kResults)**5 
            -782*sp.cosh(2*h*kResults)**4 - 741*sp.cosh(2*h*kResults)**3 - 52*sp.cosh(2*h*kResults)**2 
            + 371*sp.cosh(2*h*kResults) + 186)*9/((3*sp.cosh(2*h*kResults) + 2)*(sp.cosh(2*h*kResults) - 1)**6*128)
        B55 = (300*sp.cosh(2*h*kResults)**8 + 1579*sp.cosh(2*h*kResults)**7 + 3176*sp.cosh(2*h*kResults)**6 + 2949*sp.cosh(2*h*kResults)**5 
            + 1188*sp.cosh(2*h*kResults)**4 + 675*sp.cosh(2*h*kResults)**3 + 1326*sp.cosh(2*h*kResults)**2 + 827*sp.cosh(2*h*kResults) 
            + 130)*5/((sp.cosh(2*h*kResults) - 1)**6*(12*sp.cosh(2*h*kResults)**2 + 11*sp.cosh(2*h*kResults) + 2)*384)
        Hout = 2*a+2*(B31+B33)*kResults**2*a**3+2*(B51+B53+B55)*kResults**4*a**5
        Results['H'] = Hout


    elif H is not None and T is not None:
        omega0 = sp.sqrt(g*k*sp.tanh(k*h))
        Results['T'] = T
        omega = 2*sp.pi/T
        Results['omega'] = omega
        aguess = H/2
        
        L0 = g*T**2/2/sp.pi
        a = sp.Symbol('a')
        sigma = sp.tanh(k*h)
        B33 = (27-9*sigma**2+9*sigma**4-3*sigma**6)/64/sigma**6        
        B55 = (300*sp.cosh(2*h*k)**8 + 1579*sp.cosh(2*h*k)**7 + 3176*sp.cosh(2*h*k)**6 + 2949*sp.cosh(2*h*k)**5 
            + 1188*sp.cosh(2*h*k)**4 + 675*sp.cosh(2*h*k)**3 + 1326*sp.cosh(2*h*k)**2 + 827*sp.cosh(2*h*k) 
            + 130)*5/((sp.cosh(2*h*k) - 1)**6*(12*sp.cosh(2*h*k)**2 + 11*sp.cosh(2*h*k) + 2)*384)
        B31 = (3+8*sigma**2-9*sigma**4)/16/sigma**4
        B51 = (121*sp.cosh(2*h*k)**5 + 263*sp.cosh(2*h*k)**4 +376*sp.cosh(2*h*k)**3 
            - 1999*sp.cosh(2*h*k)**2 + 2509*sp.cosh(2*h*k) - 1108)/(192*(sp.cosh(2*h*k) - 1)**5)
        B53 = (57*sp.cosh(2*h*k)**7 + 204*sp.cosh(2*h*k)**6 - 53*sp.cosh(2*h*k)**5 
            -782*sp.cosh(2*h*k)**4 - 741*sp.cosh(2*h*k)**3 - 52*sp.cosh(2*h*k)**2 
            + 371*sp.cosh(2*h*k) + 186)*9/((3*sp.cosh(2*h*k) + 2)*(sp.cosh(2*h*k) - 1)**6*128)
        Results['H'] = H
        Hfun1 = H -(2*a+2*(B31+B33)*k**2*a**3+2*(B51+B53+B55)*k**4*a**5)
        omega2 = k**2*a**2*(2*alpha1**2+7)/4/(alpha1 -1)**2
        omega4 = (a**4*k**4*(20*alpha1**5 + 112*alpha1**4 
                - 100*alpha1**3 - 68*alpha1**2 - 211*alpha1 + 328))/(32*(alpha1 - 1)**5)
        omegaFun1 = omega0*(1+omega2+omega4) - omega
        S = sp.nsolve([Hfun1, omegaFun1], (k, a), [stokes2_k, aguess])
        kResults = abs(S[0])
        aResults = abs(S[1])

        Results['a'] = aResults
        sigma = sp.tanh(kResults*h)
        B31 = (3+8*sigma**2-9*sigma**4)/16/sigma**4
        B51 = (121*sp.cosh(2*h*kResults)**5 + 263*sp.cosh(2*h*kResults)**4 +376*sp.cosh(2*h*kResults)**3 
        - 1999*sp.cosh(2*h*kResults)**2 + 2509*sp.cosh(2*h*kResults) - 1108)/(192*(sp.cosh(2*h*kResults) - 1)**5)
        Results['aw'] = aResults*(1+kResults**2.*aResults**2*B31+kResults**4.*aResults**4*B51)
        Results['omega0'] = sp.sqrt(g*kResults*sp.tanh(kResults*h))


    elif a is not None and omega0 is not None: 
        Results['omega0'] = omega0
        kResults = sp.nsolve(omega0**2 - g*k*sp.tanh(k*h), stokes2_k)
        alpha1 = sp.cosh(2*kResults*h)   
        sigma = sp.tanh(kResults*h)
        B33 = (27-9*sigma**2+9*sigma**4-3*sigma**6)/64/sigma**6    
        B55 = (300*sp.cosh(2*h*kResults)**8 + 1579*sp.cosh(2*h*kResults)**7 + 3176*sp.cosh(2*h*kResults)**6 + 2949*sp.cosh(2*h*kResults)**5 
            + 1188*sp.cosh(2*h*kResults)**4 + 675*sp.cosh(2*h*kResults)**3 + 1326*sp.cosh(2*h*kResults)**2 + 827*sp.cosh(2*h*kResults) 
            + 130)*5/((sp.cosh(2*h*kResults) - 1)**6*(12*sp.cosh(2*h*kResults)**2 + 11*sp.cosh(2*h*kResults) + 2)*384)
        B31 = (3+8*sigma**2-9*sigma**4)/16/sigma**4
        B51 = (121*sp.cosh(2*h*kResults)**5 + 263*sp.cosh(2*h*kResults)**4 +376*sp.cosh(2*h*kResults)**3 
            - 1999*sp.cosh(2*h*kResults)**2 + 2509*sp.cosh(2*h*kResults) - 1108)/(192*(sp.cosh(2*h*kResults) - 1)**5)
        B53 = (57*sp.cosh(2*h*kResults)**7 + 204*sp.cosh(2*h*kResults)**6 - 53*sp.cosh(2*h*kResults)**5 
            -782*sp.cosh(2*h*kResults)**4 - 741*sp.cosh(2*h*kResults)**3 - 52*sp.cosh(2*h*kResults)**2 
            + 371*sp.cosh(2*h*kResults) + 186)*9/((3*sp.cosh(2*h*kResults) + 2)*(sp.cosh(2*h*kResults) - 1)**6*128)

        if modeNo ==1:
                Results['a'] = a
                Results['aw'] = a*(1+kResults**2.*a**2*B31+kResults**4.*a**4*B51)
                Hout = 2*a+2*(B31+B33)*kResults**2*a**3+2*(B51+B53+B55)*kResults**4*a**5
                Results['H'] = Hout
                omega2 = kResults**2*a**2*(2*alpha1**2+7)/4/(alpha1 -1)**2
                omega4 = (a**4*kResults**4*(20*alpha1**5 + 112*alpha1**4 
            - 100*alpha1**3 - 68*alpha1**2 - 211*alpha1 + 328))/(32*(alpha1 - 1)**5)
                omega = omega0*(1+omega2+omega4)
                Results['omega'] = float(omega)
                Results['T'] = 2*sp.pi/omega
        elif modeNo ==2:
                Results['aw'] = a # in this case, users gave aw actually
                a0 = sp.Symbol('a0')
                a0Fun = a - a0*(1+k**2.*a0**2*B31+k**4.*a0**4*B51)
                Results['a'] = sp.nsolve(a0Fun, a0, a)
                a = Results['a']
                omega2 = kResults**2*a**2*(2*alpha1**2+7)/4/(alpha1 -1)**2
                omega4 = (a**4*kResults**4*(20*alpha1**5 + 112*alpha1**4 
                    - 100*alpha1**3 - 68*alpha1**2 - 211*alpha1 + 328))/(32*(alpha1 - 1)**5)
                omega = omega0*(1+omega2+omega4)
                Results['omega'] = float(omega)
                Results['T'] = 2*sp.pi/omega
                
                Hout = 2*a+2*(B31+B33)*kResults**2*a**3+2*(B51+B53+B55)*kResults**4*a**5
                Results['H'] = Hout
                
    elif H is not None and omega0 is not None:
        Results['omega0'] = omega0
        kResults = sp.nsolve(omega0**2 - g*k*sp.tanh(k*h), stokes2_k)
        alpha1 = sp.cosh(2*kResults*h)
        a = sp.Symbol('a')
        sigma = sp.tanh(kResults*h)
        B33 = (27-9*sigma**2+9*sigma**4-3*sigma**6)/64/sigma**6     
        B55 = (300*sp.cosh(2*h*kResults)**8 + 1579*sp.cosh(2*h*kResults)**7 + 3176*sp.cosh(2*h*kResults)**6 + 2949*sp.cosh(2*h*kResults)**5 
            + 1188*sp.cosh(2*h*kResults)**4 + 675*sp.cosh(2*h*kResults)**3 + 1326*sp.cosh(2*h*kResults)**2 + 827*sp.cosh(2*h*kResults) 
            + 130)*5/((sp.cosh(2*h*kResults) - 1)**6*(12*sp.cosh(2*h*kResults)**2 + 11*sp.cosh(2*h*kResults) + 2)*384)
        B31 = (3+8*sigma**2-9*sigma**4)/16/sigma**4
        B51 = (121*sp.cosh(2*h*kResults)**5 + 263*sp.cosh(2*h*kResults)**4 +376*sp.cosh(2*h*kResults)**3 
            - 1999*sp.cosh(2*h*kResults)**2 + 2509*sp.cosh(2*h*kResults) - 1108)/(192*(sp.cosh(2*h*kResults) - 1)**5)
        B53 = (57*sp.cosh(2*h*kResults)**7 + 204*sp.cosh(2*h*kResults)**6 - 53*sp.cosh(2*h*kResults)**5 
            -782*sp.cosh(2*h*kResults)**4 - 741*sp.cosh(2*h*kResults)**3 - 52*sp.cosh(2*h*kResults)**2 
            + 371*sp.cosh(2*h*kResults) + 186)*9/((3*sp.cosh(2*h*kResults) + 2)*(sp.cosh(2*h*kResults) - 1)**6*128)              
        Results['H'] = H
        Hfun1 = H -(2*a+2*(B31+B33)*kResults**2*a**3+2*(B51+B53+B55)*kResults**4*a**5)
        aResults = sp.nsolve(Hfun1, H/2)
        omega2 = kResults**2*aResults**2*(2*alpha1**2+7)/4/(alpha1 -1)**2
        omega4 = (aResults**4*kResults**4*(20*alpha1**5 + 112*alpha1**4 
            - 100*alpha1**3 - 68*alpha1**2 - 211*alpha1 + 328))/(32*(alpha1 - 1)**5)
        omega = omega0*(1+omega2+omega4)
        Results['a'] = aResults
        Results['omega'] = omega
        Results['T'] = 2*sp.pi/omega
        sigma = sp.tanh(kResults*h)
        B31 = (3+8*sigma**2-9*sigma**4)/16/sigma**4
        B51 = (121*sp.cosh(2*h*kResults)**5 + 263*sp.cosh(2*h*kResults)**4 +376*sp.cosh(2*h*kResults)**3 
            - 1999*sp.cosh(2*h*kResults)**2 + 2509*sp.cosh(2*h*kResults) - 1108)/(192*(sp.cosh(2*h*kResults) - 1)**5)
        Results['aw'] = aResults*(1+kResults**2.*aResults**2*B31+kResults**4.*aResults**4*B51)


    Results['k'] = kResults
    Results['L'] = 2*sp.pi/kResults
    Results['c'] = Results['L']/Results['T']
    return Results


if __name__ == '__main__':
    # Example 1, 
    Results =StokesDispSolver(h=30,T=5, a=2.5, mode=1)
    print(Results)
    # Example 2, 
    Results =StokesDispSolver(h=30,T=5, a=2.5, mode=2)
    print(Results)
    # Example 3, 
    Results =StokesDispSolver(h=30,T=5, H=5, mode=1)
    print(Results)
    # Example 4, 
    Results =StokesDispSolver(h=1, omega0=1.2566, a=0.15, mode=1)
    print(Results)
    # Example 5, 
    Results =StokesDispSolver(h=1, omega0=1.2566, H=0.3, mode=1)
    print(Results)