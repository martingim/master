import sympy as sp

def StokesDispSolver(**kwargs):
    h = kwargs['h']
    a = kwargs.get('a', None)
    H = kwargs.get('H', None)
    T = kwargs.get('T', None)
    omega0 = kwargs.get('omega0', None)
    modeNo = kwargs.get('mode', None)

    Results = {}
    Results['h'] = h
    k = sp.symbols('k')
    g = 9.81
    alpha1 = sp.cosh(2 * k * h)
    sigma = sp.tanh(k * h)
    
    if a is not None and T is not None:
        omega0 = sp.sqrt(g * k * sp.tanh(k * h))
        Results['T'] = T
        omega = 2 * sp.pi / T
        Results['omega'] = omega
        L0 = g * T ** 2 / (2 * sp.pi)

        omega2 = k ** 2 * a ** 2 * (2 * alpha1 ** 2 + 7) / 4 / (alpha1 - 1) ** 2
        omega4 = (a ** 4 * k ** 4 * (20 * alpha1 ** 5 + 112 * alpha1 ** 4
                                    - 100 * alpha1 ** 3 - 68 * alpha1 ** 2 - 211 * alpha1 + 328)) / (32 * (alpha1 - 1) ** 5)
        omegaFun = omega0 * (1 + omega2 + omega4) - omega
        kResults = sp.solve(omegaFun, k)
        kResults = [k_val.evalf() for k_val in kResults]

        Results['a'] = a
        sigma = sp.tanh(kResults[0] * h)
        B31 = (3 + 8 * sigma ** 2 - 9 * sigma ** 4) / 16 / sigma ** 4
        B51 = (121 * sp.cosh(2 * h * kResults[0]) ** 5 + 263 * sp.cosh(2 * h * kResults[0]) ** 4 +
               376 * sp.cosh(2 * h * kResults[0]) ** 3 - 1999 * sp.cosh(2 * h * kResults[0]) ** 2 +
               2509 * sp.cosh(2 * h * kResults[0]) - 1108) / (192 * (sp.cosh(2 * h * kResults[0]) - 1) ** 5)
        Results['aw'] = a * (1 + kResults[0] ** 2 * a ** 2 * B31 + kResults[0] ** 4 * a ** 4 * B51)
        Results['k'] = kResults[0]
        Results['omega0'] = sp.sqrt(g * kResults[0] * sp.tanh(kResults[0] * h))
        sigma = sp.tanh(kResults[0] * h)
        B31 = (3 + 8 * sigma ** 2 - 9 * sigma ** 4) / 16 / sigma ** 4
        B33 = (27 - 9 * sigma ** 2 + 9 * sigma ** 4 - 3 * sigma ** 6) / 64 / sigma ** 6
        B51 = (121 * sp.cosh(2 * h * kResults[0]) ** 5 + 263 * sp.cosh(2 * h * kResults[0]) ** 4 +
               376 * sp.cosh(2 * h * kResults[0]) ** 3 - 1999 * sp.cosh(2 * h * kResults[0]) ** 2 +
               2509 * sp.cosh(2 * h * kResults[0]) - 1108) / (192 * (sp.cosh(2 * h * kResults[0]) - 1) ** 5)
        B53 = (57 * sp.cosh(2 * h * kResults[0]) ** 7 + 204 * sp.cosh(2 * h * kResults[0]) ** 6 - 53 * sp.cosh(2 * h * kResults[0]) ** 5 -
               782 * sp.cosh(2 * h * kResults[0]) ** 4 - 741 * sp.cosh(2 * h * kResults[0]) ** 3 - 52 * sp.cosh(2 * h * kResults[0]) ** 2 +
               371 * sp.cosh(2 * h * kResults[0]) + 186) * 9 / ((3 * sp.cosh(2 * h * kResults[0]) + 2) *
               (sp.cosh(2 * h * kResults[0]) - 1) ** 6 * 128)
        B55 = (300 * sp.cosh(2 * h * kResults[0]) ** 8 + 1579 * sp.cosh(2 * h * kResults[0]) ** 7 +
               3176 * sp.cosh(2 * h * kResults[0]) ** 6 + 2949 * sp.cosh(2 * h * kResults[0]) ** 5 +
               1188 * sp.cosh(2 * h * kResults[0]) ** 4 + 675 * sp.cosh(2 * h * kResults[0]) ** 3 +
               1326 * sp.cosh(2 * h * kResults[0]) ** 2 + 827 * sp.cosh(2 * h * kResults[0]) +
               130) * 5 / ((sp.cosh(2 * h * kResults[0]) - 1) ** 6 * (12 * sp.cosh(2 * h * kResults[0]) ** 2 +
               11 * sp.cosh(2 * h * kResults[0]) + 2) * 384)
        Hout = 2 * a + 2 * (B31 + B33) * kResults[0] ** 2 * a ** 3 + 2 * (B51 + B53 + B55) * kResults[0] ** 4 * a ** 5
        Results['H'] = Hout

    elif H is not None and T is not None:
        omega0 = sp.sqrt(g * k * sp.tanh(k * h))
        Results['T'] = T
        omega = 2 * sp.pi / T
        Results['omega'] = omega
        L0 = g * T ** 2 / (2 * sp.pi)

        omega2 = k ** 2 * a ** 2 * (2 * alpha1 ** 2 + 7) / 4 / (alpha1 - 1) ** 2
        omega4 = (a ** 4 * k ** 4 * (20 * alpha1 ** 5 + 112 * alpha1 ** 4
                                    - 100 * alpha1 ** 3 - 68 * alpha1 ** 2 - 211 * alpha1 + 328)) / (32 * (alpha1 - 1) ** 5)
        omegaFun = omega0 * (1 + omega2 + omega4) - omega
        Hfun1 = H - (2 * a + 2 * (B31 + B33) * k ** 2 * a ** 3 + 2 * (B51 + B53 + B55) * k ** 4 * a ** 5) == 0
        S = sp.solve([Hfun1, omegaFun], [k, a])
        kResults = [sp.Abs(k_val.evalf()) for k_val in S[0]]
        aResults = [a_val.evalf() for a_val in S[1]]

        Results['a'] = aResults
        sigma = sp.tanh(kResults[0] * h)
        B31 = (3 + 8 * sigma ** 2 - 9 * sigma ** 4) / 16 / sigma ** 4
        B51 = (121 * sp.cosh(2 * h * kResults[0]) ** 5 + 263 * sp.cosh(2 * h * kResults[0]) ** 4 +
               376 * sp.cosh(2 * h * kResults[0]) ** 3 - 1999 * sp.cosh(2 * h * kResults[0]) ** 2 +
               2509 * sp.cosh(2 * h * kResults[0]) - 1108) / (192 * (sp.cosh(2 * h * kResults[0]) - 1) ** 5)
        Results['aw'] = aResults * (1 + kResults[0] ** 2 * aResults ** 2 * B31 + kResults[0] ** 4 * aResults ** 4 * B51)
        Results['k'] = kResults[0]
        Results['omega0'] = sp.sqrt(g * kResults[0] * sp.tanh(kResults[0] * h))

    elif a is not None and omega0 is not None:
        Results['omega0'] = omega0
        kResults = sp.solve(omega0 ** 2 - g * k * sp.tanh(k * h), k)
        kResults = [k_val.evalf() for k_val in kResults]
        alpha1 = sp.cosh(2 * kResults[0] * h)

        omega2 = kResults[0] ** 2 * a ** 2 * (2 * alpha1 ** 2 + 7) / 4 / (alpha1 - 1) ** 2
        omega4 = (a ** 4 * kResults[0] ** 4 * (20 * alpha1 ** 5 + 112 * alpha1 ** 4
                                             - 100 * alpha1 ** 3 - 68 * alpha1 ** 2 - 211 * alpha1 + 328)) / (32 * (alpha1 - 1) ** 5
        Results['a'] = a
        Results['aw'] = a * (1 + kResults[0] ** 2 * a ** 2 * B31 + kResults[0] ** 4 * a ** 4 * B51)
        Hout = 2 * a + 2 * (B31 + B33) * kResults[0] ** 2 * a ** 3 + 2 * (B51 + B53 + B55) * kResults[0] ** 4 * a ** 5
        Results['H'] = Hout

    elif H is not None and omega0 is not None:
        Results['omega0'] = omega0
        kResults = sp.solve(omega0 ** 2 - g * k * sp.tanh(k * h), k)
        kResults = [k_val.evalf() for k_val in kResults]
        alpha1 = sp.cosh(2 * kResults[0] * h)

        omega2 = kResults[0] ** 2 * a ** 2 * (2 * alpha1 ** 2 + 7) / 4 / (alpha1 - 1) ** 2
        omega4 = (a ** 4 * kResults[0] ** 4 * (20 * alpha1 ** 5 + 112 * alpha1 ** 4
                                             - 100 * alpha1 ** 3 - 68 * alpha1 ** 2 - 211 * alpha1 + 328)) / (32 * (alpha1 - 1) ** 5)
        Hfun1 = H - (2 * a + 2 * (B31 + B33) * k ** 2 * a ** 3 + 2 * (B51 + B53 + B55) * k ** 4 * a ** 5) == 0
        S = sp.solve([Hfun1, omega0 ** 2 - g * k * sp.tanh(k * h)], [k, a])
        kResults = [sp.Abs(k_val.evalf()) for k_val in S[0]]
        aResults = [a_val.evalf() for a_val in S[1]]

        Results['a'] = aResults
        Results['aw'] = aResults * (1 + kResults[0] ** 2 * aResults ** 2 * B31 + kResults[0] ** 4 * aResults ** 4 * B51)
        Results['k'] = kResults[0]
        Results['H'] = H

    if modeNo is not None:
        Results['mode'] = modeNo
        Results['L'] = modeNo * sp.pi / kResults[0]

    return Results

# Example usage:
Results = StokesDispSolver(h=30, a=2.5, T=5, mode=1)
print(Results)
