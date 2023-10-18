import sympy as sp

def StokesDispSolver(**kwargs):
    # Extract input parameters
    h = kwargs.get('h', 1)
    T = kwargs.get('T', None)
    omega0 = kwargs.get('omega0', None)
    a = kwargs.get('a', None)
    H = kwargs.get('H', None)
    modeNo = kwargs.get('mode', 1)

    Results = {}
    Results['h'] = h
    # Define some symbolic variables
    k = sp.symbols('k')
    g = 9.81
    alpha1 = sp.cosh(2 * k * h)
    sigma = sp.tanh(k * h)

    # Check which input parameters are provided and calculate the rest
    if a is not None and T is not None:
        omega0 = sp.sqrt(g * k * sp.tanh(k * h))
        Results['T'] = T
        omega = 2 * sp.pi / T
        Results['omega'] = omega
        L0 = g * T**2 / (2 * sp.pi)

        if modeNo == 1:
            omega2 = k**2 * a**2 * (2 * alpha1**2 + 7) / 4 / (alpha1 - 1)**2
            omega4 = (a**4 * k**4 * (20 * alpha1**5 + 112 * alpha1**4 - 100 * alpha1**3 - 68 * alpha1**2 - 211 * alpha1 + 328)) / (32 * (alpha1 - 1)**5)
            omegaFun = omega0 * (1 + omega2 + omega4) - omega
            kResults = [k.evalf() for k in sp.solve(omegaFun, k, domain=sp.S.Reals)]
            Results['a'] = a
            sigma = sp.tanh(kResults[0] * h)
            B31 = (3 + 8 * sigma**2 - 9 * sigma**4) / (16 * sigma**4)
            B51 = (121 * sp.cosh(2 * h * kResults[0])**5 + 263 * sp.cosh(2 * h * kResults[0])**4 + 376 * sp.cosh(2 * h * kResults[0])**3
                - 1999 * sp.cosh(2 * h * kResults[0])**2 + 2509 * sp.cosh(2 * h * kResults[0]) - 1108) / (192 * (sp.cosh(2 * h * kResults[0]) - 1)**5)
            Results['aw'] = a * (1 + kResults[0]**2 * a**2 * B31 + kResults[0]**4 * a**4 * B51)

    elif H is not None and T is not None:
        omega0 = sp.sqrt(g * k * sp.tanh(k * h))
        Results['T'] = T
        omega = 2 * sp.pi / T
        Results['omega'] = omega
        L0 = g * T**2 / (2 * sp.pi)

        if modeNo == 1:
            B33 = (27 - 9 * sigma**2 + 9 * sigma**4 - 3 * sigma**6) / (64 * sigma**6)
            B55 = (300 * sp.cosh(2 * h * k)**8 + 1579 * sp.cosh(2 * h * k)**7 + 3176 * sp.cosh(2 * h * k)**6 + 2949 * sp.cosh(2 * h * k)**5
                + 1188 * sp.cosh(2 * h * k)**4 + 675 * sp.cosh(2 * h * k)**3 + 1326 * sp.cosh(2 * h * k)**2 + 827 * sp.cosh(2 * h * k)
                + 130) * (5 / ((sp.cosh(2 * h * k) - 1)**6 * (12 * sp.cosh(2 * h * k)**2 + 11 * sp.cosh(2 * h * k) + 2) * 384))
            B31 = (3 + 8 * sigma**2 - 9 * sigma**4) / (16 * sigma**4)
            B51 = (121 * sp.cosh(2 * h * k)**5 + 263 * sp.cosh(2 * h * k)**4 + 376 * sp.cosh(2 * h * k)**3
                - 1999 * sp.cosh(2 * h * k)**2 + 2509 * sp.cosh(2 * h * k) - 1108) / (192 * (sp.cosh(2 * h * k) - 1)**5)
            Results['H'] = H
            Hfun1 = H - (2 * a + 2 * (B31 + B33) * k**2 * a**3 + 2 * (B51 + B55) * k**4 * a**5)
            omega2 = k**2 * a**2 * (2 * alpha1**2 + 7) / 4 / (alpha1 - 1)**2
            omega4 = (a**4 * k**4 * (20 * alpha1**5 + 112 * alpha1**4 - 100 * alpha1**3 - 68 * alpha1**2 - 211 * alpha1 + 328)) / (32 * (alpha1 - 1)**5)
            omegaFun1 = omega0 * (1 + omega2 + omega4) - omega
            S = sp.solve([Hfun1, omegaFun1], [k, a], [2 * sp.pi / L0, H / 2])
            kResults = abs(S[0][0].evalf())
            aResults = S[0][1].evalf()
            Results['a'] = abs(aResults)
            sigma = sp.tanh(kResults * h)
            B31 = (3 + 8 * sigma**2 - 9 * sigma**4) / (16 * sigma**4)
            B51 = (121 * sp.cosh(2 * h * kResults)**5 + 263 * sp.cosh(2 * h * kResults)**4 + 376 * sp.cosh(2 * h * kResults)**3
                - 1999 * sp.cosh(2 * h * kResults)**2 + 2509 * sp.cosh(2 * h * kResults) - 1108) / (192 * (sp.cosh(2 * h * kResults) - 1)**5)
            Results['aw'] = aResults * (1 + kResults**2 * aResults**2 * B31 + kResults**4 * aResults**4 * B51)
            Results['omega0'] = sp.sqrt(g * kResults * sp.tanh(kResults * h))

    elif a is not None and omega0 is not None:
        Results['omega0'] = omega0
        kResults = sp.solve(omega0**2 - g * k * sp.tanh(k * h), k, domain=sp.S.Reals)
        alpha1 = sp.cosh(2 * kResults[0] * h)
        sigma = sp.tanh(kResults[0] * h)

        B55 = (300 * sp.cosh(2 * h * kResults[0])**8 + 1579 * sp.cosh(2 * h * kResults[0])**7 + 3176 * sp.cosh(2 * h * kResults[0])**6 + 2949 * sp.cosh(2 * h * kResults[0])**5
            + 1188 * sp.cosh(2 * h * kResults[0])**4 + 675 * sp.cosh(2 * h * kResults[0])**3 + 1326 * sp.cosh(2 * h * kResults[0])**2 + 827 * sp.cosh(2 * h * kResults[0])
            + 130) * (5 / ((sp.cosh(2 * h * kResults[0]) - 1)**6 * (12 * sp.cosh(2 * h * kResults[0])**2 + 11 * sp.cosh(2 * h * kResults[0]) + 2) * 384))
        B33 = (27 - 9 * sigma**2 + 9 * sigma**4 - 3 * sigma**6) / (64 * sigma**6)

        switcher = {
            1: (a, a * (1 + kResults[0]**2 * a**2 * B31 + kResults[0]**4 * a**4 * B51)),
            2: (abs(S[0][1].evalf()), a)
        }
        a, aw = switcher.get(modeNo)

        Results['a'] = a
        Results['aw'] = aw
        Results['H'] = Hout
        omega2 = kResults[0]**2 * a**2 * (2 * alpha1**2 + 7) / 4 / (alpha1 - 1)**2
        omega4 = (a**4 * kResults[0]**4 * (20 * alpha1**5 + 112 * alpha1**4 - 100 * alpha1**3 - 68 * alpha1**2 - 211 * alpha1 + 328)) / (32 * (alpha1 - 1)**5)
        omega = omega0 * (1 + omega2 + omega4)
        Results['omega'] = omega
        Results['T'] = 2 * sp.pi / omega

    elif H is not None and omega0 is not None:
        Results['omega0'] = omega0
        kResults = sp.solve(omega0**2 - g * k * sp.tanh(k * h), k, domain=sp.S.Reals)
        alpha1 = sp.cosh(2 * kResults[0] * h)
        B33 = (27 - 9 * sigma**2 + 9 * sigma**4 - 3 * sigma**6) / (64 * sigma**6)
        B55 = (300 * sp.cosh(2 * h * kResults[0])**8 + 1579 * sp.cosh(2 * h * kResults[0])**7 + 3176 * sp.cosh(2 * h * kResults[0])**6 + 2949 * sp.cosh(2 * h * kResults[0])**5
            + 1188 * sp.cosh(2 * h * kResults[0])**4 + 675 * sp.cosh(2 * h * kResults[0])**3 + 1326 * sp.cosh(2 * h * kResults[0])**2 + 827 * sp.cosh(2 * h * kResults[0])
            + 130) * (5 / ((sp.cosh(2 * h * kResults[0]) - 1)**6 * (12 * sp.cosh(2 * h * kResults[0])**2 + 11 * sp.cosh(2 * h * kResults[0]) + 2) * 384))

        Results['H'] = H
        Hfun1 = H - (2 * a + 2 * (B31 + B33) * kResults[0]**2 * a**3 + 2 * (B51 + B55) * kResults[0]**4 * a**5)
        a = sp.solve(Hfun1, a, domain=sp.S.Reals)
        a = abs(a[0])
        omega2 = kResults[0]**2 * a**2 * (2 * alpha1**2 + 7) / 4 / (alpha1 - 1)**2
        omega4 = (a**4 * kResults[0]**4 * (20 * alpha1**5 + 112 * alpha1**4 - 100 * alpha1**3 - 68 * alpha1**2 - 211 * alpha1 + 328)) / (32 * (alpha1 - 1)**5)
        omega = omega0 * (1 + omega2 + omega4)
        Results['a'] = a
        Results['omega'] = omega
        Results['T'] = 2 * sp.pi / omega
        B31 = (3 + 8 * sigma**2 - 9 * sigma**4) / (16 * sigma**4)
        B51 = (121 * sp.cosh(2 * h * kResults[0])**5 + 263 * sp.cosh(2 * h * kResults[0])**4 + 376 * sp.cosh(2 * h * kResults[0])**3
            - 1999 * sp.cosh(2 * h * kResults[0])**2 + 2509 * sp.cosh(2 * h * kResults[0]) - 1108) / (192 * (sp.cosh(2 * h * kResults[0]) - 1)**5)
        Results['aw'] = a * (1 + kResults[0]**2 * a**2 * B31 + kResults[0]**4 * a**4 * B51)

    return Results
