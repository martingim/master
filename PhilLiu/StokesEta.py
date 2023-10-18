import numpy as np

def StokesEta(k, h, a, theta):
    sigma = np.tanh(k * h)
    alpha1 = np.cosh(2 * k * h)
    
    # First and second order
    eta1 = a * np.cos(theta)
    eta2 = (k * a**2 / 4) * (3 - sigma**2) / (sigma**3) * np.cos(2 * theta)
    
    # Dingemans corrected, 3rd order
    B31 = (3 + 8 * sigma**2 - 9 * sigma**4) / (16 * sigma**4)
    B33 = (27 - 9 * sigma**2 + 9 * sigma**4 - 3 * sigma**6) / (64 * sigma**6)
    eta3 = k**2 * a**3 * (B31 * np.cos(theta) + B33 * np.cos(3 * theta))
    
    # Our derived 4th order equation
    sigma1 = 24 * (3 * alpha1 + 2) * (alpha1 - 1)**4 * np.sinh(2 * k * h)
    eta4 = (k**3 * a**4 / sigma1) * (
        (60 * alpha1**6 + 232 * alpha1**5 - 118 * alpha1**4 - 989 * alpha1**3 - 607 * alpha1**2 +
         352 * alpha1 + 260) * np.cos(2 * theta) +
        (24 * alpha1**6 + 116 * alpha1**5 + 214 * alpha1**4 + 188 * alpha1**3 + 133 * alpha1**2 +
         101 * alpha1 + 34) * np.cos(4 * theta)
    )
    
    # Our derived 5th order equation
    eta5 = (k**4 * a**5 / (192 * (alpha1 - 1)**5)) * (
        (121 * alpha1**5 + 264 * alpha1**4 + 376 * alpha1**3 - 1999 * alpha1**2 + 2509 * alpha1 - 1108) * np.cos(theta) +
        (k**4 * a**5 / (128 * (alpha1 - 1)**6 * (3 * alpha1 + 2))) * 9 * (
            57 * alpha1**7 + 204 * alpha1**6 - 53 * alpha1**5 - 782 * alpha1**4 - 741 * alpha1**3 -
            52 * alpha1**2 + 371 * alpha1 + 186
        ) * np.cos(3 * theta) +
        (k**4 * a**5 / (384 * (alpha1 - 1)**6 * (12 * alpha1**2 + 11 * alpha1 + 2))) * 5 * (
            300 * alpha1**8 + 1579 * alpha1**7 + 3176 * alpha1**6 + 2949 * alpha1**5 +
            1188 * alpha1**4 + 675 * alpha1**3 + 1326 * alpha1**2 + 827 * alpha1 + 130
        ) * np.cos(5 * theta)
    )
    
    eta = eta1 + eta2 + eta3 + eta4 + eta5
    return eta
