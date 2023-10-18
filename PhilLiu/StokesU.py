import numpy as np

def stokes_u(k, h, a, theta, z):
    g = 9.81
    omega0 = np.sqrt(g * k * np.tanh(k * h))
    sigma = np.tanh(k * h)
    alpha = np.cosh(2 * k * h)
    A11 = 1 / np.sinh(k * h)
    A22 = 3 / (8 * (np.sinh(k * h) ** 4))
    A33 = ((9 - 4 * (np.sinh(k * h) ** 2)) / (64 * (np.sinh(k * h) ** 7))
    A44 = (10 * alpha**3 - 174 * alpha**2 - 291 * alpha + 278) / (48 * (3 * alpha + 2) * (alpha - 1)**5)
    A55 = (-6 * alpha**5 + 272 * alpha**4 - 1552 * alpha**3 + 852 * alpha**2 + 2029 * alpha + 430) / (
            64 * (alpha - 1)**6 * (12 * alpha**2 + 11 * alpha + 2) * np.sinh(k * h))

    A31 = 0
    A51 = 0
    A42 = (12 * alpha**4 + 22 * alpha**3 - 84 * alpha**2 - 135 * alpha + 104) / (24 * (alpha - 1)**5)
    A53 = (8 * alpha**6 + 138 * alpha**5 + 384 * alpha**4 - 568 * alpha**3 - 2388 * alpha**2 + 237 * alpha + 974) / (
            64 * (3 * alpha + 2) * (alpha - 1)**6 * np.sinh(k * h))
    C2 = (sigma**2 - 1) / (4 * sigma)
    C4 = -9 / (4 * (alpha - 1)**3 * np.sinh(2 * h * k))

    phi1 = A11 * omega0 * a / k * np.cosh(k * (z + h)) * np.sin(theta)
    phi2 = A22 * omega0 * a**2 * np.cosh(2 * k * (z + h)) * np.sin(2 * theta)
    phi3 = A31 * omega0 * k * a**3 * np.cosh(k * (z + h)) * np.sin(theta) + A33 * omega0 * k * a**3 * np.cosh(
        3 * k * (z + h)) * np.sin(3 * theta)
    phi4 = A42 * omega0 * k**2 * a**4 * np.cosh(2 * k * (z + h)) * np.sin(2 * theta) + A44 * omega0 * k**2 * a**4 * np.cosh(
        4 * k * (z + h)) * np.sin(4 * theta)
    phi5 = A51 * omega0 * k**3 * a**5 * np.cosh(k * (z + h)) * np.sin(theta) + A53 * omega0 * k**3 * a**5 * np.cosh(
        3 * k * (z + h)) * np.sin(3 * theta) + A55 * a**5 * k**3 * omega0 * np.cosh(5 * k * (z + h)) * np.sin(5 * theta)
    phi = (phi1 + phi2 + phi3 + phi4 + phi5 + omega0**2 * a / k * sigma * (k * a * C2 + k**3 * a**3 * C4))

    u1 = A11 * omega0 * a * np.cosh(k * (z + h)) * np.cos(theta)
    u2 = 2 * A22 * omega0 * k * a**2 * np.cosh(2 * k * (z + h)) * np.cos(2 * theta)
    u3 = A31 * omega0 * k**2 * a**3 * np.cosh(k * (z + h)) * np.cos(theta) + 3 * A33 * omega0 * k**2 * a**3 * np.cosh(
        3 * k * (z + h)) * np.cos(3 * theta)
    u4 = 2 * A42 * omega0 * k**3 * a**4 * np.cosh(2 * k * (z + h)) * np.cos(2 * theta) + 4 * A44 * omega0 * k**3 * a**4 * np.cosh(
        4 * k * (z + h)) * np.cos(4 * theta)
    u5 = A51 * omega0 * k**4 * a**5 * np.cosh(k * (z + h)) * np.cos(theta) + 3 * A53 * omega0 * k**4 * a**5 * np.cosh(
        3 * k * (z + h)) * np.cos(3 * theta) + 5 * A55 * k**4 * a**5 * omega0 * np.cosh(5 * k * (z + h)) * np.cos(5 * theta)
    u1 = u1
    u2 = u2
    u3 = u3
    u4 = u4
    u5 = u5
    u = u1 + u2 + u3 + u4 + u5

    w1 = A11 * omega0 * a * np.sinh(k * (z + h)) * np.sin(theta)
    w2 = 2 * A22 * omega0 * k * a**2 * np.sinh(2 * k * (z + h)) * np.sin(2 * theta)
    w3 = A31 * omega0 * k**2 * a**3 * np.sinh(k * (z + h)) * np.sin(theta) + 3 * A33 * omega0 * k**2 * a**3 * np.sinh(
        3 * k * (z + h)) * np.sin(3 * theta)
    w4 = 2 * A42 * omega0 * k**3 * a**4 * np.sinh(2 * k * (z + h)) * np.sin(2 * theta) + 4 * A44 * omega0 * k**3 * a**4 * np.sinh(
        4 * k * (z + h)) * np.sin(4 * theta)
    w5 = A51 * omega0 * k**4 * a**5 * np.sinh(k * (z + h)) * np.sin(theta) + 3 * A53 * omega0 * k**4 * a**5 * np.sinh(
        3 * k * (z + h)) * np.sin(3 * theta) + 5 * A55 * a**5 * k**4 * omega0 * np.sinh(5 * k * (z + h)) * np.sin(5 * theta)
    w1 = w1
    w2 = w2
    w3 = w3
    w4 = w4
    w5 = w5
    w = w1 + w2 + w3 + w4 + w5

    return phi, u, u1, u2, u3, u4, u5, w, w1, w2, w3, w4, w5
