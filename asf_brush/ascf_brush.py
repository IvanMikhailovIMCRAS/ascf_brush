import numpy as np


# Function for zeros finding
def func_zeros(a, b, func, delta=1e-7):
    if func(a) * func(b) > 0:
        print(a, func(a), b, func(b))
        raise ValueError("ends have the same sign")
    while abs(b - a) > delta:
        y1 = func(a)
        y2 = func(0.5 * (a + b))
        y3 = func(b)

        if y1 * y2 < 0:
            b = 0.5 * (a + b)
        else:
            a = 0.5 * (a + b)
    return (a + b) * 0.5


def KH(level, delta=1e-7):
    f = (
        lambda x: 2 * x
        - np.sin(2 * x) * 0.5
        - np.cos(x) ** 3 * np.log(abs(np.tan(x / 2 + np.pi / 4)))
        - level
    )
    return func_zeros(delta, np.pi / 2, f)


# Secant calculation
def sec(x):
    return 1 / np.cos(x)


class BCC_model:
    """
    Brush statistical properties under good solvent conditon on BCC lattice
    """

    def __init__(self, N: int, sigma: float, eta: float, min_val:float=1e-7) -> None:
        # Initial parameters
        self.N = N
        self.sigma = sigma
        self.eta = eta
        self.min_val = min_val

        # Errors proccesing
        if N < 1:
            raise ValueError("N < 1")
        if sigma < 0.0 or sigma > 1.0:
            raise ValueError("sigma < 0.0 or sigma > 1.0")
        if eta < 1.0:
            raise ValueError("eta < 1.0")
        
        # Parameters for calculations
        self.K = np.pi * eta / 2 / N
        self.H = KH(np.pi * eta * sigma) / self.K
        self.integral = (
            lambda z: sec(self.K*self.H)
            * sec(self.K*z)**6
            * (-1 * (np.cos(self.K * z)**2 - np.cos(self.K * self.H)**2)**0.5)
            * (
                3
                * np.cos(self.K * self.H)**4
                * (
                    1
                    / 8
                    * sec(self.K*self.H) ** 4
                    * np.cos(self.K * z) ** 4
                    * (5 - 3 * (1 - np.cos(self.K*self.H) ** 2 * sec(self.K*z) ** 2))
                    + 3
                    / np.tanh((1 - np.cos(self.K*self.H) ** 2 * sec(self.K*z) ** 2)
                    ** 0.5)
                    / 8
                    / (1 - np.cos(self.K*self.H) ** 2 * sec(self.K*z) ** 2) ** 0.5
                )
                - np.cos(self.K * z) ** 4
            )
        )

    # Theoretical volume fraction profile
    def phi(self, z):
        return 1 - (np.cos(self.K * self.H) / np.cos(self.K * z)) ** 3
    
    def der_phi(self, z):
        return -3*self.K*np.cos(self.K*self.H)**3 / np.cos(self.K*z)**3 * np.tan(self.K*z)
    
    # Theoretical ends distribution
    def g(self, z):
        if z[-1] > self.H - self.min_val:
            raise ValueError("z can`t be equal to H")
        return (
            1
            / self.sigma
            / 2
            / self.N
            * np.sin(2*self.K*z)
            * (
                self.phi(z)
                / np.cos(self.K * self.H)
                / (np.cos(self.K * z)**2 - np.cos(self.K * self.H)**2)**0.5
                - self.integral(z)
            )
        )   
