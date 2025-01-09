from typing import Callable, Literal, Union

import numpy as np
from scipy import integrate


def func_zeros(
    a: float, b: float, func: Callable[[float], float], delta: float = 1e-7
) -> float:
    """The root(zero)-finding of f(x) an interval (a, b) by the bisection method

    Args:
        a (float): left border
        b (float): right border
        func (Callable): target function f(x)
        delta (float, optional): the finding accuracy of x0 (f(x0)=0). Defaults to 1e-7.

    Raises:
        ValueError: Ends have the same sign

    Returns:
        float: x0, which is the solution f(x0) = 0
    """
    while abs(b - a) > delta:
        y1 = func(a)
        y2 = func(0.5 * (a + b))
        y3 = func(b)
        if y1 * y3 > 0:
            raise ValueError(
                f"Ends have the same sign \\ func({a})={func(a)} and func({b})={func(b)}"
            )
        if y1 * y2 < 0:
            b = 0.5 * (a + b)
        else:
            a = 0.5 * (a + b)
    return (a + b) * 0.5


def KH(level: float, delta: float = 1e-7) -> float:
    """A helper function calculating K*H,
    where K = pi * eta / (2 * N), H is the brush thikness

    Args:
        level (float): Value of combination pi * eta * sigma
        delta (float, optional): Accuracy of the K*H calculation. Defaults to 1e-7.

    Returns:
        float: Value of K*H
    """
    f = (
        lambda x: 2 * x
        - np.sin(x) * np.cos(x)
        - np.cos(x) ** 3 * np.arctanh(np.sin(x))
        - level
    )
    return func_zeros(a=delta, b=np.pi / 2 - delta, func=f, delta=delta)


class BrushModel:
    """
    Brush statistical properties under good solvent conditon
    """

    def __init__(
        self,
        N: int,
        sigma: float,
        eta: float,
        model: Literal["gauss", "bcc"],
        delta: float = 1e-7,
    ) -> None:
        """_summary_

        Args:
            N (int): Degry of polymerization
            sigma (float): grafting density
            eta (float): topological ratio k(m,0)/k(m,n)
            model (Literal["gauss","bcc"]): type of the chain stretchig (Gaussian or on the BCC-lattice)
            delta (float, optional): _description_. Defaults to 1e-7.
        """
        # Initial parameters
        self.N = N
        self.sigma = sigma
        self.eta = eta
        self.model = model
        # Parameters for calculations
        self.K = 0.0
        self.KH = 0.0
        self.H = 0.0
        self.recalc()

    def recalc(self, delta: float = 1e-7) -> None:
        """Recalculation of operating parameters (use when changing input parameters)
        Args:
            delta (float, optional): Accuracy of the K*H calculation. Defaults to 1e-7.
        Raises:
            ValueError: N < 1
            ValueError: sigma < 0.0 or sigma > 1.0
            ValueError: eta < 1.0
            ValueError: Unknown type of model
        """
        # Input parameters validation
        if self.N < 1:
            raise ValueError("N < 1")
        if self.sigma < 0.0 or self.sigma > 1.0:
            raise ValueError("sigma < 0.0 or sigma > 1.0")
        if self.eta < 1.0:
            raise ValueError("eta < 1.0")
        # Operating parameters recaclulation
        self.K = 0.5 * np.pi * self.eta / self.N
        if self.model == "gauss":
            self.H = (
                (0.5 * np.pi * self.eta) ** (-2 / 3) * self.sigma ** (1 / 3) * self.N
            )
        elif self.model == "bcc":
            self.KH = KH(np.pi * self.eta * self.sigma, delta=delta)
            self.H = self.KH / self.K
        else:
            raise ValueError("Unknown type of model")

    def U(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Calculation of potential

        Args:
            z (Union[float, np.ndarray]): Distance to the gafting surface

        Returns:
            Union[float, np.ndarray]: potential U(z)
        """
        if self.model == "gauss":
            return 1.5 * self.K**2 * (self.H**2 - z**2)
        if self.model == "bcc":
            return 3.0 * np.log(np.cos(self.K * z) / np.cos(self.K * self.H))

    def phi(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Calculation of volume fraction

        Args:
            z (Union[float, np.ndarray]): Distance to the gafting surface

        Returns:
            Union[float, np.ndarray]: Volume fraction phi(z)
        """
        if self.model == "gauss":
            return 1.5 * self.K**2 * (self.H**2 - z**2)
        if self.model == "bcc":
            return 1.0 - (np.cos(self.K * self.H) / np.cos(self.K * z)) ** 3

    def g(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Calculation of the ends distribution

        Args:
            z (Union[float, np.ndarray]): Distance to the gafting surface

        Returns:
            Union[float, np.ndarray]: The ends distribution g(z)
        """
        ratio = np.cos(self.K * z) / np.cos(self.KH)
        if self.model == "gauss":
            return (
                3.0 * self.K**2 * z / (self.N * self.sigma) * np.sqrt(self.H**2 - z**2)
            )
        if self.model == "bcc":
            return (
                3.0
                * np.tan(self.K * z)
                / (2.0 * self.N * self.sigma * ratio**3)
                * (ratio * np.sqrt(ratio**2 - 1.0) + np.arccosh(ratio))
            )

    def first_moment(self, distrib_name: Literal["g", "phi"]) -> float:
        """Calculation of the first moment of a given distribution (g(z) or phi(z))

        Args:
            distrib_name (Literal["g","phi"]): Type of distribution

        Raises:
            ValueError: Unknown type of distrib_name

        Returns:
            float: First moment
        """
        if distrib_name == "g":
            f = lambda x: self.g(x) * x
            return integrate.quad(func=f, a=0.0, b=self.H)[0]
        elif distrib_name == "phi":
            f = lambda x: self.phi(x) * x
            return integrate.quad(func=f, a=0.0, b=self.H)[0] / (self.N * self.sigma)
        else:
            raise ValueError("Unknown type of distrib_name")
