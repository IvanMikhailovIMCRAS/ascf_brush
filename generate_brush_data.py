import numpy as np
from sfbox_api import Composition, Frame, Lat, Mol, Mon, Sys


def comb_brush(m: int, n: int, p: int, sigma: float) -> Frame:
    """_summary_

    Args:
        m (int): the spacer length
        n (int): the side chain length
        p (int): the number of sumermonomers (branching points)
        sigma (float): the grafting density

    Returns:
        Frame: sfbox_api framework for SCF calculation
        see: https://github.com/IvanMikhailovIMCRAS/sfbox_api/blob/main/sfbox_api/frame.py
    """

    # Input data validation
    if m < 5 or m % 2 == 0 or p < 3:
        raise ValueError(f"m={m}, p={p}, but must be (m >=5 and m%2 != 0 and p >= 3)")

    # Topological script
    comp = f"(X)1(A){m//2}[(A){n}](A){m//2}((A){m//2+1}[(A){n}](A){m//2}){p-2}(A){m//2+1}[(A){n}](A){m//2-1}(E)1"
    if n <= 0:
        comp = f"(X)1(A){m*p-2}(E)1"

    # Degree of polymerization
    N = Composition(comp).N
    theta = sigma * N

    # Framework building
    lat = Lat(
        **{
            "n_layers": p * m + 1,
            "geometry": "flat",
        }
    )

    mons = [
        Mon(**{"name": "X", "freedom": "pinned", "pinned_range": "1"}),
        Mon(**{"name": "A", "freedom": "free"}),
        Mon(**{"name": "E", "freedom": "free"}),
        Mon(**{"name": "W", "freedom": "free"}),
    ]

    mols = [
        Mol(**{"name": "Water", "composition": "(W)1", "freedom": "solvent"}),
        Mol(
            **{
                "name": "pol",
                "composition": comp,
                "freedom": "restricted",
                "theta": theta,
            }
        ),
    ]

    sys = Sys()

    frame = Frame(lat, sys, mols, mons)
    frame.text += "sys : name : overflow_protection : true"
    frame.chi = 0.0
    frame.theta = theta
    frame.run()
    return frame


if __name__ == "__main__":
    # Set structure of the grafted bottle-brush
    m = 5
    n = 75
    p = 100
    # topological ratio ( k(m,0)/k(m,n) )
    eta = (1 + n / m) ** 0.5
    print(f"eta={eta}, m={m}, n={n}, p={p}")
    # maximum grafting density for a given structure
    sigma_max = m / (n + m)

    for factor in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        sigma = sigma_max * factor
        print(f"sigma={sigma}, sigma_max={sigma_max}, factor={factor}")
        frame = comb_brush(m, n, p, sigma)
        z = frame.profile["layer"][1:-1]
        phi = frame.profile["pol"][1:-1]
        g = frame.profile["E"][1:-1]
        g = g / np.sum(g)
        data = np.vstack([z, phi, g]).T
        np.savetxt(f"data/data_m_{m}_n_{n}_p_{p}_factor_{factor}.txt", data)
