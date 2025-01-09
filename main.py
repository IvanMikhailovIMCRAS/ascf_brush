import matplotlib.pyplot as plt
import numpy as np

from ascf_brush import BrushModel

if __name__ == "__main__":
    m = 5
    n = 75
    p = 100
    N = p * (n + m)
    eta = (1.0 + n / m) ** 0.5
    factor = 0.1
    sigma = factor / eta**2

    data = np.loadtxt(f"data/data_m_{m}_n_{n}_p_{p}_factor_{factor}.txt").T
    plt.plot(data[0], data[2], label="simulation")

    BM = BrushModel(N, sigma, eta, model="gauss")
    z = np.linspace(0.0, BM.H, 100)
    plt.plot(z, BM.g(z), label="gauss")
    BM.model = "bcc"
    BM.recalc()
    z = np.linspace(0.0, BM.H, 100)
    plt.xlabel("$z/a$")
    plt.ylabel("$g(z)$")
    plt.title(f"N={N}, f={factor}, $\eta$={eta}")
    plt.plot(z, BM.g(z), label="bcc")
    plt.legend()
    plt.show()
