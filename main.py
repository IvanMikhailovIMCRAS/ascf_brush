from asf_brush import BCC_model
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


if __name__ == "__main__":
    BCC = BCC_model(1006, 0.2, 2)

    data = pd.read_csv(f"P_50_n_15_m_5_0.2.txt", delim_whitespace=True, skiprows=1)

    x = data.iloc[:, 0]
    # phi = data.iloc[:, 1]
    end = data.iloc[:, 3]

    print(BCC.H)
    z = np.linspace(0.0, BCC.H - 2e-7, 100)
    plt.plot(x, end, "-o")
    plt.plot(z, BCC.g(z))
    plt.show()

    # plt.plot(z, BCC.der_phi(z))

    # zz = (z[1:] + z[:-1])/2
    # plt.plot(zz, np.diff(BCC.phi(z))/np.diff(z))
