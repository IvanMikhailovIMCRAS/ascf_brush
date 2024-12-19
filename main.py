from asf_brush import BCC_model
import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":
    BCC = BCC_model(503*9, 0.01, 3)
    
    z = np.linspace(0.0, BCC.H-1, 200)
    print(BCC.H)
    plt.plot(z, BCC.g(z))
    plt.show()
    
    
    #     plt.plot(z, BCC.der_phi(z))
    
    # zz = (z[1:] + z[:-1])/2
    # plt.plot(zz, np.diff(BCC.phi(z))/np.diff(z))
    