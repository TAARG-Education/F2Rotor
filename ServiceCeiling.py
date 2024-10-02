


import numpy as np

def ServiceCeilingTest(Vcmax,h):
    Vcmax_theo = 0
    Vcmax_prac = 0.508 # = 100 ft/min
    index_theo = np.argmin(np.abs(Vcmax - Vcmax_theo))
    index_prac = np.argmin(np.abs(Vcmax - Vcmax_prac))

    SC_theo = h[index_theo]
    SC_prac = h[index_prac]

    return SC_theo, SC_prac

    