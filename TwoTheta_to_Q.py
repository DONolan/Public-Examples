import numpy as np

def twotheta_to_q(twotheta, wavelength):
    q=((4*((np.pi))*(np.sin(np.deg2rad(twotheta/2))))/wavelength)
    return(q)