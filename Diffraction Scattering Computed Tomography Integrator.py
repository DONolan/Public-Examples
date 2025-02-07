#General Libraries
import os
import time
import glob
import copy
import tomopy
import pybaselines
import numpy as np
import pandas as pd
import xarray as xr
import PIL
from PIL import Image
import fabio
from scipy.ndimage import median_filter
#Matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import subplots
#pyFAI
import pyFAI
from pyFAI.gui import jupyter
import fabio
from pyFAI.test.utilstest import UtilsTest
from pyFAI.calibrant import CALIBRANT_FACTORY
from pyFAI.goniometer import SingleGeometry
from pyFAI.gui.jupyter.calib import Calibration
print(f'Using pyFAI version: {pyFAI.version}')

class DSCT:
    #the init method
    def __init__(self, name):
        self.name=name
        
    #1st instance variable 
    def set_data(self, path):
        self.path = path
        
    def get_data(self):
        return([os.path.basename(x) for x in sorted(glob.glob(self.path+"*.nc"))])
        
    
    #2nd instance variable
    def config(self, calibrant, mask):
        self.calibrant = pyFAI.load(calibrant)
        self.mask = fabio.open(mask).data
        return(self.calibrant)

    #3rd instance variable
    def integrate(self, 
                  radial_range, #pass as 2-list e.g. [0,24]
                  delta_q, #stepsize e.g. 0.0025
                  data
                 ):
        
        self.radial_range=radial_range
        self.delta_q=delta_q
        q = []
        dsct = []
        phi = []
        x = []

        for i in [self.path + x for x in data]:
            with xr.open_dataset(i) as ds:
                ds.load()
            
            print(ds.attrs['save_name'])    
            phi_value = ds.attrs['mPhi']
            phi.append(phi_value)
            
            x_value = ds['mBaseX'].data
            x.append(x_value)
            
            intensities = []
            
            npt = int(np.ceil((self.radial_range[1] - self.radial_range[0]) / self.delta_q ))
            radial_range = [self.radial_range[0],self.radial_range[0]+self.delta_q*npt]
            
            for j in x_value:
                img = (ds.pe1_imgs.sel(mBaseX=j,method='nearest').astype('float32') - ds.pe1_img_dark.astype('float32')).values
                img = median_filter(img, size=2)
                
                i2d = self.calibrant.integrate2d(data=img, npt_rad=npt, 
                                                 npt_azim=360,correctSolidAngle=True, 
                                                 radial_range=self.radial_range,
                                                 mask=self.mask,dummy=np.NaN,
                                                 method='bbox',unit='q_A^-1',
                                                 safe=True, normalization_factor=1.0
                                                )
                da2d = xr.DataArray(data=i2d.intensity.astype('float32'),
                                    coords=[i2d.azimuthal.astype('float32'),i2d.radial.astype('float32')],
                                    dims=['azimuthal','radial'])
                da1d = da2d.mean(dim='azimuthal')
                q.append(da1d['radial'].data.tolist())
                intensities.append(da1d.data.tolist())

            intensities = [z for _, z in sorted(zip(x_value, intensities))]
            dsct.append(intensities)
        
        q_space = np.array(q[0], dtype='f8')
        i_obs = np.array(dsct, dtype='f8')
        phi_obs = np.array(phi, dtype='f8')
        x_obs = np.array(x[0], dtype='f8')
        
        DSCT = xr.Dataset(coords=dict(Q = ("Q", q_space), #q-space for xrd data in inverse angstroms
                                      Phi = ("Phi", phi_obs), #rotational angle 
                                      X = ("X", x_obs) #translational x position
                                     ),
                          data_vars=dict(XRD_CT = (["Phi", "X", "Q"], i_obs), #XRD-CT values
                                        ),
                          attrs=dict(Q_coordinate = "Inverse Angstroms",
                                     Phi_coordinate = "rotational degrees",
                                     X_coordinate = "translation position",
                                     XRDCT_variable = "Q-space CT values"
                                    ),
                         )
                
        return(DSCT)