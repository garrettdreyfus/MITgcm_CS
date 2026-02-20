import os, pickle, glob
from scipy.io import loadmat
import numpy as np
from xmitgcm import open_mdsdataset
from pathlib import Path
import re
from datainput import getIterNums

fname = "/data/jbod/gdf/MITgcm_CS/experiments/varysf16/imroved-coldstart-warm-sf0/results"
times = getIterNums(fname)
ds = open_mdsdataset(fname,prefix=["THETA","SALT","momKE","SHIfwFlx","VVEL","UVEL","WVEL","UVEL_inst","VVEL_inst","WVEL_inst","RHOAnoma","MXLDEPTH"],ignore_unknown_vars=True,iters=times)
ice = np.fromfile("../experiments/varysf16/imroved-coldstart-warm-sf0/input/SHELFICEtopoFile.bin",">f8")
ice = ice.reshape([ds.YC.shape[0],ds.XC.shape[0]])

bath = np.fromfile("../experiments/varysf16/imroved-coldstart-warm-sf0/input/bathyFile.bin",">f8")
bath = bath.reshape([ds.YC.shape[0],ds.XC.shape[0]])

ds['ice'] = (('YC', 'XC'), ice)
ds['bath'] = (('YC', 'XC'), bath)

ds.to_netcdf("imroved-coldstart-warm-sf0.nc")
