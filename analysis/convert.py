import os, pickle, glob
from scipy.io import loadmat
import numpy as np
from xmitgcm import open_mdsdataset
from pathlib import Path
import re
from datainput import getIterNums

fname = "/data/jbod/gdf/MITgcm_CS/experiments/width16/nobathsaltedsf10150/results"
times = getIterNums(fname)
ds = open_mdsdataset(fname,prefix=["THETA","SALT","momKE","SHIfwFlx","VVEL","UVEL","WVEL","RHOAnoma","MXLDEPTH"],ignore_unknown_vars=True,iters=times)
ice = np.fromfile("../experiments/width16/nobathsaltedsf10150/input/SHELFICEtopoFile.bin",">f8")
ice = ice.reshape([ds.YC.shape[0],ds.XC.shape[0]])

bath = np.fromfile("../experiments/width16/nobathsaltedsf10150/input/bathyFile.bin",">f8")
bath = bath.reshape([ds.YC.shape[0],ds.XC.shape[0]])

ds['ice'] = (('YC', 'XC'), ice)
ds['bath'] = (('YC', 'XC'), bath)

ds.to_netcdf("nobathsaltedsf10150.nc")
