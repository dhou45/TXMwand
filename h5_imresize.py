import os, sys
import h5py, tifffile
from pathlib import Path
import numpy as np
import scipy.ndimage as spn

#change h5 file, make a copy first!
fn ='/Users/dhou45/Desktop/Rod1_250.h5'
f = h5py.File(fn, 'r+')
imgs_all = f['/XANES_results/whiteline']

print(imgs_all.shape)
rebin_factor = 2
rebin_imgs = spn.zoom(imgs_all, (1/rebin_factor, 1/rebin_factor, 1/rebin_factor))
print(rebin_imgs.shape)


del f['/XANES_results/whiteline']
dset = f.create_dataset('/XANES_results/whiteline', data=rebin_imgs)
f.close()
print('rebin data saved')