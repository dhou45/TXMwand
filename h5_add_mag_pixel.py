import os, sys
import h5py, tifffile
from pathlib import Path
import numpy as np

#change h5 file, make a copy first!
fn ='/run/media/VTLinlab/Lin_Lab_5/dh_drop/fly_scan_id_27402.h5'
f = h5py.File(fn, 'r+')
Magnification = 331.3699296346655;
Pixel_Size = '19.615539669414886nm'
f.create_dataset('/Magnification', data=Magnification)
f.create_dataset('/Pixel Size', data=Pixel_Size)
f.close()
