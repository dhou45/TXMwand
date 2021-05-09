import os, sys
import h5py, tifffile
from pathlib import Path
import numpy as np
from shutil import copyfile

# for 3DTXM with missing angles, especially for electrode samples. 
# check the total projection number or missing angle slice number.
raw_data_top_dir = "/run/media/VTLinlab/Lin_Lab_10/MR_3DTXM/Mn_Processed"
proj_num = 377
missing_wedge_s = 200
missing_wedge_e = 300
whiteline_dir = '/whiteline_trim'
Path(raw_data_top_dir+whiteline_dir).mkdir(parents=True, exist_ok=True)

# idx start, end+1
# scan_idx = np.arange(74533, 74546,1) #C4p5V
scan_idx = np.arange(75070, 75089,1) #D1p5V

for jj in range(len(scan_idx)):
	fn_raw = os.path.join(raw_data_top_dir,'fly_scan_id_'+str(scan_idx[jj])+'.h5')
	print(fn_raw)
	fn = os.path.join(raw_data_top_dir+whiteline_dir,'fly_scan_id_'+str(scan_idx[jj])+'.h5')
	copyfile(fn_raw, fn)

	f = h5py.File(fn, 'r+')
	angle = f['/angle']
	angle_temp = angle[0:proj_num]
	angle_new = np.concatenate((angle_temp[0:missing_wedge_s], angle_temp[missing_wedge_e:proj_num]))
	del f['/angle']
	dset = f.create_dataset('/angle', data=angle_new)

	img_tomo = f['/img_tomo']
	img_tomo_temp = img_tomo[0:proj_num,:,:]
	img_tomo_new = np.concatenate((img_tomo_temp[0:missing_wedge_s,:,:], img_tomo_temp[missing_wedge_e:proj_num,:,:]), axis=0)

	del f['/img_tomo']
	dset = f.create_dataset('/img_tomo', data=img_tomo_new)
	f.close()
	print('data trimed '+str(img_tomo_new.shape))


# # single file test
# fn = "/run/media/VTLinlab/Lin_Lab_10/BNL_FXI_Nov2020/NMC-Ni-D1p5V/whiteline/fly_scan_id_74117.h5"
# f = h5py.File(fn, 'r+')
# proj_num = 429
# missing_wedge_s = 200
# missing_wedge_e = 285

# angle = f['/angle']
# print(angle.shape)
# # angle_new = []

# angle_temp = angle[0:proj_num]
# angle_new = np.concatenate((angle_temp[0:missing_wedge_s], angle_temp[missing_wedge_e:proj_num]))
# print(angle_new.shape)

# del f['/angle']
# dset = f.create_dataset('/angle', data=angle_new)

# img_tomo = f['/img_tomo']
# print(img_tomo.shape)

# img_tomo_temp = img_tomo[0:proj_num,:,:]
# img_tomo_new = np.concatenate((img_tomo_temp[0:missing_wedge_s,:,:], img_tomo_temp[missing_wedge_e:proj_num,:,:]), axis=0)
# print(img_tomo_new.shape)

# del f['/img_tomo']
# dset = f.create_dataset('/img_tomo', data=img_tomo_new)
# f.close()
# print('data trimed '+str(img_tomo_new.shape))