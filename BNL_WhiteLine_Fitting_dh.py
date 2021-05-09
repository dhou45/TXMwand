#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'ipympl')

#original
import xanes_math as xm
import xanes_analysis as xa

import h5py, tifffile
from pathlib import Path
import os, sys
import matplotlib.pyplot as plt
import numpy as np
from importlib import reload
import skimage.morphology as skm
import scipy.ndimage as spn
from numpy import savetxt
xa = reload(xa)
plt.rcParams['figure.figsize'] = [14, 20]
print('\033[04m\033[01m\033[34m     Section 1 finished     ')


# In[2]:


# edge_offset from Co
edge_offset_2_Co = 8.333 - 8.333
# estimated edge energy
edge_eng = 8.333 + edge_offset_2_Co
# end poit of the pre-edge relative to the edge_eng in keV
pre_ee = -0.05
# start point of the post-edge relative to the edge_eng in keV
post_es = 0.1
# how many times of the edge jump magnitude should be compared to the pre-edge standard deviation
edge_jump_threshold = 3
# how much should the pre_edge be offset up for validating if the post edge trend is in a reasonable range
# this is a factor to the pre-edge deviation
pre_edge_threshold = 3.5
# define an energy range for 0.5 absorption postion fitting
ep_eng_s = 8.335 + edge_offset_2_Co
ep_eng_e = 8.350 + edge_offset_2_Co
# define an energy range for whiteline peak postion fitting
#wl_eng_s = 8.345 + edge_offset_2_Co + 0.000
#wl_eng_e = 8.355 + edge_offset_2_Co - 0.000

# define an energy range for edge_pos display
ep_vmin = 8.338 + edge_offset_2_Co
ep_vmax = 8.348 + edge_offset_2_Co
# define an energy range for whiteline display
wl_vmin = 8.340 + edge_offset_2_Co
wl_vmax = 8.360 + edge_offset_2_Co
# define path and file name to save xanes analysis results; if you use the default path and name as below,
# you don't need to change anything. otherwise, give your full path and file name below.
#out_fn = os.path.join(str(Path(fn_template).parent), 'xanes_analysis_' + str(Path(fn_template).name)).format(scan_id)
#print(out_fn)
print('\033[04m\033[01m\033[34m     Section 2 finished     ')


# In[4]:


fn = []
eng = []
fn = np.array(fn)
eng = np.array(eng)


#f_path ='/run/media/VTLinlab/Elements/Gravel2/Gravel2_250/'
f_path ='/run/media/VTLinlab/Elements/Gravel2/Gravel2_150/'
f_name = "Single_3D_trial_reg_scan_id_27465-27485_2020-04-26-17-10-56.h5"
fn = f_path+f_name

#change h5 file
f = h5py.File(fn, 'r')
eng = np.array(f['/registration_results/reg_results/eng_list'])
imgs_all = f['/registration_results/reg_results/registered_xanes3D']
ny = imgs_all.shape[2]
nx = imgs_all.shape[3]
nz = imgs_all.shape[1]
print(imgs_all.shape)

print('\033[04m\033[01m\033[34m     Section 3 finished     ')


# In[6]:


imgs = np.ndarray([imgs_all.shape[0], ny, nx])
mask = np.ndarray([imgs_all.shape[1], ny, nx])

## mask with threshold_ratio
#mask_threshold_ratio = 0.25
#mask_dilation = 3
#for ii in range(imgs_all.shape[1]):
#    imgs[idx] = imgs_all[idx, ii]
#    mask_threshold = imgs[idx].min().min()+ (imgs[idx].max().max()-imgs[idx].min().min())*mask_threshold_ratio
#    mask[ii] = skm.binary_dilation((spn.gaussian_filter(imgs[idx],mask_dilation) > mask_threshold).astype(np.uint8), np.ones([mask_dilation,mask_dilation])).astype(np.uint8)[:]

## mask with fixed threshold
mask_threshold = 0.002
mask_dilation = 3

idx = int(imgs_all.shape[0]/2)
for ii in range(imgs_all.shape[1]):
    imgs[idx] = imgs_all[idx, ii]
#    mask[ii] = skm.binary_dilation((imgs[idx] > mask_threshold).astype(np.uint8), np.ones([mask_dilation,mask_dilation])).astype(np.uint8)[:]
    mask[ii] = skm.binary_dilation((spn.gaussian_filter(imgs[idx],mask_dilation) > mask_threshold).astype(np.uint8), np.ones([mask_dilation,mask_dilation])).astype(np.uint8)[:]

## plot top, middle, bottom slices of the mask * image

ii = int(0.05*imgs_all.shape[1])
fig = plt.figure()
ax3 = fig.add_subplot(1, 3, 1)
ax3.set_title('image', fontdict={'fontsize':12})
ax3.imshow(imgs_all[idx, ii, :, :])
ax4 = fig.add_subplot(1, 3, 2)
ax4.set_title('mask', fontdict={'fontsize':12})
ax4.imshow(mask[ii, :, :])
ax5 = fig.add_subplot(1, 3, 3)
ax5.set_title('image x mask', fontdict={'fontsize':12})
ax5.imshow(imgs_all[idx, ii, :, :]*mask[ii, :, :])

ii = int(0.5*imgs_all.shape[1])
fig2 = plt.figure()
ax32 = fig2.add_subplot(1, 3, 1)
ax32.set_title('image', fontdict={'fontsize':12})
ax32.imshow(imgs_all[idx, ii, :, :])
ax42 = fig2.add_subplot(1, 3, 2)
ax42.set_title('mask', fontdict={'fontsize':12})
ax42.imshow(mask[ii, :, :])
ax52 = fig2.add_subplot(1, 3, 3)
ax52.set_title('image x mask', fontdict={'fontsize':12})
ax52.imshow(imgs_all[idx, ii, :, :]*mask[ii, :, :])

ii = int(0.95*imgs_all.shape[1])
fig3 = plt.figure()
ax33 = fig3.add_subplot(1, 3, 1)
ax33.set_title('image', fontdict={'fontsize':12})
ax33.imshow(imgs_all[idx, ii, :, :])
ax43 = fig3.add_subplot(1, 3, 2)
ax43.set_title('mask', fontdict={'fontsize':12})
ax43.imshow(mask[ii, :, :])
ax53 = fig3.add_subplot(1, 3, 3)
ax53.set_title('image x mask', fontdict={'fontsize':12})
ax53.imshow(imgs_all[idx, ii, :, :]*mask[ii, :, :])

plt.show()


# In[7]:


# plot xanes for slice ii at specific point (x,y)
ii = int(0.95*imgs_all.shape[1])

for jj in range(imgs_all.shape[0]):
    imgs[jj] = imgs_all[jj, ii]
    
x=450
y=320

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
# rebin 5x5
#ax1.plot(eng, imgs[:, y:y+5, x:x+5].mean(axis=(1,2))*mask[ii, y, x])
ax1.plot(eng, imgs[:, y, x]*mask[ii, y, x])
#ax1.xlim([8.347, 8.355])
plt.show()




# whiteline fitting of imgs, norebin, poly2, with mask
#change out file name
out_fn = f_path+"Whiteline_Fitting_norebin_mask_poly2_"+f_name
poly_order = 2

#wl_eng_s = 8.344 + edge_offset_2_Co + 0.000
#wl_eng_e = 8.353 + edge_offset_2_Co - 0.000
wl_eng_s = 7.724 + edge_offset_2_Co + 0.000
wl_eng_e = 7.733 + edge_offset_2_Co - 0.000

imgs = np.ndarray([imgs_all.shape[0], ny, nx])
xanes3d = np.ndarray([imgs_all.shape[1], ny, nx])
xana = xa.xanes_analysis(imgs, eng, edge_eng, pre_ee=pre_ee, post_es=post_es, 
                         edge_jump_threshold=edge_jump_threshold, pre_edge_threshold=pre_edge_threshold)

for ii in range(imgs_all.shape[1]):
    for jj in range(imgs_all.shape[0]):
        imgs[jj] = imgs_all[jj, ii]

    xana.spectrum[:] = imgs[:]
    xana.fit_whiteline(wl_eng_s, wl_eng_e, poly_order)
    xanes3d[ii] = xana.wl_pos[:]*mask[ii]
#    xanes3d[ii] = xana.whiteline_pos[:]
    
xana.save_results(out_fn, dtype='3D_XANES', **{'whiteline': xanes3d})
print('results saved')




# whiteline fitting of imgs, norebin, poly3, with mask
#change out file name
out_fn = f_path+"Whiteline_Fitting_norebin_mask_poly3_"+f_name
poly_order = 3

#wl_eng_s = 8.344 + edge_offset_2_Co + 0.000
#wl_eng_e = 8.353 + edge_offset_2_Co - 0.000
wl_eng_s = 7.724 + edge_offset_2_Co + 0.000
wl_eng_e = 7.733 + edge_offset_2_Co - 0.000

imgs = np.ndarray([imgs_all.shape[0], ny, nx])
xanes3d = np.ndarray([imgs_all.shape[1], ny, nx])
xana = xa.xanes_analysis(imgs, eng, edge_eng, pre_ee=pre_ee, post_es=post_es, 
                         edge_jump_threshold=edge_jump_threshold, pre_edge_threshold=pre_edge_threshold)

for ii in range(imgs_all.shape[1]):
    for jj in range(imgs_all.shape[0]):
        imgs[jj] = imgs_all[jj, ii]

    xana.spectrum[:] = imgs[:]
    xana.fit_whiteline(wl_eng_s, wl_eng_e, poly_order)
    xanes3d[ii] = xana.wl_pos[:]*mask[ii]
#    xanes3d[ii] = xana.whiteline_pos[:]
    
xana.save_results(out_fn, dtype='3D_XANES', **{'whiteline': xanes3d})
print('results saved')


# In[8]:


f.close()









# whiteline fitting of imgs, norebin, poly3, no mask
#change out file name
out_fn = "/run/media/VTLinlab/Elements/Gravel3/Gravel3_250/Whiteline_Fitting_norebin_nomask_poly3_3D_trial_reg_scan_id_27528-27548.h5"
poly_order = 3

wl_eng_s = 8.344 + edge_offset_2_Co + 0.000
wl_eng_e = 8.353 + edge_offset_2_Co - 0.000

imgs = np.ndarray([imgs_all.shape[0], ny, nx])
xanes3d = np.ndarray([imgs_all.shape[1], ny, nx])
xana = xa.xanes_analysis(imgs, eng, edge_eng, pre_ee=pre_ee, post_es=post_es, 
                         edge_jump_threshold=edge_jump_threshold, pre_edge_threshold=pre_edge_threshold)


for ii in range(imgs_all.shape[1]):
    for jj in range(imgs_all.shape[0]):
        imgs[jj] = imgs_all[jj, ii]

    xana.spectrum[:] = imgs[:]
    xana.fit_whiteline(wl_eng_s, wl_eng_e, poly_order)
#    xanes3d[ii] = xana.wl_pos[:]*mask[ii]
    xanes3d[ii] = xana.wl_pos[:]
    
xana.save_results(out_fn, dtype='3D_XANES', **{'whiteline': xanes3d})
print('results saved')









# test poly4 with manully input array, DH
y1 = np.array([0.00067914, 0.00050809, 0.00060443, 0.00095565, 0.00094147, 0.00120585, 0.00168194, 0.00216041, 0.00209089, 0.00239944, 0.00218799, 0.00214025, 0.0022635,  0.00155437, 0.00143575, 0.00170442, 0.00139527, 0.00134593, 0.00192241, 0.0014839 ] )
x1 = np.array([8.339998,  8.3410015, 8.3429985, 8.344002,  8.344997,  8.346002,  8.346997 ,8.3480015, 8.348997,  8.350002,  8.350999,  8.351996,  8.353002,  8.353999, 8.354997,  8.356003,  8.357001,  8.358,     8.358998,  8.359998 ])
z = np.polyfit(x1, y1, 4)
p = np.poly1d(z)
xp = np.linspace(8.34, 8.36, 100)

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1 = plt.plot(x1, y1, '.', xp, p(xp), '-')
plt.show()

dEFM_curve = np.polyder(p)
xvals_g = np.roots(dEFM_curve)
print(xvals_g)

# choose the xvals in correct x-range
xvals_g = xvals_g[(xvals_g >= 8.347) & (xvals_g <= 8.352)]
peak_positions = xvals_g;
print(peak_positions)

# Xiao's way to choose peak position
whiteline_pos = []
whiteline_pos = np.squeeze(xp[np.argmax(np.polyval(z, xp), axis=0)])
print(whiteline_pos)




# interp and gauss filter method for polyfit.
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import numpy as np

y1 = np.array([0.00067914, 0.00050809, 0.00060443, 0.00095565, 0.00094147, 0.00120585, 0.00168194, 0.00216041, 0.00209089, 0.00239944, 0.00218799, 0.00214025, 0.0022635,  0.00155437, 0.00143575, 0.00170442, 0.00139527, 0.00134593, 0.00192241, 0.0014839 ] )
x1 = np.array([8.339998,  8.3410015, 8.3429985, 8.344002,  8.344997,  8.346002,  8.346997 ,8.3480015, 8.348997,  8.350002,  8.350999,  8.351996,  8.353002,  8.353999, 8.354997,  8.356003,  8.357001,  8.358,     8.358998,  8.359998 ])

xvals = np.linspace(x1.min(), x1.max(), 200)
yinterp = np.interp(xvals, x1, y1)

x = xvals[50:150]
y = yinterp[50:150]
z = np.polyfit(x, y, 2)
p = np.poly1d(z)
xp = np.linspace(x.min(), x.max(), 200)

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1 = plt.plot(x, y, '.', xp, p(xp), '-')


y = gaussian_filter(y, sigma = 5)
z = np.polyfit(x, y, 2)
p = np.poly1d(z)

fig3 = plt.figure()
ax13 = fig3.add_subplot(1, 1, 1)
ax13 = plt.plot(x, y, '.', xp, p(xp), '-')


x = x1[5:15]
y = y1[5:15]
z = np.polyfit(x, y, 2)
p = np.poly1d(z)
xp = np.linspace(x.min(), x.max(), 200)

fig2 = plt.figure()
ax12 = fig2.add_subplot(1, 1, 1)
ax12 = plt.plot(x, y, '.', xp, p(xp), '-')


x = x1[5:15]
y = y1[5:15]
y = gaussian_filter(y, sigma = 1)
z = np.polyfit(x, y, 2)
p = np.poly1d(z)
xp = np.linspace(x.min(), x.max(), 200)

fig4 = plt.figure()
ax14 = fig4.add_subplot(1, 1, 1)
ax14 = plt.plot(x, y, '.', xp, p(xp), '-')


plt.show()

