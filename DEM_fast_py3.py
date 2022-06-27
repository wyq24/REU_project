import numpy as np
#import math
#import matplotlib.pyplot as plt
import matplotlib
import scipy.io as io
# from dn2dem_pos import dn2dem_pos
from dn2dem_pos_selfnorm import dn2dem_pos_selfnorm
#import glob

import astropy.time as atime
from astropy.coordinates import SkyCoord
from astropy import units as u
import sunpy.map

from aiapy.calibrate import degradation
#from aiapy.calibrate.util import get_correction_table
from aiapy.calibrate import register, update_pointing
import pickle
import os
import copy
# from tqdm import tqdm
#from functools import partial
from astropy.time import Time

import warnings
from multiprocessing import Process, Manager, Pool, Queue, active_children
import time

warnings.simplefilter('ignore')
matplotlib.rcParams['font.size'] = 16

#tr_1 = int(input("Enter a tr_1:"))
#tr_2 = int(input("Enter a tr_2:"))

cur_dic = '/Volumes/Data/20170820/20220511/info/20220511_10s_long_aia.p'
in_dic_file = open(cur_dic, 'rb')
in_dic = pickle.load(in_dic_file)
in_dic_file.close()
kw_list = ['Z.131.', 'Z.171.', 'Z.211.', 'Z.94.', 'Z.335.', 'Z.193.']
def get_inp_files(tim):
    cfile_list = []
    for kw in kw_list:
        cfile_list.append(in_dic[tim][kw])
    return cfile_list
# define filed of view
#fov = [[-1040,20],[-910,150]]
#for test onlt: fov = [[-920,140],[-910,150]]
fov = [[900,-360],[1000,-220]]
start_tim = time.time()
def make_sub_map(cur_map, fov=None, ref_map=None):
    if not ref_map:
        ref_map = cur_map
    bole = SkyCoord(fov[0][0] * u.arcsec, fov[0][1] * u.arcsec, frame=ref_map.coordinate_frame)
    tori = SkyCoord(fov[1][0] * u.arcsec, fov[1][1] * u.arcsec, frame=ref_map.coordinate_frame)
    sub_map = cur_map.submap(bole, tori)
    return sub_map

#number of image files:
def file_name_time_reader(cfilename):
    tmp_date_str = os.path.basename(cfilename).split('12s.')[1].split('Z.')[0]
    tmp_date_str = tmp_date_str[:13]+':'+tmp_date_str[13:15]+':'+tmp_date_str[15:]
    return Time(tmp_date_str,format='isot')


#get files at a centain time.
#tim = 100
def file_list(tim,f_dict):
    tmp_list = []
    for kwi,ckw in enumerate(kw_list):
        tmp_list.append(f_dict[ckw][tim])
    return tmp_list

#---------------------------------------
#iter by time
def iter_by_time(ti):
    iter_st = time.time()
    cf_list = get_inp_files(ti)
    jd_array = np.zeros((6))
    for cfi, cf in enumerate(cf_list):
        jd_array[cfi] = file_name_time_reader(cf).mjd
    time_string = Time(np.mean(jd_array[:]),format='mjd').isot
    print('cur_file_list_ is:', cf_list)
    amaps = sunpy.map.Map(cf_list)
    # Get the wavelengths of the maps and get index of sort for this list of maps
    wvn0 = [m.meta['wavelnth'] for m in amaps]
    #print(wvn0)
    srt_id = sorted(range(len(wvn0)), key=wvn0.__getitem__)
    print('sorted_list is : ',srt_id)
    amaps = [amaps[i] for i in srt_id]
    #print([m.meta['wavelnth'] for m in amaps])
    channels = [94, 131, 171, 193, 211, 335] * u.angstrom
    cctime = atime.Time(time_string, scale='utc')
    nc = len(channels)
    degs = np.empty(nc)
    for i in np.arange(nc):
        degs[i] = degradation(channels[i], cctime, calibration_version=10)
    #prep aia images
    aprep = []
    for m in amaps:
        m_temp = update_pointing(m)
        aprep.append(register(m_temp))
    # Get the durations for the DN/px/s normalisation and
    # wavenlength to check the order - should already be sorted above
    wvn = [m.meta['wavelnth'] for m in aprep]
    durs = [m.meta['exptime'] for m in aprep]
    # Convert to numpy arrays as make things easier later
    durs = np.array(durs)
    wvn = np.array(wvn)
    #print(durs)
    #print(wvn)

    # pix_dict = eovsa_pixels()
    # --prepare for the calculatioin-------------------
    worder = np.argsort(wvn)
    print(worder)

    #trin = io.readsav('/Users/walter/software/external_package/demreg/python/aia_tresp_en.dat')
    trin = io.readsav('/Users/walterwei/software/external_package/demreg/python/aia_tresp_en.dat')
    tresp_logt = np.array(trin['logt'])
    nt = len(tresp_logt)
    nf = len(trin['tr'][:])
    trmatrix = np.zeros((nt, nf))
    for i in range(0, nf):
        trmatrix[:, i] = trin['tr'][i]
    gains = np.array([18.3, 17.6, 17.7, 18.3, 18.3, 17.6])
    dn2ph = gains * np.array([94, 131, 171, 193, 211, 335]) / 3397.
    rdnse = np.array([1.14, 1.18, 1.15, 1.20, 1.20, 1.18])
    temps = np.logspace(5.7, 7.6, num=42)
    # Temperature bin mid-points for DEM plotting
    mlogt = ([np.mean([(np.log10(temps[i])), np.log10((temps[i + 1]))]) \
              for i in np.arange(0, len(temps) - 1)])

    # ----------------------------------------------------------------------------
    print('iterate each pixel')
    px_loc = aprep[0].world_to_pixel(SkyCoord(fov[0][0]*u.arcsec, fov[0][1]*u.arcsec, frame=aprep[0].coordinate_frame))
    tmp_aprep = []
    for api,cap in enumerate(aprep):
        tmp_aprep.append(make_sub_map(cur_map=cap,fov=fov))
    aprep = tmp_aprep
    cmap_shape = aprep[0].data.shape
    data_cube = np.zeros((cmap_shape[0],cmap_shape[1],6))
    for mi,m in enumerate(aprep):
        data_cube[:,:,mi] = m.data
    pe_time = time.time()
    print('it take to {} to prep in iter_time'.format(pe_time-iter_st))
    print('{} * {}, {} pixels to calculate in total'.format(data_cube.shape[0],data_cube.shape[1], data_cube.shape[0]*data_cube.shape[1]))
    res_list = []
    for xx in range(cmap_shape[0]):
        for yy in range(cmap_shape[1]):
            data = data_cube[xx,yy,:]
            data = np.array(data)
            # Only doing 1 pixel so immediately in units of DN/px
            # print(data)
            cor_data = data / degs
            dn_in = cor_data / durs
            num_pix = 1
            shotnoise = (dn2ph * data * num_pix) ** 0.5 / dn2ph / num_pix / degs
            # Combine errors and put into DN/px/s
            edn_in = (rdnse ** 2 + shotnoise ** 2) ** 0.5 / durs
            # print('edn_in: ', edn_in)
            #  This function internally runs the regularisation twice, first time no weighting to
            # the constraint matrix, second time uses a smoothed form of the first solution
            dem, edem, elogt, chisq, dn_reg = dn2dem_pos_selfnorm(dn_in, edn_in, trmatrix, tresp_logt, temps)
            res_tuple = (dem, edem, elogt, chisq, dn_reg)
            cur_res = {}
            cur_res['fx'] = int(px_loc[0].value) + xx
            cur_res['x'] = xx
            cur_res['fy'] = int(px_loc[1].value) + yy
            cur_res['y'] = yy
            cur_res['res'] = res_tuple
            res_list.append(cur_res)
            #os._exit(0)
        #eo_time = time.time()
        #print('it take {0} to finish {1} column(s), which is {2} pixels'.format(eo_time - pe_time, xx+1, cmap_shape[1]*(xx+1)))
    si_time = time.time()
    print('it take {} to finish one map'.format(si_time - start_tim))
    res_info = {}
    res_info['f_list'] = cf_list
    res_info['fov'] = fov
    res_info['mtime'] = cctime.isot
    #pickle.dump(res_list, open('/Volumes/WD6T/working/20170820/dem/res_save/multi_t{0:0=4d}.p'.format(tim), 'wb'),protocol=2)
    #pickle.dump((res_list,res_info), open('/Volumes/WD6T/working/20170820/dem/res_save/multi_t{0:0=4d}.p'.format(ti), 'wb'),protocol=2)
    pickle.dump((res_list, res_info),
                open('/Volumes/Data/20170820/20220511/dem/res_save/multi_t{0:0=4d}.p'.format(ti), 'wb'))
#iter_by_time(0)

def dem_calculator(odem_res, get_em=False):
    dem_res = copy.deepcopy(odem_res)
    temps = np.logspace(5.7, 7.6, num=42)
    # Temperature bin mid-points for DEM plotting
    mlogt = ([np.mean([(np.log10(temps[i])), np.log10((temps[i + 1]))]) \
              for i in np.arange(0, len(temps) - 1)])
    cur_em = 0.0
    cur_em_wt = 0.0
    for ti, ct in enumerate(mlogt):
        # if ti!= len(mlogt)-1:
        cur_tinterval = 10 ** (mlogt[ti] + 0.046341 / 2) - 10 ** (mlogt[ti] - 0.046341 / 2)
        # else:
        #    cur_tinterval = 10 ** (mlogt[ti]+0.046341) - 10 ** (mlogt[ti])
        cur_em += dem_res['res'][0][ti] * cur_tinterval
        cur_em_wt += dem_res['res'][0][ti] * cur_tinterval * (10 ** mlogt[ti])
    mean_t = cur_em_wt / cur_em
    if np.isnan(mean_t):
        mean_t=0
    if np.isnan(cur_em):
        cur_em = 0
    if get_em:
        return cur_em
    else:
        return mean_t

def dem_map(tim, get_em=False):
    kw_list = ['dem', 'em']
    cur_fits_file = '/Volumes/WD6T/working/20170820/dem/dem_fits/t{0:0=4d}_{1}.fits'.format(tim, kw_list[get_em])
    if not os.path.exists(cur_fits_file):
        #pix_save = open('/Volumes/WD6T/working/20170820/dem/res_save/multi_t{0:0=4d}.p'.format(tim), 'rb')
        pix_save = open('/Volumes/WD6T/working/20170820/dem/res_save/multi_t{0:0=4d}.p'.format(tim), 'rb')
        dem_res_list = pickle.load(pix_save)
        pix_save.close()
        camap = sunpy.map.Map(dem_res_list[1]['f_list'][0])
        camap = update_pointing(camap)
        camap = register(camap)
        rmap = make_sub_map(camap,fov=dem_res_list[1]['fov'])
        new_data = np.zeros_like(rmap.data)
        print('shape', new_data.shape)
        print('number of res', len(dem_res_list[0]))
        for li, cres in enumerate(dem_res_list[0]):
            print(cres['y'],cres['x'])
            print(dem_calculator(cres,get_em=get_em))
            new_data[int(cres['x']),int(cres['y'])] = dem_calculator(cres,get_em=get_em)
        new_map = sunpy.map.Map(new_data, rmap.meta)
        new_map.save(cur_fits_file)
    else:
        new_map = sunpy.map.Map(cur_fits_file)
    return new_map
#
#for ttti in range(tr_1,tr_2):
    #iter_by_time(ttti)
#    if os.path.exists('/Volumes/WD6T/working/20170820/dem/res_save/multi_t{0:0=4d}.p'.format(ttti)):
#        dem_map(ttti,get_em=False)


#demap = dem_map(231)
#os._exit(0)
#demap = dem_map(0)
#demap.peek()
#os._exit(0)
#read the res to maps
#how many cores you want to use
ncpu = 4
#n_time = len(in_dic)
#inp_tim_list = np.arange(184,362).tolist()
#inp_tim_list = np.arange(50,350).tolist()
#inp_tim_list = np.arange(200,250).tolist()
inp_tim_list = [60,100,140,150,160]

cp = Pool(ncpu)
#final_resl = cp.map(cal_func, multi_input_tuple_list)
final_resl = cp.map(iter_by_time, inp_tim_list)
print("Sub-process(es) done.")

# if you want to get T/em map, please use [partial] function:
#to_get_temperature = partial(dem_map,get_em=False)
#to_get_em = partial(dem_map,get_em=True)
#cp = Pool(6)
#inp_tim_list = np.arange(n_time).tolist()
#final_resl = cp.map(cal_func, multi_input_tuple_list)
#tmap_res = cp.map(to_get_temperature, inp_tim_list)
#print("Sub-process(es) done.")



