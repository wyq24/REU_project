import plotting_tools as pot
import os
import matplotlib.pyplot as plt
from suncasa.utils import plot_mapX
#from suncasa.utils import plot_map
import matplotlib as mpl
from matplotlib.colors import SymLogNorm
import md_dspec2 as mdds
import numpy as np
#import sunpy.cm as cm
import sunpy.visualization.colormaps as cm
import sunpy.map as smap
import pickle
#import scifitting_2020 as sf
from matplotlib.colors import ListedColormap
from astropy import units as u
import scipy.ndimage.measurements
from astropy.io import fits
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import copy
from scipy.io import readsav
from astropy.time import Time
from suncasa.utils import DButil
import numpy.ma as ma
import astropy.visualization.mpl_normalize
from astropy.time import TimeDelta
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import subprocess
#2import demreg_prep as dprp
from sunpy.io import write_file
from sunpy.time import TimeRange, parse_time
import enhance
import matplotlib.dates as mdates
import matplotlib.image as mpimg
from PIL import Image
from matplotlib.ticker import ScalarFormatter
#import skp_plt as skp
import multiprocessing as mlp
import matplotlib.colors as colors
import glob
from matplotlib.dates import DateFormatter
from matplotlib import ticker
#import tmp_plt as tp
from astropy.coordinates import SkyCoord
from sunpy.physics.differential_rotation import diff_rot, solar_rotate_coordinate
#from sunpy.coordinates import Helioprojective,  propagate_with_solar_surface
from astropy.wcs import WCS
#from aiapy.calibrate import register, update_pointing
from functools import partial
from tqdm import *
#test for branching
from multiprocessing import Process, Manager, Pool, Queue, active_children
from sunpy.instr.aia import aiaprep

# import threadpoolctl
# threadpoolctl.threadpool_limits(8)
with open('/Volumes/Data/20170820/20220511/info/bmsize.p', 'rb') as fbmsize:
    bmsize = pickle.load(fbmsize, encoding='latin1')
fbmsize.close()
with open('/Volumes/Data/20170820/20220511/info/cfreqs.p', 'rb') as fcfreq:
    cfreq = pickle.load(fcfreq, encoding='latin1')
fcfreq.close()
cur_dic = '/Volumes/Data/20170820/20220511/info/20220511_10s_long_aia.p'
in_dic_file = open(cur_dic, 'rb')
in_dic_10s_la = pickle.load(in_dic_file)
in_dic_file.close()

# fnb_result = pickle.load(open('/Volumes/Data/20170906/info/spw_info.p', 'rb'),encoding='latin1')
# cfreq = fnb_result['cfreqs'][0:30] / 1.e9
fov = [[580, -253], [607, -211]]
fov2 = [[560, -263], [627, -191]]

with open('/Volumes/Data/20170820/info/eovsa_bmsize.p', 'rb') as fbmsize:
    bmsize = pickle.load(fbmsize,encoding='latin1')
fbmsize.close()
with open('/Volumes/Data/20170820/info/eovsa_cfreqs.p', 'rb') as fcfreq:
    cfreq = pickle.load(fcfreq,encoding='latin1')
fcfreq.close()
cfreq = cfreq/1.e9
#cfreq = cfreq[1:]

def diffmap(cax, st, et, kw, fov2=None, enhance_it=False, high_tres=True):
    if fov2 is None:
        fov2 = [[560, -263], [627, -191]]
    if high_tres:
        in_dic = in_dic_1s
    else:
        in_dic = in_dic_10s
    sub_map = pot.rdiff_map(kw=kw, tim=et, timeint=st, fov=fov2, high_tres=high_tres)
    # sub_map = pot.make_sub_map(cur_map=diff_map, fov=fov2)
    #if kw != 'bbso':
    if 'ha' not in kw:
    # kwargs = {'norm': SymLogNorm(linthresh=0, vmin=np.nanmin(sub_map.data), vmax=np.nanmax(sub_map.data)),
        kwargs = {
            'cmap': plt.get_cmap(name='sdoaia' + kw.replace('Z', '').replace('.', ''))}
    else:
        # kwargs = {'norm': SymLogNorm(linthresh=0, vmin=np.nanmin(sub_map.data), vmax=np.nanmax(sub_map.data))}
        kwargs = {'cmap': 'gist_gray'}
        # kwargs = {'cmap': 'jet'}
    if enhance_it:
        if 'ha' in kw:
            sub_map = smap.Map(enhance.mgn(data=sub_map.data, sigma=[21.25, 42.5, 85, 170, 340, 680]), sub_map.meta)
        else:
            sub_map = smap.Map(enhance.mgn(data=sub_map.data), sub_map.meta)
    sub_map_obj = plot_mapX.Sunmap()
    cim = sub_map_obj.imshow(axes=cax, **kwargs)
    cax.set_xlabel('')
    cax.set_ylabel('')
    return cim


def justmap(cax, tim, kw, fov=None, alp=None, updown=None, custom_color=None, inorm=None, enhance_it=False,
            start_mid_range=None, norm_exp=False, custom_expt=None, judge_material=None, custom_file=None,
            high_tres=False,custom_kernel_size=None, detrending=False,long_aia = True, blooming_block=False,prep_aia=True):
    if long_aia:
        in_dic = in_dic_10s_la
    #cur_dic = '/Volumes/WD6T/working/20170703/info/new_20170703_10s.p'
    #in_dic_file = open(cur_dic, 'rb')
    #in_dic = pickle.load(in_dic_file)
    #in_dic_file.close()
    else:
        print('using the original in_dic? Sure about that?')
        in_dic = in_dic_10s
    if kw in ['dem', 'em','br']:
        # cur_map = dprp.dem_sav_to_map(tim)
        if fov is None:
            fov = [[-1150,-40],[-910,140]]
        if kw == 'dem':
            cur_map = dem_map(tim)
        elif kw == 'em':
            cur_map = dem_map(tim=tim, get_em=True)
        elif kw == 'br':
            cur_map = get_br_on_date(tim)
    elif custom_file is not None:
        cur_map = smap.Map(custom_file)
    else:
        if detrending and 'Z.' in kw:
            if kw not in ['Z.131.','Z.193.','Z.211.']:
                cur_map = smap.Map(find_detrending_substitute(kw=kw,filename=in_dic[tim][kw]))
            else:
                cur_map = smap.Map(in_dic[tim]['det_'+kw[2:-1]])
        else:
            if not blooming_block:
                    cur_map = smap.Map(in_dic[tim][kw])
                    if prep_aia:
                        cur_map = aiaprep(cur_map)
            else:
                blooming_save = '/Volumes/Data/20170820/info/blooming_list_{}.p'.format(kw.split('.')[1])
                try:
                    with open(blooming_save, 'rb') as fff:
                        blooming_list = pickle.load(fff,encoding='latin1')
                    fff.close()
                    tmp_tim = copy.deepcopy(tim)
                    while in_dic[tmp_tim][kw] in blooming_list:
                        tmp_tim += 1
                    cur_file = in_dic[tmp_tim][kw]
                except:
                    print('no list saved. Check your input')
                    cur_file = in_dic[tim][kw]
                cur_map = smap.Map(cur_file)
    if kw == 'hmi':
        cur_map = cur_map.rotate(angle=180 * u.deg)
    #elif kw == 'bbso':
    #    cur_map = gst_image_correction_2nd_round(tim=tim)
    if fov != 'full':
        sub_map = pot.make_sub_map(cur_map=cur_map, fov=fov)
    else:
        sub_map = copy.deepcopy(cur_map)

    # cnorm=astropy.visualization.mpl_normalize.simple_norm(data=sub_map.data,stretch=inorm)
    #if norm_exp and 'Z' in kw:
        #tmp_jud = pot.aia_exposure_helper(kw=kw, tim=tim, fov=fov, custom_exp=custom_expt,
        #                                  judge_material=judge_material,high_tres=high_tres)
        #if tmp_jud is not 1:
        #    print('gocha')
        #    sub_map = tmp_jud
    '''
    if 'Z' in kw and custom_expt is not None:
        cur_exp = cur_map.meta['exptime']
        normed_data=cur_map.data*(custom_expt/cur_exp)
        cur_map=smap.Map(normed_data,cur_map.meta)
    '''
    if 'Z' in kw:
        cur_max = np.nanmax(sub_map.data)
        #kwargs = {'norm': SymLogNorm(linthresh=aia_max_list[kw] * 0.005, vmin=0.0,
        #                             vmax=aia_max_list[kw]),
        kwargs = {'norm': SymLogNorm(linthresh=cur_max * 0.005, vmin=0.0,
                                     vmax=cur_max),
        #kwargs = {
                  'cmap': plt.get_cmap(name='sdoaia' + kw.replace('Z', '').replace('.', ''))}

    elif kw == 'hmi':
        # kwargs = {'norm': SymLogNorm(linthresh=100, vmin=np.nanmin(sub_map.data), vmax=np.nanmax(sub_map.data))}
        kwargs = {'cmap': 'bwr'}
    elif 'ha' in kw:
        kwargs = {'cmap': 'gist_gray',
                  'norm': SymLogNorm(linthresh=np.nanmax(sub_map.data) * 0.7, vmin=np.nanmin(sub_map.data),
                                     vmax=np.nanmax(sub_map.data))}
    elif kw == 'em':
        kwargs = {'norm': SymLogNorm(linthresh=np.nanmax(sub_map.data) * 0.1, vmin=np.nanmin(sub_map.data),
                                     vmax=np.nanmax(sub_map.data)),
                  'cmap': 'viridis'}
    elif kw=='br':
        kwargs={}
    elif kw == 'dem':
        kwargs = {'cmap': 'viridis'}
    elif 'eu' in kw:
        kwargs = {'cmap': 'gist_gray'}
    if custom_file is not None:
        kwargs = {'cmap': 'Spectral'}
    # plot_map.imshow(sunpymap=sub_map, axes=cax, **kwargs)
    if alp is not None:
        kwargs.update({'alpha': alp})
    if updown is not None:
        kwargs.update({'origin': updown})
    if custom_color is not None:
        kwargs.update({'cmap': custom_color})
    if inorm is not None:
        kwargs.update({'norm': inorm})
    if enhance_it:
        if inorm is not None:
            Warning('!!!! you are plottinh a normlized image')
        if custom_kernel_size is not None:
            sub_map = smap.Map(enhance.mgn(data=sub_map.data, sigma=custom_kernel_size), sub_map.meta)
        elif 'ha' in kw:
            sub_map = smap.Map(enhance.mgn(data=sub_map.data, sigma=[21.25, 42.5, 85, 170, 340, 680]), sub_map.meta)
            # pool = mlp.Pool()
            # partial_enhance=partial(enhance.mgn,sub_map.data, [21.25, 42.5, 85, 170, 340, 680], 0.7, 3.2,0.7,None,3,True)
            # new_data = pool.map(partial_enhance, sub_map.data)
            # pool.close()
            # pool.join()
            # sub_map=smap.Map(enhance.mgn(data=sub_map.data),sub_map.meta)
            # sub_map = smap.Map(new_data, sub_map.meta)
        else:
            # pool = mlp.Pool()
            # new_data = pool.map(enhance.mgn, sub_map.data)
            # pool.close()
            # pool.join()
            # sub_map=smap.Map(enhance.mgn(data=sub_map.data),sub_map.meta)
            # sub_map=smap.Map(enhance.mgn(data=sub_map.data,sigma=[15.0]),sub_map.meta)
            #sub_map = smap.Map(enhance.mgn(data=sub_map.data, sigma=[2, 5, 10, 20, 40], gamma=2.2), sub_map.meta)
            #sub_map = smap.Map(enhance.mgn(data=sub_map.data, sigma=[2, 5, 10, 20], gamma=1.4), sub_map.meta)
            sub_map = smap.Map(enhance.mgn(data=sub_map.data, sigma=[1.5,5,10,40], gamma=3.0), sub_map.meta)
            # sub_map = smap.Map(new_data, sub_map.meta)
        kwargs.pop('norm', None)
    # elif 'Z' in kw:
    #     sub_map = norm_aia_map(sunpymap=sub_map)
    if start_mid_range is not None:
        tmp_sub_data = np.clip(sub_map.data, np.nanmax(sub_map.data) * start_mid_range, np.nanmax(sub_map.data))
        sub_map = smap.Map(tmp_sub_data, sub_map.meta)
    sub_map_obj = plot_mapX.Sunmap(sunmap=sub_map)
    cur_im = sub_map_obj.imshow(axes=cax, **kwargs)
    if fov != 'full':
        cax.set_xlim(fov[0][0],fov[1][0])
        cax.set_ylim(fov[0][1], fov[1][1])
    cax.set_xlabel('')
    cax.set_ylabel('')
    return cur_im

def aia_composite_map(cax,kws,tim,fov,blooming_block=False,custom_kernel_size=None,
                      enhance_it=False,high_tres=False,upper_alpha=0.6,custom_file = None,custom_norm=None,larger_fov=1):
    if high_tres:
        in_dic = in_dic_1s
    else:
        in_dic = in_dic_10s
    sub_maps = []
    for ci,ckw in enumerate(kws):
        if not blooming_block:
            cur_file = in_dic[tim][ckw]
        else:
            blooming_save = '/Volumes/Data/20170820/info/blooming_list_{}.p'.format(ckw.split('.')[1])
            try:
                with open(blooming_save, 'rb') as fff:
                    blooming_list = pickle.load(fff, encoding='latin1')
                fff.close()
                tmp_tim = copy.deepcopy(tim)
                while in_dic[tmp_tim][ckw] in blooming_list:
                    tmp_tim += 1
                cur_file = in_dic[tmp_tim][ckw]
            except:
                print('no list saved. Check your input')
                cur_file = in_dic[tim][ckw]
        if custom_file is not None:
            if ci == custom_file[0]:
                cur_file = custom_file[1]
            cur_map = smap.Map(cur_file)
        csub_map = pot.make_sub_map(cur_map=cur_map, fov=fov)
        if enhance_it:
            sub_map = smap.Map(enhance.mgn(data=csub_map.data, sigma=custom_kernel_size[ci]), csub_map.meta)
        sub_maps.append(pot.make_sub_map(cur_map=cur_map, fov=fov))
    comp_map = smap.Map(sub_maps[0],sub_maps[1],composite=True)
    for ci, ckw in enumerate(kws):
        c_setting = comp_map.get_plot_settings(ci)
        c_setting['norm'] = SymLogNorm(linthresh=np.nanmax(comp_map.get_map(0).data)*0.005, vmin=0.0,
                                       vmax=np.nanmax(comp_map.get_map(0).data))
        if custom_norm is not None:
            if ci == custom_norm[0]:
                c_setting['norm'] = custom_norm[1]
        comp_map.set_plot_settings(ci,c_setting)
    comp_map.set_alpha(1,upper_alpha)
    #ctr = comp_map.get_map(larger_fov).top_right_coord
    #cbl = comp_map.get_map(larger_fov).bottom_left_coord
    #comp_map.top_right_coord = ctr
    #comp_map.bottom_left_coord = cbl
    #sub_map_obj = plot_mapX.Sunmap(sunmap=comp_map)
    #cur_im = sub_map_obj.imshow(axes=cax, **kwargs)
    #cur_im = sub_map_obj.imshow(axes=cax)
    cur_im = comp_map.plot(axes=cax)
    return cur_im

def just_eovsa_map(tim, cax,spw, fov='full',custom_file=None,high_tres=False, custom_color=None,inorm=None):
    if high_tres:
        in_dic = in_dic_1s
    else:
        in_dic = in_dic_10s
    curfits = in_dic[tim]['radio_sbs'][spw]
    if custom_file is not None:
        curfits = '/Volumes/Data/20170820/eovsa_bkg_subtraction/slfcal/images_slfcaled/slfcaled_tb_final_XX_10s_XX_3820_3830_3730_3740_s{0:0=2d}.fits'.format(spw+1)
    cur_map = smap.Map(curfits)
    if fov != 'full':
        cur_map = pot.make_sub_map(cur_map=cur_map,fov=fov)
    kwargs={'cmap':'viridis'}
    if custom_color is not None:
        kwargs.update({'cmap': custom_color})
    if inorm is not None:
        kwargs.update({'norm': inorm})
    cur_map_obj = plot_mapX.Sunmap(sunmap=cur_map)
    cim = cur_map_obj.imshow(sunpymap=cur_map, axes=cax, **kwargs)
    #plt.colorbar(cim, cax=cax)
    if fov != 'full':
        cax.set_xlim(fov[0][0],fov[1][0])
        cax.set_ylim(fov[0][1], fov[1][1])
    return cim



def overlap_maps(kw1, kw2, tim, cur_ax, fov=None, renorm=False, over_mid=None, cust_alp=0.28):
    if kw1 == 'hmi':
        justmap(cur_ax, tim, kw1, fov=fov, enhance_it=False, custom_color=plt.cm.gray)
    else:
        justmap(cur_ax, tim, kw1, fov=fov, enhance_it=False)
    justmap(cur_ax, tim, kw2, fov=fov, enhance_it=True, alp=cust_alp, start_mid_range=over_mid, custom_expt=2.9)
    '''
    map_list=[]
    for ckw in [kw1,kw2]:
        if 'fits' in ckw:
            cur_map=smap.Map(ckw)
        else:
            cur_map=smap.Map(in_dic[tim][ckw])
        if 'magnetogram' in ckw:
            cur_map = cur_map.rotate(angle=180 * u.deg)
        if fov is not None:
            cur_map=pot.make_sub_map(cur_map=cur_map,fov=fov)
        if renorm and 'Z.' in ckw:
            tmp_data = np.clip(cur_map.data, np.median(cur_map.data)*over_mid, np.nanmax(cur_map.data))
            cur_map=smap.Map(tmp_data,cur_map.meta)
        map_list.append(cur_map)
    comp_map = smap.Map(map_list[0], composite=True, alpha=1.0)
    comp_map.add_map(map_list[1], alpha=0.25)
    kwargs={}
    cim=comp_map.plot(axes=cur_ax)
    #plot_map.imshow(sunpymap=comp_map, axes=cax, **kwargs)
    return cim
    '''


def get_minmax_aia(tim, kw, cfov=None):
    cur_map = smap.Map(in_dic[tim][kw])
    if fov is not None:
        cur_map = pot.make_sub_map(cur_map=cur_map, fov=cfov)
    return [np.nanmin(cur_map.data), np.nanmax(cur_map.data)]


def plot_goes_on_axes(trange, cur_ax, cur_color='k'):
    from matplotlib.dates import (DateFormatter, rrulewrapper, RRuleLocator, drange)
    from sunpy import lightcurve
    if 'T' in trange[0]:
        t1 = Time(trange[0], format='isot')
        t2 = Time(trange[1], format='isot')
    else:
        t1 = Time(trange[0], format='iso')
        t2 = Time(trange[1], format='iso')
    goes_lightcurve = lightcurve.GOESLightCurve.create(t1.iso, t2.iso)
    tim_list = []
    for ttti in range(len(goes_lightcurve.data)):
        tim_list.append(goes_lightcurve.data.index[ttti].isoformat())
    ctim_list = Time(tim_list, format='isot')
    tim_plt = ctim_list.plot_date
    # normed_lc=np.zeros_like(goes_lightcurve.data['xrsa'].to_numpy())
    normed_lc = [cval * (1.0 / np.nanmax(goes_lightcurve.data['xrsb'].to_numpy())) for cval in
                 goes_lightcurve.data['xrsb'].to_numpy()]
    #cur_ax.plot_date(tim_plt, normed_lc, c=cur_color, marker='', linestyle='-', label='GOES 1.0--8.0 $\AA$')
    cur_ax.plot_date(tim_plt, goes_lightcurve.data['xrsb'], c=cur_color, marker='', linestyle='-', label='GOES 1.0--8.0 $\AA$')
    cur_ax.format_coord = format_coord
    cur_ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    # cur_ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    # cur_ax.xaxis.set_tick_params(rotation=25)
    cur_ax.xaxis.set_tick_params(rotation=0)
    cur_ax.legend()


def format_coord(x, y):
    return 'x = {0}, y = {1:.3f}'.format(x, y)


def norm_aia_map(sunpymap):
    new_data = sunpymap.data
    if 'AIA' in sunpymap.meta['telescop']:
        new_data = new_data / (sunpymap.meta['exptime']) * 3
        return smap.Map(new_data, sunpymap.meta)
    else:
        return sunpymap


def tick_helper(cax, tsta):
    if tsta[0] == 0:
        cax.set_xlabel('')
        cax.set_xticks([])
    if tsta[1] == 0:
        cax.set_ylabel('')
        cax.set_yticks([])
    return cax


def plot_1():
    fig = plt.figure(figsize=(8, 8))
    for i, ckw in enumerate(['bbso', 'Z.171.']):
        cax = fig.add_subplot(3, 4, i + 1)
        justmap(cax, 100, ckw)
        # diffmap(cax,133,135,ckw)
        # pot.plot_eovsa_contourf(tim=124, ax=cax, fov=fov, abs_level=eovsa_max_list * 0.69,
        #                        dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p', abs_contourf=True)
    plt.show(aspect='auto')


def plot_spectrum(x, y, ax, timerange, hm_t=None, phase=2,cur_marker='x', label_swt=[1, 1], custom_cm=plt.cm.gist_rainbow,
                  custom_alp=1, extra_norm=None,radius=1,custom_color=None,to_sfu=False):
    if hm_t is None:
        trange = np.arange(timerange[0], timerange[1])
    else:
        trange = np.linspace(timerange[0], timerange[1], hm_t, dtype=int, endpoint=True)
    data_cube = np.zeros((30, len(trange)))
    if extra_norm is None:
        norm = mpl.colors.Normalize(vmin=timerange[0], vmax=timerange[1])
    else:
        norm = mpl.colors.Normalize(vmin=extra_norm[0], vmax=extra_norm[1])
    for ti, t in enumerate(trange):
        print(t)
        data_cube[:, ti] = sf.get_data_fits(x=x, y=y, tim=t, phase=phase, world=1, apply_ratio=True, bkg_on=False,to_sfu=
                                            to_sfu, radius=radius)
        bkg_data = sf.get_rms(tim=t,phase=phase,to_sfu=to_sfu)
        cim = ax.errorbar(cfreq, data_cube[:, ti], yerr=bkg_data, linestyle='None', marker=cur_marker,
                    # label='flare spectrum, T, delta, n_th, n_nonth, B',
                    color=custom_cm(norm(t)), alpha=custom_alp)
    ax.set_xscale("log")
    ax.set_yscale("log")
    if not to_sfu:
        ax.set_ylim(1e6, 4e7)
    ax.xaxis.set_minor_formatter(ScalarFormatter())
    if label_swt[0] == 1:
        ax.set_xlabel('Freq [GHz]')
    else:
        ax.set_xlabel('')
        ax.set_xticks([])

    if label_swt[1] == 1:
        if to_sfu:
            ax.set_ylabel('Flux Density [sfu]')
        else:
            ax.set_ylabel('${T_b}$ [k]')
    else:
        ax.set_ylabel('')
        ax.set_yticks([])
    # if pleg:
    #    add_custom_cb(target_ax=ax,cmin)
    # ax.tick_params(axis='x',direction='in')
    # ax.tick_params(axis='y',direction='in')
    # ax.ticklabel_format(axis='x', style='plain')
    return cim

def plot_spectrum_single(x, y, ax, tim, phase=2,cur_marker='x', label_swt=[1, 1], custom_cm=plt.cm.gist_rainbow,
                  custom_alp=1, extra_norm=None,radius=None,custom_color=None,to_sfu=False, bkg_subed=False):
    bkg_data = sf.get_rms(tim,to_sfu=to_sfu,phase=phase)
    tmp_data = sf.get_data_fits(x=x, y=y, tim=tim, phase=phase, world=1, apply_ratio=True, bkg_on=False,to_sfu=to_sfu,
                                radius=radius,bkg_subed=bkg_subed)
    ax.errorbar(cfreq, tmp_data, yerr=bkg_data, linestyle='None', marker=cur_marker,
                    # label='flare spectrum, T, delta, n_th, n_nonth, B',
                    color=custom_color, alpha=custom_alp)

    ax.set_xscale("log")
    ax.set_yscale("log")
    if not to_sfu:
        ax.set_ylim(1e6, 4e7)
    ax.xaxis.set_minor_formatter(ScalarFormatter())
    if label_swt[0] == 1:
        ax.set_xlabel('Freq [GHz]')
    else:
        ax.set_xlabel('')
        ax.set_xticks([])

    if label_swt[1] == 1:
        if to_sfu:
            ax.set_ylabel('Flux Density [sfu]')
        else:
            ax.set_ylabel('${T_b}$ [k]')
    else:
        ax.set_ylabel('')
        ax.set_yticks([])
    # if pleg:
    #    add_custom_cb(target_ax=ax,cmin)
    # ax.tick_params(axis='x',direction='in')
    # ax.tick_params(axis='y',direction='in')
    # ax.ticklabel_format(axis='x', style='plain')


def cmap_generator(cur_rgba):
    alpha = np.ones(256)
    # cr,cg,cb = cur_rgba[0:3]
    celement = (cur_rgba[0:3] / np.nanmax(cur_rgba[0:3]))
    lumi = np.linspace(1, 0, 256)
    ccmap = np.c_[lumi * celement[0], lumi * celement[1], lumi * celement[2], alpha]
    print(ccmap)
    res_cmap = ListedColormap(ccmap)
    return res_cmap


def plot_B_contour(tim, ax, levels, cfov=None, alp=None):
    if alp is None:
        alp = 0.5
    if cfov is None:
        cfov = fov
    # cur_map = smap.Map(in_dic[tim]['hmi'])
    #cur_map = smap.Map('/Volumes/Data/20170715/extrapolation/hmi.M_720s.20170715_192400_br_TAI.fits')
    cur_map = get_br_on_date('0820')
    #cur_map = cur_map.rotate(angle=180 * u.deg)
    sub_map = pot.make_sub_map(cur_map=cur_map, fov=cfov)
    cmap1 = plt.cm.Reds
    cmap2 = plt.cm.Blues
    norm1 = mpl.colors.Normalize(vmin=0, vmax=np.nanmax(sub_map.data))
    norm2 = mpl.colors.Normalize(vmin=np.nanmin(sub_map.data), vmax=0)
    print(np.nanmax(sub_map.data), np.nanmin(sub_map.data))
    # kwargs1 = {"levels": levels, "colors": [cmap1(norm1(cl * np.nanmax(sub_map.data))) for cl in levels]}
    kwargs1 = {'levels': [levels[1], 10000], 'colors': ['b', 'b'], 'linewidth': 0.1, 'alpha': alp}
    kwargs2 = {'levels': [-10000, levels[0]], 'colors': ['r', 'r'], 'linewidth': 0.1, 'alpha': alp}
    # kwargs2 = {"levels": levels, "colors": [cmap2(norm2(cl * np.nanmin(sub_map.data))) for cl in levels]}
    sub_map_obj = plot_mapX.Sunmap(sunmap=sub_map)
    im1 = sub_map_obj.contour( axes=ax, rot=0,label= '{} G'.format(levels[1]), **kwargs1)
    im2 = sub_map_obj.contour( axes=ax, rot=0, label= '{} G'.format(levels[0]),**kwargs2)
    #ax.legend()
    return im1,im2


def plot_aia_contour(tim, ax, ckw, levels, cco='r', cfov=None, alp=None, txt_label=None,high_tres=False):
    if alp is None:
        alp = 0.5
    if cfov is None:
        cfov = fov
    if high_tres:
        in_dic = in_dic_1s
    else:
        in_dic = in_dic_10s
    cur_map = smap.Map(in_dic[tim][ckw])
    sub_map = pot.make_sub_map(cur_map=cur_map, fov=cfov)
    cmap1 = plt.cm.Reds
    cmap2 = plt.cm.Blues
    norm1 = mpl.colors.Normalize(vmin=0, vmax=np.nanmax(sub_map.data))
    norm2 = mpl.colors.Normalize(vmin=np.nanmin(sub_map.data), vmax=0)
    print(np.nanmax(sub_map.data), np.nanmin(sub_map.data))
    # kwargs1 = {"levels": levels, "colors": [cmap1(norm1(cl * np.nanmax(sub_map.data))) for cl in levels]}
    if levels[0] < 1:
        kwargs = {'levels': [levels[0] * np.nanmax(sub_map.data), 10000000.0], 'colors': [cco, cco], 'linewidth': 0.1,
                  'alpha': alp, 'label': txt_label}
    else:
        kwargs = {'levels': [levels[0], 10000], 'colors': [cco, cco], 'linewidth': 0.1, 'alpha': alp,
                  'label': txt_label}
    # kwargs2 = {"levels": levels, "colors": [cmap2(norm2(cl * np.nanmin(sub_map.data))) for cl in levels]}
    sub_map_obj = plot_mapX.Sunmap(sunmap=sub_map)
    im = sub_map_obj.contour(axes=ax, rot=0, **kwargs)
    return im


def plot_eovsa_contourf(tim, ax, spws=range(50), fov=None, abs_contourf=None, level=None, abs_level=None, dic_file=None,
                        inp_rgba=None, cs_start=None, single_spw=False, alt_cmap=None, calp=None, reverse_zorder=False,
                        high_tres=False, bkg_sub=False):
    cmap2 = plt.cm.Spectral
    # cmap2 = plt.cm.gray
    if inp_rgba is None:
        # cmap = plt.cm.winter
        # cmap = plt.cm.gist_gray
        if alt_cmap is None:
            cmap = plt.cm.viridis
        else:
            cmap = alt_cmap
    else:
        cmap = cmap_generator(inp_rgba)
    if cs_start is None:
        norm = mpl.colors.Normalize(vmin=0, vmax=50)
    else:
        norm = mpl.colors.Normalize(vmin=cs_start[0], vmax=cs_start[1])
    norm2 = mpl.colors.Normalize(vmin=-10, vmax=40)
    if dic_file is None:
        dic_file = '/Volumes/Data/20170820/info/20170820_10s.p'
        in_dic = pickle.load(open(dic_file, 'rb'),encoding='latin1')
    if high_tres:
        in_dic = in_dic_1s
    else:
        # in_dic = pickle.load(open('/home/walter/Downloads/From_NJIT/doppler_20190604_bt_dic_subed.p', 'rb'),encoding='latin1')
        # in_dic = pickle.load(open('/home/walter/Downloads/From_NJIT/time_dics/20190610_bt_dic.p', 'rb'),encoding='latin1')
        pass
    if calp is None:
        calp = 0.6
    for i, spw in enumerate(spws):
        if reverse_zorder:
            tmp_spws = copy.deepcopy(spws)
            tmp_spws.reverse()
            spw = tmp_spws[i]
        if bkg_sub:
            if high_tres:
                print('current spw is', spw)
                cur_emap = smap.Map(
                    '/Volumes/Data/20170820/eovsa_bkg_subtraction/slfcal/images_slfcaled/subed_slfcaled_tb_final_XX_1s_XX_t{0}_s{1:0=2d}.fits'.format(tim,
                        spw+1))
            else:
                cur_emap = smap.Map(
                    '/Volumes/Data/20170820/eovsa_bkg_subtraction_10s_3730_40/slfcal/images_slfcaled/subed_slfcaled_tb_final_XX_10s_XX_t{0}_s{1:0=2d}.fits'.format(
                        tim-184, spw+1))
        else:
            cur_emap = smap.Map(in_dic[tim]['radio_sbs'][spw])
        if fov:
            cur_emap = pot.make_sub_map(cur_map=cur_emap, fov=fov)
            if abs_contourf:
                if not single_spw:
                    print(abs_level[spw])
                    cur_level = [abs_level[spw], 1.e20]
                    # cur_level=[abs_level,1.e20]
                else:
                    cur_level = [abs_level[spw] * cl for cl in level]
            else:
                if not single_spw:
                    cur_level = [np.nanmax(cur_emap.data) * level, 1.e20]
                else:
                    print(type(np.nanmax(cur_emap.data)))
                    cur_level = [np.nanmax(cur_emap.data) * cl for cl in level]
        # kwargs={"levels":cur_level,"colors":sm.to_rgba(in_efits),"linewidth":0.5}
        # kwargs = {"levels": cur_level, "colors": [cmap(norm(spw))], "alpha": 0.2,"linewidth":0.02}
        print('current_levelis : ', cur_level)
        if cs_start is not None:
            kwargs = {"levels": cur_level, "colors": [cmap(norm(spw))], "alpha": calp}
        elif not single_spw:
            kwargs = {"levels": cur_level, "colors": [cmap(norm(spw))], "alpha": calp}
        else:
            # kwargs = {"levels": cur_level, "colors": [cmap(norm(spw))], "alpha": 0.6}
            kwargs = {"levels": cur_level, "colors": [cmap(norm(cl * 50)) for cl in level], "alpha": calp}
            print(kwargs)

        # kwargss = {"levels": cur_level, "colors": [cmap2(norm2(spw))], "alpha": 0.6}
        '''
        if i not in range(10, 25):
            plot_map.contourf(sunpymap=cur_emap, axes=ax, rot=0, **kwargss)
        else:
            plot_map.contourf(sunpymap=cur_emap, axes=ax, rot=0, **kwargs)
        '''
        cur_emap_obj = plot_mapX.Sunmap(sunmap=cur_emap)
        cur_emap_obj.contourf(axes=ax, rot=0, **kwargs)


def plot_eovsa_contour(tim, ax, spws=range(50), fov=None, abs_contourf=None, level=None, abs_level=None, dic_file=None,
                       inp_rgba=None, cs_start=None, single_spw=False, contour_label=False, calp=None,
                       reverse_zorder=False, high_tres=False,bkg_sub=False):
    cmap2 = plt.cm.Spectral
    # cmap2 = plt.cm.gray
    # if inp_rgba is None:
    if inp_rgba is not None:
        # cmap = plt.cm.winter
        cmap = inp_rgba
        # cmap = plt.cm.gist_gray
        # cmap = plt.cm.Spectral
    else:
        # cmap = cmap_generator(inp_rgba)
        cmap = plt.cm.Spectral
    if cs_start is None:
        norm = mpl.colors.Normalize(vmin=0, vmax=50)
    else:
        norm = mpl.colors.Normalize(vmin=cs_start[0], vmax=cs_start[1])
    norm2 = mpl.colors.Normalize(vmin=-10, vmax=50)
    if dic_file is None:
        #dic_file = '/Volumes/WD6T/working/20170703/info/20170703_1s.p'
        #in_dic = pickle.load(open(dic_file, 'rb'),encoding='latin1')
    #     in_dic = in_dic_10s
    # if high_tres:
    #     in_dic = in_dic_1s
    # # if high_tres:
    #    in_dic = in_dic_1s
        cur_dic = '/Volumes/Data/20170820/20220511/info/20220511_10s_long_aia.p'
        in_dic_file = open(cur_dic, 'rb')
        in_dic = pickle.load(in_dic_file)
        in_dic_file.close()
    if calp is None:
        calp = 1
    for i, spw in enumerate(spws):
        if reverse_zorder:
            tmp_spws = spws
            tmp_spws.reverse()
            spw = tmp_spws[i]
        #print('cur spw is :', spw)
        if bkg_sub:
            if high_tres:
                print('eovsa file is: ', '/Volumes/Data/20170820/eovsa_bkg_subtraction/slfcal/images_slfcaled/subed_slfcaled_tb_final_XX_1s_XX_t{0}_s{1:0=2d}.fits'.format(
                        tim,
                        spw+1))
                cur_emap = smap.Map(
                    '/Volumes/Data/20170820/eovsa_bkg_subtraction/slfcal/images_slfcaled/subed_slfcaled_tb_final_XX_1s_XX_t{0}_s{1:0=2d}.fits'.format(
                        tim,
                        spw+1))
            else:
                cur_emap = smap.Map(
                    '/Volumes/Data/20170820/eovsa_bkg_subtraction_10s_3730_40/slfcal/images_slfcaled/subed_slfcaled_tb_final_XX_10s_XX_t{0}_s{1:0=2d}.fits'.format(
                        tim-184, spw+1))
                print('/Volumes/Data/20170820/eovsa_bkg_subtraction_10s_3730_40/slfcal/images_slfcaled/subed_slfcaled_tb_final_XX_10s_XX_t{0}_s{1:0=2d}.fits'.format(
                        tim-184, spw+1))
        else:
            cur_emap = smap.Map(in_dic[tim]['radio_sbs'][spw])
        if not fov:
            print('define a fov, youah')
        cur_emap = pot.make_sub_map(cur_map=cur_emap, fov=fov)
        if abs_contourf:
            if not single_spw:
                print(abs_level[spw])
                cur_level = [abs_level[spw], 1.e20]
                # cur_level=[abs_level,1.e20]
            else:
                cur_level = [abs_level[spw] * cl for cl in level]
        else:
            if not single_spw:
                cur_level = [np.nanmax(cur_emap.data) * level, 1.e20]
            else:
                print(type(np.nanmax(cur_emap.data)))
                cur_level = [np.nanmax(cur_emap.data) * cl for cl in level]
        # kwargs={"levels":cur_level,"colors":sm.to_rgba(in_efits),"linewidth":0.5}
        # kwargs = {"levels": cur_level, "colors": [cmap(norm(spw))], "alpha": 0.2,"linewidth":0.02}
        if cs_start is not None:
            kwargs = {"levels": cur_level, "colors": [cmap(norm(spw))], "alpha": calp}
        elif single_spw:
            kwargs = {"levels": cur_level, "colors": [cmap(norm(cl * 50)) for cl in level], "alpha": calp,
                      "linewidths": 1}
        else:
            kwargs = {"levels": cur_level, "colors": [cmap(norm(spw))], "alpha": calp}

        # kwargss = {"levels": cur_level, "colors": [cmap2(norm2(spw))], "alpha": 0.6}
        '''
        if i not in range(10, 25):
            plot_map.contourf(sunpymap=cur_emap, axes=ax, rot=0, **kwargss)
        else:
            plot_map.contourf(sunpymap=cur_emap, axes=ax, rot=0, **kwargs)
        '''
        cur_emap_obj = plot_mapX.Sunmap(sunmap=cur_emap)
        cim = cur_emap_obj.contour(axes=ax, rot=0, **kwargs)
        if contour_label:
            cim.clabel(fmt='%1.1e', levels=[1, 1, 1, 1, 1])
            print(cim.levels)


def eovsa_max_position(tim, weight=False, fov=None,subed=False,high_tres=False):
    res_list = []
    in_dic = in_dic_10s_la
    if high_tres:
        in_dic = in_dic_1s
    for i in range(30):
        if not subed:
            cur_file = in_dic[tim]['radio_sbs'][i]
        else:
            if high_tres:
                cur_file = '/Volumes/Data/20170820/eovsa_bkg_subtraction/slfcal/images_slfcaled/subed_slfcaled_tb_final_XX_1s_XX_t{0}_s{1:0=2d}.fits'.format(
                        tim,
                        i+1)
            else:
                cur_file='/Volumes/Data/20170820/eovsa_bkg_subtraction_10s_3730_40/slfcal/images_slfcaled/subed_slfcaled_tb_final_XX_10s_XX_t{0}_s{1:0=2d}.fits'.format(
                        tim-184, i+1)
        crmap = smap.Map(cur_file)
        if fov is not None:
            crmap = pot.make_sub_map(crmap, fov=fov)
        cdata = sf.mapdata(crmap)
        if weight:
            cpos = scipy.ndimage.measurements.center_of_mass(cdata)
        else:
            cpos = zip(*np.where(cdata == np.nanmax(cdata)))[0]
        res_list.append(pot.get_world(cpos[0], cpos[1], crmap))
    return res_list


def weighted_center(cdata):
    cdata[np.isnan(cdata)] = 0
    # in x axis:
    tt = 0
    wtt = 0
    for xi in range(cdata.shape[0]):
        cysum = np.sum(cdata[xi, :])
        tt += cysum
        wtt += cysum * xi
    x_res = wtt / tt
    tt = 0
    wtt = 0
    for yi in range(cdata.shape[1]):
        cysum = np.sum(cdata[:, yi])
        tt += cysum
        wtt += cysum * yi
    y_res = wtt / tt
    return [x_res, y_res]


def plot_rhessi_spectrum(fits_file):
    with fits.open(fits_file, mode='readonly') as hdul:
        print(len(hdul))
        # fits_out=copy.deepcopy(hdul)
        print(hdul[0].data.shape)
    hdul.close()


def add_custom_cb(target_ax, ctitle, cmin=None, cmax=None, cmap=None, label_top=True, vertical=True, insert=False,
                  cloc=0, nbins=4, time_cb=None, within_height="30%", named_axes=False, bound=None, left_tick=False,
                  cur_im=None, discrete_cb=None, tick_color=None, diy_cmap=False,spine_color='k'):
    if cur_im is not None:
        cmin = cur_im.get_array().min()
        cmax = cur_im.get_array().min()
        new_norm = cur_im.norm
        cmap = cur_im.cmap
    else:
        new_norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
    if diy_cmap:
        def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
            new_cmap = colors.LinearSegmentedColormap.from_list(
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                cmap(np.linspace(minval, maxval, n)))
            return new_cmap

        cmap = truncate_colormap(cmap, 8. / 30., 1.0)
    if discrete_cb is not None:
        cur_bounds = np.linspace(cmin, cmax, discrete_cb[0])
        cur_tick = np.linspace(cmin, cmax, discrete_cb[1], dtype=int)
    if tick_color is not None:
        cco = tick_color
    else:
        cco = 'k'
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    if insert:
        if vertical:
            if bound is None:
                addax = inset_axes(target_ax, width="1.5%", height=within_height, loc=cloc)
            else:
                addax = target_ax.inset_axes(bound)
            ccb = mpl.colorbar.ColorbarBase(addax, norm=new_norm, cmap=cmap, orientation='vertical')
            addax.set_ylabel(ctitle, color=cco)
        else:
            addax = inset_axes(target_ax, width=within_height, height="3%", loc=cloc)
            ccb = mpl.colorbar.ColorbarBase(addax, norm=new_norm, cmap=cmap, orientation='horizontal')
            addax.set_title(ctitle, color=cco)
        addax.xaxis.set_tick_params(rotation=30)
        plt.setp(addax.spines.values(), color=spine_color)
    elif named_axes:
        addax = target_ax
        if vertical:
            ccb = mpl.colorbar.ColorbarBase(addax, norm=new_norm, cmap=cmap, orientation='vertical')
        else:
            ccb = mpl.colorbar.ColorbarBase(addax, norm=new_norm, cmap=cmap, orientation='horizontal')
        if label_top:
            addax.set_title(ctitle, color=cco)
        else:
            ccb.set_label(ctitle, color=cco)
    else:
        divider = make_axes_locatable(target_ax)
        if vertical:
            addax = divider.append_axes('right', size='5%', pad=0.05)
            ccb = mpl.colorbar.ColorbarBase(addax, norm=new_norm, cmap=cmap, orientation='vertical')
        else:
            addax = divider.append_axes('bottom', size='8%', pad=0.05)
            ccb = mpl.colorbar.ColorbarBase(addax, norm=new_norm, cmap=cmap, orientation='horizontal')
        if label_top:
            addax.set_title(ctitle, color=cco)
        else:
            ccb.set_label(ctitle, color=cco)
    if vertical and insert and left_tick:
        addax.yaxis.tick_left()
        addax.yaxis.set_label_position("left")
    if time_cb is not None:
        ccb.ax.text(.5, .3, ' UT  @ 2017-07-03', horizontalalignment='center', transform=ccb.ax.transAxes, color='k')
        ticknlabel = daterize_colorbar(time_cb, nbins)
        ccb.set_ticks(ticknlabel[0])
        ccb.set_ticklabels(ticknlabel[1])
        # loc = mdates.AutoDateLocator()
        # ccb.ax.xaxis.set_major_locator(loc)
        # ccb.ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
        ccb.ax.xaxis.set_tick_params(rotation=0)
        ccb.update_ticks()
    else:
        tick_locator = ticker.MaxNLocator(nbins=nbins)
        ccb.locator = tick_locator
        ccb.update_ticks()
    '''
    ccb.ax.spines['left'].set_color(spine_color)
    ccb.ax.spines['right'].set_color(spine_color)
    ccb.ax.spines['up'].set_color(spine_color)
    ccb.ax.spines['down'].set_color(spine_color)
    addax.spines['left'].set_color(spine_color)
    addax.spines['right'].set_color(spine_color)
    addax.spines['up'].set_color(spine_color)
    addax.spines['down'].set_color(spine_color)
    '''
    plt.setp(ccb.ax.spines.values(), color=spine_color)
    plt.setp(addax.spines.values(), color=spine_color)


def plot_eovsa_max(tim, cax, cmap, minmax):
    # mp_list = eovsa_max_position(tim-100, weight=False, fov=fov)
    mp_list = eovsa_peak_list[tim - 100, :, :]
    print(type(mp_list))
    norm = mpl.colors.Normalize(vmin=minmax[0], vmax=minmax[1])
    for imp in range(mp_list.shape[0]):
        cur_world = pot.get_world(mp_list[imp, 0], mp_list[imp, 1], cur_map=smap.Map(in_dic[tim]['radio_sbs'][imp]))
        cax.plot(cur_world[0], cur_world[1], 'o', markersize='6.', color=cmap(norm(int(imp))))

#def gst_image_correction_1st_round(inp_data):

def gst_image_correction_2nd_round(tim, t_diff=-20000, t_value=12000, apply_first_round=False):
    gmap = smap.Map(in_dic[tim]['bbso'])
    gdata = copy.deepcopy(sf.mapdata(gmap))
    # gdata=gmap.data
    new_data = copy.deepcopy(gdata)
    new_data = new_data.astype(np.float64)
    print(np.nanmin(gdata))
    # fig1 = plt.figure(figsize=(8, 4))
    # ax11 = fig1.add_subplot(1, 2, 1)
    # ax21 = fig1.add_subplot(1, 2, 2)
    for yi in range(gdata.shape[1]):
        # if yi == 1154:
        # ax11.plot(np.arange(2555), gdata[:, yi])
        cdata = copy.deepcopy(gdata[:, yi])
        cur_tag = np.zeros((len(cdata) - 1))
        cswitch = 0
        for xi, cx in enumerate(cdata):
            if xi == 0: continue
            if xi < 800 - yi or xi < yi - 1400 or xi > yi + 1900 or xi > 3600 - yi: continue
            if cx - cdata[xi - 1] < t_diff and cdata[xi - 1] > t_value:
                cswitch = 1
            elif cx - cdata[xi - 1] > (-t_diff) and cx > t_value:
                cswitch = 0
            cur_tag[xi - 1] = cswitch
        # if yi == 1154:
        # ax21.plot(np.arange(2554), cur_tag)
        for xi, cx in enumerate(cdata):
            if xi == 0: continue
            if cur_tag[xi - 1] == 1:
                new_data[xi, yi] = float(cx) + 32767.
    '''
    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    ax1.imshow(gdata.transpose(), origin='lower')
    ax2.imshow(new_data.transpose(), origin='lower')
    '''
    return smap.Map(new_data.transpose(), gmap.meta)


def save_cored_ha_data_to_fitsfile():
    for idi, cdic in enumerate(in_dic):
        cbmap = gst_image_correction_2nd_round(tim=idi)
        cur_name = '/Volumes/Data/20170906/bbso/corrected/T{0:0=3d}.fits'.format(idi)
        write_file(fname=cur_name, data=cbmap.data, header=cbmap.meta, filetype='fits')


# def hessi_lightcurve(savfile,cax,trange=None,hessi_smoth=0):
def hessi_lightcurve(savfile, cax=None, trange=None, hessi_smoth=0):
    if trange is not None:
        tr_plt = [Time(in_dic[trange[0]]['time'], format='iso').isot,
                  Time(in_dic[trange[1]]['time'], format='iso').isot]
    else:
        tr_plt = [Time(in_dic[0]['time'], format='iso').isot, Time(in_dic[-1]['time'], format='iso').isot]
    hessi = readsav(savfile)
    hessi_tim = Time(list(hessi['obs_times_str']))
    hessi_tim_plt = hessi_tim.plot_date
    tidx_hessi, = np.where((hessi_tim >= tr_plt[0]) & (hessi_tim <= tr_plt[1]))
    if cax is None:
        fig = plt.figure()
        cax = fig.add_subplot(111)
    for idx, eg in enumerate(hessi['obs_energies']):
        flux_hessi = ma.masked_array(hessi['obs_data'][tidx_hessi, idx])
        if hessi_smoth > 0:
            flux_hessi = DButil.smooth(flux_hessi, 20)
            flux_hessi = flux_hessi / np.nanmax(flux_hessi)
            cax.step(hessi_tim_plt[tidx_hessi], DButil.smooth(flux_hessi, hessi_smoth),
                     label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))
        else:
            cax.step(hessi_tim_plt[tidx_hessi], flux_hessi, label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))
    cax.xaxis_date()
    # cax2.xaxis.set_major_formatter(DateFormatter("%m:%s"))
    plt.gcf().autofmt_xdate()


def plot_fit_res(tim, tim2, tim3, pfig=None):
    params_list = ['theta', 'temperature', 'delta', 'n_th', 'n_nth', 'mag']
    ylim_list = [[20, 90], [6e4, 1.5e7], [3, 9], [5e10, 5e11], [1e5, 1e11], [400, 1800]]
    title_list = ['theta', 'temperature', 'delta', 'thermal_n', 'non-thermal_n', 'B']
    res_list_file = open('/Volumes/Data/20170906/fitting_result/fit_save_all_t{}.p'.format(tim), 'rb')
    res_list = pickle.load(res_list_file)
    res_list_file.close()
    res_list_file_2 = open('/Volumes/Data/20170906/fitting_result/fit_save_all_t{}.p'.format(tim2), 'rb')
    res_list_2 = pickle.load(res_list_file_2)
    res_list_file_2.close()
    res_list_file_3 = open('/Volumes/Data/20170906/fitting_result/fit_save_all_t{}.p'.format(tim3), 'rb')
    res_list_3 = pickle.load(res_list_file_3)
    res_list_file_3.close()
    # print(res_list_2[0])
    if pfig is None:
        fig2 = plt.figure(figsize=(14, 8))
    for parmi in range(6):
        c_array = np.zeros((len(res_list)))
        e_array = np.zeros((len(res_list)))
        c_array_2 = np.zeros((len(res_list_2)))
        e_array_2 = np.zeros((len(res_list_2)))
        c_array_3 = np.zeros((len(res_list_3)))
        e_array_3 = np.zeros((len(res_list_3)))
        if pfig is None:
            ccax = fig2.add_subplot(2, 3, parmi + 1)
        else:
            ccax = pfig.add_subplot(4, 3, 7 + parmi)
        if parmi == 5:
            ext_list = []
        for pointi, res_point in enumerate(res_list):
            if parmi == 5:
                if pointi < 8:
                    ext_list.append((abs(pointi - 7) ** 1.7 / 28 + 1) * 750)
                else:
                    ext_list.append((abs(pointi - 7) ** 2 / 25 + 1) * 750)
            print(res_point)
            c_array[pointi] = res_point.params[params_list[parmi]].value
            e_array[pointi] = res_point.params[params_list[parmi]].stderr
            c_array_2[pointi] = res_list_2[pointi].params[params_list[parmi]].value
            e_array_2[pointi] = res_list_2[pointi].params[params_list[parmi]].stderr
            c_array_3[pointi] = res_list_3[pointi].params[params_list[parmi]].value
            e_array_3[pointi] = res_list_3[pointi].params[params_list[parmi]].stderr
        print(c_array)
        ccax.errorbar(range(len(res_list)), c_array, yerr=e_array, linestyle='None', marker='x')
        if pointi != 1:
            ccax.errorbar(range(len(res_list_2)), c_array_2, yerr=e_array_2, linestyle='None', marker='x', color='r')
            ccax.errorbar(range(len(res_list_3)), c_array_3, yerr=e_array_3, linestyle='None', marker='x', color='g')

        if parmi == 5:
            ccax.plot(range(len(res_list)), ext_list, linestyle='None', marker='o')
        ccax.set_ylim(ylim_list[parmi])
        ccax.set_yscale('log')
        ccax.set_title(title_list[parmi])
    if pfig is None:
        plt.show()

def hessi_maps_20170906(like_eovsa=False):
    # fits_file='/Volumes/Data/20170906/rhessi/images/hsi_imagecube_clean_20170906_1914_15tx4e.fits'
    # fits_file='/Volumes/Data/20170906/rhessi/images/hsi_imagecube_17tx2e_20170906_192000_snr15.fits'
    fits_file = '/Users/walter/hsi_imagecube_90tx1e_20170820_191900.fits'
    # fits_file='/Volumes/Data/20170906/rhessi/hsi_imagecube_66tx4e_20170906_192000.fits'
    if like_eovsa:
        ref_map = smap.Map(in_dic[0]['radio_sbs'][0])
    else:
        ref_map = smap.Map(fits_file)
    curinit = Time('2017-09-06T19:18:00', format='isot')
    timed = TimeDelta(10, format='sec')
    '''
    time_list=[[' 6-Sep-2017 19:14:44.200', ' 6-Sep-2017 19:16:30.344'],
    [' 6-Sep-2017 19:16:30.344', ' 6-Sep-2017 19:18:19.344'], [' 6-Sep-2017 19:18:19.344',
    ' 6-Sep-2017 19:19:58.344'], [' 6-Sep-2017 19:19:58.344', ' 6-Sep-2017 19:21:37.344'],
    [' 6-Sep-2017 19:21:37.344', ' 6-Sep-2017 19:23:10.344'], [' 6-Sep-2017 19:23:10.344',
    ' 6-Sep-2017 19:23:26.344'], [' 6-Sep-2017 19:23:43.160', ' 6-Sep-2017 19:25:39.421'],
    [' 6-Sep-2017 19:25:39.421', ' 6-Sep-2017 19:27:39.421'], [' 6-Sep-2017 19:28:10.180',
    ' 6-Sep-2017 19:30:10.800'], [' 6-Sep-2017 19:30:10.800', ' 6-Sep-2017 19:32:09.800'],
    [' 6-Sep-2017 19:32:11.800', ' 6-Sep-2017 19:34:11.800'], [' 6-Sep-2017 19:34:11.800',
    ' 6-Sep-2017 19:36:08.800'], [' 6-Sep-2017 19:36:08.800', ' 6-Sep-2017 19:37:57.800'],
    [' 6-Sep-2017 19:37:57.800', ' 6-Sep-2017 19:39:41.800'], [' 6-Sep-2017 19:39:41.800',
    ' 6-Sep-2017 19:40:45.800']]
    '''
    file1 = open('/Volumes/Data/IDL/time_log', 'r')
    lines = file1.readlines()
    tr_list = []
    for ii in range(21):
        tisot1 = lines[ii * 2 + 2].replace('06-Sep-17 ', '2017-09-06T').rstrip('\r\n')
        tisot2 = lines[ii * 2 + 3].replace('06-Sep-17 ', '2017-09-06T').rstrip('\r\n')
        tr_list.append([tisot1, tisot2])
    # e_list=[[6.0,12.0],[12.0,25.0],[25.0,50.0],[50.0,100.0]]
    e_list = [[6.0, 11.0], [11.0, 23.0]]
    # e_list=[[6.0,12.0],[12.0,25.0]]
    map_list = []
    with fits.open(fits_file, mode='readonly') as hdul:
        for i in range(2):
            cur_map_list = []
            for ii in range(21):
                # print(e_list[i],tr_list[ii])
                cur_meta = copy.deepcopy(ref_map.meta)
                print(cur_meta)
                if like_eovsa:
                    re_list = ['NAXIS1', 'NAXIS2', 'CDELT1', 'CDELT2', 'DATE', 'DATE-OBS']
                    cmeta = copy.deepcopy(smap.Map(fits_file).meta)
                    for rk in re_list:
                        cur_meta[rk] = cmeta[rk]
                        cur_meta['CRPIX1'] = 56.0
                        cur_meta['CRPIX2'] = 56.0
                        cur_meta['CRVAL1'] = cmeta['xcen']
                        cur_meta['CRVAL2'] = cmeta['ycen']
                else:
                    cur_meta['date_obs'] = tr_list[ii][0]
                    cur_meta['date_end'] = tr_list[ii][1]
                    cur_meta['energy_l'] = e_list[i][0]
                    cur_meta['energy_h'] = e_list[i][1]
                    cur_meta['wavelnth'] = ''
                cur_meta['wavelnth'] = ''
                cexp_time = Time(tr_list[ii][1], format='isot').mjd - Time(tr_list[ii][0], format='isot').mjd
                print(cur_meta)
                print(cexp_time)
                new_map_data = hdul[0].data[ii, i, :, :] * (0.0008680555555555555 / cexp_time)
                if ii > 14:
                    if i == 0:
                        new_map_data = new_map_data * 3.5
                        # new_map_data = new_map_data * 7
                    if i == 1:
                        new_map_data = new_map_data * 30.0
                # print(np.nanmax(hdul[0].data[ii,i,:,:]))
                # new_map=smap.Map(hdul[0].data[ii,i,:,:],cur_meta)
                new_map = smap.Map(new_map_data, cur_meta)
                # print(np.nanmax(new_map.data))
                cur_map_list.append(new_map)
            map_list.append(cur_map_list)
    hdul.close()
    return map_list

def hessi_maps(like_eovsa=False):
    # fits_file='/Volumes/Data/20170906/rhessi/images/hsi_imagecube_clean_20170906_1914_15tx4e.fits'
    # fits_file='/Volumes/Data/20170906/rhessi/images/hsi_imagecube_17tx2e_20170906_192000_snr15.fits'
    #fits_file = '/Volumes/Data/20170820/hessi/images/hsi_imagecube_90tx1e_12_25_det_1_3_6_20170820_191900.fits'
    #fits_file = '/Volumes/Data/20170820/hessi/images/hsi_imagecube_90tx1e_20170820_191900_det_1_3_6_6_12.fits'
    #fits_file = '/Volumes/Data/20170820/hessi/images/hsi_imagecube_90tx1e_20170820_191900_det_1_3_6_6_12.fits'
    fits_file = '/Volumes/Data/20170820/hessi/images/hsi_imagecube_90tx1e_20170820_191900_25_50.fits'
    # fits_file='/Volumes/Data/20170906/rhessi/hsi_imagecube_66tx4e_20170906_192000.fits'
    if like_eovsa:
        ref_map = smap.Map(in_dic_10s[0]['radio_sbs'][0])
    else:
        ref_map = smap.Map(fits_file)
    init_t = Time('2017-08-20T19:19:00', format='isot')
    time_interval = 20.0
    init_t=Time(init_t,format='isot')
    tr_list = []
    for tint in range(90):
        timed1=TimeDelta((tint*time_interval)*1.0, format='sec')
        timed2=TimeDelta(time_interval, format='sec')
        start_time=init_t+timed1
        end_time=start_time+timed2
        tr_list.append([start_time.isot, end_time.isot])
    #for ii in range(21):
    #    tisot1 = lines[ii * 2 + 2].replace('06-Sep-17 ', '2017-09-06T').rstrip('\r\n')
    #    tisot2 = lines[ii * 2 + 3].replace('06-Sep-17 ', '2017-09-06T').rstrip('\r\n')
    #    tr_list.append([tisot1, tisot2])
    # e_list=[[6.0,12.0],[12.0,25.0],[25.0,50.0],[50.0,100.0]]
    #e_list = [[6.0, 11.0], [12.0, 25.0]]
    #e_list = [[6.0, 12.0]]
    e_list = [[25.0, 50.0]]
    # e_list=[[6.0,12.0],[12.0,25.0]]
    map_list = []
    with fits.open(fits_file, mode='readonly') as hdul:
        for i in range(1):
            cur_map_list = []
            for ii in range(90):
                # print(e_list[i],tr_list[ii])
                cur_meta = copy.deepcopy(ref_map.meta)
                if like_eovsa:
                    re_list = ['NAXIS1', 'NAXIS2', 'CDELT1', 'CDELT2', 'DATE', 'DATE-OBS']
                    cmeta = copy.deepcopy(smap.Map(fits_file).meta)
                    for rk in re_list:
                        cur_meta[rk] = cmeta[rk]
                        cur_meta['CRPIX1'] = 56.0
                        cur_meta['CRPIX2'] = 56.0
                        cur_meta['CRVAL1'] = cmeta['xcen']
                        cur_meta['CRVAL2'] = cmeta['ycen']
                else:
                    cur_meta['date_obs'] = tr_list[ii][0]
                    cur_meta['date_end'] = tr_list[ii][1]
                    cur_meta['energy_l'] = e_list[i][0]
                    cur_meta['energy_h'] = e_list[i][1]
                #cexp_time = Time(tr_list[ii][1], format='isot').mjd - Time(tr_list[ii][0], format='isot').mjd
                '''
                print(cur_meta)
                print(cexp_time)
                new_map_data = hdul[0].data[ii, i, :, :] * (0.0008680555555555555 / cexp_time)
                if ii > 14:
                    if i == 0:
                        new_map_data = new_map_data * 3.5
                        # new_map_data = new_map_data * 7
                    if i == 1:
                        new_map_data = new_map_data * 30.0
                '''
                # print(np.nanmax(hdul[0].data[ii,i,:,:]))
                # new_map=smap.Map(hdul[0].data[ii,i,:,:],cur_meta)
                cur_meta.pop('wavelnth')
                new_map = smap.Map(hdul[0].data[ii, i, :, :], cur_meta)
                print(new_map.meta)
                #new_map_name = '/Volumes/Data/20170820/hessi/images/hsi_12_25_det_1_3_6_20170820_t{0:0=3d}.fits'.format(ii)
                #new_map.save(new_map_name)
                # print(np.nanmax(new_map.data))
                cur_map_list.append(new_map)
            map_list.append(cur_map_list)
    hdul.close()
    pickle.dump(map_list,open('/Volumes/Data/20170820/hessi/images/20_50_all90maps.p','wb'))
    return map_list

def hessi_map_time_list(hmaps):
    res = np.zeros(len(hmaps[0]))
    for i, hmap in enumerate(hmaps[0]):
        res[i] = (Time(hmap.meta['date_obs'], format='isot').mjd + Time(hmap.meta['date_end'], format='isot').mjd) / 2
    return res


def hessi_movie(tim, mjd_list=None, hmaps=None):
    if hmaps is None:
        hmaps = hessi_maps()
    if hmaps is not None and mjd_list is None:
        mjd_list = hessi_map_time_list(hmaps)
    plt.ioff()
    fig = plt.figure(figsize=(6, 12))
    c_residual = [abs(xm - Time(in_dic[tim]['time'], format='iso').mjd) for xm in mjd_list]
    rhessi_index = np.nanargmin(c_residual)
    print(Time(mjd_list[rhessi_index], format='mjd').isot, in_dic[tim]['time'])
    # kwlist=['Z.304.','Z.171.','Z.94.','Z.1600.']
    kwlist = ['Z.193.', 'Z.131.', 'Z.94.', 'Z.1600.']
    fov2 = [[580, -253], [607, -211]]
    # 50% and 80% contour level
    # level_list=[40,6]
    # level_list=np.array([0.3,0.5,0.9])
    level_list = np.array([0.6, 0.8, 0.9])
    hessi_peak_0 = 150.0
    hessi_peak_1 = 30.0
    # level_list=[100,30]
    cur_time = Time(in_dic[tim]['time'], format='iso').mjd
    for kwi, kw in enumerate(kwlist):
        cax = fig.add_subplot(3, 2, kwi + 1)
        print(in_dic[tim][kw])
        justmap(cax=cax, tim=tim, kw=kw, fov=fov2)
        if kwi == 3:
            plot_eovsa_contourf(tim=tim, ax=cax, fov=fov2, spws=range(9, 30),
                                abs_contourf=False, level=0.95,
                                alt_cmap=plt.cm.Spectral)
            cax.text(.75, .8, 'EOVSA\nabove 8GHz', horizontalalignment='center', transform=cax.transAxes, color='y')
            add_custom_cb(target_ax=cax, cmin=3, cmax=18, cmap=plt.cm.Spectral,
                          ctitle='GHz')
        # im1 = pot.plot_rhessi_contour(erange=0, trange=rhessi_index, ax=cax, fov=fov2, levels=[0.6, 0.92], inp_color='r',maplist=hmaps)
        # im2 = pot.plot_rhessi_contour(erange=1, trange=rhessi_index, ax=cax, fov=fov2, levels=[0.6, 0.92], inp_color='g',maplist=hmaps)
        # im3 = pot.plot_rhessi_contour(erange=2, trange=tim-100, ax=cax, fov=fov2, levels=[0.25, 0.8], inp_color='b',maplist=hmaps)
        # im1 = pot.plot_rhessi_contour(erange=0, trange=rhessi_index, ax=cax, fov=fov2, abs_contour=True,levels=[level_list[0]*0.4,level_list[0]*0.92], inp_color='r',
        # maplist=hmaps,plot_centroid=True)
        print('for debug!!!!!', np.nanmax(hmaps[1][rhessi_index].data))
        im1 = pot.plot_rhessi_contour(erange=0, trange=rhessi_index, ax=cax, fov=fov2,
                                      abs_contour=False,
                                      levels=level_list,
                                      inp_color='r',
                                      maplist=hmaps, plot_centroid=False)
        im1 = pot.plot_rhessi_contour(erange=1, trange=rhessi_index, ax=cax, fov=fov2,
                                      abs_contour=False,
                                      levels=level_list,
                                      inp_color='g',
                                      maplist=hmaps, plot_centroid=False)
        # im2 = pot.plot_rhessi_contour(erange=1, trange=rhessi_index, ax=cax, fov=fov2, abs_contour=True,levels=[level_list[1]*0.4,level_list[1]*0.92], inp_color='g',
        #                              maplist=hmaps,plot_centroid=True)
        # im3 = pot.plot_rhessi_contour(erange=2, trange=rhessi_index, ax=cax, fov=fov2, abs_contour=True,levels=[level_list[2][0],level_list[2][1]], inp_color='b',
        #                              maplist=hmaps)

        legend_elements = [Line2D([0], [0], color='r', lw=0.2, label='6~12KeV,  40%N92%'),
                           Line2D([0], [0], color='g', lw=0.2, label='12~25KeV, 40%N92%')]
        if kwi == 0:
            cleg = cax.legend(handles=legend_elements, loc='best', fontsize='5')
    with open('/Volumes/Data/20170906/tpower/from_big_fov.p', 'rb') as dspec_save:
        aspec = pickle.load(dspec_save)
    dspec_save.close()
    dsax = fig.add_subplot(3, 1, 3)
    mdds.plt_dspec(specdata=aspec, timerange=[79, 248], pax=dsax, time_line=[Time(in_dic[tim]['time'], format='iso')],
                   rhessisave='/Volumes/Data/20170906/rhessi/lc/corrected_hessi.sav')
    # cur_name='/Volumes/Data/20170906/plotting/hessi/abs_sec_part_hessi_p40_92_6_11_23_T{0:0=3d}.png'.format(tim)
    cur_name = '/Volumes/Data/20170906/plotting/hessi/selected/hessi_p40_6_11_23_T{0:0=3d}.png'.format(tim)
    fig.suptitle(in_dic[tim]['time'])

    plt.savefig(cur_name)
    # plt.show()


def sdo_downloader(t_begin, t_end, t_end_hmi):
    print('use this GD form: 2010/12/06 00:00:00')
    from sunpy.net.vso import VSOClient
    client = VSOClient()
    wave_list = ['94', '131', '171', '193', '211', '304', '335', '1600', '1700']
    for cwave in wave_list:
        query_response = client.query_legacy(tstart=t_begin, tend=t_end, instrument='AIA',
                                             wave=cwave, sample=3600)
        if len(query_response) > 1:
            query_response = query_response[0]
        cur_folder = t_begin.replace('/', '_').replace(' ', 'T').replace(':', '_')
        #results = client.get(query_response, path='/Volumes/Data/20170906/aia_segment/' + cur_folder, site='rob')
        results = client.get(query_response, path='/Volumes/WD6T/working/20170703/aia_helioviewer_movie/previous/' + cur_folder, site='rob')
        files = results.wait()
    client = VSOClient()
    query_response = client.query_legacy(tstart=t_begin, tend=t_end_hmi, instrument='HMI',
                                         physobs='los_magnetic_field', sample=3600)
    if len(query_response) > 1:
        query_response = query_response[0]
    #results = client.get(query_response, path='/Volumes/Data/20170906/aia_segment/' + cur_folder, site='rob')
    results = client.get(query_response, path='/Volumes/WD6T/working/20170703/aia_helioviewer_movie/previous/' + cur_folder, site='rob')
    files = results.wait()


def tmp_dl_sdo():
    cur_time = Time('2017-09-06 12:05:00.000', format='iso')
    bi_timed = TimeDelta(600, format='sec')
    si_timed = TimeDelta(20, format='sec')
    hmi_timed = TimeDelta(730, format='sec')
    for i in range(23):
        sdo_downloader(t_begin=cur_time.iso.replace('-', '/'), t_end=(cur_time + si_timed).iso.replace('-', '/'),
                       t_end_hmi=(cur_time + hmi_timed).iso.replace('-', '/'))
        cur_time = cur_time + bi_timed


def multi_rhessi_contour_in_single_ax(tim):
    hmaps = hessi_maps()
    mjd_list = hessi_map_time_list(hmaps)
    fig = plt.figure(figsize=(12, 6))
    c_residual = [abs(xm - Time(in_dic[tim]['time'], format='iso').mjd) for xm in mjd_list]
    rhessi_index = np.nanargmin(c_residual)
    # fov2 = [[580, -253], [607, -211]]
    fov2 = [[560, -273], [627, -191]]
    cax1 = fig.add_subplot(2, 2, 1)
    cax3 = fig.add_subplot(2, 2, 3)
    justmap(cax=cax1, tim=tim, kw='Z.1600.', fov=fov2)
    justmap(cax=cax3, tim=tim, kw='Z.94.', fov=fov2)
    # im1 = pot.plot_rhessi_contour(erange=0, trange=rhessi_index, ax=cax, fov=fov2, abs_contour=True,
    #                              levels=[level_list[0] * 0.4, level_list[0] * 0.92], inp_color='r',
    #                              maplist=hmaps)
    print(rhessi_index)
    diff_map_1 = smap.Map((hmaps[0][rhessi_index].data - hmaps[0][rhessi_index - 1].data), hmaps[0][rhessi_index].meta)
    diff_map_2 = smap.Map((hmaps[1][rhessi_index].data - hmaps[1][rhessi_index - 1].data), hmaps[1][rhessi_index].meta)
    sub_map_1 = pot.make_sub_map(cur_map=diff_map_1, fov=fov2)
    sub_map_2 = pot.make_sub_map(cur_map=diff_map_2, fov=fov2)
    kwargs = {'cmap': 'jet'}
    cax2 = fig.add_subplot(2, 2, 2)
    cax4 = fig.add_subplot(2, 2, 4)
    sc_map_obj_1 = plot_mapX.Sunmap(sunmap=sub_map_1)
    sc_map_obj_2 = plot_mapX.Sunmap(sunmap=sub_map_2)
    sc_map_obj_1.imshow(axes=cax2, **kwargs)
    sc_map_obj_2.imshow(axes=cax4, **kwargs)
    plt.show()


def magnetogram_footpoint(tim):
    fig = plt.figure(figsize=(8, 8))
    anno_gap = 2
    fov2 = [[560, -263], [627, -191]]
    # fov2 = [[575, -253], [625, -190]]
    # point_list=[[582,-223],[585,-231],[587,-240],[597,-215],[587,-228]]
    # point_list=[[587,-240],[582,-223],[585,-231],[597,-215],[587,-228]]
    # ------------------------
    # point_list=[[588,-246],[582,-223],[585,-231],[597,-215],[587,-228],[591,-226],[588,-235]]
    point_list = [[587, -240], [582, -223], [585, -231], [597, -215], [587, -228], [591, -226], [588, -235]]
    # point_list=[[590,-249],[582,-223],[585,-231],[597,-215],[587,-228],[591,-226],[588,-235]]
    # ----------------------
    # point_list=[[582,-223],[587,-240],[588,-246],[590,-249],[585,-231],[597,-215],[617,-202],[589,-234],[587,-228],[591,-226]]
    # point_list=[[590,-249],[588,-246],[587,-240],[582,-223],[585,-231],[597,-215],[617,-202],[589,-234],[587,-228],[591,-226]]
    # color_list=['coral','blue','crimson','darkblue']
    # label_list=['P1','P2','P3','P4','N1','N2','N3','L1','L2','L3']
    label_list = ['P1', 'P2', 'N1', 'N2', 'L1', 'L2', 'L3']
    # kw_st=['Z.1600.','Z.1600.','Z.1600.','Z.1600.','Z.1600.','Z.1600.','Z.1600.','Z.94.','Z.94.','Z.94.']
    kw_st = ['Z.1600.', 'Z.1600.', 'Z.1600.', 'Z.1600.', 'Z.94.', 'Z.94.', 'Z.94.']
    # color_list=[plt.cm.Reds(0.5),plt.cm.Reds(0.6),plt.cm.Reds(0.7),plt.cm.Reds(0.9),plt.cm.Blues(0.5),plt.cm.Blues(0.7),plt.cm.Blues(0.9),plt.cm.Purples(0.5),plt.cm.Purples(0.7),plt.cm.Purples(0.9)]
    color_list = [plt.cm.Reds(0.5), plt.cm.Reds(0.9), plt.cm.Blues(0.5), plt.cm.Blues(0.9), plt.cm.Purples(0.9),
                  plt.cm.Purples(0.75), plt.cm.Purples(0.4)]
    cax1 = fig.add_subplot(3, 3, 2)
    cax3 = fig.add_subplot(3, 3, 1)
    cax2 = plt.subplot2grid((3, 3), (1, 1), colspan=2, rowspan=1)
    cax4 = plt.subplot2grid((3, 3), (2, 1), colspan=2, rowspan=1)
    # cax2=fig.add_subplot(4,2,2)
    im_list = []
    cax5 = fig.add_subplot(3, 3, 3)
    # cax6 = fig.add_subplot(4, 2, 6)
    # cax8 = fig.add_subplot(4, 2, 8)
    # cax4 = fig.add_subplot(4, 2, 4)
    justmap(cax=cax1, tim=tim, kw='hmi', fov=fov2)
    justmap(cax=cax3, tim=tim, kw='Z.1600.', fov=fov2)
    # justmap(cax=cax3,tim=tim,kw='hmi',fov=fov2,alp=0.25)
    plot_B_contour(tim, cax3, [-100, 100], cfov=fov2, alp=0.4)
    justmap(cax=cax5, tim=tim, kw='Z.94.', fov=fov2)
    '''
    for ccpi,ccp in enumerate(point_list):
        for ccax in [cax1,cax3]:
            if ccpi<7:
                ccax.plot(ccp[0], ccp[1], 's', markersize='6.', color=color_list[ccpi], markerfacecolor='none')
        if ccpi >=7:
            cax5.plot(ccp[0], ccp[1], 's', markersize='6.', color=color_list[ccpi], markerfacecolor='none')
        print(kw_st[ccpi])
        if ccpi<4:
            cim,=pot.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p', x=ccp[0], y=ccp[1],
                          timerange=[88, 198], kw_list=[kw_st[ccpi]], sample_range=2, normq=True,
                          cur_color=color_list[ccpi],pax=cax2,ltex=label_list[ccpi])
        elif ccpi<7:
            cim,=pot.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p', x=ccp[0], y=ccp[1],
                          timerange=[88, 198], kw_list=[kw_st[ccpi]], sample_range=2, normq=True,
                          cur_color=color_list[ccpi],pax=cax6,ltex=label_list[ccpi])
        else:
            cim,=pot.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p', x=ccp[0], y=ccp[1],
                          timerange=[88, 198], kw_list=[kw_st[ccpi]], sample_range=2, normq=True,
                          cur_color=color_list[ccpi],pax=cax8,ltex=label_list[ccpi])
    '''
    for ccpi, ccp in enumerate(point_list):
        for ccax in [cax1, cax3]:
            if ccpi < 4:
                ccax.plot(ccp[0], ccp[1], 's', markersize='6.', color=color_list[ccpi], markerfacecolor='none')
                ccax.text(ccp[0] + anno_gap, ccp[1] + anno_gap, label_list[ccpi], fontsize=10, color=color_list[ccpi])
        if ccpi >= 4:
            cax5.plot(ccp[0], ccp[1], 's', markersize='6.', color=color_list[ccpi], markerfacecolor='none')
            cax5.text(ccp[0] + anno_gap, ccp[1] + anno_gap, label_list[ccpi], fontsize=10, color=color_list[ccpi])
        print(kw_st[ccpi])
        cim, = pot.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p', x=ccp[0], y=ccp[1],
                              timerange=[88, 208], kw_list=[kw_st[ccpi]], sample_range=2, normq=False,
                              cur_color=color_list[ccpi], pax=cax2, ltex=label_list[ccpi])
    cax2.xaxis.set_major_locator(plt.MaxNLocator(5))
    locs = cax2.get_xticks()
    '''
    cim, = pot.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p', x=596, y=-222,
                          timerange=[88, 198], kw_list=['Z.1600.'], sample_range=134, normq=True,
                          cur_color='k', pax=cax2, ltex='AIA 1600 A all fov')
    '''
    with open('/Volumes/Data/20170906/tpower/from_big_fov.p', 'rb') as dspec_save:
        aspec = pickle.load(dspec_save)
    dspec_save.close()
    mdds.plt_dspec(specdata=aspec, timerange=[88, 208], pax=cax4,
                   rhessisave='/Volumes/Data/20170906/rhessi/lc/corrected_hessi.sav')
    cax4.set_xticks([])
    cax4.set_xticklabels([])

    # justmap(cax=cax5, tim=tim, kw='Z.94.', fov=fov2)
    # cax5.plot(589,-234 , 's', markersize='6.', color='y', markerfacecolor='none')
    # cim, = pot.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p', x=587, y=-228,
    # cim, = pot.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p', x=589, y=-234,
    #                                            timerange=[88, 198], kw_list=['Z.94.'], sample_range=3, normq=True,
    #                      cur_color='y', pax=cax2, ltex='AIA 94 A')
    cax2.legend()
    # cax6.legend()
    # cax8.legend()
    cax2.set_title('light curve of AIA 1600' + r'$\AA$' + 'and 94' + r'$\AA$')
    # cax6.set_title('light curve of AIA 1600 A N1~N3')
    # cax8.set_title('light curve of AIA 1600 A L1~L3')
    cur_time = Time(in_dic[tim]['time'], format='iso')
    # cax2.axvline(cur_time.plot_date, linewidth=2, color='b')
    cur_time_1 = Time(in_dic[157]['time'], format='iso')
    cur_time_2 = Time(in_dic[166]['time'], format='iso')
    cax2.axvline(cur_time_1.plot_date, linewidth=2, color='b')
    cax2.axvline(cur_time_2.plot_date, linewidth=2, color='b')
    # cax2.axvspan(cur_time_1, cur_time_2, color='red', alpha=0.5)
    # cax6.axvline(cur_time.plot_date, linewidth=2, color='b')
    # cax8.axvline(cur_time.plot_date, linewidth=2, color='b')
    plt.show()


# def custom_lightcurve_for_dspec():

def ttmm(maps_save):
    for i in range(2):
        for ii, cur_map in enumerate(maps_save[i]):
            curname = '/Volumes/Data/20170906/rhessi/images/single_files/like_eovsa/hessi_image_{0}_{1:0=2d}.fits'.format(
                i, ii)
            subprocess.call('cp {0} {1}'.format(in_dic[0]['radio_sbs'][0], curname), shell=True)
            with fits.open(curname, mode='update') as chdul:
                re_list = ['NAXIS1', 'NAXIS2', 'CDELT1', 'CDELT2', 'DATE', 'DATE-OBS']
                cmeta = copy.deepcopy(maps_save[i][ii].meta)
                for rk in re_list:
                    chdul[0].header.set(rk, cmeta[rk])
                chdul[0].header.set('CRPIX1', 56.0)
                chdul[0].header.set('CRPIX2', 56.0)
                chdul[0].header.set('CRVAL1', cmeta['xcen'])
                chdul[0].header.set('CRVAL2', cmeta['ycen'])
                odi = maps_save[i][ii].data.shape
                chdul[0].data = maps_save[i][ii].data.reshape((1, 1, odi[0], odi[1]))
                chdul.flush()
                chdul.close()

    '''
    #maps_save = pickle.load(open('/Volumes/Data/20170906/rhessi/images/single_files/maps.p','rb'),encoding='latin1')
    res_out = np.zeros((2,len(maps_save[0]), 2))
    for i in range(2):
        for ii,cur_map in enumerate(maps_save[i]):
            #cur_map.meta['wavelnth']=cur_map.meta['wavelnth'][0]
            curname='/Volumes/Data/20170906/rhessi/images/single_files/like_eovsa/hessi_image_{0}_{1:0=2d}.fits'.format(i,ii)
            cur_map.save(curname)
    '''


def rhessi_cent(hmaps):
    rhessi_max = (([48.8, 48.9], [52.47, 52.90]),
                  ([49.45, 49.56], [53.53, 54.42], [54.96, 49.67]),
                  ([49.57, 49.42], [53.53, 54.42]),
                  ([49.53, 49.36], [52.60, 53.74]),
                  ([48.84, 51.00], [53.02, 53.54]),
                  ([51.93, 50.36], [53.84, 55.47]))
    rh_cen = pickle.load(open('/Volumes/Data/20170906/rhessi/images/hessi_pimfit_peak.p', 'rb'),encoding='latin1')
    res = np.zeros_like(rh_cen)
    for i in range(2):
        for ii in range(rh_cen.shape[1]):
            if i == 0 and ii >= 15:
                print(rhessi_max[ii - 15][0])
                w_coor = hmaps[i][ii].pixel_to_world(rhessi_max[ii - 15][0][0] * u.pix,
                                                     rhessi_max[ii - 15][0][1] * u.pix)
            else:
                w_coor = hmaps[i][ii].pixel_to_world(rh_cen[i][ii][0] * u.pix, rh_cen[i][ii][1] * u.pix)
            res[i][ii][0] = w_coor.Tx.value
            res[i][ii][1] = w_coor.Ty.value
    print(res)
    return res

    # since 15
    rhessi_max = (([48.8, 48.9], [52.47, 52.90]),
                  ([49.45, 49.56], [53.53, 54.42], [54.96, 49.67]),
                  ([49.57, 49.42]),
                  ([49.53, 49.36], [52.60, 53.74]),
                  ([48.84, 51.00], [53.02, 53.54]),
                  ([51.93, 50.36], [53.84, 55.47]))


def multi_wl_aia_map(tim, ax, kw_list=None):
    fov2 = [[560, -263], [627, -191]]
    if kw_list is None:
        kw_list = ['Z.1600.', 'Z.94.']
    cmap_list = ['Greens', 'Blues']
    alp_list = [1, 0.53]
    norm_list_arr = [[]]
    norm_list = []
    # for ii in range(len(kw_list)):
    # cur_norm=SymLogNorm(linthresh=0, vmin=norm_list_arr[ii][0], vmax=norm_list_arr[ii][1])
    # norm_list.append(cur_norm)
    for i, ckw in enumerate(kw_list):
        if i != 1: continue
        # justmap(ax,tim,fov=fov2,alp=0.4,custom_color=cmap_list[i],inorm=norm_list[i])
        justmap(ax, tim, kw=ckw, fov=fov2, alp=alp_list[0], custom_color=cmap_list[i])
        # justmap(ax,tim,ckw,fov=fov2,alp=alp_list[i])


def lightcurve_over_ribbon(kw, yrange):
    #kwlist = ['north', 'south', 'east']
    save_file = '/Volumes/Data/20170820/stack/long/100_300_1600_ribbon_{0}_wrap.p'.format(kw)
    opfile = open(save_file, 'rb')
    cdspec = pickle.load(opfile)
    opfile.close()
    world_range = np.zeros((2))
    world_range[0] = np.nanargmin([abs(yrange[0] - cpy) for cpy in cdspec['y']])
    world_range[1] = np.nanargmin([abs(yrange[1] - cpy) for cpy in cdspec['y']])
    res = np.zeros((len(cdspec['x'])-1))
    for ri in range(len(res)):
        res[ri] = np.nanmax(cdspec['dspec'][int(world_range[0]):int(world_range[1]), ri])
    #fig, ax = plt.subplots()
    # ax.plot(np.arange(len(res)),res)
    #ax.plot(np.arange(len(res)), res)
    #plt.show()
    return (res, cdspec['x'][:-1])


def stp_maximum(kw,yrange,tar_tim = None):
    save_file = '/Volumes/Data/20170820/stack/long/100_300_1600_ribbon_{0}_wrap.p'.format(kw)
    opfile = open(save_file, 'rb')
    cdspec = pickle.load(opfile)
    opfile.close()
    res = np.zeros((len(cdspec['x'])-1, 2))
    world_range = np.zeros((2))
    world_range[0] = np.nanargmin([abs(yrange[0] - cpy) for cpy in cdspec['y']])
    world_range[1] = np.nanargmin([abs(yrange[1] - cpy) for cpy in cdspec['y']])
    res = np.zeros((len(cdspec['x'])-1, 2))
    for ri in range(len(res)):
        #res[ri] = np.nanmax(cdspec['dspec'][int(world_range[0]):int(world_range[1]), ri])
        tmp_arg = np.nanargmax(cdspec['dspec'][int(world_range[0]):int(world_range[1]), ri])
        res[ri,:] = [cdspec['cutslit']['cutslit']['xcen'][tmp_arg + int(world_range[0])],cdspec['cutslit']['cutslit']['ycen'][tmp_arg +
                                                                                                        int(world_range[0])]]
    if tar_tim is not None:
        res_list = []
        for ctar_tim in tar_tim:
            tar = Time(in_dic_10s[ctar_tim]['time'],format='iso').plot_date
            tar_index = np.nanargmin([abs(tar - cxt) for cxt in cdspec['x']])
            res_list.append(res[tar_index,:])
        return res_list
    return (res, cdspec['x'][:-1])

def func_lightcurve_ribbon(witch_one='upper'):
    if witch_one == 'upper':
        yrange = [23.5, 40]
    else:
        yrange = [20.5, 23.5]
    # kwlist=['north','south','east']
    # save_file='/Volumes/Data/20170906/plotting/stack_plot/wrap_save/dspec_1600_{}_0_208.p'.format(kwlist[kwindex])
    save_file = '/Volumes/Data/20170906/plotting/stack_plot/wrap_save/wide_north_detailed_0_200.p'
    opfile = open(save_file, 'rb')
    cdspec = pickle.load(opfile)
    opfile.close()
    world_range = np.zeros((2))
    world_range[0] = np.nanargmin([abs(yrange[0] - cpy) for cpy in cdspec['y']])
    world_range[1] = np.nanargmin([abs(yrange[1] - cpy) for cpy in cdspec['y']])
    res = np.zeros((len(cdspec['x']) - 1))
    if witch_one == 'upper':
        for ri in range(len(res)):
            if ri < 89:
                low_b = world_range[0]
            elif ri >= 89 and ri < 159:
                tmp_pix = (36 - 23.8) / 70 * (ri - 89) + 23.5
                low_b = np.nanargmin([abs(tmp_pix - cpy) for cpy in cdspec['y']])
            else:
                low_b = np.nanargmin([abs(36 - cpy) for cpy in cdspec['y']])
            #res[ri] = np.nanmax(cdspec['dspec'][int(low_b):int(world_range[1]), ri])
            res[ri] = np.nanmax(cdspec['dspec'][int(low_b):int(world_range[1]), ri])
    else:
        for ri in range(len(res)):
            if ri < 89:
                up_b = world_range[1]
            elif ri >= 89 and ri < 159:
                tmp_pix = (36 - 23.8) / 70 * (ri - 89) + 23.5
                up_b = np.nanargmin([abs(tmp_pix - cpy) for cpy in cdspec['y']])
            else:
                up_b = np.nanargmin([abs(36 - cpy) for cpy in cdspec['y']])
            res[ri] = np.nanmax(cdspec['dspec'][int(world_range[0]):int(up_b), ri])
    # fig,ax=plt.subplots()
    # ax.plot(np.arange(len(res)),res)
    # ax.plot(np.arange(len(res)),res)
    # plt.show()
    return res


def aia_lightcurve_from_stp(kw, yrange, cur_marker = '-',timerange=None, normq=False, pax=None, ltex='', cur_color='r', sep_north=None):
    if timerange is not None:
        tr_plt = Time(timerange)
    cur_res = lightcurve_over_ribbon(kw=kw, yrange=yrange)
    cur_data = cur_res[0]
    #cur_date = Time(list(cur_res[1] / 3600. / 24.),format='mjd')
    cur_date = Time(list(cur_res[1]),format='plot_date')
    cur_max = np.nanmax(cur_data)
    dates_index = np.arange(0, len(cur_data))
    if timerange is not None:
        tr_plt = Time(timerange)
        dates_index, = np.where((cur_date >= tr_plt[0]) & (cur_date <= tr_plt[1]))
    if normq:
        cur_data = [cx / cur_max for cx in cur_data]
    if pax is None:
        fig, pax = plt.subplots()
    #cim = pax.plot_date(dates, cur_data[timerange[0]:timerange[1]], '-', color=cur_color, label=ltex)
    cim = pax.plot_date(cur_date[dates_index[0]:dates_index[-1]].plot_date, cur_data[dates_index[0]:dates_index[-1]],
            cur_marker, color=cur_color, label=ltex, linewidth=2)
    #cim = pax.plot(cur_date[dates_index[0]:dates_index[-1]].plot_date, cur_data[dates_index[0]:dates_index[-1]], '-', color=cur_color, label=ltex)
    # pax.set_ylim([0,1.1])
    pax.set_xlabel('time')
    formatter = mpl.dates.DateFormatter('%H:%M:%S')
    pax.xaxis.set_major_formatter(formatter)
    pax.fmt_xdata = mpl.dates.DateFormatter('%H:%M:%S')
    #pax.set_xticklabels(rotation=45)
    pax.xaxis.set_tick_params(rotation=30)
    return cim


def aia_movie(ti):
    # plt.ioff()
    fig = plt.figure(figsize=(10, 10))
    kwlist = ['Z.94.', 'Z.131.', 'Z.171.', 'Z.193.', 'Z.211.', 'Z.304.', 'Z.335.', 'bbso', 'Z.1600.']
    tick_list = [[1, 1], [1, 0], [1, 0], [0, 1], [0, 0], [0, 0], [0, 1], [0, 0], [0, 0]]
    # for ti,cdic in enumerate(in_dic):
    cdic = in_dic[ti]
    for ki, ckw in enumerate(kwlist):
        cax = fig.add_subplot(4, 3, ki + 1)
        justmap(cax=cax, kw=ckw, fov=[[568, -255], [618, -205]], tim=ti, enhance_it=True)

        if 'Z' in ckw:
            cur_text = 'SDO/AIA' + ckw.replace('Z', '').replace('.', '')
        else:
            cur_text = 'H-Alpha line core'
        cax.text(.4, .9, cur_text, horizontalalignment='center', transform=cax.transAxes,
                 color='y')
        # if ki==0:
        # plot_eovsa_shift_arrow(cax=cax,start_time=152,end_time=157,fov=[[568, -255], [618, -205]])
        cax.xaxis.tick_top()
        cax.set_xlabel('solar X arcsec')
        cax.set_ylabel('solar Y arcsec')
        cax.xaxis.set_label_position('top')
        tick_helper(cax=cax, tsta=tick_list[ki])
    with open('/Volumes/Data/20170906/tpower/from_big_fov.p', 'rb') as dspec_save:
        aspec = pickle.load(dspec_save)
    dspec_save.close()
    dsax = fig.add_subplot(4, 1, 4)
    mdds.plt_dspec(specdata=aspec, timerange=[0, 297], pax=dsax,
                   time_line=[Time(cdic['time'], format='iso')],
                   rhessisave='/Volumes/Data/20170906/rhessi/lc/corrected_hessi.sav')
    # cur_time = Time(cdic['time'], format='iso')
    # cax.axvline(cur_time.plot_date, linewidth=2, color='b')
    plt.subplots_adjust(wspace=0, hspace=0, left=0.21, right=0.79)
    curname = '/Volumes/Data/20170906/plotting/justmap/all_channels/T{0:0=3d}_all_cha.png'.format(ti)
    fig.suptitle(Time(cdic['time'], format='iso').isot)
    # plt.savefig(curname)
    # plt.clf()
    plt.show()
    # return


def mp_aia_movie(fff):
    # pool=mlp.Pool()
    larr = np.linspace(0, 240, 9, dtype=int)
    for iii in range(larr[fff], larr[fff + 1]):
        print('plotting the {}th time'.format(iii))
        aia_movie(iii)
        # p = mlp.Process(target=aia_movie, args=(iii,))
        # p.start()

        # pool.apply_async(aia_movie,(iii,))
    # pool.close()
    # pool.join()
    # pool.map(worker,range(296))
    # pool.close()
    # pool.join()


# def worker(args):

def shift_maps(tim, kw):
    curmap = smap.Map(in_dic[tim][kw])
    print(curmap.meta)
    new_meta = curmap.meta
    # for T150 19:28:20
    new_meta['crval1'] = 1.0
    new_meta['crval2'] = -2.0
    new_data = np.log(curmap.data)
    curname = '/Volumes/Data/20170906/gx_simulator/shifted_T{0:0=3d}_1600.fits'.format(tim)
    write_file(fname=curname, data=new_data, header=new_meta, filetype='fits')


def plot_gx_as_png(tim):
    # fov: 654,-171   571,-260
    png_fov = [[558, -266], [654, -171]]
    fig = plt.figure(figsize=(6, 6.5))
    # img1=mpimg.imread('/Volumes/Data/20170906/gx_simulator/img_gx.png')
    # img2=mpimg.imread('/Volumes/Data/20170906/gx_simulator/img_eovsa.png')
    cax1 = fig.add_subplot(2, 2, 1)
    cax2 = fig.add_subplot(2, 2, 2)
    cax3 = fig.add_subplot(2, 2, 3)
    cax4 = fig.add_subplot(2, 2, 4)
    img3 = Image.open('/Volumes/Data/20170906/gx_simulator/img_gx.png')
    img1 = Image.open('/Volumes/Data/20170906/gx_simulator/notick/cartoon_20170906_before_nt.png')
    img2 = Image.open('/Volumes/Data/20170906/gx_simulator/notick/cartoon_20170906_after_nt.png')
    justmap(cax=cax4, kw='Z.94.', fov=png_fov, tim=tim, enhance_it=True)
    cax4.xaxis.set_major_locator(plt.MaxNLocator(3))
    cax4.yaxis.set_major_locator(plt.MaxNLocator(4))
    cxt = cax4.get_xticks()
    cyt = cax4.get_yticks()
    '''
    img2=Image.open('/Volumes/Data/20170906/gx_simulator/img_eovsa.png')
    img2.thumbnail((np.asarray(img1).shape[0], np.asarray(img1).shape[1]), Image.ANTIALIAS)
    arr_img2 = np.asarray(img2)
    white = np.sum(arr_img2[:, :, :3], axis=2)
    white_mask = np.where(white == 255 * 3, 1, 0)
    alpha = np.where(white_mask, 0, arr_img2[:, :, -1])
    tmp_arr = copy.deepcopy(arr_img2)
    tmp_arr[:, :, -1] = alpha
    trans_img = Image.fromarray(np.uint8(tmp_arr))
    '''
    # justmap(cax=cax2,kw='Z.1600.',fov=png_fov,tim=tim)
    # text_img = Image.new('RGBA', (img1.shape[0], img1.shape[1]), (0, 0, 0, 0))
    # text_img.paste(img1, (0, 0))
    # text_img.paste(img2, (0, 0), mask=img2)
    # cax1.imshow(img1)
    plot_eovsa_contourf(tim=tim, ax=cax4, fov=png_fov, abs_contourf=False, spws=range(9, 30),
                        level=0.95, alt_cmap=plt.cm.Spectral)
    cax1.imshow(img1)
    cax2.imshow(img2)
    cax3.imshow(img3)
    # --tick control-------------
    cxt = cax4.get_xticks()
    cyt = cax4.get_yticks()
    cxt = cxt[1:-1]
    cyt = cyt[1:-1]
    # cax1.set_xticks(cxt,cxtnl)
    # cax1.set_yticks(cyt,cytnl)
    for ccax in [cax1, cax2, cax3]:
        cur_xticks = np.zeros_like(cxt)
        cur_yticks = np.zeros_like(cyt)
        for cii, cticks in enumerate(cxt):
            cur_xticks[cii] = (cticks - cax4.get_xlim()[0]) / (cax4.get_xlim()[1] - cax4.get_xlim()[0]) * (
                        ccax.get_xlim()[1] - ccax.get_xlim()[0]) + ccax.get_xlim()[0]
        for cii, cticks in enumerate(cyt):
            cur_yticks[cii] = (cticks - cax4.get_ylim()[0]) / (cax4.get_ylim()[1] - cax4.get_ylim()[0]) * (
                        ccax.get_ylim()[1] - ccax.get_ylim()[0]) + ccax.get_ylim()[0]
        ccax.set_xticks(cur_xticks)
        ccax.set_yticks(cur_yticks)
        ccax.set_xticklabels([int(ttt) for ttt in cxt])
        ccax.set_yticklabels([int(ttt) for ttt in cyt])
    # end of tick control---------------------------

    # cax1.imshow(trans_img,alpha=.89, interpolation='bilinear')
    # cax3.set_title('19:00:00')
    cax1.set_title('pre-reconnection')
    cax3.text(.5, .9, '19:00:00', horizontalalignment='center', transform=cax3.transAxes, color='y')
    cax4.text(.5, .9, '19:28:20', horizontalalignment='center', transform=cax4.transAxes, color='y')
    # cax4.set_title('19:28:20')
    cax2.set_title('post-reconnection')
    tick_helper(cax1, [0, 1])
    tick_helper(cax2, [0, 0])
    tick_helper(cax4, [1, 0])
    plt.subplots_adjust(left=0.12, bottom=0.16, right=0.90, top=0.88, wspace=0.02, hspace=0)
    plt.show()
    # return cax4.get_xlim()


def plot_eovsa_shift_arrow(cax, start_time, end_time, spw_range=[8, 30], fov=None):
    with open('/Volumes/Data/20170906/info/ori_pimfit_out.p', 'rb') as pimfit_peak:
        eovsa_peak_list = copy.deepcopy(pickle.load(pimfit_peak))
    pimfit_peak.close()
    norm = mpl.colors.Normalize(vmin=0, vmax=30)
    cmap = plt.cm.viridis
    if fov is None:
        ref_map = smap.Map(in_dic[start_time]['radio_sbs'][15])
    else:
        ref_map = pot.make_sub_map(cur_map=smap.Map(in_dic[start_time]['radio_sbs'][15]), fov=fov)
    for spwi in range(spw_range[0], spw_range[1]):
        print(spwi)
        cur_pix_s = eovsa_peak_list[start_time - 100][spwi]['outputs'][0]['XX']['results']['component0']['pixelcoords']
        cur_pix_e = eovsa_peak_list[end_time - 100][spwi]['outputs'][0]['XX']['results']['component0']['pixelcoords']
        cur_world_s = pot.get_world(cur_pix_s[0], cur_pix_s[1], cur_map=smap.Map(in_dic[start_time]['radio_sbs'][15]))
        cur_world_e = pot.get_world(cur_pix_e[0], cur_pix_e[1], cur_map=smap.Map(in_dic[end_time]['radio_sbs'][15]))
        new_pix_s = pot.get_pixel(cur_map=ref_map, x=cur_world_s[0], y=cur_pix_s[1])
        new_pix_e = pot.get_pixel(cur_map=ref_map, x=cur_world_e[0], y=cur_pix_e[1])
        kwargs = {"color": cmap(norm(spwi))}
        print(cur_world_e[0] - cur_world_s[0])
        print(cur_world_e[1] - cur_world_s[1])
        # cax.arrow(x=new_pix_s[0],y=new_pix_s[1],dx=new_pix_e[0]-new_pix_s[0],dy=new_pix_e[1]-new_pix_s[1],**kwargs)
        cax.arrow(x=cur_world_s[0], y=cur_world_s[1], dx=cur_world_e[0] - cur_world_s[0],
                  dy=cur_world_e[1] - cur_world_s[1], **kwargs)
        # cax.plot(cur_world_s[0], cur_world_e[1], 's', markersize='6.', color='k', markerfacecolor='k')


def plot_arrow_ax(cax, y, dy, cx=None, cdx=None, pdate=True, tim=None, acolor=None, text=None, cwidth=0.05):
    if pdate:
        cx = Time(in_dic[tim]['time'], format='iso').plot_date
        cdx = 0
        anno_gap = [0, -(dy) / 2]
    else:
        anno_gap = [-(cdx) / 2, -(dy) / 2]
    kwargs = {"color": acolor, "width": cwidth, "length_includes_head": True, "head_length": abs(dy) * 0.3}
    kwargst = {"color": acolor}
    cim = cax.arrow(x=cx, y=y, dx=cdx, dy=dy, **kwargs)
    cax.text(cx + anno_gap[0], y + anno_gap[1], text, fontsize=10, **kwargst)
    return cim

def axes_helper(cax,ctitle=None,ti_di=['out','out'],right_y=False,top_x=False,no_xtick=False,no_ytick=False,xlabel=None,ylabel=None,cur_in_dic=None,
                index_letter=None,ind_let_loc=None,ind_let_color='r',ind_let_size=10,legend_loc=None,num_of_xticks=None,num_of_yticks=None,xspan=None,add_text=None,ind_let_bkg=None,high_tres = False):
    font = {'weight': 'bold',
            'size': 12}
    mpl.rc('font', **font)
    in_dic = cur_in_dic
    if ctitle is not None:
        cax.set_title(ctitle)
    cax.tick_params(axis='x', direction=ti_di[0])
    cax.tick_params(axis='y', direction=ti_di[1])
    cax.tick_params(axis='both', bottom=True, top=True, left=True, right=True, width=1.6,labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    if right_y:
        #cax.yaxis.tick_right()
        cax.yaxis.set_label_position("right")
        cax.tick_params(axis='both', labelleft=False, labelright=True)
    if top_x:
        cax.xaxis.set_label_position("top")
        cax.tick_params(axis='both', labelbottom=False, labeltop=True)
        #cax.xaxis.tick_top()
    if no_xtick:
        #cax.set_xticks([])
        cax.set_xlabel('')
        #cax.xaxis.set_ticklabels([])
        cax.tick_params(axis='both', labelbottom=False, labeltop=False)
        #cax.xaxis.set_visible(False)
    if no_ytick:
        #cax.set_yticks([])
        cax.set_ylabel('')
        cax.tick_params(axis='both', labelleft=False, labelright=False)
        #cax.yaxis.set_ticklabels([])
        #cax.yaxis.set_visible(False)
    if xlabel is not None:
        cax.set_xlabel(xlabel)
    if ylabel is not None:
        cax.set_ylabel(ylabel)
    if index_letter is not None:
        #font_size_list=['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large']
        #textkw={'size':font_size_list[ind_let_size]}
        if ind_let_loc is None:
            ind_let_loc=[0.15,0.9]
        if ind_let_bkg is not None:
            props = dict(boxstyle='round', facecolor=ind_let_bkg, alpha=0.8)
            cax.text(ind_let_loc[0], ind_let_loc[1], '({0})'.format(index_letter), horizontalalignment='center', transform=cax.transAxes, color=ind_let_color,fontsize=ind_let_size,bbox=props)
        else:
            cax.text(ind_let_loc[0], ind_let_loc[1], '({0})'.format(index_letter), horizontalalignment='center',
             transform=cax.transAxes, color=ind_let_color, fontsize=ind_let_size)
    if legend_loc is not None:
        cax.legend(loc=legend_loc)
    if num_of_xticks is not None:
        cax.xaxis.set_major_locator(plt.MaxNLocator(num_of_xticks))
    if num_of_yticks is not None:
        cax.yaxis.set_major_locator(plt.MaxNLocator(num_of_yticks))
    if xspan is not None:
        if xspan [2] == 1:
            sxv = Time(in_dic[xspan[0][0]]['time'], format='iso').plot_date
            exv = Time(in_dic[xspan[0][1]]['time'], format='iso').plot_date
        else:
            sxv=xspan[0][0]
            exv = xspan[0][1]
        cax.axvspan(sxv,exv,alpha=0.4,color=xspan[1])
    if add_text is not None:
        cax.text(add_text[1], add_text[2], add_text[0], transform=cax.transAxes,
                 color=add_text[3])
    font = {'weight': 'bold',
            'size': 12}
    mpl.rc('font', **font)
    mpl.rcParams.update({'font.size': 12})

def axes_helper_bkp(cax, ctitle=None, ti_di=['in', 'in'], right_y=False, top_x=False, no_xtick=False, no_ytick=False,
                xlabel=None, ylabel=None,
                index_letter=None, ind_let_loc=None, ind_let_color='r', ind_let_size=10, legend_loc=None,
                num_of_xticks=None, num_of_yticks=None, xspan=None, add_text=None, data_1s=True, ind_let_bkg=None):
    if data_1s:
        in_dic = in_dic_1s
    else:
        cur_dic = '/Volumes/WD6T/working/20170703/info/20170703_1s.p'
        in_dic_file = open(cur_dic, 'rb')
        in_dic = pickle.load(in_dic_file)
        in_dic_file.close()
    if ctitle is not None:
        cax.set_title(ctitle)
    cax.tick_params(axis='x', direction=ti_di[0])
    cax.tick_params(axis='y', direction=ti_di[1])
    if right_y:
        cax.yaxis.tick_right()
        cax.yaxis.set_label_position("right")
    if top_x:
        cax.xaxis.set_label_position("top")
        cax.xaxis.tick_top()
    if no_xtick:
        cax.set_xticks([])
        cax.set_xlabel('')
        cax.xaxis.set_visible(False)
    if no_ytick:
        cax.set_yticks([])
        cax.set_ylabel('')
        cax.yaxis.set_visible(False)
    if xlabel is not None:
        cax.set_xlabel(xlabel)
    if ylabel is not None:
        cax.set_ylabel(ylabel)
    if index_letter is not None:
        # font_size_list=['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large']
        # textkw={'size':font_size_list[ind_let_size]}
        if ind_let_loc is None:
            ind_let_loc = [0.15, 0.9]
            if ind_let_bkg is not None:
                props = dict(boxstyle='round', facecolor=ind_let_bkg, alpha=0.8)
                cax.text(ind_let_loc[0], ind_let_loc[1], '({0})'.format(index_letter), horizontalalignment='center',
                         transform=cax.transAxes, color=ind_let_color, fontsize=ind_let_size, bbox=props)
            else:
                cax.text(ind_let_loc[0], ind_let_loc[1], '({0})'.format(index_letter), horizontalalignment='center',
                         transform=cax.transAxes, color=ind_let_color, fontsize=ind_let_size)
    if legend_loc is not None:
        cax.legend(loc=legend_loc)
    if num_of_xticks is not None:
        cax.xaxis.set_major_locator(plt.MaxNLocator(num_of_xticks))
    if num_of_yticks is not None:
        cax.yaxis.set_major_locator(plt.MaxNLocator(num_of_yticks))
    if xspan is not None:
        if xspan[2] == 1:
            sxv = Time(in_dic[xspan[0][0]]['time'], format='iso').plot_date
            exv = Time(in_dic[xspan[0][1]]['time'], format='iso').plot_date
        else:
            sxv = xspan[0][0]
            exv = xspan[0][1]
        cax.axvspan(sxv, exv, alpha=0.4, color=xspan[1])
    if add_text is not None:
        cax.text(add_text[1], add_text[2], add_text[0], transform=cax.transAxes,
                 color=add_text[3])


def derivative_lightcurve(ori_lightcurve, cur_color='r', pax=None, ltex=''):
    der_res = np.zeros_like(ori_lightcurve[1])
    for cl in range(len(ori_lightcurve[1])):
        if cl == 0:
            der_res[cl] = 0
            continue
        else:
            der_res[cl] = ori_lightcurve[1][cl] - ori_lightcurve[1][cl - 1]
    cim = pax.plot_date(ori_lightcurve[0], der_res, '-', color=cur_color, label=ltex)
    return cim


def get_all_points_to_be_fitted(ctim):
    cur_eovsa_map = smap.Map(in_dic[804]['radio_sbs'][9])
    cur_data = sf.mapdata(cur_eovsa_map)
    fig, cax = plt.subplots()
    # plt.imshow(cur_data>0.5*eovsa_max_list[9],origin='lower')
    # res=np.zeros_like(cur_data)
    # res=sf.get_data_block(x_size=15,y_size=21,initx=123,inity=119,time_arr=[124])
    # res=sf.get_data_block(x_size=13,y_size=19,initx=123,inity=120,time_arr=[124])
    # print(res[:,:,9,0].shape)
    # plt.imshow(cur_data[119:140,123:138])
    # plt.imshow(res[:,:,9,0],origin='lower')

    # plt.imshow(cur_data==10735983.0,origin='lower')
    # plot_eovsa_contourf(tim=124,spws=[9],abs_contourf=True,abs_level=eovsa_max_list,level=[0.45],ax=cax, single_spw=True)
    # plot_map.imshow(sunpymap=cur_eovsa_map,axes=cax)
    # plt.imshow(cur_data>0.5*eovsa_max_list[9])
    plt.show()
    pix_list = []
    res_list = []
    data_list = []
    world_list = []
    cur_fov = [[570, -255], [610 - 210]]
    data_arr = np.zeros((256, 256, 30))
    for xx in range(570, 610):
        for yy in range(-270, -210):
            print(xx, yy)
            c_pix = pot.get_pixel(cur_map=cur_eovsa_map, x=xx, y=yy)
            if c_pix not in pix_list:
                cdata = sf.get_data_fits(x=xx, y=yy, tim=ctim, world=1, phase=13)
                data_arr[c_pix[0], c_pix[1], :] = cdata
                data_list.append(cdata)
                world_list.append([xx, yy])
                pix_list.append(c_pix)
    print(pix_list)
    final_data = (data_arr, data_list, pix_list, world_list)
    pickle.dump(final_data, open('/Volumes/Data/20170906/fitting_result/map/T{0:0=3d}_data.p'.format(ctim), 'wb'))
    # return (data_arr,data_list,pix_list,world_list)


def save_enhanceed_halpha_to_fits(ct):
    sub_map = smap.Map('/Volumes/Data/20170906/bbso/corrected/T{0:0=3d}.fits'.format(ct))
    ori_name = '/Volumes/Data/20170906/bbso/corrected/T{0:0=3d}.fits'.format(ct)
    cur_name = ori_name.replace('/corrected/', '/corrected_enhanced/')
    if not os.path.exists(cur_name):
        sub_map = smap.Map(enhance.mgn(data=sub_map.data, sigma=[21.25, 42.5, 85, 170, 340, 680]), sub_map.meta)
        print(cur_name)
        write_file(cur_name, sub_map.data, header=sub_map.meta, filetype='fits')


def multi_proc(ncpu):
    import decay_phase
    pool = mlp.Pool(ncpu)
    # res=pool.map(decay_phase.plot_bt_z,range(240))
    res = pool.map(get_all_points_to_be_fitted, range(100, 102))
    # res=pool.map(save_enhanceed_halpha_to_fits,range(296))
    # res=pool.map(eovsa_integ_lc,range(30))
    # res=pool.map(dprp.save_tmap_2_fits,range(210))
    # res=pool.map(single_plot,range(200))

    pool.close()
    pool.join()


def eovsa_integ_lc(cspw):
    cname = '/Volumes/Data/20170906/plotting/eovsa_integrated_lc/spw_{0:0=2d}_all_positive.p'.format(cspw)
    if os.path.exists(cname):
        cflie = open(cname, 'rb')
        res = pickle.load(cflie)
        cflie.close()
    else:
        print('make spw', cspw)
        res = np.zeros((len(in_dic)))
        for ti, cdic in enumerate(in_dic):
            semap = smap.Map(cdic['radio_sbs'][cspw])
            tmp_array = semap.data
            tmp_array[tmp_array < 0] = 0
            res[ti] = np.nansum(tmp_array)
        cfile = open(cname, 'wb')
        pickle.dump(res, cfile)
        cfile.close()
    return res


def customized_stp_on_axes(cax, dspec, min_max=None,min_max_factor=None, plot_log=False, kw=None, yrange=None, custom_cmap=None):
    if type(dspec) == str:
        with open(dspec, 'rb') as cfile:
            dspec = pickle.load(cfile)
    if min_max is not None:
        cvmin = min_max[0]
        cvmax = min_max[1]
    elif min_max_factor is not None:
        cvmin = np.nanmin(dspec['dspec'])*min_max_factor[0]
        cvmax = np.nanmax(dspec['dspec'])*min_max_factor[1]
    else:
        cvmin = np.nanmin(dspec['dspec'])
        cvmax = np.nanmax(dspec['dspec'])
    if plot_log:
        norm = colors.LogNorm(vmin=cvmin, vmax=cvmax)
    else:
        norm = colors.Normalize(vmin = cvmin, vmax=cvmax)
    kwargs = {'norm': norm}
    if custom_cmap is not None:
        kwargs['cmap'] = custom_cmap
    if kw is not None:
        kwargs['cmap'] = plt.get_cmap(name='sdoaia' + kw.replace('Z', '').replace('.', ''))
    # if uni_cm:
    #    dspec['args']['norm'] = norm
    print(dspec['y'].shape, dspec['dspec'].shape)
    if yrange is not None:
        im = cax.pcolormesh(dspec['x'], dspec['y'][yrange[0]:yrange[1]], dspec['dspec'][yrange[0]:yrange[1], :],
                            **kwargs)
    else:
        im = cax.pcolormesh(dspec['x'], dspec['y'], dspec['dspec'], **kwargs)
    date_format = mpl.dates.DateFormatter('%H:%M:%S')
    cax.xaxis_date()
    cax.xaxis.set_major_formatter(date_format)
    cax.set_xlabel('Time[UT]')
    cax.set_ylabel('Distance[arcsec]')

    return im


def parameters_on_axes(cax, fitsavefile, dspec, uni_cm=True, kw=None, yrange=None, pfig=None):
    if type(dspec) == str:
        with open(dspec, 'rb') as cfile:
            dspec = pickle.load(cfile)
    if type(fitsavefile) == str:
        with open(fitsavefile, 'rb') as fsfile:
            res_list = pickle.load(fsfile)
    title_list = ['theta', 'temperature', 'delta', 'thermal_n', 'non-thermal_n', 'B']
    params_list = ['theta', 'temperature', 'delta', 'n_th', 'n_nth', 'mag']
    for parmi in range(6):
        c_array = np.zeros((len(res_list)))
        e_array = np.zeros((len(res_list)))
        if parmi == 5:
            ext_list = []
        for pointi, res_point in enumerate(res_list):
            cefile = glob.glob(
                '/Volumes/Data/20170906/fitting_result/emcee/time_sequense/600_230_vari_bpar/world_emcee_res_at_*_T{0:0=3d}.p'.format(
                    pointi + 110))[0]
            with open(cefile, 'rb') as emcee_res_file:
                emcee_res = pickle.load(emcee_res_file)
            if parmi == 5:
                if pointi < 8:
                    ext_list.append((abs(pointi - 7) ** 1.7 / 28 + 1) * 750)
                else:
                    ext_list.append((abs(pointi - 7) ** 2 / 25 + 1) * 750)
            print(res_point)
            c_array[pointi] = res_point.params[params_list[parmi]].value
            # e_array[pointi] = res_point.params[params_list[parmi]].stderr
            e_array[pointi] = emcee_res.params[params_list[parmi]].stderr
    print(e_array)
    axe = cax.twinx()
    print(c_array.shape)
    im = axe.errorbar(dspec['x'][10:59], c_array, yerr=e_array, linestyle='None', marker='x', color='r')
    axe.set_ylim(300, 800)
    axe.xaxis_date()
    date_format = mpl.dates.DateFormatter('%H:%M:%S')
    axe.xaxis.set_major_formatter(date_format)

    return im


def rhessi_eovsa_lc(cax, timerange, temp_curve=True):
    tidx = range(timerange[0], timerange[1])
    # fig=plt.figure()
    # cax=fig.add_subplot(111)
    rhe_yy = np.zeros((len(tidx)))
    rhessisave = '/Volumes/Data/20170906/rhessi/lc/corrected_hessi.sav'
    rhe_sav = readsav(rhessisave)
    rhessi_mjd_list = Time([str(ll) for ll in rhe_sav['obs_times_str']], format='isot').mjd
    tim_plt = []
    for ti in range(len(tidx)):
        ctim = Time(in_dic[tidx[ti]]['time'], format='iso')
        rhe_yy[ti] = rhe_sav['obs_data'][:, 1][np.argmin(a=[abs(x - ctim.mjd) for x in rhessi_mjd_list])] / np.nanmax(
            rhe_sav['obs_data'][:, 1]) * 20
        tim_plt.append(ctim.plot_date)
    cax.plot_date(tim_plt, rhe_yy, linestyle='-', marker='', color='k')
    if temp_curve:
        with open('/Volumes/Data/20170906/dem/demreg_py/n_pix_hotter_than_20MK.p', 'rb') as ct_file:
            temp_cv = pickle.load(ct_file)[timerange[0]:timerange[1]]
        cax.plot_date(tim_plt, [ccc / np.nanmax(temp_cv) * 20 for ccc in temp_cv], linestyle='-', marker='', color='b')
    cmap = plt.cm.gist_rainbow
    norm = mpl.colors.Normalize(vmin=10, vmax=25)
    # for spwi in range(10,25):
    for spwi in range(16, 17):
        continue
        ceovsa_lc_file = open(
            '/Volumes/Data/20170906/plotting/eovsa_integrated_lc/spw_{0:0=2d}_all_positive.p'.format(spwi), 'rb')
        ceovsa_lc = pickle.load(ceovsa_lc_file)
        ceovsa_lc_file.close()
        ceovsa_lc = [cd / np.nanmax(ceovsa_lc) * 20 for cd in ceovsa_lc]
        cax.plot_date(tim_plt, ceovsa_lc[99:182], linestyle=':', color=cmap(norm(spwi)), marker='')
    return cax


# def multi_traj_on_axes():
def check_aia_exptime():
    y_res = np.zeros((len(in_dic)))
    for ti in range(len(in_dic)):
        cur_map = smap.Map(in_dic[ti]['Z.193.'])
        y_res[ti] = cur_map.meta['exptime']
    fig, ax = plt.subplots()
    ax.plot(np.arange(len(in_dic)), y_res, linestyle='', marker='.')
    plt.show()
    return y_res

def number_pixel_hotter_than(cax=None, thre=None, exp_time_save=None, xrange=None):
    if cax is None:
        fig, cax = plt.subplots()
    if exp_time_save is None:
        exp_time_save = check_aia_exptime()
    temp_curve = np.zeros((211))
    init_val = 0
    for ti, cdic in enumerate(in_dic):
        if ti > 210: continue
        cur_tmap = dprp.dem_sav_to_map(sav_file=ti)
        if exp_time_save[ti] > 1.9 or ti in [141, 153, 158, 159, 160, 165, 177]:
            temp_curve[ti] = init_val
        else:
            temp_curve[ti] = (cur_tmap.data > thre).sum()
            init_val = temp_curve[ti]
    cax.plot(np.arange(211), temp_curve)
    plt.show()
    return temp_curve


def single_plot(tim):
    plt.ioff()
    fig = plt.figure(figsize=(5, 5))
    cax = fig.add_subplot(111)
    fov2 = [[560, -263], [627, -191]]
    justmap(cax, tim, 'Z.1600.', fov=fov2)
    plot_eovsa_contourf(tim=tim, ax=cax, fov=fov2, spws=range(9, 30),
                        # abs_contourf=False, level=0.9,
                        abs_contourf=True, abs_level=eovsa_max_list * 0.5,
                        alt_cmap=plt.cm.Spectral, cs_start=[9, 30])
    cax.set_title(in_dic[tim]['time'])
    cname = '/Volumes/Data/20170906/plotting/justmap/halpha_eovsa_abs/T{0:0=3d}.png'.format(tim)
    plt.savefig(cname)
    # plt.show()
    # plt.close('all')
    fig.clf()


def add_patch(cax, cfov, ccc='r'):
    import matplotlib.patches as patches
    cax.add_patch(
        patches.Rectangle(
            xy=(cfov[0][0], cfov[0][1]),  # point of origin.
            width=cfov[1][0] - cfov[0][0],
            height=cfov[1][1] - cfov[0][1],
            linewidth=1,
            color=ccc,
            fill=False
        )
    )


def daterize_colorbar(trange, nbins):
    from datetime import datetime
    from matplotlib.dates import (DateFormatter, rrulewrapper, RRuleLocator, drange)
    if not isinstance(trange[0], str):
        start_t = mdates.date2num(datetime.strptime(str(in_dic[trange[0]]['time']).split('.')[0], "%Y-%m-%d %H:%M:%S"))
        end_t = mdates.date2num(datetime.strptime(str(in_dic[trange[1]]['time']).split('.')[0], "%Y-%m-%d %H:%M:%S"))
    else:
        start_t = mdates.date2num(datetime.strptime(trange[0], "%Y-%m-%d %H:%M:%S"))
        end_t = mdates.date2num(datetime.strptime(trange[1], "%Y-%m-%d %H:%M:%S"))
    seq = np.linspace(start_t, end_t, nbins)
    label_list = []
    for i in range(nbins):
        ct = mdates.num2date(seq[i])
        label_list.append(ct.strftime('%H:%M:%S'))
    return (seq, label_list)

    loc = mdates.AutoDateLocator()
    return loc


def plot_fermi_dynamic_spec_in_ax(cax, spec_dict, timerange, dmin_max=None, erg_max=None, nde=2):
    fermi_tim = Time(spec_dict['tim'], format='mjd')
    fermi_tim_plt = fermi_tim.plot_date
    tr_plt = Time(timerange)
    tidx_spec, = np.where((fermi_tim >= tr_plt[0]) & (fermi_tim <= tr_plt[1]))
    if not dmin_max:
        dmin_max = [np.nanmin(spec_dict['spec'][:, :]),
                    np.nanmax(spec_dict['spec'][:, :])]
    (ntim, nerg) = spec_dict['spec'].shape
    tim = spec_dict['tim']
    if erg_max is not None:
        emax_idx = min(range(len(spec_dict['erg'])), key=lambda i: abs(spec_dict['erg'][i] - erg_max))
    erg = spec_dict['erg'][0:emax_idx + 1]
    spec_tim = Time(spec_dict['tim'], format='mjd')
    timstrr = spec_tim[tidx_spec[0]:tidx_spec[-1]].plot_date
    plt.ion()
    # fig = plt.figure(figsize=(12, 7), dpi=100)
    # gs1 = gridspec.GridSpec(3, 1)
    # gs1.update(left=0.08, right=0.32, wspace=0.05)
    # gs2 = gridspec.GridSpec(2, 2)
    # gs2.update(left=0.38, right=0.98, hspace=0.02, wspace=0.02)
    # if npol > 1:
    #    spec_1 = np.absolute(spec[0, 0, :, :])
    #    spec_2 = np.absolute(spec[1, 0, :, :])
    # else:
    #    spec_1 = np.absolute(spec[0, 0, :, :])
    #    spec_2 = np.zeros_like(spec_1)
    # polstr = ['XX', 'YY']
    # plot begin
    # ax1 = plt.subplot(gs1[0])
    cax.pcolormesh(timstrr, erg, np.swapaxes(spec_dict['spec'][tidx_spec[0]:tidx_spec[-1], 0:emax_idx + 1], 0, 1),
                   cmap='jet', vmin=dmin_max[0],
                   vmax=dmin_max[1])
    # cax.set_xlim(timstrr[tidx[0]], timstrr[tidx[-1]])
    cax.xaxis_date()
    cax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
    # ax1.set_xticklabels(['']*10)
    # cax.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])
    cax.set_ylabel('KeV Det_No.{}'.format(nde), fontsize=10)
    # cax.set_title(observatory + ' ' + datstr + ' ' + polstr[0] + ' & ' + polstr[1], fontsize=12)
    cax.set_autoscale_on(False)
    plt.xticks(rotation=15)
    # cax.add_patch(patches.Rectangle((bt, bfreqghz), et - bt, efreqghz - bfreqghz, ec='w', fill=False))
    # cax.plot([(bt + et) / 2.], [(bfreqghz + efreqghz) / 2.], '*w', ms=12)
    for tick in cax.get_xticklabels():
        tick.set_fontsize(8)
    for tick in cax.get_yticklabels():
        tick.set_fontsize(8)


def aia_files_before_after_timepoint(tim, ckw):
    aia_time = tp.read_time(fits_file=in_dic_1s[tim][ckw], instrument='aia').jd
    ctime = Time(in_dic_1s[tim]['time'], format='iso').jd
    if aia_time < ctime:
        t_plus = 1
        while in_dic_1s[tim][ckw] == in_dic_1s[tim + t_plus][ckw]:
            t_plus += 1
        return [tim, tim + t_plus]
    else:
        t_plus = 1
        while in_dic_1s[tim][ckw] == in_dic_1s[tim - t_plus][ckw]:
            t_plus += 1
        return [tim - t_plus, tim]


def convert_fermi_spectrum_to_lc_rebin_energy(fermi_save, erange):
    with open(fermi_save, 'rb') as fdspec_file:
        fdspec = pickle.load(fdspec_file)
    fdspec_file.close()
    sindex = (np.abs(erange[0] - fdspec['erg'])).argmin()
    eindex = (np.abs(erange[1] - fdspec['erg'])).argmin()
    res_lc = np.nansum(fdspec['spec'][:, sindex:eindex], axis=1)
    return (res_lc, fdspec['tim'])


def plot_parameters_map_on_axes(cax, tim, cparam, fov=None, alp=None):
    with open('/Volumes/WD6T/working/20170703/fitting_result/fit/dictionary/fitting_dict_t{0:0=4d}.p'.format(tim),
              'rb') as fpmap:
        params_dic = pickle.load(fpmap)
    fpmap.close()
    ylim_list = {'theta': [15, 75], 'temperature': [2.e6, 1.1e7], 'n_th': [7.4, 10.4], 'mag': [300, 1200],
                 'chisqr': [3572643.0, 357264350.0]}
    if cparam == 'chisqr':
        param_map = params_dic[cparam][1]
    else:
        param_map = params_dic[cparam][2]
    if fov is not None:
        param_map = pot.make_sub_map(cur_map=param_map, fov=fov)
    kwargs = {'cmap': 'viridis', 'vmin': ylim_list[cparam][0], 'vmax': ylim_list[cparam][1]}
    if alp is not None:
        kwargs.update({'alpha': alp})
    param_map_obj = plot_mapX.Sunmap(sunmap=param_map)
    param_map_obj.imshow(axes=cax, **kwargs)


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


def dem_map_old(tim, get_em=False,high_tres=False):
    if high_tres:
        in_dic = in_dic_1s
    else:
        in_dic = in_dic_10s
    kw_list = ['dem', 'em']
    cur_fits_file = '/Volumes/Data/20170820/dem/dem_fits/t{0:0=4d}_{1}.fits'.format(tim, kw_list[get_em])
    if not os.path.exists(cur_fits_file):
        pix_save = open('/Volumes/Data/20170820/dem/res_save/multi_t{0:0=4d}.p'.format(tim), 'rb')
        dem_res_list = pickle.load(pix_save)
        pix_save.close()
        rmap = smap.Map(in_dic[231]['radio_sbs'][0])
        pi = 0
        new_data = np.zeros_like(rmap.data)
        for xi in range(rmap.data.shape[0]):
            for yi in range(rmap.data.shape[1]):
                new_data[yi, xi] = dem_calculator(odem_res=dem_res_list[pi], get_em=get_em)
                pi += 1
        new_map = smap.Map(new_data, rmap.meta)
        new_map.save(cur_fits_file)
    else:
        new_map = smap.Map(cur_fits_file)
    return new_map

def dem_map(tim, get_em=False, high_tres=False):
    if high_tres:
        in_dic = in_dic_1s
    else:
        in_dic = in_dic_10s
    kw_list = ['dem', 'em']
    cur_fits_file = '/Volumes/Data/20170820/dem/dem_fits/t{0:0=4d}_{1}.fits'.format(tim, kw_list[get_em])
    if not os.path.exists(cur_fits_file):
        #pix_save = open('/Volumes/Data/20170820/dem/res_save/multi_t{0:0=4d}.p'.format(tim), 'rb')
        pix_save = open('/Volumes/Data/20170820/dem/res_save/multi_t{0:0=4d}.p'.format(tim), 'rb')
        dem_res_list = pickle.load(pix_save)
        pix_save.close()
        #camap = smap.Map(dem_res_list[1]['f_list'][0])
        camap = smap.Map(in_dic_10s[tim]['Z.94.'])
        #camap = update_pointing(camap)
        #camap = register(camap)
        rmap = pot.make_sub_map(camap,fov=dem_res_list[1]['fov'])
        new_data = np.zeros_like(rmap.data)
        print('shape', new_data.shape)
        print('number of res', len(dem_res_list[0]))
        for li, cres in enumerate(dem_res_list[0]):
            print(cres['y'],cres['x'])
            print(dem_calculator(cres,get_em=get_em))
            new_data[int(cres['y']),int(cres['x'])] = dem_calculator(cres,get_em=get_em)
        new_map = smap.Map(new_data, rmap.meta)
        new_map.save(cur_fits_file)
    else:
        new_map = smap.Map(cur_fits_file)
    return new_map

def dem_result(wx,wy,tim,plot_it=True):
    font = {'weight': 'bold',
            'size': 15}
    mpl.rc('font', **font)
    dem_file = '/Volumes/Data/20170820/dem/res_save/multi_t{0:0=4d}.p'.format(tim)
    cur_fits_file = '/Volumes/Data/20170820/dem/dem_fits/t{0:0=4d}_{1}.fits'.format(tim, 'dem')
    dopen = open(dem_file, 'rb')
    dem_sav = pickle.load(dopen)
    dopen.close()
    cur_tmap = smap.Map(cur_fits_file)
    pix_coord = pot.get_pixel(cur_tmap,wx,wy)
    sorted_dem_sav = sorted(dem_sav[0], key=lambda x:abs(x['x']-pix_coord[1]) + abs(x['y']-pix_coord[0]))
    #sorted_dem_sav = sorted(dem_sav, key=lambda x:abs(x['wy'].value) )
    print(sorted_dem_sav[0]['res'][0])
    if plot_it:
        fig = plt.figure(figsize=(6,6))
        cax = fig.add_subplot(111)
        #temps = np.logspace(5.7, 7.6, num=42)
        temps = np.logspace(5.7, 7.6, num=42)
        mlogt = ([np.mean([(np.log10(temps[i])), np.log10((temps[i + 1]))]) \
                  for i in np.arange(0, len(temps) - 1)])
        #cax.errorbar(temps,sorted_dem_sav[0]['res'][0],xerr=sorted_dem_sav[0]['res'][2],yerr=sorted_dem_sav[0]['res'][1],fmt='or',ecolor='gray', elinewidth=3, capsize=0)
        cax.errorbar(mlogt,sorted_dem_sav[0]['res'][0],xerr=sorted_dem_sav[0]['res'][2],yerr=sorted_dem_sav[0]['res'][1])
        cax.set_xlabel('$\mathrm{\log_{10}T\;[K]}$')
        cax.set_ylabel('$\mathrm{DEM\;[cm^{-5}\;K^{-1}]}$')
        #cax.set_ylim([1e20, 4e22])
        #cax.set_xlim([5.7, 7.6])
        #plt.rcParams.update({'font.size': 16})
        cax.set_yscale('log')
        #plt.savefig('demregpy_aiapxl.png', bbox_inches='tight')
        plt.show()
    return (sorted_dem_sav[0],dem_calculator(sorted_dem_sav[0],get_em=False))

def get_br_on_date(cdate='0820'):
    #e.g. date:     '0702'
    if cdate == '0820':
        #new_map = smap.Map('/Volumes/WD6T/working/20170703/magnetogram/0702/drot_20170702_16_00_Br.fits')
        new_map = smap.Map('/Volumes/Data/20170820/hmi/br.fits')
    else:
        #save_file = '/Volumes/WD6T/working/20170703/magnetogram/{0}/bptr_2017{1}_1630.sav'.format(cdate,cdate)
        save_file = '/Volumes/Data/20170820/hmi/bptr.save'
        #c_mag_folder = '/Volumes/WD6T/working/20170703/magnetogram/{0}/tmp/2017-{1}-{2}/'.format(cdate,cdate[0:2],cdate[2:4])
        #c_mag_fits = glob.glob(c_mag_folder+'hmi.M_720s.*_TAI.magnetogram.fits')[0]
        c_mag_fits = '/Volumes/Data/20170820/hmi/hmi.sharp_720s_nrt.5630.20170820_193600_TAI.magnetogram.fits'
        c_mag_map = smap.Map(c_mag_fits)
        b_save = readsav(save_file)
        br = b_save['bptr'][2,:,:]
        new_map = smap.Map(br, c_mag_map.meta)
        new_map = new_map.rotate(angle=180 * u.deg)
    return new_map

def read_imfit_res_position(tim,ref_map=None, ret_fr_avg=None,subed=False,phase=2):
    if phase==2:
        if subed:
            with open('/Volumes/Data/20170820/fitting_result/imfit/subed_10s_imfit_t{0:0=4d}.p'.format(tim), 'rb') as imfit_save:
                imf = pickle.load(imfit_save,encoding='latin1')
            imfit_save.close()
        else:
            with open('/Volumes/Data/20170820/fitting_result/imfit/10s_imfit_t{0:0=4d}.p'.format(tim),'rb') as imfit_save:
                imf = pickle.load(imfit_save,encoding='latin1')
            imfit_save.close()
    if phase == 1:
        if subed:
            with open('/Volumes/Data/20170820/fitting_result/imfit/2022/subed_1s_imfit_t{0:0=4d}.p'.format(tim),
                      'rb') as imfit_save:
                imf = pickle.load(imfit_save,encoding='latin1')
            imfit_save.close()
        else:
            with open('/Volumes/Data/20170820/fitting_result/imfit/1s_imfit_t{0:0=4d}.p'.format(tim),'rb') as imfit_save:
                imf = pickle.load(imfit_save,encoding='latin1')
            imfit_save.close()
    if ref_map is None:
        if phase==2:
            ref_map = smap.Map(in_dic_10s[tim]['radio_sbs'][10])
        elif phase==1:
            ref_map = smap.Map(in_dic_1s[tim]['radio_sbs'][10])
    wx = np.zeros((30))
    wxe= np.zeros((30))
    wy = np.zeros((30))
    wye= np.zeros((30))
    for i in range(29):
        if imf[i]['converged'][0]:
            cpx = imf[i]['results']['component0']['pixelcoords'][0]
            cpy = imf[i]['results']['component0']['pixelcoords'][1]
            [wx[i],wy[i]] = pot.get_world(cpx,cpy,cur_map=ref_map)
            wxe[i] = max(0.6-0.025*i, 0.05)
            wye[i] = max(0.6-0.025*i, 0.05)
    if ret_fr_avg is not None:
        return [np.nanmean(wx[ret_fr_avg[0]:ret_fr_avg[1]]), np.nanmean(wy[ret_fr_avg[0]:ret_fr_avg[1]])]
    else:
        return (wx,wy,wxe,wye)


def add_patch_to_existing_legend(legend,label,color='k',marker=None,finished=False):
    from matplotlib.patches import Patch
    ax = legend.axes

    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([0], [0], marker=marker, color=color,
                          markerfacecolor='g', markersize=15))
    labels.append(label)
    legend._legend_box = None
    legend._init_legend_box(handles, labels)
    legend._set_loc(legend._loc)
    legend.set_title(legend.get_title().get_text())

def find_detrending_substitute(filename, kw):
    target_path = '/Volumes/Data/20170820/aia_detrending/{}'.format(kw.split('.')[1])
    if kw=='Z.131.':
        target_path = '/Volumes/Data/20170820/aia_detrending/131_non_blooming'
    #print(target_path)
    print(target_path, os.path.basename(filename))
    res = tp.makelist(tdir=target_path, keyword1=os.path.basename(filename))[0]
    return res

def differential_rotate_coord(wxy,delta_sec, ref_map):
    start_coord = SkyCoord(wxy[0]* u.arcsec, wxy[1]* u.arcsec, frame=ref_map.coordinate_frame)
    dt = TimeDelta(delta_sec * u.second)
    future_date = Time(ref_map.meta['t_obs'],format='isot') + dt
    rotated_coord = solar_rotate_coordinate(start_coord, new_observer_time=future_date)
    return (rotated_coord.Tx, rotated_coord.Ty)
    #coord = SkyCoord([start_coord.Tx, rotated_coord.Tx],
    #                 [start_coord.Ty, rotated_coord.Ty],
    #                 frame=aia_map.coordinate_frame)

def diff_rotate_map(kw,tim,delta_t_sec):
    aiamap = smap.Map(in_dic_10s[tim][kw])
    in_time = aiamap.date

    out_time = in_time + delta_t_sec * u.second
    out_frame = Helioprojective(observer='earth', obstime=out_time,
                                rsun=aiamap.coordinate_frame.rsun)

    out_center = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=out_frame)
    header = smap.make_fitswcs_header(aiamap.data.shape,
                                           out_center,
                                           scale=u.Quantity(aiamap.scale))
    out_wcs = WCS(header)
    print('here')
    with propagate_with_solar_surface():
        out_warp = aiamap.reproject_to(out_wcs)
    return out_warp

def plot_velocity_on_axes(save_file,start_tim_string,cax=None):
    from PIL import Image
    import cv2
    st = Time(start_tim_string,format='iso')
    with open(save_file,'rb') as sf:
        save_res = pickle.load(sf,encoding='latin1')
    sf.close()
    time_list = []
    for xi, cx in enumerate(save_res[0]['x']):
        td = TimeDelta(cx, format='sec')
        ctime = st+td
        time_list.append(ctime.mjd)
    #time_list_i = [time_list]
    #time_list_o = Image.fromarray(np.array(time_list_i).astype('float'))
    #time_list_o = time_list_o.resize((len(save_res[0]['cutslit']['speed']),1),Image.NEAREST)
    #time_list_o = np.array(time_list_o)
    time_list_o = np.linspace(time_list[0],time_list[-1],len(save_res[0]['cutslit']['speed']))
    plt_d = Time(time_list_o,format='mjd').plot_date
    cax.plot_date(plt_d,save_res[0]['cutslit']['speed'], 'r', label='Speed')
    cax.set_ylabel('Velocity $\mathrm{ARCSEC\;[S^{-1}]}$')

def total_non_thermal_energy_calculator(iii):
    #pf3 = '/Volumes/Data/20170820/fitting_result/emcee/guide_field_fix_delta_nnth_area/pwl_world_emcee_res_at_-946_111_T215.p'
    pf3 = glob.glob('/Volumes/Data/20170820/fitting_result/emcee/20220503/pwl_world_emcee_res_at_*_T{0}.p'.format(iii))[0]
    r3 = []
    dc3={}
    with open(pf3, 'rb') as emcee_save:
        emcee_res_3 = pickle.load(emcee_save, encoding='latin1')
    emcee_save.close()
    erg_bin = np.zeros(10000)
    nod = emcee_res_3.flatchain['n_nth']
    for ni, cnd in enumerate(tqdm(nod)):
        delta = emcee_res_3.flatchain['delta'][ni]
        vo_n = 4 / 3 * np.pi * ((emcee_res_3.flatchain['area'][ni] / 2) ** 3)
        to_n = cnd * vo_n
        for ci in range(24, 10000):
            erg_bin[ci] = 1.e9 * (ci ** (-delta)) * ci
        erg_bin *= to_n / 1.e9
        r3.append(np.sum(erg_bin[24:10000]))
    dc3['res'] = r3
    dc3['idx'] = iii
    csave_file = '/Volumes/Data/20170820/fitting_result/emcee/20220503/tne_res_{}.p'.format(iii)
    pickle.dump(dc3,open(csave_file,'wb'))
    return 0

def alfven_velocity_calculator(bg,ne):
    #https://en.wikipedia.org/wiki/Alfvén_wave
    #bg in Gauss, ne in cm-3
    me = 9.11e-31
    mp = 1.67e-27
    va = 2.18e11* (me/mp)**-0.5 *ne**-0.5 * bg
    print(va)
    #return va in km/s
    return va/100./1000.

def mp_wrapper():
    inp_list = np.arange(213,247).tolist()
    ncpu=8
    cp = Pool(ncpu)
    final_resl = cp.map(total_non_thermal_energy_calculator, inp_list)

def gradient_with_error(inp_arr, step, inp_err):
    res = np.ones_like(inp_arr)
    err = np.ones_like(inp_arr)
    for i in range(inp_arr.shape[0]):
        if i==0: continue
        res[i] = (inp_arr[i]-inp_arr[i-1])/step
        upl = (inp_arr[i]+inp_err[i]/2. - (inp_arr[i-1]-inp_err[i-1]/2.))/step
        lol = (inp_arr[i]-inp_err[i]/2. - (inp_arr[i-1]+inp_err[i-1]/2.))/step
        err[i] = abs(upl-lol)
    res[0] = res[1]
    err[0] = err[1]
    return(res, err)



'''
def main():
    mp_wrapper()


if __name__ == '__main__':
    main()
'''