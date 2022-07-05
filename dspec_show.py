import plotting_tools as pt
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import SymLogNorm
import numpy as np
#import sunpy.cm as cm
#import sunpy.map as smap
import pickle
from astropy.time import Time
import copy

from astropy.io import fits
from scipy.signal import savgol_filter
import warnings
from make_dict_list import makelist
import sunpy.map as smap
import warnings


def tmp_spectrum_creator():
    import math
    warnings.filterwarnings('ignore')
    with open('/Volumes/Data/20170820/20220511/info/bmsize.p', 'rb') as fbmsize:
        bmsize = pickle.load(fbmsize, encoding='latin1')
    fbmsize.close()
    with open('/Volumes/Data/20170820/20220511/info/cfreqs.p', 'rb') as fcfreq:
        cfreq = pickle.load(fcfreq, encoding='latin1')
    fcfreq.close()
    #cfov = [[865,-408],[1070,-174]]
    #cfov = [[913,-288],[961,-232]]
    fov_list = [[[862, -374], [1066 ,-176]], [[909, -349], [989, -223]], [[909, -345], [972, -226]],
                [[910, -298], [957, -235]], [[910, -295], [952, -243]]]
    fov_list = [[[822, -414], [1106 ,-136]], [[879, -379], [1019, -193]], [[909, -345], [972, -226]],
                [[910, -298], [957, -235]], [[910, -295], [952, -243]]]
    p_list = makelist(tdir='/Volumes/Data/20170820/20220511/eovsa/eovsa_full/slfcal/images_slfcaled/',keyword1='_t41_'
                      , keyword2='fits')
    p_list.sort()
    b_list = makelist(tdir='/Volumes/Data/20170820/20220511/eovsa/eovsa_full/slfcal/images_slfcaled/',keyword1='_t146_'
                      , keyword2='fits')
    b_list.sort()
    print(p_list)
    print(b_list)
    p_arr = np.zeros((50))
    b_arr = np.zeros((50))
    for spwi in range(50):
        findex = min(int(math.floor(spwi/5)),4)
        cur_sub_map = pt.make_sub_map(cur_map=smap.Map(p_list[spwi]), fov=fov_list[findex])
        new_sum_data = cur_sub_map.data
        old_sum = np.nansum(new_sum_data.clip(min=0.0))
        cur_sfu = pt.test_convertion(tb=old_sum, pix=2.0, bmmin=bmsize[spwi], freq=cfreq[spwi], switch='tb2sfu')
        p_arr[spwi] = cur_sfu
    for spwi in range(50):
        cur_sub_map = pt.make_sub_map(cur_map=smap.Map(b_list[spwi]), fov=fov_list[min(int(math.floor(spwi/5)),4)])
        new_sum_data = cur_sub_map.data
        old_sum = np.nansum(new_sum_data.clip(min=0.0))
        cur_sfu = pt.test_convertion(tb=old_sum, pix=2.0, bmmin=bmsize[spwi], freq=cfreq[spwi], switch='tb2sfu')
        b_arr[spwi] = cur_sfu
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6, 4))
    axs[0].plot(cfreq, b_arr, label='bkg')
    axs[0].plot(cfreq, p_arr , label='Peak')
    axs[1].plot(cfreq, np.subtract(p_arr, b_arr),
                label='BKGsubed_Peak')
    axs[0].legend()
    axs[1].legend()
    axs[1].set_ylim([-1, 90])
    axs[0].set_ylim([-1, 150])
    plt.show()
    return (p_arr,b_arr,np.subtract(p_arr, b_arr))
def all_eovsa_image(tim):
    fig, axs = plt.subplots(nrows=5, ncols=7, sharex=True, sharey=True, figsize=(12, 8))
    axs = [eaxes for erow in axs for eaxes in erow]
    cfov = [[-1050,0],[-850,200]]
    for i in range(len(cfreq)):
        cax = axs[i]
        ep.justmap(cax = cax, tim = tim, kw='radio', fov=cfov,custom_file=in_dic_10s[tim]['radio_sbs'][i])
        plt.show()
def all_sub_dspec():
    '''
    fig=plt.figure(figsize=(8,6))
    cax1=fig.add_subplot(3,4,1)
    cax2=fig.add_subplot(3,4,5)
    cax3=fig.add_subplot(3,4,9)
    cax4 = plt.subplot2grid((3, 4), (0, 2), colspan=2, rowspan=1)
    cax5 = plt.subplot2grid((3, 4), (1, 2), colspan=2, rowspan=1)
    cax6 = plt.subplot2grid((3, 4), (2, 2), colspan=2, rowspan=1)
    '''
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8, 6))
    dmin=3.
    dmax=100.
    #cdspec1=pickle.load(open('/Volumes/Data/20170715/dspec/200_400_770_120_850_40.p','rb'),encoding='latin1')
    #cdspec2=pickle.load(open('/Volumes/Data/20170715/dspec/200_400_730_200_810_120.p','rb'),encoding='latin1')
    #cdspec1=pickle.load(open('/Volumes/Data/20170715/dspec/bkg_1s_770_120_850_40.p','rb'),encoding='latin1')
    #cdspec2=pickle.load(open('/Volumes/Data/20170715/dspec/bkg_1s_730_200_810_120.p','rb'),encoding='latin1')
    #cdspec1=spatial_dspec_minus(phase=1,source=1)
    #cdspec2=pickle.load(open('/Volumes/Data/20170715/dspec/all_1s_770_120_850_40.p','rb'),encoding='latin1')
    #cdspec2=pickle.load(open('/Volumes/Data/20170715/dspec/full_fov_1s.p','rb'),encoding='latin1')
    #mdds.plt_dspec(specdata='/Volumes/Data/20170715/eovsa/msdata/IDB20170715_concat.ms.dspec.npz',plot_goes_hessi=False,pax=axs[0],dmin=dmin,dmax=dmax)
    cdspec1=pickle.load(open('/Volumes/Data/20170715/dspec/10s_all_770_120_850_40.p','rb'),encoding='latin1')
    #cdspec1=pickle.load(open('/Volumes/Data/20170715/dspec/10s_all_730_200_810_120.p','rb'),encoding='latin1')
    cdspec2=spatial_dspec_minus(phase=2,source=1)
    #cdspec2=spatial_dspec_minus(phase=2,source=2)

    mdds.plt_dspec(specdata=cdspec1,plot_goes_hessi=False,pax=axs[0],dmax=dmax,dmin=dmin)
    mdds.plt_dspec(specdata=cdspec2,plot_goes_hessi=False,pax=axs[1],dmax=dmax,dmin=dmin)
    plt.show()
    #return (cdspec1,cdspec2)
def spatial_dspec_minus(phase=1,source=1,whole_fov=False,do_plot=True):
    with open('/Volumes/Data/20170820/20220511/info/bmsize.p', 'rb') as fbmsize:
        bmsize = pickle.load(fbmsize, encoding='latin1')
    fbmsize.close()
    with open('/Volumes/Data/20170820/20220511/info/cfreqs.p', 'rb') as fcfreq:
        cfreq = pickle.load(fcfreq, encoding='latin1')
    fcfreq.close()
    cur_dic = '/Volumes/Data/20170820/20220511/info/20220511_10s_long_aia.p'
    in_dic_file = open(cur_dic, 'rb')
    in_dic = pickle.load(in_dic_file)
    in_dic_file.close()
    if True:
        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6, 4))

        with open('/Volumes/Data/20170820/20220511/eovsa/image_dspec/10s_dspec_865_408_1070_174.p', 'rb') as image_dspec_save:
            image_dspec = pickle.load(image_dspec_save)
        image_dspec_save.close()
        axs[0].plot(cfreq,image_dspec['spec'][0,0,:,178-34],label='bkg')
        axs[0].plot(cfreq,image_dspec['spec'][0,0,:,73-34],label='Peak')
        axs[1].plot(cfreq,np.subtract(image_dspec['spec'][0,0,:,73-34],image_dspec['spec'][0,0,:,178-34]),label='BKGsubed_Peak')
        axs[0].legend()
        axs[1].legend()
        axs[1].set_ylim([-1,90])
        axs[0].set_ylim([-1,90])
        #ep.axes_helper(cax=axs[1], ylabel='FD[sfu]', xlabel='Freq_Hz')
        #ep.axes_helper(cax=axs[0], ylabel='FD[sfu]')
        fig.suptitle('image integrated FD')

        plt.show()
    #return new_sig
    #return np.subtract(np.mean(image_dspec['spec'][0,0,:,80:88],axis=1),np.mean(image_dspec_bkg['spec'][0,0,:,0:3],axis=1))
    #return np.subtract(np.max(image_dspec['spec'][0,0,:,1:5],axis=1),np.mean(image_dspec_bkg['spec'][0,0,:,0:3],axis=1))
    #return np.subtract(image_dspec['spec'][0,0,:,81],np.min(image_dspec_bkg['spec'][0,0,:,0:3],axis=1))
    #return np.subtract(image_dspec['spec'][0,0,:,49],image_dspec_bkg['spec'][0,0,:,0])
    return np.subtract(image_dspec['spec'][0,0,:,73-34],image_dspec['spec'][0,0,:,178-34])

def cali_tp_spectrum():
    cfits='/Volumes/Data/20170820/20220511/eovsa/EOVSA_TPall_20220511.fts'
    hdul=fits.open(cfits,mode='readonly')
    dfreq=   copy.deepcopy(hdul[1].data)
    ds_data= copy.deepcopy(hdul[0].data)
    ds_time= copy.deepcopy(hdul[2].data)
    hdul.close()
    #timerange1=['2017-07-03 16:02:01','2017-07-03 16:02:31']
    #timerange1=['2017-07-03 16:12:01','2017-07-03 16:12:31']
    #timerange1=['2017-07-03 16:12:01','2017-07-03 16:12:31']
    #timerange1=['2017-07-03 16:12:01','2017-07-03 16:12:31']
    timerange2=['2022-05-11 18:42:10','2022-05-11 18:42:20']
    #timerange2=['2017-07-03 16:13:41','2017-07-03 16:13:59']
    #timerange1=['2017-07-15 18:52:00','2017-07-15 19:10:00']
    #-------for 1s in_dic----------
    #timerange2=['2017-07-03 16:13:46','2017-07-03 16:13:50']
    timerange1=['2022-05-11 18:59:40','2022-05-11 18:59:50']
    #timerange2=['2017-07-03 16:13:26','2017-07-03 16:13:33']
    #timerange2=['2017-07-15 19:31:07','2017-07-15 19:31:17']
    #timerange2=['2017-07-15 19:50:00','2017-07-15 20:00:00']
    tr_plt1 = Time(timerange1)
    tr_plt2 = Time(timerange2)
    mjd_list=[]
    for cds_t in ds_time:
        tc1=cds_t[0]
        tc2=cds_t[1]/1000./8640.*0.1
        mjd_list.append(tc1*1.0+tc2*1.0)
    spec_tim=Time(mjd_list, format='mjd')
    tidx_spec1, = np.where((spec_tim >= tr_plt1[0]) & (spec_tim <= tr_plt1[1]))
    tidx_spec2, = np.where((spec_tim >= tr_plt2[0]) & (spec_tim <= tr_plt2[1]))
    print(tidx_spec1)
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6, 4))
    #axs.plot(np.arange(134),np.median(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1))
    #axs[0].plot(np.arange(134),np.mean(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1),label='background')
    axs[0].plot(dfreq.astype(float)[:],np.mean(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1),label='background')
    #axs[0].plot(np.arange(134),np.median(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),label='flare peak')
    #axs[0].plot(dfreq.astype(float)[:],np.max(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),label='flare peak')
    axs[0].plot(dfreq.astype(float)[:],np.mean(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),label='flare peak')
    #axs[0].plot(np.arange(134),np.mean(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),label='flare peak')
    #axs[1].plot(np.arange(134),np.median(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1)-np.median(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1),label='bkg_subtracted_flare_peak')
    subed_sig=np.subtract(np.max(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),np.mean(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1))
    #axs[1].plot(dfreq.astype(float)[:],np.subtract(np.max(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),np.mean(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1)),label='bkg_subtracted_flare_peak')
    axs[1].plot(dfreq.astype(float)[:],np.subtract(np.mean(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),np.mean(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1)),label='bkg_subtracted_flare_peak')
    #axs[1].plot(dfreq.astype(float)[:],np.subtract(np.median(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),np.median(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1)),label='bkg_subtracted_flare_peak')
    axs[0].legend()
    axs[1].legend()
    #axs[1].set_ylim([-1, 40])
    #axs[0].set_ylim([-1, 7])
    #ep.axes_helper(cax=axs[0],ylabel='FD[sfu]')
    #ep.axes_helper(cax=axs[1],ylabel='FD[sfu]',xlabel='Freq_GHz')
    #print(np.median(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1))
    #axs.set_ylim([55.,560.])
    fig.suptitle('EOVSA_TP')
    plt.show()
    return (np.subtract(np.max(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),np.mean(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1)),dfreq.astype(float)[:])
    #return (np.subtract(np.max(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),np.median(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1)))
'''
def tp_fits_to_spectrum():
    cfits='/Volumes/WD6T/working/20170703/eovsa/tpall/EOVSA_TPall_20170703.fts'
    hdul=fits.open(cfits,mode='readonly')
    dfreq=   copy.deepcopy(hdul[1].data)
    ds_data= copy.deepcopy(hdul[0].data)
    ds_time= copy.deepcopy(hdul[2].data)
    hdul.close()
    timerange1=['2017-07-03 16:00:01','2017-07-03 16:00:29']
    #timerange1=['2017-07-15 18:52:00','2017-07-15 19:10:00']
    timerange2=['2017-07-03 16:13:10','2017-07-03 16:13:20']
    #timerange2=['2017-07-03 16:13:26','2017-07-03 16:13:33']
    #timerange2=['2017-07-15 19:31:07','2017-07-15 19:31:17']
    #timerange2=['2017-07-15 19:50:00','2017-07-15 20:00:00']
    tr_plt1 = Time(timerange1)
    tr_plt2 = Time(timerange2)
    mjd_list=[]
    for cds_t in ds_time:
        tc1=cds_t[0]
        tc2=cds_t[1]/1000./8640.*0.1
        mjd_list.append(tc1*1.0+tc2*1.0)
    spec_tim=Time(mjd_list, format='mjd')
    time_interval = 1.0
    for tint in range(180):
        init_t = Time(init_time, format='isot')
        timed1 = TimeDelta((tint * time_interval) * 1.0, format='sec')
        timed2 = TimeDelta(time_interval, format='sec')
        start_time = init_t + timed1
        end_time = start_time + timed2


    tidx_spec1, = np.where((spec_tim >= tr_plt1[0]) & (spec_tim <= tr_plt1[1]))
    tidx_spec2, = np.where((spec_tim >= tr_plt2[0]) & (spec_tim <= tr_plt2[1]))
    print(tidx_spec1)
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6, 4))
    #axs.plot(np.arange(134),np.median(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1))
    #axs[0].plot(np.arange(134),np.mean(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1),label='background')
    axs[0].plot(dfreq.astype(float)[:],np.median(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1),label='background')
    #axs[0].plot(np.arange(134),np.median(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),label='flare peak')
    axs[0].plot(dfreq.astype(float)[:],np.mean(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),label='flare peak')
    #axs[0].plot(np.arange(134),np.mean(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),label='flare peak')
    #axs[1].plot(np.arange(134),np.median(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1)-np.median(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1),label='bkg_subtracted_flare_peak')
    subed_sig=np.subtract(np.max(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),np.mean(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1))
    axs[1].plot(dfreq.astype(float)[:],np.subtract(np.mean(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),np.mean(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1)),label='bkg_subtracted_flare_peak')
    #axs[1].plot(dfreq.astype(float)[:],np.subtract(np.median(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),np.median(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1)),label='bkg_subtracted_flare_peak')
    axs[0].legend()
    axs[1].legend()
    axs[1].set_ylim([-1, 8])
    #axs[0].set_ylim([-1, 7])
    ep.axes_helper(cax=axs[0],ylabel='FD[sfu]')
    ep.axes_helper(cax=axs[1],ylabel='FD[sfu]',xlabel='Freq_GHz')
    #print(np.median(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1))
    #axs.set_ylim([55.,560.])
    fig.suptitle('EOVSA_TP')
    plt.show()
    return (np.subtract(np.mean(ds_data[:,tidx_spec2[0]:tidx_spec2[1]],axis=1),np.mean(ds_data[:,tidx_spec1[0]:tidx_spec1[1]],axis=1)),dfreq.astype(float)[:])
'''

def cal_ratio():
    with open('/Volumes/Data/20170820/20220511/info/cfreqs.p', 'rb') as fcfreq:
        cfreq = pickle.load(fcfreq, encoding='latin1')
    fcfreq.close()
    if True:
        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6, 4))

        with open('/Volumes/Data/20170820/20220511/eovsa/image_dspec/10s_dspec_865_408_1070_174.p', 'rb') as image_dspec_save:
            image_dspec = pickle.load(image_dspec_save)
        image_dspec_save.close()
        img_res_all = tmp_spectrum_creator()
        img_res = img_res_all[2]
        axs[0].plot(cfreq,img_res_all[1],label='bkg')
        axs[0].plot(cfreq,img_res_all[0],label='Peak')
        axs[1].plot(cfreq,img_res,label='BKGsubed_Peak')
        axs[0].legend()
        axs[1].legend()
        axs[1].set_ylim([-1,100])
        axs[0].set_ylim([-1,100])
        #ep.axes_helper(cax=axs[1], ylabel='FD[sfu]', xlabel='Freq_Hz')
        #ep.axes_helper(cax=axs[0], ylabel='FD[sfu]')
        fig.suptitle('image integrated FD')

        plt.show()

    tp_res=cali_tp_spectrum()
    #img_res=spatial_dspec_minus(do_plot=True)
    print('img integ: ',img_res )
    print('tp: ',tp_res )
    #img_res=spatial_dspec_minus(phase=1,source='whole_fov',do_plot=True)
    #img_res=spatial_dspec_minus(phase=1,source='whole_fov')
    freq_index=[]
    cfreq=cfreq/1.e9
    for cur_cfreq in cfreq:
        freq_index.append(min(range(len(tp_res[1])), key=lambda ii: abs(tp_res[1][ii] - cur_cfreq)))
    print('freq_index is : ', freq_index, len(freq_index))
    tp_plot=tp_res[0][freq_index]
    ratio=np.zeros((50))
    ratio[0] = 1.0
    for iii in range(50):
        if iii == 0: continue
        #ratio[iii]=tp_plot[iii]/img_res['spec'][0,0,iii,309]
        #ratio[iii]=tp_res[iii]/img_res[iii]
        ratio[iii]=tp_plot[iii]/img_res[iii]
        #ratio[iii] = tp_plot[iii] / img_res['spec'][0, 0, iii, 253-210]
    sratio= savgol_filter(ratio, 11, 3)
    fig, axs = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=(4, 4))
    axs.plot(cfreq,ratio,label='original')
    axs.plot(cfreq,sratio,label='smoothed')
    axs.legend()
    axs.set_ylim([0.0,3.0])
    #ep.axes_helper(cax=axs, ylabel='EOVSA_tp/img_integed', xlabel='Freq_GHz')
    print(sratio)
    print(ratio)
    plt.show()
def verify_std(ctim=0):
    cmap = plt.cm.Spectral
    cnorm = mpl.colors.Normalize(vmin=0, vmax=660)
    fig, axs = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(6, 6))
    cur_dic_1s = '/Volumes/WD6T/working/20170703/info/20170703_1s.p'
    in_dic_file_1s = open(cur_dic_1s, 'rb')
    in_dic_1s = pickle.load(in_dic_file_1s)
    in_dic_file_1s.close()
    for ti,cur_dic in enumerate(in_dic_1s):
        if ti%10 != 0: continue
        cstd = sf.get_rms(phase=1,tim=ti)
        axs.plot(cfreq,cstd,color=cmap(cnorm(ti)))
    ep.axes_helper(cax=axs, ylabel='STD Tb [K]', xlabel='Freq [Hz]')
    plt.show()

def dip_during_the_peak(tim=307):
    fig, axs = plt.subplots(nrows=5, ncols=1, sharex=True, figsize=(8, 10))
    pt.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p',timerange=[200,400],x=798,y=-98,kw_list=np.arange(30),radio=True,sample_range=2,pax=axs[0])
    pt.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p',timerange=[200,400],x=802,y=-80,kw_list=np.arange(30),radio=True,sample_range=2,pax=axs[1])
    pt.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p',timerange=[200,400],x=823,y=-104,kw_list=np.arange(30),radio=True,sample_range=2,pax=axs[2])
    pt.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p',timerange=[200,400],x=820,y=-86,kw_list=np.arange(30),radio=True,sample_range=2,pax=axs[3])
    pt.lightcurve(dic_file='/Volumes/WD6T/working/20170703/info/20170703_1s.p',timerange=[200,400],x=820,y=-86,kw_list=np.arange(30),radio=True,sample_range=2,pax=axs[4],dspec_radio_file='/Volumes/Data/20170715/dspec/S1_1s_include_neg.p')
    for uci,ucax in enumerate(axs):
        if uci <=3:
            ep.axes_helper(cax=ucax, ti_di=['in', 'in'],ylabel='${T_b}$ [k]')
            ucax.set_ylim([0.0, 2.8e7])
        else:
            ep.axes_helper(cax=ucax, ti_di=['in', 'in'], ylabel='Flux Density [sfu]')

    fig2=plt.figure(figsize=(7,7))
    cur_dic_1s = '/Volumes/WD6T/working/20170703/info/20170703_1s.p'
    in_dic_file_1s = open(cur_dic_1s, 'rb')
    in_dic_1s = pickle.load(in_dic_file_1s)
    in_dic_file_1s.close()
    kwlist = ['Z.94.', 'Z.131.', 'Z.171.', 'Z.193.', 'Z.211.', \
              'Z.304.', 'Z.335.', 'Z.1600.']
    for ik, ckw in enumerate(kwlist):
        cax = fig2.add_subplot(3,3,ik+1)
        cst=ep.aia_files_before_after_timepoint(tim=tim,ckw=ckw)[0]
        cet=ep.aia_files_before_after_timepoint(tim=tim,ckw=ckw)[1]
        cst_str = tp.read_time(fits_file=in_dic_1s[cst][ckw],instrument='aia').iso.split(' ')[1]
        cet_str = tp.read_time(fits_file=in_dic_1s[cet][ckw],instrument='aia').iso.split(' ')[1]
        ep.diffmap(cax=cax,st=ep.aia_files_before_after_timepoint(tim=tim,ckw=ckw)[0],et=ep.aia_files_before_after_timepoint(tim=tim,ckw=ckw)[1],kw=ckw,fov2=[[775, -150], [875, -50]],high_tres=True,enhance_it=True)
        ep.axes_helper(cax=cax, ti_di=['in', 'in'],add_text=[ckw.replace('Z', '').replace('.', '')+'\n'+cst_str+'~\n'+cet_str,0.1,0.05,'w'])
    cax = fig2.add_subplot(3,3,9)
    ep.justmap(cax=cax,tim=tim,kw='Z.1600.',fov=[[775, -150], [875, -50]],high_tres=True)
    ep.plot_B_contour(tim=tim,ax=cax,levels=[-100,100],cfov=[[775, -150], [875, -50]])
    plt.show()
def dip_ratio():
    tim=307
    tim1=302
    tim2=314
    dspec_radio_file='/Volumes/Data/20170715/dspec/S1_1s_include_neg.p'
    if dspec_radio_file is not None:
        with open(dspec_radio_file, 'rb') as rdspec:
            rds = pickle.load(rdspec)
        rdspec.close()
        rds_data = rds['spec'][0, 0, :, :]
    ratio1=np.zeros((30))
    ratio2=np.zeros((30))
    fig, ax = plt.subplots()
    for spwi in range(30):
        ratio1[spwi]= rds_data[spwi,tim]/rds_data[spwi,tim1]
        ratio2[spwi]= rds_data[spwi,tim]/rds_data[spwi,tim2]
    ax.plot(cfreq,ratio1,linestyle='dotted',label='ratio to 1st peak')
    ax.plot(cfreq,ratio2,linestyle='dashed',label='ratio to 2nd peak')
    ep.axes_helper(cax=ax,ti_di=['in', 'in'],ylabel='Ratio _FluxDensity',xlabel='Frequency [Hz]')
    ax.legend()
    plt.show()
def main():
    #spatial_dspec_minus()
    #cali_tp_spectrum()
    #tmp_spectrum_creator()
    cal_ratio()

if __name__ == '__main__':
    main()
