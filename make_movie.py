import ex_plot_1 as ep
import plotting_tools as pt
import matplotlib.pyplot as plt
#from md_utils import plot_map
import matplotlib as mpl
from matplotlib.colors import SymLogNorm
#import md_dspec2 as mdds
import numpy as np
#import sunpy.cm as cm
#import sunpy.map as smap
import pickle
from astropy.time import Time
from matplotlib.lines import Line2D

import lightcurves_two_axes as lcta
import warnings

import copy

def all_aia_movie_contour(tim,plot_rhessi=None,high_tres=False, inp_fig=None):
    warnings.filterwarnings('ignore')
    with open('/Volumes/Data/20170820/20220511/info/bmsize.p', 'rb') as fbmsize:
        bmsize = pickle.load(fbmsize, encoding='latin1')
    fbmsize.close()
    with open('/Volumes/Data/20170820/20220511/info/eovsa_max_ref_105.p', 'rb') as contour_ref:
        max_ref = pickle.load(contour_ref, encoding='latin1')
    contour_ref.close()
    max_ref = np.asarray(max_ref)
    with open('/Volumes/Data/20170820/20220511/info/cfreqs.p', 'rb') as fcfreq:
        cfreq = pickle.load(fcfreq, encoding='latin1')
    fcfreq.close()
    cur_dic = '/Volumes/Data/20170820/20220511/info/20220511_10s_long_aia.p'
    in_dic_file = open(cur_dic, 'rb')
    in_dic = pickle.load(in_dic_file)
    in_dic_file.close()
    #plt.ioff()
    if inp_fig is None:
        fig = plt.figure(figsize=(8, 10))
    else:
        fig = inp_fig

    #cfov = [[800,-450],[1100,-150]]
    cfov = [[800,-450],[1220,150]]

    kwlist = ['Z.94.', 'Z.131.', 'Z.171.', 'Z.304.']
    aia_axs_list = [1,2,3,4,5,6,7,8,13]
    eovsa_spw_axs_list = [0,3,6,9,12,15,19,23]
    eovsa_spwlist = [0,1,2,3,4,5,6,7,8,9,10,15,16,19,20,21,22]

    #with open('/Volumes/WD6T/working/20170703/info/aia_exposure.p', 'rb') as femax:
    #    aia_exp = pickle.load(femax)
    #femax.close()
    map_axes_list=[]
    for ki, ckw in enumerate(kwlist):
        #if ki!=0: continue
        cax = fig.add_subplot(3, 2, aia_axs_list[ki])
        map_axes_list.append(cax)
        #if ki!=5:
        #    ep.justmap(cax=cax, kw=ckw, fov=cfov, tim=tim, enhance_it=True,judge_material=aia_exp)
        #else:
        ep.justmap(cax=cax, kw=ckw, fov=cfov, tim=tim, enhance_it=True,custom_kernel_size=[3,15],custom_color='gist_gray',prep_aia=False)
        #if ki == 1 and plot_rhessi:

        cl_list = [0.2,0.4,0.6,0.8,0.9]
    #if tim > 33 and tim<209:
    if tim < 0:
        ep.plot_eovsa_contour(tim=tim, ax=cax, fov=cfov, spws=range(40),
                          abs_contourf=True,  abs_level=max_ref * 0.25,
                            #abs_contourf=False, level=0.3,
                          inp_rgba=plt.cm.Spectral,  calp=0.9)
        ep.add_custom_cb(target_ax=cax, cmin=cfreq[0], cmax=cfreq[-1], ctitle='GHz', cmap=plt.cm.Spectral,
                     insert=True, cloc=7, nbins=4, left_tick=True, within_height=0.8)
        cax.text(.5, .8, 'all spws @ 30%', horizontalalignment='center', transform=cax.transAxes, color='w', fontsize=12.0)
    lax2 = plt.subplot2grid((6, 2), (4, 0), colspan=2, rowspan=1)
    lax3 = plt.subplot2grid((6, 2), (5, 0), colspan=2, rowspan=1, sharex=lax2)
    lcta.lightcurves(specfile='/Volumes/Data/20170820/20220511/eovsa/image_dspec/fov_list_10s_dspec.p',
                     goes=True,only_dspec=False,
                     # timerange=['2017-07-15 19:25:00.000', '2017-07-15 19:55:00.000'], custom_ax=[lax1, lax2, lax3])
                     timerange=['2022-05-11 18:20:00.000', '2022-05-11 19:50:00.000'],
                     custom_ax=[lax2, lax3], vmin=0.1, vmax=300.)
    # lcsf.lightcurves(specfile='/Volumes/Data/20170820/eovsa/msdata/IDB20170820_1950_2050.ms.dspec.npz',
    #                      hessifile='/Volumes/Data/20170820/hessi/hessi.sav', goes=True,
    #                      only_dspec=False,
    #                      # timerange=['2017-07-15 19:25:00.000', '2017-07-15 19:55:00.000'], custom_ax=[lax1, lax2, lax3])
    #                      timerange=['2017-08-20 19:00:00.000', '2017-08-20 20:10:00.000'],
    #                      custom_ax=[lax2, lax3, lax4])
    #for axi, iax in enumerate(map_axes_list):
    #        ep.axes_helper(cax=iax,
    #                   add_text=['SDO/AIA' + kwlist[axi].replace('Z', '').replace('.', '') + r'$\AA$', 0.2, 0.1, 'white'])
    ep.axes_helper(cax=map_axes_list[0], xlabel='Solar X [arcsec]', ylabel='Solar Y [arcsec]', ti_di=['out', 'out'],index_letter='A',cur_in_dic=in_dic,
                       ind_let_size=14, ind_let_color='white', no_xtick=False, no_ytick=False,
                top_x=True,add_text=['SDO/AIA' + kwlist[0].replace('Z', '').replace('.', '') + r'$\AA$', 0.2, 0.1, 'white'])
    ep.axes_helper(cax=map_axes_list[1], xlabel='Solar X [arcsec]', ylabel='Solar Y [arcsec]', ti_di=['out', 'out'],index_letter='B',cur_in_dic=in_dic,
                       ind_let_size=14, ind_let_color='white', no_xtick=False, no_ytick=False,right_y=True,
                top_x=True,add_text=['SDO/AIA' + kwlist[1].replace('Z', '').replace('.', '') + r'$\AA$', 0.2, 0.1, 'white'])
    ep.axes_helper(cax=map_axes_list[2], xlabel='Solar X [arcsec]', ylabel='Solar Y [arcsec]', ti_di=['out', 'out'],index_letter='C',cur_in_dic=in_dic,
                       ind_let_size=14, ind_let_color='white', no_xtick=True, no_ytick=False
                ,add_text=['SDO/AIA' + kwlist[2].replace('Z', '').replace('.', '') + r'$\AA$', 0.2, 0.1, 'white'])
    ep.axes_helper(cax=map_axes_list[3], xlabel='Solar X [arcsec]', ylabel='Solar Y [arcsec]', ti_di=['out', 'out'],index_letter='D',cur_in_dic=in_dic,
                       ind_let_size=14, ind_let_color='white', no_xtick=True, no_ytick=False
                ,add_text=['SDO/AIA' + kwlist[3].replace('Z', '').replace('.', '') + r'$\AA$', 0.2, 0.1, 'white'],right_y=True)
    ep.axes_helper(cax=lax2, ti_di=['out', 'out'],cur_in_dic=in_dic,
                   ind_let_size=14, ind_let_color='w', no_xtick=False, no_ytick=False, index_letter='E',xspan=[[tim,tim+1],'b',1],ind_let_bkg='k')
    ep.axes_helper(cax=lax3, ti_di=['out', 'out'],cur_in_dic=in_dic,
                   ind_let_size=14, ind_let_color='w', no_xtick=False, no_ytick=False, index_letter='F',xspan=[[tim,tim+1],'b',1],ind_let_bkg='k')
    fig.suptitle(in_dic[tim]['time'])
    #plt.subplots_adjust(hspace=0.0, wspace=0.0, left=0.1, bottom=0.11, right=0.9, top=0.88)
    # plt.savefig(
    #     '/Volumes/Data/20170820/20220511/movie/10s_larger_aia/aia_eovsa_contour_abs_t{0:0=4d}.png'.format(tim))
    # plt.clf()
    plt.show()
#def only_one_wavelength(tim):

def single_aia_wl(tim, rhessi_maps_1=None,rhessi_maps_2=None,rhessi_maps_3=None,inp_fig=None,patch_fov=None):
    plt.ioff()
    if inp_fig is None:
        fig = plt.figure(figsize=(16, 10))
    else:
        fig = inp_fig
    #small fov figures
    cfov = [[-1050, 0], [-900, 150]]
    #smallest fov:
    #ssfov = [[-955,95],[-930,120]]
    #ssfov = [[-955,70],[-930,95]]
    ssfov=[[-1050, 0], [-900, 150]]
    #big fov figures
    #cfov = [[-1200, -100], [-800, 300]]
    eufov = [[600, 0], [1000, 400]]
    in_dic = in_dic_10s
    #kwlist = ['Z.131.', 'Z.131.']
    kwlist = ['Z.131.', 'Z.131.']
    #kwlist = ['dem', 'em']
    map_axes_list=[]
    cax1 = plt.subplot2grid((3, 4), (0, 0), colspan=2, rowspan=2)
    cax2 = plt.subplot2grid((3, 4), (0, 2), colspan=2, rowspan=2)
    #ep.justmap(cax=cax1, kw=kwlist[0], fov=cfov, tim=tim, enhance_it=True,norm_exp=True,custom_color='gist_gray',high_tres=False)
    #ep.justmap(cax=cax1, kw=kwlist[0], fov=ssfov, tim=tim, enhance_it=True,norm_exp=True,custom_color='gist_gray',high_tres=False)
    if kwlist[0] == 'dem':
        cnorm1 = mpl.colors.Normalize(vmin=5.e5, vmax=3.5e+7)
        im1 = ep.justmap(cax=cax1, kw=kwlist[0], fov=ssfov, inorm=cnorm1,tim=tim,custom_color='gist_gray',high_tres=False)
    else:
        im1 = ep.justmap(cax=cax1, kw=kwlist[0], fov=ssfov, tim=tim, custom_color='gist_gray', high_tres=False)
    #cb1 = plt.colorbar(im1,cax=cax1)
    cnorm2 = SymLogNorm(linthresh=2.2e27 , vmin=2e27,
               vmax=3.55e30)
    if kwlist[0] == 'em':
        cnorm2 = SymLogNorm(linthresh=2.2e27, vmin=2e27,
                            vmax=3.55e30)
        im2 = ep.justmap(cax=cax2, kw=kwlist[1], fov=ssfov, tim=tim, inorm=cnorm2,custom_color='gist_gray',high_tres=False)
    else:
        im2 = ep.justmap(cax=cax2, kw=kwlist[1], fov=ssfov, tim=tim, custom_color='gist_gray', high_tres=False)
    #cb2 = plt.colorbar(im2,cax=cax2)
    #ep.justmap(cax=cax2, kw=kwlist[1], fov=cfov, tim=tim, enhance_it=True,norm_exp=True,custom_color='gist_gray',high_tres=False)im2 =
    #ep.justmap(cax=cax2, kw=kwlist[1], fov=ssfov, tim=tim, enhance_it=True,norm_exp=True,custom_color='gist_gray',high_tres=False)
    if patch_fov is not None:
        for pfi, pf in enumerate(patch_fov):
            ep.add_patch(cax = cax1,cfov = pf)
    if tim >= 184 and tim < 362:
        #ep.plot_eovsa_contour(tim=tim, ax=cax1, fov=cfov, spws=range(31),
        #                      abs_contourf=False, level=0.95,
        #                      inp_rgba=plt.cm.Spectral, cs_start=[0, 31], calp=0.9, high_tres=False)
        #ep.plot_eovsa_contour(tim=tim, ax=cax1, fov=ssfov, spws=range(31),
        #                      abs_contourf=False, level=0.9,
        #                      inp_rgba=plt.cm.Spectral, cs_start=[0, 31], calp=0.9, high_tres=False)
        ep.plot_eovsa_contour(tim=tim, ax=cax2, fov=cfov, spws=range(31),
                              abs_contourf=False, level=0.2,
                              inp_rgba=plt.cm.Spectral, cs_start=[0, 31], calp=0.2, high_tres=False)
        #ep.plot_eovsa_contour(tim=tim, ax=cax1, fov=cfov, spws=range(31),
        #                      abs_contourf=False, level=0.8,
        #                      inp_rgba=plt.cm.Spectral, cs_start=[0, 31], calp=0.7, high_tres=False)
        ep.add_custom_cb(target_ax=cax2, cmin=cfreq[0], cmax=cfreq[30], ctitle='GHz', cmap=plt.cm.Spectral,
                         insert=True, cloc=7, nbins=4, left_tick=True, within_height=1.8)
        ep.add_custom_cb(target_ax=cax1, cmin=5e5, cmax=3.5e7, ctitle='K', cmap=plt.cm.gist_gray,
                         insert=True, cloc=3, nbins=4, left_tick=False, within_height=1.8,spine_color='w')
        ep.add_custom_cb(target_ax=cax2, cmin=2.e27, cmax=3.55e30, ctitle='EM', cmap=plt.cm.gist_gray,
                         insert=True, cloc=2, nbins=4, left_tick=False, within_height=1.8,spine_color='w')

    if tim >= 114 and tim <294:
        hessi_index = int(np.floor((tim-114)/2))
        cur_hessi_map_1 = rhessi_maps_1[hessi_index]
        cur_hessi_map_2 = rhessi_maps_2[hessi_index]
        cur_hessi_map_3 = rhessi_maps_3[hessi_index]
        #him1=pt.plot_rhessi_contour_file(cur_file=cur_hessi_map_1,levels = [0.6,0.9],fov = cfov, ax=cax1,inp_color='tab:green')
        #him2=pt.plot_rhessi_contour_file(cur_file=cur_hessi_map_2,levels = [0.6,0.9],fov = cfov, ax=cax1,inp_color='tab:purple')
        him11=pt.plot_rhessi_contour_file(cur_file=cur_hessi_map_1,levels = [0.5],fov = ssfov, ax=cax1, calpha=0.95,inp_color='tab:green')
        him22=pt.plot_rhessi_contour_file(cur_file=cur_hessi_map_2,levels = [0.5],fov = ssfov, ax=cax1, calpha=0.95,inp_color='tab:purple')
        him23=pt.plot_rhessi_contour_file(cur_file=cur_hessi_map_3,levels = [0.5],fov = ssfov, ax=cax1, calpha=0.95,inp_color='tab:red')
        #legend_elements = [Line2D([0], [0], color='tab:green', lw=0.6, label='6$-$12 keV'),
        #                   Line2D([0], [0], color='tab:red', lw=0.6, label='12$-$25 keV')]
        legend_elements = [Line2D([0], [0], color='tab:green', lw=0.6, label='6$-$12 keV'),
                           Line2D([0], [0], color='tab:purple', lw=0.6, label='12$-$25 keV'),
                            Line2D([0], [0], color='tab:red', lw=0.6, label='25$-$50 keV')]
        cleg = cax1.legend(handles=legend_elements, loc=8, fontsize='12')

    cax1.text(.7, .09, 'SDO/AIA 131'+r'$\AA$'+'\nEOVSA spws @ 80%', horizontalalignment='center', transform=cax1.transAxes, color='w',
                 fontsize=12.0)
    cax2.text(.7, .09, 'SDO/AIA 131'+r'$\AA$'+'\nEOVSA spws @ 50%', horizontalalignment='center', transform=cax2.transAxes, color='w',
                 fontsize=12.0)
    lax1 = plt.subplot2grid((6, 4), (4, 0), colspan=4, rowspan=1)
    lax2 = plt.subplot2grid((6, 4), (5, 0), colspan=4, rowspan=1, sharex=lax1)
    lcta.lightcurves(specfile='/Volumes/Data/20170820/eovsa/msdata/IDB20170820_1950_2050.ms.dspec.npz',
                     hessifile='/Volumes/Data/20170820/hessi/hessi.sav', goes=True,
                     only_dspec=False,
                     # timerange=['2017-07-15 19:25:00.000', '2017-07-15 19:55:00.000'], custom_ax=[lax1, lax2, lax3])
                     timerange=['2017-08-20 19:00:00.000', '2017-08-20 20:10:00.000'],
                     custom_ax=[lax1, lax2])
    cur_time_text = Time(in_dic_10s[tim]['time'], format='iso').isot.split('T')[-1].split('.')[0]
    cax1.text(.2, .8, cur_time_text, horizontalalignment='center', transform=cax1.transAxes,
             color='gold')
    ep.axes_helper(cax=cax1, index_letter='A',xlabel='Solar X [arcsec]', ylabel='Solar Y [arcsec]', ti_di=['out', 'out'],
                   ind_let_size=14, ind_let_color='white', no_xtick=False, no_ytick=False,
                   top_x=True)
    ep.axes_helper(cax=cax2, index_letter='B',xlabel='Solar X [arcsec]', ylabel='Solar Y [arcsec]', ti_di=['out', 'out'],
                   ind_let_size=14, ind_let_color='white', no_xtick=False, no_ytick=False,right_y=True,
                   top_x=True)
    ep.axes_helper(cax=lax1,ti_di=['out', 'out'],index_letter='C',
                   ind_let_size=14, ind_let_color='w', no_xtick=False, no_ytick=False,
                   xspan=[[tim, tim + 1], 'b', 1], ind_let_bkg='k', high_tres=False)
    ep.axes_helper(cax=lax2,ti_di=['out', 'out'],index_letter='D',
                   ind_let_size=14, ind_let_color='w', no_xtick=False, no_ytick=False,
                   xspan=[[tim, tim + 1], 'b', 1], ind_let_bkg='k', high_tres=False)
    #plt.savefig(
    #    '/Volumes/Data/20170820/movie/10s/hessi_eovsa_lower_source/aia_131_cartoon_t{0:0=4d}.png'.format(tim))
    #plt.clf()
    plt.show()

def aia_wl_only(tim, inp_fig=None,patch_fov=None, high_tres=False):
    #plt.ioff()
    if inp_fig is None:
        fig = plt.figure(figsize=(12, 10))
    else:
        fig = inp_fig
    #small fov figures
    cfov = [[-1050, 0], [-900, 150]]
    ckw='Z.171.'
    #smallest fov:
    #ssfov = [[-955,95],[-930,120]]
    #ssfov = [[-955,70],[-930,95]]
    #ssfov=[[-1050, 0], [-900, 150]]
    #ssfov=[[-1100, -50], [-900, 150]]
    ssfov=[[-1200, -100], [-900, 200]]
    #ssfov=[[-960, 55], [-900, 130]]
    #big fov figures
    #cfov = [[-1200, -100], [-800, 300]]
    eufov = [[600, 0], [1000, 400]]
    in_dic = in_dic_10s_la
    if high_tres:
        in_dic=in_dic_1s
    #kwlist = ['Z.131.', 'Z.131.']
    kwlist = [ckw, ckw]
    #kwlist = ['Z.171.', 'Z.171.']
    #kwlist = ['dem', 'em']
    map_axes_list=[]
    cax1 = plt.subplot2grid((3, 4), (0, 0), colspan=2, rowspan=2)
    cax2 = plt.subplot2grid((3, 4), (0, 2), colspan=2, rowspan=2)
    #cnorm1 = mpl.colors.Normalize(vmin=-2000., vmax=9000.)
    #cnorm2 = mpl.colors.Normalize(vmin=-4000., vmax=25000.)
    cnorm1 = mpl.colors.LogNorm(vmin=600, vmax=30000)
    cnorm2 = mpl.colors.LogNorm(vmin=100, vmax=10000)
    cnorm_171_hc = mpl.colors.Normalize(vmin=150, vmax=6000)
    #cnorm_211_hc = mpl.colors.Normalize(vmin=30, vmax=3000)
    #cnorm_193_hc = mpl.colors.Normalize(vmin=100, vmax=6000)
    #im1 = ep.justmap(cax=cax1, kw=kwlist[0], fov=ssfov, tim=tim, inorm=cnorm1,  custom_color='gist_gray', high_tres=high_tres,detrending=False,long_aia=True)
    im1 = ep.justmap(cax=cax1, kw=kwlist[0], fov=ssfov, tim=tim, inorm=cnorm_171_hc,  custom_color='gist_gray', high_tres=high_tres,detrending=False,long_aia=True, enhance_it=True,  custom_kernel_size=[11])
    im2 = ep.justmap(cax=cax2, kw=kwlist[1], fov=ssfov, tim=tim, inorm=cnorm_171_hc,  custom_color='gist_gray', high_tres=high_tres,detrending=False,long_aia=True, enhance_it=False, custom_kernel_size=[11])
    #im2 = ep.justmap(cax=cax2, kw=kwlist[1], fov=ssfov, tim=tim, inorm=cnorm2,  custom_color='gist_gray', high_tres=high_tres,detrending=False,long_aia=True)
    try:
        pass
        #ep.plot_eovsa_contourf(tim=tim, ax=cax2, fov=ssfov, spws=range(1, 29),
        #                   abs_contourf=False, level=0.93,
        #                   alt_cmap=plt.cm.Spectral, cs_start=[1, 31], calp=0.1, high_tres=high_tres,
        #                   bkg_sub=False, single_spw=False)
    except:
        print('no bkg-subed images')
    #ep.plot_eovsa_contourf(tim=tim, ax=cax1, fov=ssfov, spws=range(1, 29),
    #                       abs_contourf=False, level=0.9,
    #                       alt_cmap=plt.cm.Spectral, cs_start=[1, 31], calp=0.3, high_tres=True,
    #                       bkg_sub=False, single_spw=False)
    ep.add_custom_cb(target_ax=cax2, cmin=cfreq[1], cmax=cfreq[30], ctitle='GHz', cmap=plt.cm.Spectral,
                     insert=True, cloc='center left', nbins=4, left_tick=True, within_height=1.8,spine_color='w')
    #ep.plot_B_contour(tim = tim, ax = cax1,levels=[-1000,-100],cfov = ssfov)
    #ep.plot_B_contour(tim = tim, ax = cax2,levels=[-1000,-100],cfov = ssfov)

    #------hessi images---------------
    '''
    hessi_file4 = '/Volumes/Data/20170820/hessi/images/2022/hsi_image_20170820_193820_194000_det_1_3_eg_6_12_i100.fits'
    hessi_file3 = '/Volumes/Data/20170820/hessi/images/2022/hsi_image_20170820_193523_193711_det_1_3_eg_6_12_i100.fits'
    hessi_file2 = '/Volumes/Data/20170820/hessi/images/2022/hsi_image_20170820_193245_193425_det_1_3_eg_6_12_i100.fits'
    hessi_file1 = '/Volumes/Data/20170820/hessi/images/2022/hsi_image_20170820_192651_192809_det_1_3_eg_6_12_i100.fits'
    hessi_file24 = '/Volumes/Data/20170820/hessi/images/2022/hsi_image_20170820_193820_194000_det_1_3_eg_12_25_i100.fits'
    hessi_file23 = '/Volumes/Data/20170820/hessi/images/2022/hsi_image_20170820_193523_193711_det_1_3_eg_12_25_i100.fits'
    hessi_file22 = '/Volumes/Data/20170820/hessi/images/2022/hsi_image_20170820_193245_193425_det_1_3_eg_12_25_i100.fits'
    hessi_file21 = '/Volumes/Data/20170820/hessi/images/2022/hsi_image_20170820_192651_192809_det_1_3_eg_12_25_i100.fits'
    him1 = pt.plot_rhessi_contour_file(cur_file=hessi_file4,   levels=[0.4], fov=ssfov, ax=cax2, inp_color='tab:green',line_width=0.9)
    him21 = pt.plot_rhessi_contour_file(cur_file=hessi_file24, levels=[0.4], fov=ssfov, ax=cax2, inp_color='tab:red'  ,line_width=0.9)
    him11 = pt.plot_rhessi_contour_file(cur_file=hessi_file4,   levels=[0.4], fov=ssfov, ax=cax1, inp_color='tab:green',line_width=0.9)
    him221 = pt.plot_rhessi_contour_file(cur_file=hessi_file24, levels=[0.4], fov=ssfov, ax=cax1, inp_color='tab:red'  ,line_width=0.9)
    '''
    #legend_elements = [Line2D([0], [0], color='tab:green',lw=0.6, label='@6$-$12 keV\n@40%'),
    #                   Line2D([0], [0], color='tab:red',  lw=0.6, label='@12$-$25 keV\n@40%'),
    #                   Line2D([0], [0], color='tab:blue', lw=0.6, label='@1600A\n@20%')]
    legend_elements = [Line2D([0], [0], color='tab:blue', lw=0.6, label='@1600A\n@20%')]
    cleg = cax1.legend(handles=legend_elements, loc='lower left', fontsize='8')
    #------hessi images---------------
    #-----plot_aia_contour------------------------
    ep.plot_aia_contour(tim=tim, ax=cax1, ckw='Z.1600.', levels=[0.2], cco='b', cfov=ssfov, alp=None, txt_label=None, high_tres=high_tres)
    #ep.plot_aia_contour(tim=tim, ax=cax2, ckw='Z.1600.', levels=[0.2], cco='b', cfov=ssfov, alp=None, txt_label=None, high_tres=high_tres)

    if patch_fov is not None:
        for pfi, pf in enumerate(patch_fov):
            ep.add_patch(cax = cax1,cfov = pf)
    #--------------add eovsa contour------------
    '''
    if tim >= 184 and tim < 362:
        ep.plot_eovsa_contour(tim=tim, ax=cax2, fov=cfov, spws=range(31),
                              abs_contourf=False, level=0.2,
                              inp_rgba=plt.cm.Spectral, cs_start=[0, 31], calp=0.7, high_tres=False)
        ep.plot_eovsa_contour(tim=tim, ax=cax1, fov=cfov, spws=range(31),
                              abs_contourf=False, level=0.8,
                              inp_rgba=plt.cm.Spectral, cs_start=[0, 31], calp=0.7, high_tres=False)
        ep.add_custom_cb(target_ax=cax2, cmin=cfreq[0], cmax=cfreq[30], ctitle='GHz', cmap=plt.cm.Spectral,
                      insert=True, cloc=7, nbins=4, left_tick=True, within_height=1.8)
    '''

        #ep.add_custom_cb(target_ax=cax1, cmin=5e5, cmax=3.5e7, ctitle='K', cmap=plt.cm.gist_gray,
        #                 insert=True, cloc=3, nbins=4, left_tick=False, within_height=1.8,spine_color='w')
        #ep.add_custom_cb(target_ax=cax2, cmin=2.e27, cmax=3.55e30, ctitle='EM', cmap=plt.cm.gist_gray,
        #                 insert=True, cloc=2, nbins=4, left_tick=False, within_height=1.8,spine_color='w')
    #cax1.text(.7, .09, 'SDO/AIA 171'+r'$\AA$'+'\nEOVSA spws @ 80%', horizontalalignment='center', transform=cax1.transAxes, color='w',
    #             fontsize=12.0)
    cax1.text(.7, .09, 'SDO/AIA 193'+r'$\AA$', horizontalalignment='center', transform=cax1.transAxes, color='w',
                 fontsize=12.0)
    cax2.text(.7, .09, 'SDO/AIA 193'+r'$\AA$'+'\nsubed EOVSA spws @ 90%', horizontalalignment='center', transform=cax2.transAxes, color='w',
                 fontsize=12.0)
    lax1 = plt.subplot2grid((6, 4), (4, 0), colspan=4, rowspan=1)
    lax2 = plt.subplot2grid((6, 4), (5, 0), colspan=4, rowspan=1, sharex=lax1)
    lcta.lightcurves(specfile='/Volumes/Data/20170820/eovsa/msdata/IDB20170820_1950_2050.ms.dspec.npz',
                     hessifile='/Volumes/Data/20170820/hessi/hessi.sav', goes=True,
                     only_dspec=False,
                     # timerange=['2017-07-15 19:25:00.000', '2017-07-15 19:55:00.000'], custom_ax=[lax1, lax2, lax3])
                     timerange=['2017-08-20 19:00:00.000', '2017-08-20 20:30:00.000'],
                     custom_ax=[lax1, lax2])
    cur_time_text = Time(in_dic[tim]['time'], format='iso').isot.split('T')[-1].split('.')[0]
    cax1.text(.2, .8, cur_time_text, horizontalalignment='center', transform=cax1.transAxes,
             color='gold')
    ep.axes_helper(cax=cax1, index_letter='A',xlabel='Solar X [arcsec]', ylabel='Solar Y [arcsec]', ti_di=['out', 'out'],
                   ind_let_size=14, ind_let_color='white', no_xtick=False, no_ytick=False,
                   top_x=True)
    ep.axes_helper(cax=cax2, index_letter='B',xlabel='Solar X [arcsec]', ylabel='Solar Y [arcsec]', ti_di=['out', 'out'],
                   ind_let_size=14, ind_let_color='white', no_xtick=False, no_ytick=False,right_y=True,
                   top_x=True)
    ep.axes_helper(cax=lax1,ti_di=['out', 'out'],index_letter='C',
                   ind_let_size=14, ind_let_color='w', no_xtick=False, no_ytick=False,
                   xspan=[[tim, tim + 1], 'b', 1], ind_let_bkg='k', high_tres=high_tres)
    ep.axes_helper(cax=lax2,ti_di=['out', 'out'],index_letter='D',
                   ind_let_size=14, ind_let_color='w', no_xtick=False, no_ytick=False,
                   xspan=[[tim, tim + 1], 'b', 1], ind_let_bkg='k', high_tres=high_tres)
    #plt.savefig(
    #    '/Volumes/Data/20170820/movie/10s/hessi_eovsa_lower_source/aia_131_cartoon_t{0:0=4d}.png'.format(tim))
    #    '/Volumes/Data/20170820/movie/1s/subed_n_131_193/1s_131_subed_t{0:0=4d}.png'.format(tim))
    #    '/Volumes/Data/20170820/movie/10s/171_show_rc_inflow_mw/10s_171_rc_inflow_subed_t{0:0=4d}.png'.format(tim))
    #    '/Volumes/Data/20170820/movie/10s/193_show_rc_inflow/10s_193_rc_inflow_subed_t{0:0=4d}.png'.format(tim))
    #    '/Volumes/Data/20170820/movie/10s/171_show_fr_upflow/10s_171_mfr_upflow_subed_t{0:0=4d}.png'.format(tim))
    #plt.clf()
    plt.show()

def movie_wrapper(tr):##
    #for big movie!!!!!
    #cfig = plt.figure(figsize=(16, 10))
    cfig = plt.figure(figsize=(12, 8))
    hm_file1 = '/Volumes/Data/20170820/hessi/images/6_12_all90maps.p'
    hm_file2 = '/Volumes/Data/20170820/hessi/images/12_25_all90maps.p'
    hm_file3 = '/Volumes/Data/20170820/hessi/images/25_50_all90maps.p'
    hml1_o = open(hm_file1, 'rb')
    hml2_o = open(hm_file2, 'rb')
    hml3_o = open(hm_file3, 'rb')
    hml1 = pickle.load(hml1_o)
    hml1_o.close()
    hml2 = pickle.load(hml2_o)
    hml2_o.close()
    hml3 = pickle.load(hml3_o)
    hml3_o.close()
    for ii in range(tr[0],tr[1]):
        #if ii!=214: continue
        if ii!=225: continue
        #single_aia_wl(tim=ii,inp_fig=cfig,rhessi_maps_1=hml1[0],rhessi_maps_2=hml2[0],rhessi_maps_3=hml3[0],patch_fov=[[-965,80],[-951,94]])
        #single_aia_wl(tim=ii,inp_fig=cfig,rhessi_maps_1=hml1[0],rhessi_maps_2=hml2[0],rhessi_maps_3=hml3[0],patch_fov=[[[-945,18],[-903,60]],[[-968,52],[-926,94]],[[-960,94],[-926,136]]])
        single_aia_wl(tim=ii,inp_fig=cfig,rhessi_maps_1=hml1[0],rhessi_maps_2=hml2[0],rhessi_maps_3=hml3[0])
        #all_aia_movie_contour(tim=ii,high_tres=False,plot_rhessi=False,inp_fig=cfig)
        #cfig.clear()

def hi_low_filter_test(tim):
    from itertools import chain
    #fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(10, 10), sharex=True, sharey=True)
    #core_list = np.linspace(1,100,16,dtype=int)
    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(10, 10), sharex=True, sharey=True)
    #core_list = np.linspace(1,10,9,dtype=int)
    core_list = np.linspace(0.1,8,9)
    faxs = list(chain.from_iterable(axs))
    for i in range(9):
        ep.justmap(cax=faxs[i], kw='Z.131.', fov=[[812,-395],[1012,-195]], tim=tim, custom_color='gist_gray',
                    enhance_it=True, custom_kernel_size=[core_list[i]],prep_aia=False)
        faxs[i].set_title('kernel size: {}'.format(core_list[i]))
    plt.show()




def movie_wrapper_all(tr):##
    #for big movie!!!!!
    #cfig = plt.figure(figsize=(16, 10))
    cfig = plt.figure(figsize=(8, 10))
    for ii in range(tr[0],tr[1]):
        #if ii!=200: continue
        #if ii<231: continue
        print(ii)
        all_aia_movie_contour(tim=ii,inp_fig=cfig)
        #aia_wl_only(tim=ii, inp_fig=cfig, high_tres=False)
        cfig.clear()

def main():
    all_aia_movie_contour(120)
    #movie_wrapper_all([1,300])

if __name__ == '__main__':
    main()