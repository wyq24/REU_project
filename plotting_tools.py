import pickle
import numpy as np
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
import sunpy.map as smap
import warnings
import math

def image_to_dynamicspec(fov=None, dic_file=None, spec_method=None, timerange=None,source_size=None,pix_size=2.0,high_tres = False):
    fov = [[930,-307],[942,-287]]
    fov = [[930,-307],[942,-287]]
    fov_list = [[[862,-374],[1066,-176]],[[909,-349],[989,-223]],[[909,-345],[972,-226]],[[910,-298],[957,-235]],[[910,-295],[952,-243]]]
    warnings.filterwarnings('ignore')
    spec_dict = {}
    with open('/Volumes/Data/20170820/20220511/info/bmsize.p', 'rb') as fbmsize:
        bmsize = pickle.load(fbmsize, encoding='latin1')
    fbmsize.close()
    with open('/Volumes/Data/20170820/20220511/info/cfreqs.p', 'rb') as fcfreq:
        cfreq = pickle.load(fcfreq, encoding='latin1')
    fcfreq.close()
    #cur_dic = '/Volumes/Data/20170820/20220511/info/20220511_10s_long_aia.p'
    cur_dic = '/Volumes/Data/20170820/20220511/info/20220511_10s_long.p'
    in_dic_file = open(cur_dic, 'rb')
    in_dic = pickle.load(in_dic_file)
    in_dic_file.close()
    # if fov is 'custom':
    #     fov_list = [[[865,-408],[1070,-174]]]
    #     fov = fov_list[0]
        #fov_list = [[[897,-11],[1018,98]],[[901,4],[1009,92]],[[918,16],[ 986,72]],[[922,15],[980,71]],[[926,31],[964,71]],[[929,35 ],[957,76]]]
    if not timerange:
        timerange = [0, len(in_dic)-1]
    # spec_res = np.zeros((1, 1, 30, timerange[1]-timerange[0],nfov))
    spec_res = np.zeros((1, 1, 50, timerange[1] - timerange[0]))
    time_res = np.zeros((timerange[1] - timerange[0]))
    #sbeam = 35.
    for i in range(timerange[0], timerange[1]):
        print(i)
        time_res[i - timerange[0]] = Time(in_dic[i]['time'], format='iso').mjd * 3600. * 24
        for ii in range(50):
            cfov = fov_list[min(int(math.floor(ii / 5)), 4)]
            #cfov = fov
            #cfov = fov_list[0]
            cur_sub_map = make_sub_map(cur_map=smap.Map(in_dic[i]['radio_sbs'][ii]), fov=cfov)
            # if fov=='custom':
            #     print('use custom fov')
            #     #cfov = fov_list[int(np.floor(ii/5))]
            #     # cfov = fov_list[0]
            #     # fov = cfov
            #     cfov = fov_list[min(int(math.floor(ii/5)),4)]
            #     #print(in_dic[i]['radio_sbs'][ii])
            #     cur_sub_map = make_sub_map(cur_map=smap.Map(in_dic[i]['radio_sbs'][ii]), fov=cfov)
            # elif fov:
            #     cur_sub_map = make_sub_map(cur_map=smap.Map(in_dic[i]['radio_sbs'][ii]), fov=fov)
            # else:
            #     cur_sub_map =smap.Map(in_dic[i]['radio_sbs'][ii])
            new_sum_data=cur_sub_map.data
            old_sum = np.nansum(new_sum_data.clip(min=0.0))
            #old_sum = np.nansum(new_sum_data)
            #old_sum = np.nansum(new_sum_data)
            #cur_sfu=mw.new_2021_tb2sfu(tb=old_sum,pixel_size=2.0,bmsize=bmsize[ii],freq=cfreq[ii])
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!we /2. here
            cur_sfu=test_convertion(tb=old_sum,pix=pix_size,bmmin=bmsize[ii],freq=cfreq[ii],switch='tb2sfu')
            #cur_sfu=mw.tb2sfu(tb=old_sum,size=source_size,freq=cfreq[ii]*1.e9)
            print('current sfu:  ',cur_sfu)
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!we /2. here
            spec_res[0, 0, ii, i - timerange[0]] = cur_sfu
    spec_dict['spec'] = spec_res
    spec_dict['tim'] = time_res
    spec_dict['freq'] = cfreq
    spec_dict['bl'] = ''
    #cur_name = '/Volumes/Data/20170820/eovsa/image_dspec/10s_dspec_{0}_{1}_{2}_{3}.p'.format(abs(cfov[0][0]),abs(cfov[0][1]),abs(cfov[1][0]),abs(cfov[1][1]))
    #cur_name = '/Volumes/Data/20170820/eovsa/image_dspec/part/10s_dspec_{0}_{1}_{2}_{3}.p'.format(abs(fov[0][0]),abs(fov[0][1]),abs(fov[1][0]),abs(fov[1][1]))
    #cur_name = '/Volumes/Data/20170820/20220511/eovsa/image_dspec/subed_10s_dspec_{0}_{1}_{2}_{3}.p'.format(abs(cfov[0][0]),abs(cfov[0][1]),abs(cfov[1][0]),abs(cfov[1][1]))
    #cur_name = '/Volumes/Data/20170820/20220511/eovsa/image_dspec/submap/lower_fov_list_10s_{0}_{1}_{2}_{3}_dspec.p'.format(abs(cfov[0][0]),abs(cfov[0][1]),abs(cfov[1][0]),abs(cfov[1][1]))
    cur_name = '/Volumes/Data/20170820/20220511/eovsa/image_dspec/fov_list_10s_dspec.p'
    #cur_name = '/Volumes/Data/20170820/eovsa/image_dspec/bkg_10s_dspec_{0}_{1}_{2}_{3}.p'.format(abs(cfov[0][0]),abs(cfov[0][1]),abs(cfov[1][0]),abs(cfov[1][1]))
    pickle.dump(spec_dict, open(cur_name,'wb'))
    return spec_dict

def make_sub_map(cur_map, fov=None, ref_map=None):
    if not ref_map:
        ref_map = cur_map
    bole = SkyCoord(fov[0][0] * u.arcsec, fov[0][1] * u.arcsec, frame=ref_map.coordinate_frame)
    tori = SkyCoord(fov[1][0] * u.arcsec, fov[1][1] * u.arcsec, frame=ref_map.coordinate_frame)
    sub_map = cur_map.submap(bole, tori)
    return sub_map

def test_convertion(tb,freq,pix,switch='tb2sfu',bmmin=None):
    constk=1.38e-23
    wl = 3.0e8/freq
    sa = (pix*4.848e-6)**2
    if switch == 'tb2sfu':
        #!!!!!!!!!!!!!!!!!!!!!!!!make it /2.     ------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        #return (2*constk*tb)/(wl**2)*sa*1.e22/2.
        return (2*constk*tb)/(wl**2)*sa*1.e22
    elif switch=='tb2jyb':
        beamarea = bmmin**2.0 * np.pi / (4. * np.log(2.))
        return (2*constk*tb)/(wl**2)*sa*1.e22 * 10000./(pix**2.0)*beamarea
    elif switch == 'sfu2tb':
        return tb/1.e22/sa*(wl**2.0)/constk/2.0

def lightcurve(timerange, x, y, kw_list, world=True, radio=False, sample_range=1, sample_method='mean',
    #def lightcurve(dic_file, timerange, x, y, kw_list, world=True, radio=False, sample_range=1, sample_method='mean',
               cur_color='r', normq=False, pax=None, ltex='',ret_data=False,dspec_radio_file=None,plot_lc=True,dspec_norm=None
               ,high_tres=False,derivative=False,cur_legend=None,blooming_cut_off=True,cmarker='-',displayed_trange=None):
    #in_dic = pickle.load(open(dic_file, 'rb'),encoding='latin1')
    if high_tres:
        in_dic = in_dic_1s
    else:
        in_dic = in_dic_10s
    aia_file_list = []
    for tmpt in range(timerange[0],timerange[1]):
        if radio:
            aia_file_list.append(in_dic[tmpt]['radio_sbs'][kw_list[0]])
        else:
            aia_file_list.append(in_dic[tmpt][kw_list[0]])
    final_aia_list = pd.Series(aia_file_list).drop_duplicates().tolist()
    blooming_list = []
    if blooming_cut_off:
        print('cur off the bloomings')
        blooming_save = '/Volumes/Data/20170820/info/blooming_list_{}.p'.format(kw_list[0].split('.')[1])
        try:
            with open(blooming_save, 'rb') as fff:
                blooming_list = pickle.load(fff,encoding='latin1')
            fff.close()
            print(blooming_list)
        except:
            print('no list saved. Check your input')
    new_tmp_list = []
    for cfaf in final_aia_list:
        if cfaf in blooming_list:
            #print(cfaf)
            print('catch ya !', cfaf)
            #final_aia_list.remove(cfaf)
        else:
            print('???????????',cfaf)
            new_tmp_list.append(cfaf)
    final_aia_list = new_tmp_list
    #print(final_aia_list)
    final_n_elements = len(final_aia_list)
    dates = np.zeros((final_n_elements))
    lc_arr = np.zeros((len(kw_list), final_n_elements))
    #dates = np.zeros((timerange[1] - timerange[0]))
    #lc_arr = np.zeros((len(kw_list), timerange[1] - timerange[0]))
    if pax == None:
        fig = plt.figure(figsize=(5, 5))
        pax = fig.add_subplot(111)
    if radio:
        cur_name = 'lc_eovsa_{0:0=3d}_{1:0=3d}.p'.format(int(abs(x)),int(abs(y)))
    else:
        cur_name='lc_aia_{0}_T_{1:0=3d}_{2:0=3d}_radii_{3:0=3d}_{4:0=3d}_{5}.p'.format(kw_list[0].replace('Z','').replace('.',''),int(abs(x)),int(abs(y)),int(sample_range[0]),int(sample_range[1]),sample_method)
    cur_save='/Volumes/Data/20170820/aia_lc/'+cur_name
    if os.path.exists(cur_save):
        lc_arr=pickle.load(open(cur_save,'rb'),encoding='latin1')
        #for it, t in enumerate(range(timerange[0], timerange[1])):
        for it in range(final_n_elements):
            #dates[it] = mpl.dates.date2num(parse_time(in_dic[t]['time']))
            if radio:
                dates[it] = mpl.dates.date2num(parse_time(tp.read_time(final_aia_list[it],instrument='EOVSA').iso))
            else:
                dates[it] = mpl.dates.date2num(parse_time(tp.read_time(final_aia_list[it]).iso))
    else:
        for it, cur_ff in enumerate(final_aia_list):
        #for it, t in enumerate(range(timerange[0], timerange[1])):
            if radio:
                dates[it] = mpl.dates.date2num(parse_time(tp.read_time(final_aia_list[it], instrument='EOVSA').iso))
            else:
                dates[it] = mpl.dates.date2num(parse_time(tp.read_time(final_aia_list[it]).iso))
            for ik, kw in enumerate(kw_list):
                if radio:
                    cur_map = smap.Map(in_dic_10s[it+timerange[0]]['radio_sbs'][kw])
                else:
                    cur_map = smap.Map(cur_ff)
                    print(cur_map.meta['EXPTIME'])
                    #tmp_jud = aia_exposure_helper(kw=kw, tim=cur_ff, high_tres=high_tres)
                    #if tmp_jud is not 1:
                    #    print('gocha')
                    #    cur_map = tmp_jud
                coor = SkyCoord(x * u.arcsec, y * u.arcsec, frame=cur_map.coordinate_frame)
                pix_coor = cur_map.world_to_pixel(coor)
                px = int(pix_coor.x.value)
                py = int(pix_coor.y.value)
                cur_data = sf.mapdata(smap=cur_map)
                # cur_data=cur_map.data
                if radio:
                    act_data = cur_data[px - sample_range[0] + 1:px + sample_range[0], py - sample_range[1] + 1:py + sample_range[1]]
                else:
                    act_data = cur_data[px - sample_range[0] + 1:px + sample_range[0],
                               #py - sample_range + 1:py + sample_range] / (cur_map.meta['exptime'] / 2)
                                py - sample_range[1] + 1:py + sample_range[1]]
                #print('using {}'.format(cur_ff))
                print('exp time is {}'.format(cur_map.meta['EXPTIME']))
                #print('lc is {}'.format(np.nanmean(act_data)))
                #fig, axdss = plt.subplots()
                #axdss.imshow(act_data,origin='lower')
                #plt.show()
                if sample_method == 'max':
                    lc_arr[ik, it] = np.nanmax(act_data)
                elif sample_method == 'min':
                    lc_arr[ik, it] = np.nanmin(act_data)
                else:
                    lc_arr[ik, it] = np.nanmean(act_data)
        pickle.dump(lc_arr,open(cur_save,'wb'))
    print(lc_arr)
    if radio:
        cur_cmap = plt.cm.Spectral
        cnorm = mpl.colors.Normalize(vmin=0, vmax=30)
        if dspec_radio_file is not None:
            with open(dspec_radio_file, 'rb') as rdspec:
                rds = pickle.load(rdspec,encoding='latin1')
            rdspec.close()
            rds_data = rds['spec'][0, 0, :, timerange[0]: timerange[1]]
    for ik, kw in enumerate(kw_list):
        if normq:
            cur_od = lc_arr[ik, :] / np.nanmax(lc_arr[ik, :])
            if derivative:
                cur_od = np.gradient(cur_od)
            # cur_ax.plot_date(dates, cur_od, '-',linestyle='None',color=cur_color)
            if radio:
                if not plot_lc: continue
                if dspec_radio_file is None:
                    cim = pax.plot_date(dates, cur_od, cmarker, color=cur_cmap(cnorm(ik)), label=ltex)
                else:
                    cur_od = rds_data[ik,:]
                    cim = pax.plot_date(dates, cur_od, cmarker, color=cur_cmap(cnorm(ik)), label=ltex)
            else:
                cim = pax.plot_date(dates, cur_od, cmarker, color=cur_color, label=ltex)
        else:
            if derivative:
                lc_arr[ik,:] = np.gradient(lc_arr[ik,:])
            # cur_ax.plot_date(dates, lc_arr[ik,:], '-', label='', color='red', lw=2, linestyle='None')
            # pax.plot_date(dates, lc_arr[ik,:], '-',linestyle='None',color=cur_color)
            if radio:
                if not plot_lc: continue
                if dspec_radio_file is None:
                    cim = pax.plot_date(dates, lc_arr[ik, :], cmarker, color=cur_cmap(cnorm(ik)), label=ltex)
                else:
                    cur_od = rds_data[ik, :]
                    cim = pax.plot_date(dates, cur_od, cmarker, color=cur_cmap(cnorm(ik)), label=ltex)
            else:
                cim = pax.plot_date(dates, lc_arr[ik, :], cmarker, color=cur_color, label=ltex)
        # pax.set_ylim([0,1.1])
        if displayed_trange is not None:
            tim_plt0 = Time(in_dic_10s[displayed_trange[0]]['time'], format='iso').plot_date
            tim_plt1 = Time(in_dic_10s[displayed_trange[1]]['time'], format='iso').plot_date
            pax.set_xlim([tim_plt0,tim_plt1])
        pax.set_xlabel('UT Time')
        formatter = mpl.dates.DateFormatter('%H:%M:%S')
        pax.xaxis.set_major_formatter(formatter)
        pax.fmt_xdata = mpl.dates.DateFormatter('%H:%M:%S')
    if not plot_lc and radio:
        print('plot dspec instead of lc')
        if dspec_norm is None:
            #dspec_norm=[np.nanmin(lc_arr),np.nanmax(lc_arr)]
            print(np.nanmin(lc_arr),np.nanmax(lc_arr))
            dspec_norm=[1.e4,1.e7]
            #dlog_norm = SymLogNorm(linthresh=0.0,vmin=np.nanmin(lc_arr),vmax=np.nanmax(lc_arr))
            dlog_norm = colors.LogNorm(vmin=dspec_norm[0], vmax = dspec_norm[1])
        #cim = pax.pcolormesh(dates, cfreq, lc_arr, cmap='jet', vmin=dspec_norm[0], vmax=dspec_norm[1])
        cim = pax.pcolormesh(dates, cfreq, lc_arr, cmap='jet',norm=dlog_norm)
        pax.set_xlabel('UT Time')
        pax.set_ylabel('Freq [GHz]')
        formatter = mpl.dates.DateFormatter('%H:%M:%S')
        pax.xaxis.set_major_formatter(formatter)
        pax.fmt_xdata = mpl.dates.DateFormatter('%H:%M:%S')
        #pcb = plt.colorbar(cim, ax=pax)
    if plot_lc:
        #pax.set_yscale("log")
        #pax.set_ylabel('${T_b}$ [k]')
        if normq:
            pax.set_ylabel('Normalized Intensity')
        else:
            pax.set_ylabel('Intensity')
        if derivative:
            pax.set_ylabel(pax.get_ylabel()+'\nderivative')
    pax.xaxis.set_tick_params(rotation=30)
    if cur_legend is not None:
       ep.add_patch_to_existing_legend(legend=cur_legend, label=ltex, marker='_')
    if pax is None:
        plt.show()
    if ret_data:
        if normq:
            return (dates, cur_od)
        else:
            return (dates, lc_arr)
    return cim

def aia_max(aia_exp_file):
    aia_exp = pickle.load(open(aia_exp_file, 'rb'))
    aia_max_dic = {}
    for cur_key in aia_exp.keys():
        aia_max_dic[cur_key] = 0.0
        for i, cfile in enumerate(aia_exp[cur_key]['file_list']):
            print(i)
            print(cfile)
            cxdata=np.nanmax(smap.Map(cfile).data)
            if cxdata > aia_max_dic[cur_key]:
                aia_max_dic[cur_key] = cxdata
    print(aia_max_dic)
    pickle.dump(aia_max_dic, open('/Volumes/Data/20170820/20220511/info/aia_max.p','wb'))
    return aia_max_dic


def aia_exp(dic_file):
    in_dic = pickle.load(open(dic_file, 'rb'))
    aia_max_dic = {}
    for cur_key in in_dic[0].keys():
        aia_max_dic[cur_key] = np.zeros(len(in_dic))
    for i, dic in enumerate(in_dic):
        print(i)
        for ckey in dic.keys():
            if 'Z' in ckey:
                # if np.nanmax(smap.Map(dic[ckey]).data) > aia_max_dic[ckey]:
                #    aia_max_dic[ckey] = np.nanmax(smap.Map(dic[ckey]).data)
                aia_max_dic[ckey][i]=smap.Map(dic[ckey]).meta['exptime']
    pickle.dump(aia_max_dic, open('/Volumes/Data/20170820/20220511/info/aia_exp.p','wb'))
    return aia_max_dic

