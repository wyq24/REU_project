from astropy.time import Time
import os
import pickle
from astropy.time import TimeDelta
import numpy as np
from astropy.io import fits
def make_time_dic_10s(workdir=None, radio_dir=None, kw2=None, start_timeindex=None, end_timeindex=None,longer_sdo_duration=True,radio_start_time=None,radio_end_time=None):
    keyword1_list = ['Z.94.', 'Z.131.', 'Z.171.', 'Z.193.', 'Z.211.', \
                     'Z.304.', 'Z.335.', 'Z.1600.', 'Z.1700.']
     #'Z.304.', 'Z.335.', 'bbso', 'Z.1700.', 'Z.1600.', 'hmi']
    file_list = []
    dic_list = []
    #if SDO coverage is longer than EOVSA:
    if longer_sdo_duration:
        #init_time = '2017-07-03T16:12:40.002'
        #init_time = '2017-07-03T16:00:00.002'
        #init_time = '2017-07-03T15:56:30.002'
        init_time = '2022-05-11T18:30:00.002'
        time_interval = 10.0
        tranges = []
        tname = []
        # for tint in range(298):
        # for tint in range(298):
        for tint in range(end_timeindex):
            init_t = Time(init_time, format='isot')
            timed1 = TimeDelta((tint * time_interval) * 1.0, format='sec')
            start_time = init_t + timed1
            tname.append(start_time.jd)
            '''
            timed2 = TimeDelta(time_interval, format='sec')
            start_time = init_t + timed1
            end_time = start_time + timed2
            tmp_trange = '{0}~{1}'.format(start_time.iso.replace('-', '/').replace(' ', '/'),
                                          end_time.iso.replace('-', '/').replace(' ', '/'))
            tranges.append(tmp_trange)
            tname.append(tmp_trange.split(':')[3] + tmp_trange.split(':')[4].split('.')[0])
            '''
    # #make ha file lists----------
    # hfile_list=[]
    # ha_kw_list = ['195eu','304eu']
    # for hkwi, hkw in enumerate(ha_kw_list):
    #     source_dir = '/Volumes/Data/20170820/euvi/'
    #     cur_list = makelist(tdir=source_dir,keyword1=hkw,keyword2='R.fts')
    #     hfile_list.append(cur_list)
    # print('ha list is:', hfile_list)
    # #make da file lists----------
    # dafile_list=[]
    # da_kw_list = ['det_131','det_193','det_211']
    # for dkwi, dkw in enumerate(da_kw_list):
    #     source_dir = '/Volumes/Data/20170820/aia_detrending/{}_non_blooming/'.format(dkw[4:])
    #     cur_list = makelist(tdir=source_dir,keyword1='detrending',keyword2='.fits')
    #     dafile_list.append(cur_list)
    # print('da list is:', dafile_list)
    #------make aia file list --------------------
    for kindex, kw1 in enumerate(keyword1_list):
        cur_list = makelist(tdir=workdir, keyword1=kw1, keyword2=kw2)
        file_list.append(cur_list)

    #-------------------aia time list---------------------
    print('make image list')
    time_list = []
    for kindex, flist in enumerate(file_list):
        cur_number = len(flist)
        jd_list = np.zeros((cur_number))
        for cn in range(cur_number):
            if 'aia.lev1' in flist[0]:
                #print(flist[cn])
                jd_list[cn] = read_time(fits_file=flist[cn], instrument='aia').jd
            elif 'hmi' in flist[0]:
                jd_list[cn] = read_time(fits_file=flist[cn], instrument='hmi').jd
            else:
                #print(flist[cn])
                jd_list[cn] = read_time(fits_file=flist[cn]).jd
            # jd_list[cn] = Time(smap.Map(flist[cn]).date).jd
        time_list.append(jd_list)
    #---------------------gst ha time list--------------------
    # print('make ha lists')
    # htime_list = []
    # for hkindex, hflist in enumerate(hfile_list):
    #     hcur_number = len(hflist)
    #     hjd_list = np.zeros((hcur_number))
    #     for hcn in range(hcur_number):
    #         if 'aia.lev1' in hflist[0]:
    #             hjd_list[hcn] = read_time(fits_file=hflist[hcn], instrument='aia').jd
    #         elif 'hmi' in hflist[0]:
    #             hjd_list[hcn] = read_time(fits_file=hflist[hcn], instrument='hmi').jd
    #         elif 'eu_R' in hflist[0]:
    #             hjd_list[hcn] = read_time(fits_file=hflist[hcn], instrument='stereo').jd
    #         else:
    #             #print(hflist[hcn])
    #             hjd_list[hcn] = read_time(fits_file=hflist[hcn]).jd
    #         # jd_list[cn] = Time(smap.Map(flist[cn]).date).jd
    #     htime_list.append(hjd_list)
    # #------make detrending_aia image list---------------
    # print('make detrending da lists')
    # datime_list = []
    # for dakindex, daflist in enumerate(dafile_list):
    #     dacur_number = len(daflist)
    #     dajd_list = np.zeros((dacur_number))
    #     for dacn in range(dacur_number):
    #         if 'aia.lev1' in daflist[0]:
    #             dajd_list[dacn] = read_time(fits_file=daflist[dacn], instrument='aia').jd
    #         elif 'hmi' in daflist[0]:
    #             dajd_list[dacn] = read_time(fits_file=daflist[dacn], instrument='hmi').jd
    #         elif 'eu_R' in daflist[0]:
    #             dajd_list[dacn] = read_time(fits_file=daflist[dacn], instrument='stereo').jd
    #         else:
    #             print(daflist[dacn])
    #             dajd_list[dacn] = read_time(fits_file=daflist[dacn]).jd
    #         # jd_list[cn] = Time(smap.Map(flist[cn]).date).jd
    #     datime_list.append(dajd_list)
    print('make radio list')
    radio_time_list=[]
    radio_sequence=[]
    for i in range(radio_start_time,radio_end_time):
        # radio_list = makelist(tdir=radio_dir, keyword1='slf_final_XX_t' + str(i) + '_T', keyword2=kw2)
        if i-radio_start_time < 2: continue
        radio_list = []
        for spwi in range(50):
            #print(i, spwi)
            radio_list.append(makelist(tdir=radio_dir, keyword1='slfcaled_tb_final_XX_10s_XX_t{0}_'.format(i-radio_start_time),
                                       keyword2='s{0:0=2d}.fits'.format(spwi + 1))[0])
            # radio_list.append(makelist(tdir=radio_dir, keyword1='image_subed_t{0:0=3d}_'.format(i), keyword2='_s{0:0=2d}.fits'.format(spwi+1))[0])
            # radio_list.append(makelist(tdir=radio_dir, keyword1='_t{0:0=3d}_'.format(i), keyword2='_s{0:0=2d}.fits'.format(spwi+1))[0])
            # radio_list.append(makelist(tdir=radio_dir, keyword1='_t{0:0=3d}_'.format(i), keyword2='_s{0:0=2d}.fits'.format(spwi+1))[0])
        # print('_t{0:0=3d}_s'.format(i),'while end at ',end_timeindex)
        ref_fits = radio_list[0]
        ref_time = read_time(fits_file=ref_fits, instrument='EOVSA').jd
        radio_time_list.append(ref_time)
        radio_sequence.append(radio_list)
        # if i < 20:
        #    timed = TimeDelta(700, format='sec')
        #    print(type(ref_time))
        #    new_ref_time=Time(ref_time,format='jd') + timed
        #    ref_time=new_ref_time.jd
    for cjdi,cur_jd in enumerate(tname):
        time_dic = {}
        time_dic['time'] = str(Time(cur_jd, format='jd').iso)
        for tindex, tlist in enumerate(time_list):
            diff_list = []
            for itime in tlist:
                diff_list.append(abs(itime - cur_jd))
            #print(diff_list)
            cur_min = diff_list.index(min(diff_list))
            time_dic[keyword1_list[tindex]] = file_list[tindex][cur_min]
        # for htindex, htlist in enumerate(htime_list):
        #     hdiff_list = []
        #     for hitime in htlist:
        #         hdiff_list.append(abs(hitime - cur_jd))
        #     #print(hdiff_list)
        #     hcur_min = hdiff_list.index(min(hdiff_list))
        #     time_dic[ha_kw_list[htindex]] = hfile_list[htindex][hcur_min]
        # for datindex, datlist in enumerate(datime_list):
        #     dadiff_list = []
        #     for daitime in datlist:
        #         dadiff_list.append(abs(daitime - cur_jd))
        #     #print(dadiff_list)
        #     dacur_min = dadiff_list.index(min(dadiff_list))
        #     time_dic[da_kw_list[datindex]] = dafile_list[datindex][dacur_min]
        #     print(dafile_list[datindex][dacur_min])
        #if cjdi in range(radio_start_time,radio_end_time):
            #notice!!!! only for 20220511
        if cjdi in range(radio_start_time+2, radio_end_time):
            time_dic['radio_sbs'] = radio_sequence[cjdi-(radio_start_time+2)]
        else:
            time_dic['radio_sbs'] ='None'
        dic_list.append(time_dic)

    pickle.dump(dic_list, open('/Volumes/Data/20170820/20220511/info/20220511_10s_long_aia.p', 'wb'))
    # pickle.dump(dic_list, open('/home/walter/Downloads/From_NJIT/20191004_1s.p', 'wb'))
    return dic_list

def makelist(tdir='', keyword1='', keyword2='', exclude=None):
    li = []
    # for root, dirs, files in os.walk(tdir):
    root = os.getcwd()
    files = os.listdir(tdir)
    for file in files:
        if exclude is None:
            if keyword1 in file and keyword2 in file:
                li.append(os.path.join(tdir, file))
        else:
            if keyword1 in file and keyword2 in file and exclude not in file:
                li.append(os.path.join(tdir, file))

    return li

def read_time(fits_file=None, instrument=None):
    hdulist = fits.open(fits_file,mode='readonly')
    if instrument == 'aia':
        tmp_list = hdulist[1].header
        #result = tmp_list['T_OBS']
        result = tmp_list['DATE-OBS']
    elif instrument == 'hmi':
        tmp_list = hdulist[1].header
        result = tmp_list['DATE-OBS']
    elif instrument == 'EOVSA':
        tmp_list = hdulist[0].header
        result = tmp_list['DATE-OBS']
    elif instrument == 'rhessi':
        tmp_list = hdulist[0].header
        result = tmp_list['date_obs']
    elif instrument == 'stereo':
        tmp_list = hdulist[0].header
        result = tmp_list['DATE-OBS']
    else:
        tmp_list = hdulist[0].header
        #result = tmp_list['DATE-OBS']
        result = tmp_list['T_OBS']
    hdulist.close()
    return Time(result)