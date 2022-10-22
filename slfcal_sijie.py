import os
import pickle
from glob import glob

import matplotlib.colors as colors
import numpy as np
from matplotlib import gridspec as gridspec
from matplotlib import pyplot as plt
from suncasa import dspec as ds
from suncasa.utils import helioimage2fits as hf
from suncasa.utils import mod_slftbs as mods
from suncasa.utils import qlookplot as ql
from sunpy import map as smap

try:
    input = raw_input
except NameError:
    pass

'''
Example script for self-calibrating EOVSA flare data
'''
# History:
#   2019-May-15 BC
#       Created a new example script based on B. Chen's practice for self-calibrating
#       the 2017 Aug 21 20:20 UT flare data. Made it available for EOVSA tutorial at
#       RHESSI XVIII Workshop (http://rhessi18.umn.edu/)

#   2022-Feb Sijie Yu
#       Fix some bugs
#       Make the script available for 50 spws

# =========== task handlers =============
dodspec = 0
dofullsun = 0  # initial full-sun imaging
doqlookplot = 0
domasks = 0  # get masks
doslfcal = 0  # main cycle of doing selfcalibration
doapply = 0  # apply the results
doclean_slfcaled = 1  # perform clean for self-calibrated data

# ============ declaring the working directories ============
# workdir = os.getcwd()+'/' #main working directory. Using current directory in this example
workdir = '/Volumes/WD6T/working/20220511/slfcal_sijie/'  # main working directory.
slfcaldir = os.path.join(workdir, 'slfcal/')  # place to put all selfcalibration products
imagedir = os.path.join(slfcaldir, 'images/')  # place to put all selfcalibration images
maskdir = os.path.join(slfcaldir, 'masks/')  # place to put clean masks
imagedir_slfcaled = os.path.join(slfcaldir, 'images_slfcaled/')  # place to put final self-calibrated images
caltbdir = os.path.join(slfcaldir, 'caltbs/')  # place to put calibration tables
# make these directories if they do not already exist
dirs = [workdir, slfcaldir, imagedir, maskdir, imagedir_slfcaled, caltbdir]
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)

os.chdir(workdir)
# ============ Split a short time for self-calibration ===========
# input visibility
ms_in = os.path.join(workdir, 'msdata/IDB20220511_1830-1850.cal.ms')
#ms_in = os.path.join(workdir, 'msdata/IDB20220511_1830-1850.ms')
# clearcal(ms_in)
# delmod(ms_in)
if dodspec:
    specfile = ms_in + '.dspec.npz'
    if os.path.exists(specfile):
        d = ds.Dspec()
        d.read(specfile)
    else:
        d = ds.Dspec(ms_in, specfile=specfile, usetbtool=False, domedian=True)
    d.plot(vmin=None, norm=colors.LogNorm(vmax=200, vmin=0.5), pol='XX')
    print(d.data.shape)  # npol, nbl, nfreq, ntime
# output, selfcaled, visibility
ms_slfcaled = os.path.join(workdir, os.path.basename(ms_in).replace('.ms', '.ms.slfcaled'))
# intermediate small visibility for selfcalbration
# selected time range for generating self-calibration solutions

# trange_peak = '2022/03/30/17:31:07.5~2022/03/30/17:31:08.5' # t0
# trange_peak = '2022/03/30/17:50:25.5~2022/03/30/17:50:30.5'  # t1

# trange = trange_peak
tranges = ['2022/05/11/18:42:12~2022/05/11/18:42:16']
subms_a = []
slfcalms_a = []
slfcaledms_a = []
tr_strs = []
for t, trange in enumerate(tranges):
    slfcalms_ = slfcaldir + 'slfcalms.t{0:d}.XX.slfcal'.format(t)
    slfcaledms_ = slfcaldir + 'slfcalms.t{0:d}.XX.slfcaled'.format(t)
    if not os.path.exists(slfcalms_):
        print('no files found, splitting')
        split(vis=ms_in, outputvis=slfcalms_, datacolumn='data', timerange=trange, correlation='XX', width=30)
        split(vis=ms_in, outputvis=slfcaledms_, datacolumn='data', timerange=trange, correlation='XX', width=30)
    slfcalms_a.append(slfcalms_)
    slfcaledms_a.append(slfcaledms_)

# ============ Prior definitions for spectral windows, antennas, pixel numbers =========
spws = [str(s) for s in range(0, 50)]
antennas = '0~12'
npix = 512
nround = 3  # number of slfcal cycles
# parameters specific to the event (found from step 1)
phasecenter = 'J2000 0.8382380452843443rad 0.30992505094358624rad'
xycen = [940, -260]
xran = xycen[0]+np.array([-1,1])*150
yran = xycen[1]+np.array([-1,1])*150

spwrans_mask = ['0~1', '2~3', '4~6', '7~10', '11~19', '20~29', '30~39', '40~49']
# spwrans_mask = ['3~6', '7~10', '11~19', '20~29', '30~39', '40~49']
# convert to a list of spws
spwrans_mask_list = []
for m in spwrans_mask:
    if '~' in m:
        spwrans_mask_list.append([str(i) for i in (np.arange(int(m.split('~')[0]), int(m.split('~')[1]) + 1))])
    else:
        spwrans_mask_list.append([m])
# print(spwrans_mask_list)
# spwrans_mask_list = [[str(i) for i in (np.arange(int(m.split('~')[0]),int(m.split('~')[1])))] for m in spwrans_mask]

tb.open(slfcalms_a[0] + '/SPECTRAL_WINDOW')
reffreqs = tb.getcol('REF_FREQUENCY')
bdwds = tb.getcol('TOTAL_BANDWIDTH')
cfreqs = reffreqs + bdwds / 2.
tb.close()
# starting beam size at 3.4 GHz in arcsec
sbeam = 40.
# setting restoring beam size (not very useful for selfcal anyway, but just to see the results)
bms = sbeam * cfreqs[5] / cfreqs
bms[bms < 6.] = 6.

# =========== Step 1, doing a full-Sun image to find out phasecenter and appropriate field of view =========
if dofullsun:
    # initial mfs clean to find out the image phase center
    im_init = 'fullsun_init'
    os.system('rm -rf ' + im_init + '*')
    tclean(vis=slfcalms_a[0],
           antenna='0~12',
           imagename=im_init,
           spw='2',
           specmode='mfs',
           timerange=tranges[0],
           imsize=[npix],
           cell=['5arcsec'],
           niter=1000,
           gain=0.05,
           stokes='XX',
           restoringbeam=['50arcsec'],
           interactive=False,
           pbcor=True)

    hf.imreg(vis=slfcalms_a[0], imagefile=im_init + '.image.pbcor', fitsfile=im_init + '.fits',
             timerange=trange, usephacenter=False, verbose=True)
    clnjunks = ['.flux', '.mask', '.model', '.psf', '.residual', '.sumwt', '.pb', '.image']
    for clnjunk in clnjunks:
        if os.path.exists(im_init + clnjunk):
            os.system('rm -rf ' + im_init + clnjunk)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    eomap = smap.Map(im_init + '.fits')
    # eomap.data=eomap.data.reshape((npix,npix))
    eomap.plot_settings['cmap'] = plt.get_cmap('jet')
    eomap.plot(axes=ax)
    eomap.draw_limb()
    plt.show()
    imview(im_init + '.image.pbcor')

if doqlookplot:
    specfile = ms_in + '.dspec.npz'
    ql.qlookplot(vis=ms_in, plotaia=True,
                 specfile=specfile,
                 overwrite=True,
                 pbcor=True, toTb=True,
                 dnorm=colors.LogNorm(vmax=2000, vmin=0.5),
                 aiafits='/Volumes/WD6T/working/20220511/aia/aia.lev1_euv_12s.2022-05-11T184221Z.171.image_lev1.fits',
                 # spw=spw,
                 spw='2',
                 # spw='30~39',
                 # spw='40~49',
                 # xycen=[940, -260],
                 xycen=[0, 0],
                 imsize=[512, 512], cell=['5.0arcsec'],
                 timerange='2022/05/11/18:42:12~2022/05/11/18:42:16',
                 usemsphacenter=False,
                 # aiafits='',
                 # aiawave='131',
                 # amin=50,
                 # robust = -0.5,
                 #restoringbeam=['8arcsec'],  # spw='30~39'
                 # restoringbeam=['5arcsec'], #spw='40~49'
                 #restoringbeam=['30arcsec'], # spw='6~10'
                 #restoringbeam=['70arcsec'], # spw='0~1'
                 restoringbeam=['45arcsec'], # spw='2~5'
                 icmap='jet',
                 docompress=True, stokes='XX')

# =========== Step 2 (optional), generate masks =========
# if skipped, will not use any masks
masks_a = []
if domasks:
    for t, trange in enumerate(tranges):
        slfcalms = slfcalms_a[t]
        clearcal(slfcalms)
        delmod(slfcalms)
        antennas = '0~12'
        pol = 'XX'
        imgprefix = maskdir + 'slf_t{0:d}'.format(t)

        # step 1: set up the clean masks
        img_init = imgprefix + '_init_ar_'
        os.system('rm -rf ' + img_init + '*')
        masks = []
        imnames = []
        for spwran in spwrans_mask:
            imname = img_init + spwran.replace('~', '-')
            mname = imname + '.mask'
            if os.path.exists(mname):
                masks.append(mname)
            else:
                spwran_idx = np.int(np.nanmean(np.array(spwran.split('~')).astype(np.float)))
                try:
                    tclean(vis=slfcalms,
                           antenna='0~12',
                           imagename=imname,
                           spw=spwran,
                           specmode='mfs',
                           timerange=trange,
                           imsize=[512],
                           # cell=['2arcsec'],
                           cell=['5arcsec'],
                           niter=1000,
                           gain=0.05,
                           stokes='XX',
                           restoringbeam='{:.1f}arcsec'.format(bms[spwran_idx]),
                           phasecenter='',
                           # phasecenter=phasecenter,
                           weighting='briggs',
                           robust=1.0,
                           interactive=True,
                           datacolumn='data',
                           pbcor=True,
                           savemodel='modelcolumn')
                    imnames.append(imname + '.image')
                    masks.append(mname)
                    clnjunks = ['.flux', '.psf', 'pb', '.residual']
                    for clnjunk in clnjunks:
                        if os.path.exists(imname + clnjunk):
                            os.system('rm -rf ' + imname + clnjunk)
                except:
                    print('error in cleaning spw: ' + spwran)
        masks_a.append(masks)

    pickle.dump(masks_a, open(slfcaldir + 'masks_a.p', 'wb'))

# =========== Step 3, main step of selfcalibration =========
if doslfcal:
    maskfile = os.path.join(slfcaldir, 'masks_a.p')
    if os.path.exists(maskfile):
        masks_a = pickle.load(open(maskfile, 'rb'))
    else:
        print('masks do not exist. Use default mask')
        masks_a = []
    os.system('rm -rf ' + imagedir + '*')
    os.system('rm -rf ' + caltbdir + '*')
    # first step: make a mock caltable for the entire database
    slftbs_a = []
    for t, trange in enumerate(tranges):
        print('Processing ' + str(t + 1) + ' in ' + str(len(tranges)) + ' times: ' + trange)
        slftbs = []
        slfcalms = slfcalms_a[t]
        slfcaledms = slfcaledms_a[t]
        if masks_a:
            masks = masks_a[t]
        else:
            masks = []
        calprefix = caltbdir + 'slf_t{0:d}'.format(t)
        imgprefix = imagedir + 'slf_t{0:d}'.format(t)

        strtmp = [m.replace(':', '') for m in trange.split('~')]
        timestr = 't' + strtmp[0] + '-' + strtmp[1]
        refantenna = '0'
        # number of iterations for each round
        niters = [100, 300, 500]
        # roburst value for weighting the baselines
        robusts = [1.0, 0.5, 0.0]
        # apply calibration tables? Set to true for most cases
        doapplycal = [1, 1, 1]
        # modes for calibration, 'p' for phase-only, 'a' for amplitude only, 'ap' for both
        calmodes = ['p', 'p', 'a']
        # setting uvranges for model image (optional, not used here)
        uvranges = ['', '', '']
        for n in range(nround):
            slfcal_tb_g = calprefix + '.G' + str(n)
            fig = plt.figure(figsize=(11, 7.))
            gs = gridspec.GridSpec(5, 10)
            for s, sp in enumerate(spws):
                print('processing spw: ' + sp)
                cfreq = cfreqs[int(sp)]

                slfcal_img = imgprefix + '.spw' + sp.zfill(2) + '.slfcal' + str(n)
                # only the first round uses nearby spws for getting initial model
                if n < 1:
                    ## v1
                    # sp_ = int(sp)
                    # if sp_ <= 2:
                    #     spbg = 0
                    #     sped = 2
                    # elif sp_ <= 10:
                    #     spbg = max(sp_ - 1, 2)
                    #     sped = sp_ + 1
                    # elif sp_ <= 30:
                    #     spbg = max(sp_ - 2, 2)
                    #     sped = sp_ + 2
                    # else:
                    #     # spbg = max(sp_ - 3, 28)
                    #     spbg = 28
                    #     sped = min(sp_ + 3, 49)

                    # ## v2
                    # sp_ = int(sp)
                    # if sp_ <= 1:
                    #     spbg = 0
                    #     sped = 1
                    # elif sp_ <= 10:
                    #     spbg = max(sp_ - 1, 2)
                    #     sped = sp_ + 1
                    # elif sp_ <= 30:
                    #     spbg = max(sp_ - 2, 2)
                    #     sped = sp_ + 2
                    # else:
                    #     spbg = max(sp_ - 3, 2)
                    #     sped = min(sp_ + 3, 49)

                    # ## v3
                    # sp_ = int(sp)
                    # if sp_ <= 1:
                    #     spbg = 0
                    #     sped = 1
                    # elif sp_ <= 10:
                    #     spbg = max(sp_ - 1, 2)
                    #     sped = sp_ + 1
                    # elif sp_ <= 20:
                    #     spbg = max(sp_ - 2, 2)
                    #     sped = sp_ + 2
                    # elif sp_ <= 30:
                    #     spbg = max(sp_ - 3, 2)
                    #     sped = sp_ + 3
                    # else:
                    #     spbg = max(sp_ - 4, 2)
                    #     sped = min(sp_ + 4, 49)

                    ## v4
                    sp_ = int(sp)
                    if sp_ <= 1:
                        spbg = 0
                        sped = 3
                    elif sp_ <= 10:
                        spbg = max(sp_ - 1, 2)
                        sped = sp_ + 1
                    elif sp_ <= 20:
                        spbg = max(sp_ - 2, 2)
                        sped = sp_ + 2
                    elif sp_ <= 30:
                        spbg = max(sp_ - 3, 2)
                        sped = sp_ + 3
                    else:
                        #spbg = 28
                        spbg = max(sp_ - 4, 2)
                        sped = min(sp_ + 4, 49)

                    spwran = str(spbg) + '~' + str(sped)
                    print('processling spw {0:s} using spw {1:s} as model'.format(sp, spwran))

                    if 'spwrans_mask_list' in vars():
                        for m, spwran_mask in enumerate(spwrans_mask_list):
                            if sp in spwran_mask:
                                mask = masks[m]
                                print('using mask {0:s}'.format(mask))
                                findmask = True
                        if not findmask:
                            mask = ''
                            print('mask not found. Do use any masks')
                else:
                    spwran = sp
                    if 'spwrans_mask_list' in vars():
                        for m, spwran_mask in enumerate(spwrans_mask_list):
                            if sp in spwran_mask:
                                mask = masks[m]
                                print('using mask {0:s}'.format(mask))
                                findmask = True
                        if not findmask:
                            mask = ''
                            print('mask not found. Do use any masks')
                # try:
                tclean(vis=slfcalms,
                       antenna=antennas,
                       imagename=slfcal_img,
                       uvrange=uvranges[n],
                       spw=spwran,
                       specmode='mfs',
                       timerange=trange,
                       imsize=[npix],
                       cell=['2.0arcsec'],
                       niter=niters[n],
                       gain=0.05,
                       stokes='XX',  # use pol XX image as the model
                       weighting='briggs',
                       robust=robusts[n],
                       phasecenter=phasecenter,
                       mask=mask,
                       restoringbeam=['{:.1f}arcsec'.format(bms[s])],
                       pbcor=False,
                       interactive=False,
                       savemodel='modelcolumn')
                if os.path.exists(slfcal_img + '.image'):
                    fitsfile = slfcal_img + '.fits'
                    hf.imreg(vis=slfcalms, imagefile=slfcal_img + '.image', fitsfile=fitsfile,
                             timerange=trange, usephacenter=False, toTb=True, verbose=False, overwrite=True)
                clnjunks = ['.mask', '.flux', '.model', '.psf', '.residual', '.image', '.pb', '.image.pbcor', '.sumwt']
                for clnjunk in clnjunks:
                    if os.path.exists(slfcal_img + clnjunk):
                        os.system('rm -rf ' + slfcal_img + clnjunk)
                ax = fig.add_subplot(gs[s])
                eomap = smap.Map(fitsfile)
                eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                eomap.plot(axes=ax)
                eomap.draw_limb()
                # eomap.draw_grid()
                ax.set_title(' ')
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                ax.set_xlim(xran)
                ax.set_ylim(yran)
                #os.system('rm -f ' + fitsfile)

                # except:
                # print('error in cleaning spw: '+sp)
                # print('using nearby spws for initial model')
                # sp_e=int(sp)+2
                # sp_i=int(sp)-2
                # if sp_i < 0:
                #   sp_i = 0
                # if sp_e > 49:
                #    sp_e = 49
                # sp_=str(sp_i)+'~'+str(sp_e)
                # try:
                #    tget(tclean)
                #    spw=sp_
                #    print('using spw {0:s} as model'.format(sp_))
                #    tclean()
                # except:
                #    print('still not successful. abort...')
                #    break

                gaincal(vis=slfcalms, refant=refantenna, antenna=antennas, caltable=slfcal_tb_g, spw=sp, uvrange='', \
                        gaintable=[], selectdata=True, timerange=trange, solint='inf', gaintype='G',
                        calmode=calmodes[n], \
                        combine='', minblperant=4, minsnr=2, append=True)
                if not os.path.exists(slfcal_tb_g):
                    print('No solution found in spw: ' + sp)
            figname = imagedir + 'slf_t{0:d}_n{1:d}.png'.format(t, n)
            plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
            plt.savefig(figname)
            # time.sleep(10)
            plt.close()

            if os.path.exists(slfcal_tb_g):
                slftbs.append(slfcal_tb_g)
                slftb = [slfcal_tb_g]
                # os.chdir(slfcaldir)
            #     if calmodes[n] == 'p':
            #         plotcal(caltable=slfcal_tb_g,antenna='1~12',xaxis='freq',yaxis='phase',\
            #                 subplot=431,plotrange=[-1,-1,-180,180],iteration='antenna',figfile=slfcal_tb_g+'.png',showgui=False)
            #     if calmodes[n] == 'a':
            #         plotcal(caltable=slfcal_tb_g,antenna='1~12',xaxis='freq',yaxis='amp',\
            #                 subplot=431,plotrange=[-1,-1,0,2.],iteration='antenna',figfile=slfcal_tb_g+'.png',showgui=False)
            #     os.chdir(workdir)

            if doapplycal[n]:
                clearcal(slfcalms)
                delmod(slfcalms)
                applycal(vis=slfcalms, gaintable=slftb, spw=','.join(spws), selectdata=True, \
                         antenna=antennas, interp='nearest', flagbackup=False, applymode='calonly', calwt=False)

            if n < nround - 1:
                # prompt=input('Continuing to selfcal?')
                prompt = 'y'
                if prompt.lower() == 'n':
                    if os.path.exists(slfcaledms):
                        os.system('rm -rf ' + slfcaledms)
                    split(slfcalms, slfcaledms, datacolumn='corrected')
                    print('Final calibrated ms is {0:s}'.format(slfcaledms))
                    break
                if prompt.lower() == 'y':
                    slfcalms_ = slfcalms + str(n)
                    if os.path.exists(slfcalms_):
                        os.system('rm -rf ' + slfcalms_)
                    split(slfcalms, slfcalms_, datacolumn='corrected')
                    slfcalms = slfcalms_
            else:
                if os.path.exists(slfcaledms):
                    os.system('rm -rf ' + slfcaledms)
                split(slfcalms, slfcaledms, datacolumn='corrected')
                print('Final calibrated ms is {0:s}'.format(slfcaledms))
        slftbs_a.append(slftbs)

# =========== Step 4: Apply self-calibration tables =========
if doapply:
    os.chdir(workdir)
    clearcal(ms_in)

    slftbs_comb = []
    for n in range(nround):
        tbin = glob(caltbdir + 'slf_t*.G{0}'.format(n))
        mods.concat(tb_in=tbin, tb_out=caltbdir + 'slf_comb.G{0}'.format(n))
        slftbs_comb.append(caltbdir + 'slf_comb.G{0}'.format(n))

    mod_slftbs = False
    ## note, this is for stokes XX only
    if mod_slftbs:
        amp_bound = [0.4, 2.5]
        for slftb2bmod in slftbs_comb[3:]:
            tb.open(slftb2bmod, nomodify=False)
            # tb.open(slftb2bmod)
            cparam = tb.getcol('CPARAM')
            flag = tb.getcol('FLAG')
            for row in range(tb.nrows()):
                cpar = abs(cparam[0, 0, row])
                if cpar > amp_bound[1] or cpar < amp_bound[0]:
                    # print cpar, flag[0,0,row]
                    print(slftb2bmod, cpar)
                    flag[0, 0, row] = True
            tb.putcol('FLAG', flag)
            tb.close()

    # clearcal(slfcalms)
    # applycal(vis=slfcalms, gaintable=slftbs, spw=','.join(spws), selectdata=True, \
    #          antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)
    applycal(vis=ms_in, gaintable=slftbs_comb, spw=','.join(spws), selectdata=True, \
             antenna=antennas, interp='linear', flagbackup=False, applymode='calonly', calwt=False)
    if os.path.exists(ms_slfcaled):
        os.system('rm -rf ' + ms_slfcaled)
    split(ms_in, ms_slfcaled, datacolumn='corrected')
    # split(ms_slfcaled, os.path.join(workdir, os.path.basename(ms_in).replace('cal.1s', 'slfcaled.30chan.12s')),
    #       datacolumn='data', width=30, timebin='12s')

    # ## IDB20220118_1719-1809XXYY.cal.1s.ms is used as the ms_in for subvs2.
    # outputvis = os.path.join(workdir, 'IDB20220118_1719-1809XXYY.slfcaled.subvs.1s.ms')
    # from suncasa.tasks.task_subvs2 import subvs2
    # subvs2(vis=ms_slfcaled, outputvis=outputvis, timerange='2022/03/30/17:22:00~2022/03/30/17:57:00', spw='',
    #        subtime1='2022/03/30/17:14:00~2022/03/30/17:19:00', subtime2='2022/03/30/17:59:00~2022/03/30/18:04:00')

# =========== Step 5: Generate final self-calibrated images (optional) =========
if doclean_slfcaled:
    pol = 'XX'
    vis = ms_slfcaled
    tb.open(vis + '/SPECTRAL_WINDOW')
    reffreqs = tb.getcol('REF_FREQUENCY')
    bdwds = tb.getcol('TOTAL_BANDWIDTH')
    cfreqs = reffreqs + bdwds / 2.
    tb.close()
    sbeam = 40.
    npix_final = 256
    for t, trange in enumerate(tranges):
        print('Processing ' + str(t + 1) + ' in ' + str(len(tranges)) + ' times: ' + trange)
        img_final = imagedir_slfcaled + '/slf_final_{0}_t{1:d}'.format(pol, t)
        fitsfiles = []
        for s, sp in enumerate(spws):
            cfreq = cfreqs[int(sp)]
            bm = max(sbeam * cfreqs[4] / cfreq, 6.)
            imname = img_final + '_s' + sp.zfill(2)
            fitsfile = imname + '.fits'
            if not os.path.exists(fitsfile):
                print('cleaning spw {0:s} with beam size {1:.1f}"'.format(sp, bm))
                try:
                    tclean(vis=vis,
                           antenna=antennas,
                           imagename=imname,
                           spw=sp,
                           specmode='mfs',
                           timerange=trange,
                           imsize=[npix_final],
                           cell=['2arcsec'],
                           niter=1000,
                           gain=0.05,
                           stokes=pol,
                           weighting='briggs',
                           robust=2.0,
                           restoringbeam=[str(bm) + 'arcsec'],
                           phasecenter=phasecenter,
                           mask='',
                           pbcor=True,
                           interactive=False)
                except:
                    print('cleaning spw ' + sp + ' unsuccessful. Proceed to next spw')
                    continue
                if os.path.exists(imname + '.image.pbcor'):
                    imn = imname + '.image.pbcor'
                    hf.imreg(vis=vis, imagefile=imn, fitsfile=fitsfile,
                             timerange=trange, usephacenter=False, toTb=True, docompress=True, verbose=False)
                fitsfiles.append(fitsfile)
                junks = ['.flux', '.model', '.psf', '.residual', '.mask', '.image', '.pb', '.image.pbcor', '.sumwt']
                for junk in junks:
                    if os.path.exists(imname + junk):
                        os.system('rm -rf ' + imname + junk)
            else:
                print('fits file ' + fitsfile + ' already exists, skip clean...')
                fitsfiles.append(fitsfile)

        fig = plt.figure(figsize=(11, 7.))
        gs = gridspec.GridSpec(5, 10)
        for s, sp in enumerate(spws):
            cfreq = cfreqs[int(sp)]
            ax = fig.add_subplot(gs[s])
            eomap = smap.Map(fitsfiles[s])
            eomap.plot_settings['cmap'] = plt.get_cmap('jet')
            eomap.plot(axes=ax)
            eomap.draw_limb()
            ax.set_title(' ')
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            plt.text(0.98, 0.85, '{0:.1f} GHz'.format(cfreq / 1e9), transform=ax.transAxes, ha='right', color='w',
                     fontweight='bold')
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        plt.show()
        plt.savefig(img_final + '.png')
        for l in fitsfiles:
            os.system('rm -rf {}'.format(l))
