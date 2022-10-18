import sys
sys.path.append('/Users/walterwei/opt/miniconda3/envs/scasa_env/lib/python3.6/site-packages')
from suncasa.utils import helioimage2fits as hf
import os
import numpy as np
import pickle
import suncasa.utils.mod_slftbs as mods
from matplotlib import gridspec as gridspec
from sunpy import map as smap
from matplotlib import pyplot as plt
import time
from suncasa.utils import plot_map
import shutil
from astropy.time import Time,TimeDelta
# from casatools import table as tbtool
# from casatools import ms as mstool
# from casatools import quanta as qatool
# from casatools import image as iatool
#
# tb = tbtool()
# ms = mstool()
# qa = qatool()
# ia = iatool()
# from casatasks import tclean
# from casatasks import split
# from casatasks import gaincal
# from casatasks import applycal
# from casatasks import delmod
# from casatasks import clearcal
# from suncasa.suncasatasks import ptclean6 as ptclean

#from casatools import ms as mstool
#ebmsize=[40.89866667, 35.        , 30.52746815, 27.06907583, 24.30586371,
#       22.0676259 , 20.20685112, 18.63547995, 17.29086809, 16.12723449,
#       15.11034483, 14.21408712, 13.41819773, 12.70671085, 12.06687648,
#       11.48838951, 10.96283059, 10.48325359, 10.04387688, 10.        ,
#       10.        , 10.        , 10.        , 10.        , 10.        ,
#       10.        , 10.        , 10.        , 10.        , 10.        ]
#task handles
dofullsun=1
domasks=0
doslfcal=0
doapply=0
dofinalclean=0 #final clean of the slfcaled data for the selected slfcal time only
doclean_slfcaled=0 #final clean of all the slfcaled data

#workdir='/Volumes/Data/20170906/msdata/2021_new_calibration/'
#workdir='/Volumes/WD6T/working/eovsa_events/20170703/'
#workdir='/Volumes/Data/20170820/20220511/eovsa/eovsa_full/'
#workdir='/Volumes/Data/20170820/20220511/eovsa/1800_2000/data/'
#workdir='/Volumes/Data/20170820/20220511/eovsa/1800_2000/full_disk_slfcal/'
workdir='/Volumes/Data/20170820/20220511/eovsa/1800_2000/channels_test/'
#workdir='/Volumes/WD6T/working/eovsa_full_disk/'
imagedir=workdir+'slfcal/images/'
maskdir=workdir+'slfcal/masks/'
imagedir_final=workdir+'slfcal/images_final/'
imagedir_slfcaled=workdir+'slfcal/images_slfcaled/'
slfcaldir=workdir+'slfcal/'
caltbdir=slfcaldir+'caltbs/'
if not os.path.exists(imagedir):
    os.makedirs(imagedir)
if not os.path.exists(maskdir):
    os.makedirs(maskdir)
if not os.path.exists(imagedir_final):
    os.makedirs(imagedir_final)
if not os.path.exists(imagedir_slfcaled):
    os.makedirs(imagedir_slfcaled)
if not os.path.exists(caltbdir):
    os.makedirs(caltbdir)
#refms = workdir+'msdata/IDB20170906191320-20170906194320.corrected.ms'
#refms = workdir+'msdata/IDB20170703_concat.ms'
#refms = workdir+'msdata/IDB20210508_1800_1920.ms'
refms = workdir+'msdata/IDB20220511_1800-2000.ms'
refms_slfcal_XXYY = refms + '.XXYY.slfcal'
#refms_slfcal_XX = refms + '.XX.slfcal'
refms_slfcaled_XXYY = refms + '.XXYY.slfcaled'
#refms_slfcaled_XX = refms + '.XX.slfcaled'
refms_slfcal = refms + '.slfcal'
refms_slfcaled = refms + '.slfcaled'
if not os.path.exists(refms_slfcal):
    os.system('cp -r '+refms+' '+refms_slfcal)
if not os.path.exists(refms_slfcal_XXYY):
#if not os.path.exists(refms_slfcal_XX):
    split(vis=refms,outputvis=refms_slfcal_XXYY,datacolumn='data',correlation='XX,YY')
    #split(vis=refms, outputvis=refms_slfcal_XX, datacolumn='data', correlation='XX')
# selected times for generating self-calibration solutions
#tranges=['2017/09/09/21:04:56~2017/09/09/21:05:00']
#tranges=['2017/09/06/19:24:00~2017/09/06/19:24:10']
#tranges=['2017/07/15/19:32:10~2017/07/15/19:32:16']
#tranges=['2021/05/08/18:38:18~2021/05/08/18:38:24']
#tranges=['2022/05/11/18:42:13~2022/05/11/18:42:17']
tranges=['2022/05/11/18:42:12~2022/05/11/18:42:18']
#tranges=['2022/05/11/19:29:00~2022/05/11/19:29:30']
#tranges=['2021/05/08/18:43:25~2021/05/08/18:43:35']
start_time = Time('2022-05-11 18:42:13',format='iso')
end_time = Time('2022-05-11 18:42:17',format='iso')
#find the phasecenter
midtime_mjd = (start_time.mjd + end_time.mjd) / 2.
eph = hf.read_horizons(t0=Time(midtime_mjd, format='mjd'))
tb.open(refms + '/FIELD')
phadir = tb.getcol('PHASE_DIR').flatten()
tb.close()
#ms.open(refms)
#metadata = ms.metadata()
#observatory = metadata.observatorynames()[0]
#spwInfo = ms.getspectralwindowinfo()
#nspwall = len(spwInfo)
#print(len(cfreqs))
#ms.close()
ra0 = phadir[0]
dec0 = phadir[1]
#xycen=[912.0,-295.0]
#xycen=[920.0,-200.0]
xycen=[0.0,0.0]
#xycen=[912.0,-295.0]
x0 = np.radians(xycen[0] / 3600.)
y0 = np.radians(xycen[1] / 3600.)
p0 = np.radians(eph['p0'][0])  # p angle in radians
raoff = -((x0) * np.cos(p0) - y0 * np.sin(p0)) / np.cos(eph['dec'][0])
decoff = (x0) * np.sin(p0) + y0 * np.cos(p0)
newra = ra0 + raoff
newdec = dec0 + decoff
phasecenter = 'J2000 ' + str(newra) + 'rad ' + str(newdec) + 'rad'
#cphasecenter = 'J2000 ' + str(ra0) + 'rad ' + str(dec0) + 'rad'
print('Phasecenter: ',phasecenter)

tb.open(refms_slfcal + '/SPECTRAL_WINDOW')
reffreqs = tb.getcol('REF_FREQUENCY')
bdwds = tb.getcol('TOTAL_BANDWIDTH')
cfreqs = reffreqs + bdwds / 2.
tb.close()
print('cfreqs are: ', cfreqs)
#sbeam = 40.
#ebmsize=np.zeros_like(cfreqs)
#for iii in range(len(cfreqs)):
    #cfreq = cfreqs[iii]
    #ebmsize[iii] = max(sbeam * cfreqs[1] / cfreq, 6.)
    #cur_spw = str(ebmsize[iii])
#print('bmsizes are: ', ebmsize)
sbeam = 40.
bmsize = sbeam * 1.6e9 / np.asarray(cfreqs, dtype=float)
bmsize[bmsize < 4.0] = 4.0
print(bmsize.tolist())
#restoringbms = mstools.get_bmsize(cfreqs, refbmsize=sbeam, reffreq=1.6, minbmsize=4.0)
ebmsize = ['{:.1f}arcsec'.format(crbm) for crbm in bmsize]
print(ebmsize)
n_spw_per_mask = 5
# found from full-sun image in viewer
#phacenter='J2000 11h11m17 05d10m28'
#xran=[850,1050]
#yran=[-250,-50]
#phacenter='J2000 5h46m20 2d40m44'
#phacenter='J2000 1.8218694030439997rad 0.400123935501163rad'
#phacenter='J2000 00h01m04.733 0d00m42.734'
#phacenter='J2000 1.7918694030439997rad 0.400123935501163rad'
#phasecenter='J2000 1.7868374584594264rad 0.4002520172843901rad'
#phacenter='J2000 2.00406343122rad 0.373837961939rad'
#cphacenter='J2000 0.000rad 0.000rad'
#xran=[500,700]
#yran=[-350,-150]
#xran=[620,1220]
#yran=[100,-500]

#spws=[str(s+1) for s in range(len(cfreqs))]
spws=[str(s) for s in range(len(cfreqs))]
antennas='!4'
nround=3 #number of slfcal cycles

subms_a=[]
slfcalms_a=[]
slfcaledms_a=[]

for t,trange in enumerate(tranges):
    subms_ =      workdir+'slfcal/IDB20220511_1842_11_17.ms.t{0:d}.slfcal'.format(t)
    slfcalms_ =   workdir+'slfcal/IDB20220511_1842_11_17.ms.t{0:d}.XXYY.slfcal'.format(t)
    slfcaledms_ = workdir+'slfcal/IDB20220511_1842_11_17.ms.t{0:d}.XXYY.slfcaled'.format(t)
    if not os.path.exists(slfcalms_):
        split(vis=refms,outputvis=slfcalms_,datacolumn='data',timerange=trange,correlation='XX,YY')
    if not os.path.exists(subms_):
        split(vis=refms,outputvis=subms_,datacolumn='data',timerange=trange,correlation='')
    slfcalms_a.append(slfcalms_)
    subms_a.append(subms_)
    slfcaledms_a.append(slfcaledms_)
print('HERE WE ARE!')
# =========== Step 1, doing a full-Sun image to find out phasecenter and appropriate field of view =========
if dofullsun:
    #initial mfs clean to find out the image phase center
    im_init='fullsun_init_0_15'
    os.system('rm -rf '+im_init+'*')
    tclean(vis=slfcalms_a[0],
            antenna=antennas,
            imagename=im_init,
            spw='2~5',
            specmode='mfs',
            timerange=trange,
            imsize=[512],
            cell=['5arcsec'],
            niter=1000,
            gain=0.05,
            stokes='I',
            restoringbeam=['30arcsec'],
            interactive=False,
            pbcor=True)

    hf.imreg(vis=slfcalms_a[0],imagefile=im_init+'.image.pbcor',fitsfile=im_init+'.fits',
             timerange=trange,usephacenter=False,verbose=True)
    clnjunks = ['.flux', '.mask', '.model', '.psf', '.residual','.sumwt','.pb','.image']
    for clnjunk in clnjunks:
        if os.path.exists(im_init + clnjunk):
            os.system('rm -rf '+im_init + clnjunk)

    from sunpy import map as smap
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    eomap=smap.Map(im_init+'.fits')
    #eomap.data=eomap.data.reshape((npix,npix))
    eomap.plot_settings['cmap'] = plt.get_cmap('jet')
    eomap.plot(axes = ax)
    eomap.draw_limb()
    plt.show()
    viewer(im_init+'.image.pbcor')
masks_a=[]
if domasks:
    for t,trange in enumerate(tranges):
        slfcalms=slfcalms_a[t]
        clearcal(slfcalms)
        delmod(slfcalms)
        #subms = subms_a[t]
        #clearcal(subms)
        #delmod(subms)
        antennas=antennas,
        pol='XX'
        calprefix=caltbdir+'slf_t{0:d}'.format(t)
        imgprefix=maskdir+'slf_t{0:d}'.format(t)

        # step 1: set up the clean masks
        img_init=imgprefix+'_init_ar_'
        os.system('rm -rf '+img_init+'*')
        spwrans=[]
        bmrans=[]
        spwrans = ['0~1','2~4', '5~9', '10~14', '15~19', '20~24', '25~29', '30~34', '35~39', '40~44', '45~49']
        bmrans = ['50.7arcsec','25.3arcsec', '18.1arcsec', '12.4arcsec', '9.4arcsec', '7.6arcsec', '6.4arcsec', '5.5arcsec', '4.8arcsec', '4.3arcsec', '4.0arcsec']
        # for nmi in range(int(len(cfreqs)/n_spw_per_mask)):
        #     #spwrans.append('{0}~{1}'.format(nmi*n_spw_per_mask+1, (nmi+1)*n_spw_per_mask))
        #     spwrans.append('{0}~{1}'.format(nmi*n_spw_per_mask, (nmi+1)*n_spw_per_mask-1))
        #     #bmrans.append(str(ebmsize[nmi*n_spw_per_mask+1]))
        #     bmrans.append(str(ebmsize[nmi * n_spw_per_mask]))
        #spwrans[0] = '4'
        print('ranges are:' , spwrans, bmrans)
        #spwrans=['1~5','6~12','13~20','21~30']
        #bmrans=['41','19','12','10']
        mask=[]
        imnames=[]
        #for spwri,spwran in enumerate(spwrans):
            #imname=img_init+spwran.replace('~','-')
        #for cspw in spws:
        for spwri,spwran in enumerate(spwrans):
            #imname=img_init+cspw
            imname = img_init + spwran.replace('~', '-')
            try:
                tclean(vis=slfcalms,
                #tclean(vis=subms,
                        antenna='!4',
                        imagename=imname,
                        spw=spwran,
                        #spw=cspw,
                        specmode='mfs',
                        timerange='',
                        #imagermode='csclean',
                        #psfmode='clark',
                        #imsize=[256],
                        imsize=[512],
                        #cell=['1.2arcsec'],
                        cell=['5arcsec'],
                        niter=1000,
                        gain=0.05,
                        stokes='XX',
                        restoringbeam=['30arcsec'],
                        #restoringbeam=[bmrans[spwri]+'arcsec'],
                        #restoringbeam=[bmrans[spwri]],
                        #restoringbeam=['{}arcsec'.format(int(ebmsize[int(cspw)-int(1)]))],
                        #restoringbeam=['{}arcsec'.format(int(bmrans[spwri]))],
                        phasecenter=phasecenter,
                        weighting='briggs',
                        robust=1.0,
                        interactive=True,
                        datacolumn='data',
                        pbcor=True,
                        savemodel='modelcolumn')
                        #usescratch=False)
                imnames.append(imname+'.image')
                mask.append(imname+'.mask')
                clnjunks = ['.flux', '.model', '.psf', '.residual']
                for clnjunk in clnjunks:
                    if os.path.exists(imname + clnjunk):
                        shutil.rmtree(imname + clnjunk)
            except:
                print('error in cleaning spw: '+spwran)
                #print('error in cleaning spw: '+cspw)
        masks_a.append(mask)

    pickle.dump(masks_a,open(workdir+'slfcal/masks_a.p','wb'))
    #viewer(imnames)

    #os.system('rm -rf '+calprefix+'*')
    #os.system('rm -rf '+imgprefix+'.spw*')
if doslfcal:
    if not os.path.exists(workdir+'slfcal/masks_a.p') or not ('masks_a' in vars()):
        print ('masks do not exist. Please run dopartsun first!')
        masks_a=[]
    if os.path.exists(workdir+'slfcal/masks_a.p'):
        masks_a=pickle.load(open(workdir+'slfcal/masks_a.p','rb'))
    os.system('rm -rf '+imagedir+'*')
    os.system('rm -rf '+caltbdir+'*')
    #first step: make a mock caltable for the entire database
    slftbs_a=[]
    for t,trange in enumerate(tranges):
        slftbs=[]
        slfcalms=slfcalms_a[t]
        slfcaledms=slfcaledms_a[t]
        print('current time is ',t)
        subms=subms_a[t]
        masks=masks_a[t]
        calprefix=caltbdir+'slf_t{0:d}'.format(t)
        imgprefix=imagedir+'slf_t{0:d}'.format(t)
        #tb.open(slfcalms+'/SPECTRAL_WINDOW')
        #reffreqs=tb.getcol('REF_FREQUENCY')
        #bdwds=tb.getcol('TOTAL_BANDWIDTH')
        #cfreqs=reffreqs+bdwds/2.
        #tb.close()
        #sbeam=40.
        strtmp=[m.replace(':','') for m in trange.split('~')]
        timestr='t'+strtmp[0]+'-'+strtmp[1]
        refantenna='0'
        niters=[150,300,500]
        robusts=[1.0,0.5,0.0]
        doapplycal=[1,1,1]
        calmodes=['p','p','a']
        uvranges=['','','']
        spwrans=[]
        bmrans=[]
        spwrans = ['0~1','2~4', '5~9', '10~14', '15~19', '20~24', '25~29', '30~34', '35~39', '40~44', '45~49']
        bmrans = ['50.7arcsec','25.3arcsec', '18.1arcsec', '12.4arcsec', '9.4arcsec', '7.6arcsec', '6.4arcsec', '5.5arcsec', '4.8arcsec', '4.3arcsec', '4.0arcsec']
        # for nmi in range(int(len(cfreqs)/n_spw_per_mask)):
        #     #spwrans.append('{0}~{1}'.format(nmi*n_spw_per_mask+1, (nmi+1)*n_spw_per_mask))
        #     spwrans.append('{0}~{1}'.format(nmi*n_spw_per_mask, (nmi+1)*n_spw_per_mask-1))
        #     #bmrans.append(str(ebmsize[nmi*n_spw_per_mask+1]))
        #     bmrans.append(str(ebmsize[nmi * n_spw_per_mask]))
        #spwrans[0] = '2~5'
        #spwrans=['1~5','6~12','13~20','21~30']
        for n in range(nround):
            slfcal_tb_g= calprefix+'.G'+str(n)
            fig = plt.figure(figsize=(8.4, 7.))
            gs = gridspec.GridSpec(5, 6)
            #for s,sp in enumerate(spws[:-1]):
            for s, sp in enumerate(spws):
                #if s!=0 and n!=0: continue
                print ('processing spw: '+sp)
                cfreq=cfreqs[s]
                #bm=max(sbeam*cfreqs[1]/cfreq,10.)
                #bm=ebmsize[s]
                #bm = max(sbeam * cfreqs[1] / cfreq, 6.)
                bm = ebmsize[s]
                slfcal_img = imgprefix+'.spw'+sp.zfill(2)+'.slfcal'+str(n)
                #!!!!!!!!!!!!!!!!!!4~49!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1spw!!!!!!!!!!!!!!!!!
                if n == 0:
                    #if s < 35:
                        #spbg=max(int(sp)-2,1)
                    #    spbg=max(int(sp)-2,0)
                    #    sped=min(int(sp)+2,len(cfreqs)-1)
                    if s < 2:
                    #if s < 3:
                        spwran='0~1'
                    else:
                        #spbg=max(int(sp)-2,1)
                        #spbg=max(int(sp)-2,0)
                        spbg = max(int(sp) - 3, 2)
                        #spbg=max(int(sp)-2,4)
                        sped=min(int(sp)+3,len(cfreqs)-1)
                        spwran=str(spbg)+'~'+str(sped)
                    #spwran=str(spbg)+'~'+str(sped)
                    print('The current spw is: ', spwran)
                    if s <5:
                        mask = masks[1]
                        if s < 2:
                            mask = masks[0]
                    else:
                        mask = masks[int(s/5)+1]
                    #mask = masks[int(s / 5)]
                    '''
                    if int(sp) < 5:
                        mask = masks[0]
                    if int(sp) > 5 and int(sp) <= 12:
                        mask = masks[1]
                    if int(sp) > 12 and int(sp) <= 20:
                        mask = masks[2]
                    if int(sp) > 20 and int(sp) <= 30:
                        mask = masks[3]
                    '''
                    print ('using spw {0:s} as model'.format(spwran))
                    print ('using mask {0:s}'.format(mask))
                else:
                    if s<0:
                    #if s<3:
                        spwran = '3'
                    else:
                        spwran = sp
                    #if int(spwran) < 2:
                    #    spwran = '2'
                    if s <5:
                        mask = masks[1]
                        if s < 2:
                            mask = masks[0]
                    else:
                        mask = masks[int(s/5)+1]
                    #mask = masks[int(s / 5)]
                    '''
                    if int(sp) < 5:
                        mask = masks[0]
                    if int(sp) > 5 and int(sp) <= 12:
                        mask = masks[1]
                    if int(sp) > 12 and int(sp) <= 20:
                        mask = masks[2]
                    if int(sp) > 20 and int(sp) <= 30:
                        mask = masks[3]
                    #mask = masks[0]
                    '''
                #try:
                if True:
                    '''
                    tclean(vis=slfcalms,
                            antenna=antennas,
                            imagename=slfcal_img,
                            uvrange=uvranges[n],
                            #spw=sp,
                            spw=spwran,
                            specmode='mfs',
                            timerange=trange,
                            deconvolver='clark_exp',
                            imsize=[256,256],
                            #imsize=[512,512],
                            cell=['2arcsec'],
                            niter=niters[n],
                            gain=0.05,
                            #stokes='I',
                            stokes='XX',
                            weighting='briggs',
                            robust=robusts[n],
                            phasecenter=phacenter,
                            #mask='box [ [ 75pix , 90pix] , [205pix, 165pix ] ]',
                            mask=masks[s],
                            restoringbeam=['{}arcsec'.format(int(bm))],
                            #restoringbeam=[str(bm)+'arcsec'],
                            pbcor=False,
                            interactive=False)
                    '''
                    if n==0 and s<=11:
                        ntry=0
                        while True:
                            ntry+=1
                            tclean(vis=slfcalms,
                            antenna=antennas,
                            imagename=slfcal_img,
                            uvrange=uvranges[n],
                            #spw=sp,
                            spw=spwran,
                            #spw='0',
                            specmode='mfs',
                            timerange=trange,
                            #imagermode='csclean',
                            #psfmode='clark',
                            imsize=[512],
                            #imsize=[512,512],
                            cell=['5arcsec'],
                            niter=niters[n],
                            gain=0.05,
                            #stokes='I',
                            stokes='XX',
                            weighting='briggs',
                            robust=robusts[n],
                            phasecenter=phasecenter,
                            #mask='box [ [ 75pix , 90pix] , [205pix, 165pix ] ]',
                            mask=mask,
                            #restoringbeam=['{}arcsec'.format(int(ebmsize[int(s) - int(1)]))],
                            #restoringbeam=[str(bm)+'arcsec'],
                            restoringbeam=[bm],
                            pbcor=False,
                            interactive=False,
                            #usescratch=True,
                            savemodel='modelcolumn')
                            cur_header_refcal = imhead(slfcal_img+'.image')['refval']
                            print('have already try {} times'.format(ntry))
                            if np.max(abs(cur_header_refcal)) < 1.e11 and cur_header_refcal[-1] > 1.e2:
                               break
                            if n==0:
                                os.system('rm -rf ' + slfcal_img + '*')
                                #clearcal(slfcalms)
                                #os.system('rm -rf ' + caltbdir + '*')
                                time.sleep(1.0)
                                # if not os.path.exists(caltbdir):
                                #     os.makedirs(caltbdir)
                                # if not os.path.exists(imagedir):
                                #     os.makedirs(imagedir)
                        if os.path.exists(slfcal_img+'.image'):
                                    fitsfile=slfcal_img+'.fits'
                                    hf.imreg(vis=slfcalms,imagefile=slfcal_img+'.image',fitsfile=fitsfile,
                                     timerange=trange,usephacenter=False,toTb=True,verbose=False,overwrite=True)
                        print('NAILED IT')
                        #clnjunks = ['.mask','.flux', '.model', '.psf', '.residual', '.image']
                        clnjunks = ['.mask','.flux', '.model', '.psf', '.residual']
                        for clnjunk in clnjunks:
                            if os.path.exists(slfcal_img + clnjunk):
                                shutil.rmtree(slfcal_img + clnjunk)
                    else:
                        tclean(vis=slfcalms,
                            antenna=antennas,
                            imagename=slfcal_img,
                            uvrange=uvranges[n],
                            #spw=sp,
                            spw=spwran,
                            #spw='0',
                            specmode='mfs',
                            timerange=trange,
                            #imagermode='csclean',
                            #psfmode='clark',
                            imsize=[512],
                            #imsize=[512,512],
                            cell=['5arcsec'],
                            niter=niters[n],
                            gain=0.05,
                            #stokes='I',
                            stokes='XX',
                            weighting='briggs',
                            robust=robusts[n],
                            phasecenter=phasecenter,
                            #mask='box [ [ 75pix , 90pix] , [205pix, 165pix ] ]',
                            mask=mask,
                            #restoringbeam=['{}arcsec'.format(int(ebmsize[int(s) - int(1)]))],
                            #restoringbeam=[str(bm)+'arcsec'],
                            restoringbeam=[bm],
                            pbcor=False,
                            interactive=False,
                            #usescratch=True,
                            savemodel='modelcolumn')
                        if os.path.exists(slfcal_img+'.image'):
                            fitsfile=slfcal_img+'.fits'
                            hf.imreg(vis=slfcalms,imagefile=slfcal_img+'.image',fitsfile=fitsfile,
                                     timerange=trange,usephacenter=False,toTb=True,verbose=False,overwrite=True)
                        print('NAILED IT')
                        #clnjunks = ['.mask','.flux', '.model', '.psf', '.residual', '.image']
                        clnjunks = ['.mask','.flux', '.model', '.psf', '.residual']
                        for clnjunk in clnjunks:
                            if os.path.exists(slfcal_img + clnjunk):
                                shutil.rmtree(slfcal_img + clnjunk)
                    '''
                    ax = fig.add_subplot(gs[s])
                    eomap=smap.Map(fitsfile)
                    eomap.data=eomap.data.reshape((256,256))
                    eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                    eomap.plot()
                    eomap.draw_limb()
                    eomap.draw_grid()
                    ax.set_title(' ')
                    ax.set_xlim(xran)
                    ax.set_ylim(yran)
                    '''

                else:
                #except:
                    print('error in cleaning spw!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!: '+sp)
                    print('using nearby spws for initial model')
                    sp_e=int(sp)+2
                    sp_i=int(sp)-2
                    if sp_i < 0:
                        sp_i = 0
                    if sp_e > 30:
                        sp_e = 30
                    sp_=str(sp_i)+'~'+str(sp_e)
                    try:
                        tget(clean)
                        spw=sp_
                        clean()
                    except:
                        print('still not successful. abort...')
                        break

                # copy model from xx to yy
                #pdb.set_trace()
                mods.cpxx2yy(slfcalms)
                gaincal(vis=slfcalms, refant=refantenna,antenna=antennas,caltable=slfcal_tb_g,spw=sp, uvrange='',\
                        gaintable=[],selectdata=True,timerange=trange,solint='inf',gaintype='G',calmode=calmodes[n],\
                        combine='',minblperant=4,minsnr=2,append=True)
                if not os.path.exists(slfcal_tb_g):
                    #listcal(vis=slfcalms, caltable=slfcal_table)
                    print('No solution found in spw: '+sp)
            figname=imagedir+'slf_t{0:d}_n{1:d}.png'.format(t,n)
            plt.savefig(figname)
            time.sleep(1)
            plt.close()

            if os.path.exists(slfcal_tb_g):
                slftbs.append(slfcal_tb_g)
                slftb=[slfcal_tb_g]
                os.chdir(slfcaldir)
                if calmodes[n] == 'p':
                    plotcal(caltable=slfcal_tb_g,antenna='1~12',xaxis='freq',yaxis='phase',\
                            subplot=431,plotrange=[-1,-1,-180,180],iteration='antenna',figfile=slfcal_tb_g+'.png',showgui=False)
                if calmodes[n] == 'a':
                    plotcal(caltable=slfcal_tb_g,antenna='1~12',xaxis='freq',yaxis='amp',\
                            subplot=431,plotrange=[-1,-1,0,2.],iteration='antenna',figfile=slfcal_tb_g+'.png',showgui=False)
                os.chdir(workdir)
            if doapplycal[n]:
                clearcal(slfcalms)
                delmod(slfcalms)
                applycal(vis=slfcalms, gaintable=slftb, spw=','.join(spws), selectdata=True, \
                         antenna=antennas, interp='nearest', flagbackup=False, applymode='calonly', calwt=False)
            if n < nround - 1:
                #prompt = raw_input('Continuing to selfcal?')
                # prompt='y'
                #if prompt.lower() == 'n':
                #    if os.path.exists(slfcaledms):
                #        os.system('rm -rf ' + slfcaledms)
                #    split(slfcalms, slfcaledms, datacolumn='corrected')
                #    print 'Final calibrated ms is {0:s}'.format(slfcaledms)
                #    break
                #if prompt.lower() == 'y':
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

if doapply:
    import glob
    slftbs_comb=[]
    for n in range(nround):
        tbin=glob.glob(caltbdir+'slf_t*.G{0}'.format(n))
        mods.concat(tb_in=tbin,tb_out=caltbdir+'slf_comb.G{0}'.format(n))
        slftbs_comb.append(caltbdir+'slf_comb.G{0}'.format(n))
    '''
    pha1='/Volumes/WD6T/working/eovsa_full_disk/full_disk_model/20170715_1.pha'
    pha2='/Volumes/WD6T/working/eovsa_full_disk/full_disk_model/20170715_2.pha'
    amp3='/Volumes/WD6T/working/eovsa_full_disk/full_disk_model/20170715_3.amp'
    slftbs_comb.append(pha1)
    slftbs_comb.append(pha2)
    slftbs_comb.append(amp3)
    '''
    os.chdir(slfcaldir)
    #plotcal(caltable=slftbs_comb[0], antenna='1~12',spw='1~10',xaxis='time',yaxis='phase',subplot=341,
    #        iteration='antenna',plotrange=[-1,-1,-100,100])
    os.chdir(workdir)
    clearcal(refms_slfcal_XXYY)
    clearcal(refms_slfcal)
    applycal(vis=refms_slfcal_XXYY,gaintable=slftbs_comb,spw=','.join(spws),selectdata=True,\
             antenna=antennas,interp='linear',flagbackup=False,applymode='calonly',calwt=False)
    applycal(vis=refms_slfcal,gaintable=slftbs_comb,spw=','.join(spws),selectdata=True,\
             antenna=antennas,interp='linear',flagbackup=False,applymode='calonly',calwt=False)
    if os.path.exists(refms_slfcaled_XXYY):
        os.system('rm -rf '+refms_slfcaled_XXYY)
    if os.path.exists(refms_slfcaled):
        os.system('rm -rf '+refms_slfcaled)
    split(refms_slfcal_XXYY,refms_slfcaled_XXYY,datacolumn='corrected')
    split(refms_slfcal,refms_slfcaled,datacolumn='corrected')

if dofinalclean:
    import glob
    pol='X'
    if os.path.exists(workdir+'slfcal/masks_a.p'):
        masks_a=pickle.load(open(workdir+'slfcal/masks_a.p','rb'))
    for t,trange in enumerate(tranges):
        trange=tranges[t]
        img_final=imagedir_final+'/slf_final_{0}_t{1:d}'.format(pol,t)
        #masks=masks_a[t]
        slfcaledms=slfcaledms_a[t]
        spws=[str(s+1) for s in range(30)]
        tb.open(slfcaledms+'/SPECTRAL_WINDOW')
        reffreqs=tb.getcol('REF_FREQUENCY')
        bdwds=tb.getcol('TOTAL_BANDWIDTH')
        cfreqs=reffreqs+bdwds/2.
        tb.close()
        sbeam=35.
        fig = plt.figure(figsize=(12,10))
        gs = gridspec.GridSpec(6, 5)
        for s,sp in enumerate(spws):
            cfreq=cfreqs[s]
            bm=max(sbeam*cfreqs[1]/cfreq,2.)
            imname=img_final+'_s'+sp.zfill(2)
            fitsfile=imname+'.fits'
            if int(sp) < 5:
                mask = masks[0]
            if int(sp) > 5 and int(sp) <= 12:
                mask = masks[1]
            if int(sp) > 12 and int(sp) <= 20:
                mask = masks[2]
            if int(sp) > 20 and int(sp) <= 30:
                mask = masks[3]
            #print 'using mask {0:s}'.format(mask)
            if not os.path.exists(fitsfile):
                print('cleaning spw {0:s} with beam size {1:.1f}"'.format(sp,bm))
                try:
                    clean(vis=slfcaledms,
                            antenna=antennas,
                            imagename=imname,
                            spw=sp,
                            #mode='channel',
                            mode='mfs',
                            timerange=trange,
                            imagermode='csclean',
                            psfmode='clark',
                            imsize=[256,256],
                            cell=['2arcsec'],
                            niter=1000,
                            gain=0.05,
                            stokes=pol,
                            weighting='briggs',
                            robust=0.0,
                            restoringbeam=[str(bm)+'arcsec'],
                            phasecenter=phacenter,
                            #mask=mask,
                            mask='',
                            pbcor=True,
                            interactive=False,
                            usescratch=False)
                except:
                    print('cleaning spw '+sp+' unsuccessful. Proceed to next spw')
                    continue
                if os.path.exists(imname+'.image'):
                    hf.imreg(vis=slfcaledms,imagefile=imname+'.image',fitsfile=fitsfile,
                             timerange=trange,usephacenter=False,toTb=True,verbose=False)
                junks=['.flux','.model','.psf','.residual','.mask','.image']
                for junk in junks:
                    if os.path.exists(imname+junk):
                        os.system('rm -rf '+imname+junk)
            else:
                ax = fig.add_subplot(gs[s])
                eomap=smap.Map(fitsfile)
                eomap.data=eomap.data.reshape((256,256))
                eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                eomap.plot()
                eomap.draw_limb()
                eomap.draw_grid()
                ax.set_title(' ')
                ax.set_xlim(xran)
                ax.set_ylim(yran)

    plt.show()

if doclean_slfcaled:
    import glob
    pol='XX'
    if os.path.exists(workdir+'slfcal/masks_a.p'):
        masks_a=pickle.load(open(workdir+'slfcal/masks_a.p','rb'))
    for t,trange in enumerate(tranges):
        trange=tranges[t]
        img_final=imagedir_slfcaled+'/slf_final_{0}_t{1:d}'.format(pol,t)
        masks=masks_a[t]
        spws=[str(s+1) for s in range(30)]
        tb.open(refms_slfcaled_XXYY+'/SPECTRAL_WINDOW')
        reffreqs=tb.getcol('REF_FREQUENCY')
        bdwds=tb.getcol('TOTAL_BANDWIDTH')
        cfreqs=reffreqs+bdwds/2.
        tb.close()
        sbeam=35.
        from matplotlib import gridspec as gridspec
        from sunpy import map as smap
        from matplotlib import pyplot as plt
        fig = plt.figure(figsize=(12,10))
        gs = gridspec.GridSpec(6, 5)
        for s,sp in enumerate(spws):
            cfreq=cfreqs[int(sp)]
            #bm=max(sbeam*cfreqs[1]/cfreq,10.)
            bm=int(ebmsize[int(s) - int(1)])
            imname=img_final+'_s{0:0=2d}'.format(s+1)
            fitsfile=imname+'.fits'
            if int(sp) < 5:
                mask = masks[0]
            if int(sp) > 5 and int(sp) <= 12:
                mask = masks[1]
            if int(sp) > 12 and int(sp) <= 20:
                mask = masks[2]
            if int(sp) > 20 and int(sp) <= 30:
                mask = masks[3]
            #print 'using mask {0:s}'.format(mask)
            if not os.path.exists(fitsfile):
                print('cleaning spw {0:s} with beam size {1:.1f}"'.format(sp,bm))
                try:
                    clean(vis=refms_slfcaled_XXYY,
                            antenna=antennas,
                            imagename=imname,
                            spw=sp,
                            #mode='channel',
                            mode='mfs',
                            timerange=trange,
                            imagermode='csclean',
                            psfmode='clark',
                            imsize=[256,256],
                            cell=['2arcsec'],
                            niter=1000,
                            gain=0.05,
                            stokes=pol,
                            weighting='briggs',
                            robust=0.0,
                            restoringbeam=[str(bm)+'arcsec'],
                            phasecenter=phacenter,
                            #mask=mask,
                            mask='',
                            pbcor=True,
                            interactive=False,
                            usescratch=False)
                except:
                    print('cleaning spw '+sp+' unsuccessful. Proceed to next spw')
                    continue
                if os.path.exists(imname+'.image'):
                    hf.imreg(vis=refms_slfcaled_XXYY,imagefile=imname+'.image',fitsfile=fitsfile,
                             timerange=trange,usephacenter=False,toTb=True,verbose=False)
                junks=['.flux','.model','.psf','.residual','.mask','.image']
                for junk in junks:
                    if os.path.exists(imname+junk):
                        os.system('rm -rf '+imname+junk)
            else:
                ax = fig.add_subplot(gs[s])
                eomap=smap.Map(fitsfile)
                eomap.data=eomap.data.reshape((256,256))
                eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                eomap.plot()
                eomap.draw_limb()
                eomap.draw_grid()
                ax.set_title(' ')
                ax.set_xlim(xran)
                ax.set_ylim(yran)

    plt.show()
