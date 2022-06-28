from suncasa.utils import helioimage2fits as hf
import os
import numpy as np

from astropy.time import Time,TimeDelta
from casatools import table as tbtool
from casatools import ms as mstool
from casatools import quanta as qatool
from casatools import image as iatool

tb = tbtool()
ms = mstool()
qa = qatool()
ia = iatool()
from casatasks import tclean
from casatasks import split
from suncasa.suncasatasks import ptclean6 as ptclean

#ebmsize=[40.89866667, 35.        , 30.52746815, 27.06907583, 24.30586371,
#       22.0676259 , 20.20685112, 18.63547995, 17.29086809, 16.12723449,
#       15.11034483, 14.21408712, 13.41819773, 12.70671085, 12.06687648,
#       11.48838951, 10.96283059, 10.48325359, 10.04387688, 10.        ,
#       10.        , 10.        , 10.        , 10.        , 10.        ,
#       10.        , 10.        , 10.        , 10.        , 10.        ]
#task handles
dofullsun=0
domasks=0
doslfcal=0
doapply=0
dofinalclean=0 #final clean of the slfcaled data for the selected slfcal time only
doclean_slfcaled=1 #final clean of all the slfcaled data

#workdir='/Volumes/Data/20170906/msdata/2021_new_calibration/'
#workdir='/Volumes/Data/20170715/eovsa_large_fov/'
#workdir='/Volumes/WD6T/working/20170703/eovsa/'
workdir='/Volumes/Data/20170820/20220511/eovsa/eovsa_full/'
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
#refms = workdir+'msdata/IDB20170703_concat.ms'
refms = workdir+'msdata/IDB20220511_1800-2000.ms.XXYY.slfcaled'
#refms = workdir+'msdata_stable/IDB20170910T154625-161625.ms.corrected'
refms_slfcal_XXYY = refms + '.XXYY.slfcal'
#refms_slfcal_XX = refms + '.XX.slfcal'
refms_slfcaled_XXYY = refms + '.XXYY.slfcaled'
#refms_slfcaled_XX = refms + '.XX.slfcaled'
refms_slfcal = refms + '.slfcal'
refms_slfcaled = refms + '.slfcaled'
#refms_slfcaled = '/Volumes/Data/20170906/msdata/2021_new_calibration/msdata/IDB20170906191320-20170906194320.corrected.ms.XXYY.slfcaled'
#refms_slfcaled = '/Volumes/WD6T/working/bkg_sub_eovsa_20170715/msdata/bkg_subed_IDB20170715_concat.ms'
#refms_slfcaled = '/Volumes/WD6T/working/bkg_sub_eovsa_20170715/msdata/IDB20170715_concat.ms'
#refms_slfcaled = '/Volumes/WD6T/working/20170703/eovsa_10s/msdata/IDB20170703_concat.ms.XXYY.slfcaled'
refms_slfcaled = '/Volumes/Data/20170820/20220511/eovsa/eovsa_full/msdata/IDB20220511_1800-2000.ms.XXYY.slfcaled'

init_time='2022-05-11T18:36:00.002'
time_interval=10.0
tranges=[]
tname=[]

start_time = Time('2022-05-11 18:40:36',format='iso')
end_time = Time('2022-05-11 18:40:40',format='iso')
#find the phasecenter
midtime_mjd = (start_time.mjd + end_time.mjd) / 2.
eph = hf.read_horizons(t0=Time(midtime_mjd, format='mjd'))
tb.open(refms + '/FIELD')
phadir = tb.getcol('PHASE_DIR').flatten()
tb.close()

ra0 = phadir[0]
dec0 = phadir[1]
xycen=[912.0,-295.0]
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

tb.open(refms + '/SPECTRAL_WINDOW')
reffreqs = tb.getcol('REF_FREQUENCY')
bdwds = tb.getcol('TOTAL_BANDWIDTH')
cfreqs = reffreqs + bdwds / 2.
tb.close()
print('cfreqs are: ', cfreqs)

sbeam = 40.
bmsize = sbeam * 1.6e9 / np.asarray(cfreqs, dtype=float)
bmsize[bmsize < 4.0] = 4.0
#restoringbms = mstools.get_bmsize(cfreqs, refbmsize=sbeam, reffreq=1.6, minbmsize=4.0)
ebmsize = ['{:.1f}arcsec'.format(crbm) for crbm in bmsize]

#sbeam = 40.
#ebmsize=np.zeros_like(cfreqs)
#for iii in range(len(cfreqs)):
#    cfreq = cfreqs[iii]
#    ebmsize[iii] = max(sbeam * cfreqs[1] / cfreq, 6.)
    #cur_spw = str(ebmsize[iii])
print('bmsizes are: ', ebmsize)
for tint in range(300):
    init_t=Time(init_time,format='isot')
    timed1=TimeDelta((tint*time_interval)*1.0, format='sec')
    timed2=TimeDelta(time_interval, format='sec')
    start_time=init_t+timed1
    end_time=start_time+timed2
    tmp_trange = '{0}~{1}'.format(start_time.iso.replace('-', '/').replace(' ', '/'),end_time.iso.replace('-', '/').replace(' ', '/'))
    tranges.append(tmp_trange)
    tname.append(tmp_trange.split(':')[3]+tmp_trange.split(':')[4].split('.')[0])

#phacenter='J2000 11h00m49 06d14m38'
#xran=[500,700]
#yran=[-350,-150]
#phacenter='J2000 2.00406343122rad 0.373837961939rad'
'''
midtime_mjd = (start_time.mjd + end_time.mjd) / 2.
eph = hf.read_horizons(t0=Time(midtime_mjd, format='mjd'))
tb.open(refms_slfcaled + '/FIELD')
phadir = tb.getcol('PHASE_DIR').flatten()
tb.close()
ra0 = phadir[0]
dec0 = phadir[1]
xycen=[955.0,50.0]
x0 = np.radians(xycen[0] / 3600.)
y0 = np.radians(xycen[1] / 3600.)
p0 = np.radians(eph['p0'][0])  # p angle in radians
raoff = -((x0) * np.cos(p0) - y0 * np.sin(p0)) / np.cos(eph['dec'][0])
decoff = (x0) * np.sin(p0) + y0 * np.cos(p0)
newra = ra0 + raoff
newdec = dec0 + decoff
phasecenter = 'J2000 ' + str(newra) + 'rad ' + str(newdec) + 'rad'
print('Phasecenter: ',phasecenter)
'''
#phasecenter='J2000 1.7868374584594264rad 0.4002520172843901rad'
#phacenter='J2000 1.7918694030439997rad 0.400123935501163rad'
#xran=[617,867]
#yran=[-271,-21]
spws=[str(s) for s in range(len(cfreqs))]
antennas='!4;!9'
nround=3 #number of slfcal cycles

subms_a=[]
slfcalms_a=[]
slfcaledms_a=[]

#for t in range(298):
#for t in range(180):
#for t in range(180):
for t in range(300):
    trange=tranges[t]
    slfcaledms_ = workdir+'slfcal/IDB20220511.ms.t{0:d}.XXYY.slfcaled'.format(t)
    #slfcaledms_ = workdir+'slfcal/IDB20170703.ms.t{0:d}.XX'.format(t)
    #if t >= 40 and t<=44:
    #if t >= 90 and t<=94:
    #if t >= 0 and t<=2:
    #if t>=40 and t<45:
    #if t>=144:
    #if t>=150:
    if not os.path.exists(slfcaledms_) and t in [36,141]:
        print(slfcaledms_+' no file found')
        split(vis=refms_slfcaled,outputvis=slfcaledms_,datacolumn='data',timerange=trange,correlation='')
    slfcaledms_a.append(slfcaledms_)

if doclean_slfcaled:
    import glob
    pol='XX'
    #if os.path.exists(workdir+'slfcal/masks_a.p'):
    #    masks_a=pickle.load(open(workdir+'slfcal/masks_a.p','rb'))
    for t,trange in enumerate(tranges):
        if not t in [36,141]: continue
        #if t<40 or t>44: continue
        #if t<90 or t>94: continue
        #if t<0 or t>2: continue
        #if t< 40 or t >=45: continue
        #if t >=20: continue
        #if t < 144: continue
        trange=tranges[t]
        slfcaledms = slfcaledms_a[t]
        #img_final=imagedir_slfcaled+'/slf_final_{0}_t{1:d}'.format(pol,t)
        #img_final=imagedir_slfcaled+'slfcaled_tb_final_XX_10s_{0}_t{1:d}'.format(pol,t)
        img_final='slfcaled_tb_final_XX_10s_{0}_t{1:d}'.format(pol,t)
        #masks=masks_a[t]
        #spws=[str(s+1) for s in range(30)]
        #tb.open(slfcaledms+'/SPECTRAL_WINDOW')
        #tb.open(refms_slfcaled_XXYY+'/SPECTRAL_WINDOW')
        #reffreqs=tb.getcol('REF_FREQUENCY')
        #bdwds=tb.getcol('TOTAL_BANDWIDTH')
        #cfreqs=reffreqs+bdwds/2.
        #tb.close()
        #sbeam=35.
        from matplotlib import gridspec as gridspec
        from sunpy import map as smap
        from matplotlib import pyplot as plt
        #fig = plt.figure(figsize=(12,10))
        gs = gridspec.GridSpec(6, 5)
        for s,sp in enumerate(spws):
            #if s<4:
            #    sp = '4'
            #if s!=1: continue
            cfreq=cfreqs[s]
            #bm = max(sbeam * cfreqs[1] / cfreq, 6.)
            bm = ebmsize[s]
            '''
            if s < 20:
                spbg = max(int(s) - 4, 1)
                sped = min(int(s) + 4, 30)
            else:
                spbg = max(int(s) - 2, 1)
                sped = min(int(s) + 2, 30)
            '''
            spbg = max(int(s) - 1, 1)
            sped = min(int(s) + 1, 30)
            spwran = str(spbg) + '~' + str(sped)
            #bm=max(sbeam*cfreqs[1]/cfreq,10.)
            #bm=int(ebmsize[s])
            imname=img_final+'_s{0:0=2d}'.format(s+1)
            fitsfile=imname+'.fits'
            print(fitsfile)
            '''
            if int(sp) < 5:
                mask = masks[0]
            if int(sp) > 5 and int(sp) <= 12:
                mask = masks[1]
            if int(sp) > 12 and int(sp) <= 20:
                mask = masks[2]
            if int(sp) > 20 and int(sp) <= 30:
                mask = masks[3]
            #print 'using mask {0:s}'.format(mask)
            '''
            if not os.path.exists(fitsfile) and t<180:
                #print 'cleaning spw {0:s} with beam size {1:.1f}"'.format(sp,bm)
                print('not existing, cleaning spw {0} with beam size {1}"'.format(sp,bm))
                # try:
                #     tclean(vis=slfcaledms,
                #     #tclean(vis=refms_slfcaled,
                #             antenna=antennas,
                #             imagename=imname,
                #             spw=sp,
                #             #spw='5',
                #             #spw=spwran,
                #             #mode='channel',
                #             specmode='mfs',
                #             timerange=trange,
                #             #imagermode='csclean',
                #             #psfmode='clark',
                #             imsize=[256],
                #             cell=['2arcsec'],
                #             niter=1000,
                #             gain=0.1,
                #             stokes='XX',
                #             weighting='briggs',
                #             #robust=2.0,
                #             robust=0.0,
                #             #restoringbeam=[str(bm)+'arcsec'],
                #             restoringbeam=[bm],
                #             #restoringbeam=['30arcsec'],
                #             phasecenter=phasecenter,
                #             #mask=mask,
                #             mask='',
                #             pbcor=True,
                #             interactive=False)
                #             #usescratch=False)

                res = ptclean(vis=slfcaledms,
                              imageprefix=imagedir_slfcaled,
                              imagesuffix=imname,
                              antenna='!4;!9',
                              timerange=trange,
                              twidth=12,
                              spw=sp,
                              restoringbeam=[bm],
                              ncpu=6,
                              niter=1000,
                              gain=0.1,
                              imsize=[256],
                              cell=['2arcsec'],
                              stokes='XX',
                              doreg=True,
                              usephacenter=True,
                              phasecenter=phasecenter,
                              docompress=False,
                              toTb=True,
                              pbcor=True,
                              weighting='briggs',
                              robust=0.0)
                # except:
                #     print('cleaning spw '+sp+' unsuccessful. Proceed to next spw')
                #     continue
                if os.path.exists(imname+'.image'):
                    hf.imreg(vis=slfcaledms,imagefile=imname+'.image',fitsfile=fitsfile,
                             timerange=trange,toTb=True,verbose=False,)
                             #timerange = trange, usephacenter = False, toTb = False, verbose = False)
                junks=['.flux','.model','.psf','.residual','.mask','.image','.sumwt','.pb','.image.pbcor']
                for junk in junks:
                    if os.path.exists(imname+junk):
                        os.system('rm -rf '+imname+junk)
            else:
                print(fitsfile+' has already been cleaned')
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

    #plt.show()
