import sys
sys.path.append('/Users/walterwei/opt/miniconda3/envs/scasa_env/lib/python3.6/site-packages')
from suncasa.utils import helioimage2fits as hf
import os
import numpy as np
import pickle
import pdb
import suncasa.utils.mod_slftbs as mods
from matplotlib import gridspec as gridspec
from sunpy import map as smap
from matplotlib import pyplot as plt
import time
from suncasa.utils import plot_map
import shutil
from astropy.time import Time,TimeDelta


#task handles
dofullsun=0
domasks=0
doslfcal=0
doapply=1
dofinalclean=0 #final clean of the slfcaled data for the selected slfcal time only
doclean_slfcaled=0 #final clean of all the slfcaled data


if doapply:
    workdir = '/Volumes/Data/20170820/20220511/eovsa/eovsa_full/'
    slfcaldir = workdir + 'slfcal/'
    if not os.path.exists(slfcaldir):
        os.makedirs(slfcaldir)
    import glob
    slftbs_comb=[]
    caltbdir = '/Volumes/Data/20170820/20220511/eovsa/eovsa_data/slfcal/caltbs/'
    refms = workdir + 'msdata/IDB20220511_1800-2000.ms'
    refms_slfcal_XXYY = refms + '.XXYY.slfcal'
    # refms_slfcal_XX = refms + '.XX.slfcal'
    refms_slfcaled_XXYY = refms + '.XXYY.slfcaled'
    # refms_slfcaled_XX = refms + '.XX.slfcaled'
    refms_slfcal = refms + '.slfcal'
    refms_slfcaled = refms + '.slfcaled'
    tb.open(refms + '/SPECTRAL_WINDOW')
    reffreqs = tb.getcol('REF_FREQUENCY')
    bdwds = tb.getcol('TOTAL_BANDWIDTH')
    cfreqs = reffreqs + bdwds / 2.
    tb.close()
    spws = [str(s) for s in range(len(cfreqs))]
    antennas = '0~12'
    for n in range(3):
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
    #clearcal(refms_slfcal_XXYY)
    clearcal(refms)
    #applycal(vis=refms_slfcal_XXYY,gaintable=slftbs_comb,spw=','.join(spws),selectdata=True,\
    #         antenna=antennas,interp='linear',flagbackup=False,applymode='calonly',calwt=False)
    applycal(vis=refms,gaintable=slftbs_comb,spw=','.join(spws),selectdata=True,\
             antenna=antennas,interp='linear',flagbackup=False,applymode='calonly',calwt=False)
    if os.path.exists(refms_slfcaled_XXYY):
        os.system('rm -rf '+refms_slfcaled_XXYY)
    if os.path.exists(refms_slfcaled):
        os.system('rm -rf '+refms_slfcaled)
    #split(refms_slfcal_XXYY,refms_slfcaled_XXYY,datacolumn='corrected')
    split(refms_slfcal,refms_slfcaled,datacolumn='corrected')


