import time, datetime
t0 = time.time()
t_elapsed = 0
t0 = time.time()
import suncasa
print("building config file...")
import os
homedir = os.path.expanduser('~/.casa')
if not os.path.exists(homedir):
    os.system("mkdir {}".format(homedir))
with (open(os.path.join(homedir,'config.py'), 'w')) as fi:
    fi.write("import sys \n")
    fi.write("import os \n")
    fi.write("import sysconfig \n")
    fi.write("import casatools \n")
    fi.write("import casadata \n")
    fi.write("import time \n")
    fi.write("logfile='casalog-{}.log'.format(time.strftime('%Y%m%d-%H',time.localtime())) \n")
    fi.write("telemetry_enabled = False \n")
    fi.write("crashreporter_enabled = True \n")
    fi.write("tb=casatools.table() \n")
    fi.write("ospathsep = os.path.sep \n")
    fi.write("libpath = sysconfig.get_paths()['purelib'] \n")
    fi.write("obsdict = {'MJD': 57447.0, 'Name': 'EOVSA', 'Type': 'WGS84', 'Long': -118.287, \n")
    fi.write("            'Lat': 37.2332, 'Height': 1207.13, 'X': 0.0, 'Y': 0.0, 'Z': 0.0,  \n")
    fi.write("            'Source': 'Dale Gary'} \n")
    fi.write("obstable = os.path.join(casadata.datapath,'geodetic','Observatories') \n")
    fi.write("tb.open(obstable, nomodify=True) \n")
    fi.write("if 'EOVSA' not in tb.getcol('Name'): \n")
    fi.write("    print('Adding EOVSA to the Observatories') \n")
    fi.write("    tb.close() \n")
    fi.write("    tb.open(obstable, nomodify=False) \n")
    fi.write("    nrows = tb.nrows() \n")
    fi.write("    tb.addrows(1) \n")
    fi.write("    for k in obsdict.keys(): \n")
    fi.write("        tb.putcell(k, nrows, obsdict[k])     \n")
    fi.write("tb.close() \n")
    fi.close()

print('Config file is written successfully!')

t_exec = time.time()-t0
t_elapsed += t_exec
print('Wall time of this cell: {}. Elapsed time since the beginning {}:'.format(datetime.timedelta(seconds=t_exec),datetime.timedelta(seconds=t_elapsed)))

import time, datetime
import matplotlib.pyplot as plt
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/eovsa_full/msdata/IDB20220511_1800-2000.ms.slfcaled'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/1700_1810/IDB20220511_1812-1814.ms'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/1800_2000/playground/IDB20220511_1830-1930.ms'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/1800_2000/IDB20220511_1800-2000.ms'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/background/IDB20220511_1805-1815.ms.XXYY.slfcaled'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/bkg_subed/playground/IDB20220511_1830_1840_sub_1805_1808.ms.XXYY.slfcaled'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/eovsa_data/msdata/IDB20220511_1830-1850.ms'
#specfile = '/Volumes/Data/20170820/20220511/eovsa/background/IDB20220511_1805-1815.ms.XXYY.slfcaled.dspec.npz'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/play_ground/original/IDB20220511_1830-1850.ms'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/1800_2000/playground/1842_ori/IDB20220511_1842-1843.ms'
visibility_data = '/Volumes/WD6T/working/20220511/slfcal_sijie/IDB20220511_1830-1850.cal.ms.slfcaled'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/1800_2000/full_disk_slfcal/playground/IDB20220511_1842_11_17.ms.t0.XXYY.slfcaled'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/1800_2000/data/playground/IDB20220511_1842_11_17.ms.t0.XXYY.slfcaled'
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/1800_2000/playground/1810/IDB20220511_1810-1811.ms.slfcaled'
#specfile = '/Volumes/Data/20170820/20220511/eovsa/play_ground/original/IDB20220511_1830-1850.ms.dspec.npz'
t0 = time.time()
from suncasa.utils import qlookplot as ql
from suncasa.utils.mstools import time2filename
#visibility_data = '/Volumes/Data/20170820/20220511/eovsa/eovsa_data/msdata/IDB20220511_1830-1850.ms.slfcal'
## (Optional) Supply the npz file of the dynamic spectrum from previous step.
## If not provided, the program will generate a new one from the visibility.
## set the time interval
#timerange = '18:40:25~18:40:45'
#timerange = '18:30:00~18:30:20'
#timerange = '18:45:00~18:45:20'
timerange = '18:42:11~18:42:17'
#timerange = '18:50:00~18:50:10'
#timerange = '18:40:36~18:40:40'
#timerange = ''
#t0 = time.time()
from suncasa import dspec as ds
import time
# The example below shows the cross-power spectrogram from a baseline selected using the parameter "bl".
# bl = '4&9' means selecting a baseline from Antenna ID 4 (Antenna Name "eo05") correlating with Antenna ID 9 (Antenna Name "eo10") - c.f., listobs outputs.
# you can also use the "bl" parameter to select multiple baseline(s), i.e., bl='0&2;4&9;8&11'.
#specfile = visibility_data + '.dspec.npz'
specfile = '/Volumes/WD6T/working/20220511/slfcal_sijie/msdata/IDB20220511_1830-1850.cal.ms.dspec.npz'
# d = ds.Dspec(visibility_data, bl='5&10', specfile=specfile)
# d.plot(vmin=0.0, vmax=50, pol='XX')
# print(d.data.shape) # npol, nbl, nfreq, ntime
aiafits='/Volumes/WD6T/working/20220511/aia/aia.lev1_euv_12s.2022-05-11T184221Z.171.image_lev1.fits'
#t_exec = time.time()-t0
#t_elapsed += t_exec
#print('Wall time of this cell: {}. Elapsed time since the beginning {}:'.format(datetime.timedelta(seconds=t_exec),datetime.timedelta(seconds=t_elapsed)))

## select (almost) all spectral windows from spw id #0 to #47
#spw = ['{}'.format(l) for l in range(48)]
spw = ['{}'.format(l) for l in range(43)]
outfits = time2filename(visibility_data,timerange=timerange)+'.outim.image.allbd.fits'
## select stokes XX
stokes = 'XX'
## image center for clean in solar X-Y in arcsec
#xycen=[912.0,-295.0]
#xycen=[912.0,-195.0]
xycen=[0.0,0.0]
#xycen=[920.0,-200.0]
## pixel scale
#cell=['2.0arcsec']
cell=['5.0arcsec']
#cell=['1.2arcsec']
## number of pixels in X and Y. If only one value is provided, NX = NY
#imsize=[128]
imsize=[512]
## field of view of the zoomed-in panels in unit of arcsec
fov = [2460,2460]
## turn off AIA image plotting, default is True
plotaia = True
## AIA passband in Å. The options are [171,131,304,335,211,193,94,1600,1700]
aiawave = 171
## Choose the coloar map for AIA images. If not provided, the program will use default AIA colormap.
acmap = 'gray_r'

ql.qlookplot(vis=visibility_data, specfile=specfile, timerange=timerange,aiafits=aiafits,
             spw=spw, stokes=stokes, plotaia=plotaia, aiawave=aiawave,
             restoringbeam=['30arcsec'], robust = 0.5, acmap=acmap,
             imsize=imsize,cell=cell,xycen=xycen,fov=fov,
             outfits=outfits,overwrite=False,clevels=[0.3,1,2])

t_exec = time.time()-t0
t_elapsed += t_exec
print('Wall time of this cell: {}. Elapsed time since the beginning {}:'.format(datetime.timedelta(seconds=t_exec),datetime.timedelta(seconds=t_elapsed)))