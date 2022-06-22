from suncasa.utils import stackplotX as stp
import glob
import os.path

# listing and sorting all the aia fits files
trange=['190000','190500']
aiafiles = glob.glob('/Volumes/Data/20170820/20220511/aia131/aia.lev1_euv_12s.2022-05-11T*Z.131.image_lev1.fits')
aiafiles.sort()
sub_aiafiles = [af for af in aiafiles if int(os.path.basename(af)[28:34]) >= int(trange[0]) and
                int(os.path.basename(af)[28:34]) <= int(trange[1])]
#print(sub_aiafiles)
cur_st = stp.Stackplot()
mapcubefile = '/Volumes/Data/20170820/20220511/st_obj/aia_131_2.mapcube'
binpix = 1
fov=[800,1100,-400,-100]
cur_st.make_mapseq(sub_aiafiles, outfile=mapcubefile, tosave=False, binpix=binpix,fov=fov,
                   superpixel=True,aia_prep=True)
# save the mapseqs to a file
cur_st.mapseq_tofile(outfile=mapcubefile)

# make running diff imgs
cur_st.mapseq_mkdiff(dt=24.)

#plot the runnning diff imgs
cur_st.plot_mapseq(diff=True, vmin=-500, vmax=500)
# draw the slit on the img

#plot_the stackplot(time-distance plot)
cur_st.plot_stackplot(diff=True, vmin=-500, vmax=500)

#click the trajectory button