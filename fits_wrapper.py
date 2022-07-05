import glob
import os
from astropy.io import fits
import pickle
import numpy as np
def fits_wrap_spwX(inp_time, high_tres=False, custom_fits_list=None):
    with open('/Volumes/Data/20170820/20220511/info/cfreqs.p', 'rb') as fcfreq:
        cfreq = pickle.load(fcfreq, encoding='latin1')
    fcfreq.close()
    cur_dic = '/Volumes/Data/20170820/20220511/info/20220511_10s_long_aia.p'
    in_dic_file = open(cur_dic, 'rb')
    in_dic = pickle.load(in_dic_file)
    in_dic_file.close()
    fitsfiles = in_dic[inp_time]['radio_sbs']
    if custom_fits_list is not None:
        custom_list = glob.glob(
            #'/Volumes/Data/20170820/20220511/eovsa/eovsa_full/slfcal/images_slfcaled/20220511T*.500slfcaled_tb_final_XX_10s_XX_t{0}_s*.fits'.format(
                '/Volumes/Data/20170820/20220511/eovsa/bkg_subed/slfcal/images_slfcaled/20220511T*.500slfcaled_tb_final_XX_10s_XX_t{0}_s*.fits'.format(
            inp_time))
        print(custom_list)
        fitsfiles = custom_list
        print(fitsfiles)
    #outfitsfile = '/Volumes/Data/20170820/20220511/eovsa/allbd_wrapped/10s/T{0:0=4d}_allbd.fits'.format(inp_time)
    outfitsfile = '/Volumes/Data/20170820/20220511/eovsa/allbd_wrapped/bkg_subed/10s/T{0:0=4d}_allbd_40spws.fits'.format(inp_time)
    # if custom_fits_list is not None:
    #     #outfitsfile = '/Volumes/Data/20170820/eovsa/allbd_wrapped/10s/subed_T_3820_3830_3730_3740_allbd.fits'
    #     #outfitsfile = '/Volumes/Data/20170820/eovsa/allbd_wrapped/10s/seperated_subed_T{0:0=4d}_allbd.fits'.format(inp_time)
    #     outfitsfile = '/Volumes/Data/20170820/eovsa/allbd_wrapped/seperated_1s_subed_T{0:0=4d}_allbd.fits'.format(
    #     inp_time)
    fitsfiles = sorted(fitsfiles)[0:40]
    nband = len(fitsfiles)
    fits_exist = []
    idx_fits_exist = []
    for sidx, fitsf in enumerate(fitsfiles):
        if os.path.exists(fitsf):
            fits_exist.append(fitsf)
            idx_fits_exist.append(sidx)
    if len(fits_exist) == 0: raise ValueError('None of the input fitsfiles exists!')
    os.system('cp {} {}'.format(fits_exist[0], outfitsfile))
    hdu0 = fits.open(outfitsfile, mode='update')
    header = hdu0[0].header
    npol, nbd, ny, nx = int(header['NAXIS4']), nband, int(header['NAXIS2']), int(header['NAXIS1'])
    data = np.zeros((npol, nbd, ny, nx))
    cfreqs = []
    for sidx, fitsf in enumerate(fits_exist):
        hdu = fits.open(fitsf)
        cfreqs.append(hdu[0].header['CRVAL3'])
        for pidx in range(npol):
            if len(hdu[0].data.shape) == 2:
                data[pidx, idx_fits_exist[sidx], :, :] = hdu[0].data
            else:
                data[pidx, idx_fits_exist[sidx], :, :] = hdu[0].data[pidx, 0, :, :]
    df = np.nanmean(np.diff(cfreqs) / np.diff(idx_fits_exist))  ## in case some of the band is missing
    header['cdelt3'] = df
    header['NAXIS3'] = nband
    header['NAXIS'] = 4
    header['CRVAL3'] = header['CRVAL3'] - df * idx_fits_exist[0]
    if os.path.exists(outfitsfile):
        os.system('rm -rf {}'.format(outfitsfile))
    fits.writeto(outfitsfile, data, header)
    print('wrapped fits writed to ' + outfitsfile)
    return
def main():
    fits_wrap_spwX(inp_time=109)

if __name__ == '__main__':
    main()