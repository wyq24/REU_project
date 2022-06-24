from sunpy.net import Fido
from sunpy.net import attrs as a
from astropy.time import Time
from astropy import units as u
import os

outdir = '/Volumes/Data/20170820/20220511/aia/'
tst = Time('2022-05-11 18:30:00',format='iso')
tet = Time('2022-05-11 19:00:00',format='iso')
wavelength = [94.,131.,171.,193.,211.,335.,304.]
for widx, wave in enumerate(wavelength):
    wave1 = wave - 3.0
    wave2 = wave + 3.0
    qr = Fido.search(a.Time(tst.iso, tet.iso),
                     a.Instrument.aia,
                     a.Wavelength(wave1 * u.AA, wave2 * u.AA))
    res = Fido.fetch(qr)
    for ll in res:
        vsonamestrs = ll.split('_')
        if vsonamestrs[2].startswith('1600') or vsonamestrs[2].startswith('1700'):
            product = 'aia.lev1_uv_24s'
        else:
            product = 'aia.lev1_euv_12s'
        jsocnamestr = product + '.' + '{}-{}-{}{}{}Z.'.format(vsonamestrs[3], vsonamestrs[4], vsonamestrs[5],
                                                              vsonamestrs[6],
                                                              vsonamestrs[7]).upper() + vsonamestrs[2][
                                                                                        :-1] + '.image_lev1.fits'
        print(ll, jsocnamestr)
        os.system('mv {} {}/{}'.format(ll, outdir, jsocnamestr))