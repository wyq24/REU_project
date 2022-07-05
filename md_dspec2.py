import matplotlib.gridspec as gridspec
import pickle
import numpy as np
import os
import datetime
import struct
from scipy.io.idl import readsav
from datetime import datetime
# import md_quanta as qa
# import quanta as qa
import matplotlib.dates as mdates
from matplotlib.dates import date2num, AutoDateFormatter, AutoDateLocator
from matplotlib.colors import SymLogNorm
import matplotlib as mpl
import numpy.ma as ma
#import ex_plot_1 as ep
import plotting_tools as pt
import warnings

font = {'weight': 'bold',
        'size': 12}
mpl.rc('font', **font)
'''
def get_dspec(vis=None, savespec=True, specfile=None, bl='', uvrange='', field='', scan='', datacolumn='data',
              domedian=False, timeran=None, spw=None, timebin='0s', verbose=False):
    from split_cli import split_cli as split

    msfile = vis
    if not spw:
        spw = ''
    if not timeran:
        timeran = ''
    if not bl:
        bl = ''
    if domedian:
        if not uvrange:
            uvrange = '0.2~0.8km'
    else:
        uvrange = ''
    # Open the ms and plot dynamic spectrum
    if verbose:
        print 'Splitting selected data...'
    vis_spl = './tmpms.splitted'
    if os.path.exists(vis_spl):
        os.system('rm -rf ' + vis_spl)

    split(vis=msfile, outputvis=vis_spl, timerange=timeran, antenna=bl, field=field, scan=scan, spw=spw,
          uvrange=uvrange, timebin=timebin, datacolumn=datacolumn)
    ms.open(vis_spl, nomodify=False)
    if verbose:
        print 'Regridding into a single spectral window...'
        # print 'Reading data spw by spw'
    ms.cvel(outframe='LSRK', mode='frequency', interp='nearest')
    ms.selectinit(datadescid=0, reset=True)
    data = ms.getdata(['amplitude', 'time', 'axis_info'], ifraxis=True)
    ms.close()
    os.system('rm -rf ' + vis_spl)
    specamp = data['amplitude']
    (npol, nfreq, nbl, ntim) = specamp.shape
    if verbose:
        print 'npol, nfreq, nbl, ntime:', data['amplitude'].shape
    spec = np.swapaxes(specamp, 2, 1)
    print(spec.shape)
    freq = data['axis_info']['freq_axis']['chan_freq'].reshape(nfreq)
    tim = data['time']

    if domedian:
        if verbose:
            print('doing median of all the baselines')
        # mask zero values before median
        spec_masked = np.ma.masked_where(spec < 1e-9, spec)
        spec_med = np.ma.filled(np.ma.median(spec_masked, axis=1), fill_value=0.)
        nbl = 1
        ospec = spec_med.reshape((npol, nbl, nfreq, ntim))
    else:
        ospec = spec
    # Save the dynamic spectral data
    if savespec:
        if not specfile:
            specfile = msfile + '.dspec.npz'
        if os.path.exists(specfile):
            os.system('rm -rf ' + specfile)
        np.savez(specfile, spec=ospec, tim=tim, freq=freq,
                 timeran=timeran, spw=spw, bl=bl, uvrange=uvrange)
        if verbose:
            print 'Median dynamic spectrum saved as: ' + specfile

    return {'spec': ospec, 'tim': tim, 'freq': freq, 'timeran': timeran, 'spw': spw, 'bl': bl, 'uvrange': uvrange}
'''


def plt_dspec(specdata, pol='I', dmin=None, dmax=None,over_all_tr=None,
              timerange=None, freqrange=None, timestr=True,
              movie=False, framedur=60., dtframe=10.,
              goessav=None, goes_trange=None,
              savepng=True, savepdf=False, pax=None, time_line=None, rhessisave=None,leg_loc='best',custom_cm='jet',customized_light_curve=False):
    """
    timerange: format: ['2012/03/10/18:00:00','2012/03/10/19:00:00']
    freqrange: format: [1000.,1500.] in MHz
    movie: do a movie of dynamic spectrum?
    framedur: time range of each frame
    dtframe: time difference of consecutive frames
    goessav: provide an IDL save file from the sswidl GOES widget output
    goes_trange: plot only the specified time range for goes
    timestr: display time as strings on X-axis -- currently the times do not update themselves when zooming in
    """
    # Set up variables 
    import matplotlib.pyplot as plt
    import numpy
    from numpy import log10
    from astropy.time import Time
    warnings.filterwarnings('ignore')
    if pol != 'RR' and pol != 'LL' and pol != 'RRLL' and pol != 'I' and pol != 'V' and pol != 'IV':
        print ("Please enter 'RR', 'LL', 'RRLL', 'I', 'V', 'IV' for pol")
        return 0

    if type(specdata) is str:
        if '.npz' in specdata:
            specdata = np.load(specdata)
        elif '.p' in specdata:
            open_dspec = open(specdata, 'rb')
            specdata = pickle.load(open_dspec,encoding='latin1')
            open_dspec.close()
            bl = ''
    elif type(specdata) is dict:
        specdata = specdata
    try:
        (npol, nbl, nfreq, ntim) = specdata['spec'].shape
        spec = specdata['spec']
        tim = specdata['tim']
        tim_ = Time(tim / 3600. / 24., format='mjd')
        #print(tim)
        tim_plt = tim_.plot_date
        freq = specdata['freq']
        print('here')
        if not 'bl' in vars():
            bl = specdata['bl']
            print(bl == '')
    except:
        print('format of specdata not recognized. Check your input')
        return -1

    if timerange:
        '''
        if type(timerange[0]) is str:
            #timerange = [qa.convert(qa.quantity(t), 's')['value'] for t in timerange]
            timerange_ = qa.convert_s(timerange)
            print(timerange_)
        tidx = np.where((tim >= float(timerange_[0])) & (tim <= float(timerange_[1])))[0]
        print(tidx)
        '''
        tidx = range(timerange[0], timerange[1])
    else:
        tidx = range(ntim)
    if freqrange:
        fidx = np.where((freq >= freqrange[0] * 1e6) & (freq <= freqrange[1] * 1e6))[0]
    else:
        fidx = range(nfreq)

    # setup plot parameters
    print ('ploting dynamic spectrum...')
    spec_med = np.median(np.absolute(spec))
    if not dmin:
        #dmin = spec_med / 20.
        dmin = np.min(np.absolute(spec))
    if not dmax:
        #dmax = spec_med * 2.
        dmax = np.max(np.absolute(spec))
    # do the plot
    for b in range(nbl):
        if pol != 'RRLL' and pol != 'IV':
            if pol == 'RR':
                spec_plt = spec[0, b, :, :]
            elif pol == 'LL':
                spec_plt = spec[1, b, :, :]
            elif pol == 'I':
                # spec_plt = (spec[0, b, :, :] + spec[1, b, :, :]) / 2.
                spec_plt = (spec[0, b, :, :] + spec[0, b, :, :]) / 2.
            elif pol == 'V':
                spec_plt = (spec[0, b, :, :] - spec[1, b, :, :]) / 2.
            if movie:
                '''
                    f = plt.figure(figsize=(16, 8), dpi=100)
                    if goessav:
                        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
                        gs.update(left=0.06, right=0.97, top=0.95, bottom=0.06)
                        ax1 = f.add_subplot(gs[0])
                        ax2 = f.add_subplot(gs[1])
                        if os.path.exists(goessav):
                            goes = readsav(goessav)
                            # IDL anytim 0 sec correspond to 1979 Jan 01, convert to mjd time
                            #anytimbase = qa.convert(qa.quantity('1979/01/01/00:00:00'), 's')['value']
                            anytimbase = qa.convert_s(['1979/01/01/00:00:00'])[0]
                            mjdbase = goes['utbase'] + anytimbase
                            ts = goes['tarray'] + mjdbase
                            lc0 = goes['yclean'][0, :]
                            lc1 = goes['yclean'][1, :]
                    else:
                        ax1 = f.add_subplot(211)
                    tstart = tim[tidx[0]]
                    tend = tim[tidx[-1]]
                    #tstartstr = qa.time(qa.quantity(tstart, 's'))[0]
                    #tendstr = qa.time(qa.quantity(tend, 's'))[0]
                    #nfrm = int((tend - tstart) / dtframe) + 1
                    #print ('Movie mode set. ' + str(nfrm) + ' frames to plot from ' + tstartstr + ' to ' + tendstr)
                    nfrm = int((tend - tstart) / dtframe) + 1
                    for i in range(nfrm):
                        if (i != 0) and (i % 10 == 0):
                            print (str(i) + ' frames done')
                        timeran = [tstart + i * dtframe, tstart + i * dtframe + framedur]
                        tidx1 = np.where((tim >= timeran[0]) & (tim <= timeran[1]))[0]
                        tim1 = tim_[tidx1]
                        freq1 = freq[fidx] / 1e9
                        spec_plt1 = spec_plt[fidx, :][:, tidx1]
                        ax1.pcolormesh(tim1.plot_date, freq1, spec_plt1, cmap='jet', vmin=dmin, vmax=dmax)
                        ax1.set_xlim(tim1[0].plot_date, tim1[-1].plot_date)
                        ax1.set_ylim(freq1[0], freq1[-1])
                        ax1.set_ylabel('Frequency (GHz)')
                        ax1.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + pol)
                        if timestr:
                            # date_format = mdates.DateFormatter('%H:%M:%S.%f')
                            # ax1.xaxis_date()
                            # ax1.xaxis.set_major_formatter(date_format)
                            locator = AutoDateLocator()
                            ax1.xaxis.set_major_locator(locator)
                            ax1.xaxis.set_major_formatter(AutoDateFormatter(locator))
                        ax1.set_autoscale_on(False)
                        if goessav:
                            if goes_trange:
                                if type(goes_trange[0]) is str:
                                    #goes_trange = [qa.convert(qa.quantity(t), 's')['value'] for t in goes_trange]
                                    goes_trange = qa.convert_s(goes_trange)
                                    idx = np.where((ts >= goes_trange[0]) & (ts <= goes_trange[1]))[0]
                            else:
                                idx = range(len(ts))
                            ts_plt = ts[idx]
                            lc0_plt = lc0[idx]
                            #utbase = qa.convert(qa.quantity('0001/01/01/00:00:00'), 'd')['value'] + 1
                            utbase = qa.convert_s(['0001/01/01/00:00:00'])[0]/24./3600.+1
                            ts_plt_d = ts_plt / 3600. / 24. - utbase
                            ax2.plot_date(ts_plt_d, lc0_plt, 'b-')
                            ax2.axvspan(tim1[0].mjd - utbase, tim1[-1].mjd - utbase, color='red',
                                        alpha=0.5)
                            ax2.set_yscale('log')
                            ax2.set_title('GOES 1-8 A')
    
                        tstartstr_ = tim1[0].datetime.strftime('%Y-%m-%dT%H%M%S.%f')[:-3]
                        tendstr_ = tim1[1].datetime.strftime('%H%M%S.%f')[:-3]
                        timstr = tstartstr_ + '-' + tendstr_
                        figfile = 'dspec_t' + timstr + '.png'
                        if not os.path.isdir('dspec'):
                            os.makedirs('dspec')
                        f.savefig('dspec/' + figfile)
                        plt.cla()
                        '''
                print('not avaliable now')
            else:
                if pax is None:
                    f = plt.figure(figsize=(10, 4), dpi=100)
                    ax = f.add_subplot(111)
                else:
                    ax = pax
                #freqghz = freq / 1.e9*10.0
                freqghz = freq
                # ax.ticklabel_format(style='sci',axis='y')
                if type(specdata) is not str:
                    kwargs = {'norm': SymLogNorm(linthresh=specdata['spec'][0, 0, 24, 124] * 2, vmin=dmin, vmax=dmax)}
                else:
                    kwargs = {}
                # ax.pcolormesh(tim_plt, freqghz, spec_plt, cmap='jet', vmin=dmin, vmax=dmax,**kwargs)
                print('here??????????????')
                pim = ax.pcolormesh(tim_plt, freqghz, spec_plt, cmap=custom_cm, vmin=dmin, vmax=dmax,zorder=1)
                tr_plt = Time(over_all_tr)
                #ax.set_xlim(tim_plt[tidx[0]], tim_plt[tidx[-1]])
                ax.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])
                ax.set_xlim(tr_plt.plot_date)
                formatter = mpl.dates.DateFormatter('%H:%M:%S')
                ax.xaxis.set_major_formatter(formatter)
                ax.fmt_xdata = formatter
                # ax.xaxis.set_major_locator(plt.MaxNLocator(5))
                # ax.colorbar()
                # plt.colorbar(pim,ax=ax,label='Amplitude [arb.units]')
                if time_line is not None:
                    if len(time_line) >= 1:
                        for cur_time in time_line:
                            if type(cur_time) == 'int':
                                cur_dic = '/Volumes/Data/20170906/info/all_time_20170906_10s.p'
                                in_dic_file = open(cur_dic, 'rb')
                                in_dic = pickle.load(in_dic_file,encoding='latin1')
                                in_dic_file.close()
                                cur_time = Time(in_dic[cur_time]['time'], format='iso')
                                ax.axvline(cur_time.plot_date, linewidth=2, color='b')
                            else:
                                ax.axvline(cur_time.plot_date, linewidth=2, color='b')

                #from sunpy import lightcurve
                from sunpy.time import TimeRange, parse_time
                import sunpy.timeseries
                import sunpy.data.sample
                goes = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)
                gax = ax.twinx()
                gax.ticklabel_format(axis='y', style='sci')
                t1 = tim_[tidx[0]]
                t2 = tim_[tidx[-1]]
                tr = TimeRange(t1.iso, t2.iso)
                #goes = lightcurve.GOESLightCurve.create(tr,overwrite=True)
                #goes = lightcurve.GOESLightCurve.create(tr)
                mjd_list = Time([str(ll) for ll in np.array(goes.data.index)], format='isot').mjd
                new_yy = np.zeros((len(tidx)))
                new_zz = np.zeros((len(tidx)))
                #if eovsa_lc:
                #    new_ee = np.zeros((len(tidx)))
                for ti in range(len(tidx)):
                    ctim = tim_[tidx[ti]].mjd
                    # new_yy[ti] = np.array(goes.data['xrsb'])[np.argmin(a=[abs(x-ctim) for x in mjd_list])]*1.e6
                    # new_zz[ti] = np.array(goes.data['xrsa'])[np.argmin(a=[abs(x-ctim) for x in mjd_list])]*1.e6
                    new_yy[ti] = np.array(goes.data['xrsb'])[
                                     np.argmin(a=[abs(x - ctim) for x in mjd_list])] / np.nanmax(
                        np.array(goes.data['xrsb']))
                    new_zz[ti] = np.array(goes.data['xrsa'])[
                                     np.argmin(a=[abs(x - ctim) for x in mjd_list])] / np.nanmax(
                        np.array(goes.data['xrsa']))
                if rhessisave is not None:
                    if '.sav' in rhessisave:
                        hessi = readsav(rhessisave)
                        hessi_tim = Time(list(hessi['obs_times_str']))
                        hessi_tim_plt = hessi_tim.plot_date
                        #tidx_hessi, = np.where((hessi_tim >= tim_[0]) & (hessi_tim <= tim_[1]))
                        tidx_hessi, = np.where((hessi_tim_plt >= tim_plt[0]) & (hessi_tim_plt <= tim_plt[-1]))
                        print(hessi_tim[471].isot)
                        print(tidx_hessi)
                        for idx, eg in enumerate(hessi['obs_energies']):
                            flux_hessi = ma.masked_array(hessi['obs_data'][tidx_hessi, idx])
                            #if hessi_smoth > 0:
                            if rhessisave is None:
                                pass
                            #    flux_hessi = DButil.smooth(flux_hessi, 20)
                            #    flux_hessi = flux_hessi / np.nanmax(flux_hessi)
                            #    ax.step(hessi_tim_plt[tidx_hessi], DButil.smooth(flux_hessi, hessi_smoth),
                            #            label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))
                            else:
                                if idx <= 3:
                                    print(flux_hessi)
                                    gax.step(hessi_tim_plt[tidx_hessi], flux_hessi,
                                            label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]),zorder=0)
                                else:
                                    gax.step(hessi_tim_plt[tidx_hessi], flux_hessi)
                    elif '.p' in rhessisave:
                        hessi_lc_comp = pickle.load(open('/Volumes/Data/20170820/hessi/lc_fitting/lc_components.p','rb'),encoding='latin1')
                        hessi_tim_plt = hessi_lc_comp[0]
                        tidx_hessi, = np.where((hessi_tim_plt >= tim_plt[0]) & (hessi_tim_plt <= tim_plt[-1]))
                        cur_hessi_max = 0.0
                        for hlci in range(5):
                            if np.nanmax(hessi_lc_comp[hlci + 1][tidx_hessi]) > cur_hessi_max:
                                cur_hessi_max = np.nanmax(hessi_lc_comp[hlci + 1][tidx_hessi])
                        for hlci in range(5):
                            if hlci==0:
                                cur_label = 'Total'
                            else:
                                #cur_label = 'Peak_{}'.format(hlci)
                                cur_label = 'Peaks'
                            if np.nanmax( hessi_lc_comp[hlci+1][tidx_hessi]) > cur_hessi_max:
                                cur_hessi_max = hessi_lc_comp[hlci+1][tidx_hessi]
                            if hlci in [0,4]:
                                gax.plot(hessi_tim_plt[tidx_hessi], hessi_lc_comp[hlci+1][tidx_hessi]/float(cur_hessi_max), label = cur_label,
                                         linewidth = 4)
                            else:
                                gax.plot(hessi_tim_plt[tidx_hessi],
                                         hessi_lc_comp[hlci + 1][tidx_hessi] / float(cur_hessi_max),linewidth=4)
                            #gax.plot(hessi_tim_plt[tidx_hessi], hessi_lc_comp[hlci + 1][tidx_hessi], label=cur_label)
                    #if ylog:
                    #    gax.set_yscale("log", nonposy='clip')
                    #else:
                    gax.set_yscale("linear", nonposy='clip')
                    gleg = gax.legend(fontsize='8')
                    # ax.legend(prop={'size': 7})
                    gax.set_ylabel('Count Rate\n [s$^{-1}$ detector$^{-1}$ ]')
                    #gax.set_xlim(tim_plt[tidx[0]], tim_plt[tidx[-1]])
                    gax.set_xlim(tr_plt.plot_date)
                    #ax.xaxis.set_visible(False)
                    '''
                    rhe_sav = readsav(rhessisave)
                    rhessi_mjd_list = Time([str(ll) for ll in rhe_sav['obs_times_str']], format='isot').mjd
                    rhe_yy = np.zeros((len(tidx)))
                    rhe_zz = np.zeros((len(tidx)))
                    for ti in range(len(tidx)):
                        ctim = tim_[tidx[ti]].mjd
                        # new_yy[ti] = np.array(goes.data['xrsb'])[np.argmin(a=[abs(x-ctim) for x in mjd_list])]*1.e6
                        # new_zz[ti] = np.array(goes.data['xrsa'])[np.argmin(a=[abs(x-ctim) for x in mjd_list])]*1.e6
                        rhe_yy[ti] = rhe_sav['obs_data'][:, 1][
                                         np.argmin(a=[abs(x - ctim) for x in rhessi_mjd_list])] / np.nanmax(
                            rhe_sav['obs_data'][:, 1])
                        #rhe_zz[ti] = rhe_sav['obs_data'][:, 2][
                        #                 np.argmin(a=[abs(x - ctim) for x in rhessi_mjd_list])] / np.nanmax(
                        #    rhe_sav['obs_data'][:, 2])
                    gax.plot(tim_plt[tidx[0]:tidx[-1] + 1], rhe_yy, c='green', label='RHESSI 6$-$13 keV')
                    '''
                if customized_light_curve:
                    #aia_lc = ep.aia_lightcurve_from_stp('n',yrange = [0,20],timerange=timerange,pax = gax,
                    #                                    normq=True, ltex = 'LC of the N ribbon',cur_marker='-.')
                    aia_94_lc = pt.lightcurve(timerange=[184,350],x=-940, y=105, kw_list=['Z.94.'],sample_range=[25,25],
                                              derivative=True, pax=gax,normq=False,ltex='AIA94 LC',cur_legend = gleg,
                                              blooming_cut_off=False,cur_color='g')
                    aia_131_lc = pt.lightcurve(timerange=[184, 350], x=-940, y=105, kw_list=['Z.131.'], sample_range=[25, 25],
                                              derivative=True, pax=gax, normq=False, ltex='AIA131 LC', cur_legend=gleg,
                                              blooming_cut_off=True,cmarker='-')
                    aia_193_lc = pt.lightcurve(timerange=[184, 350], x=-940, y=105, kw_list=['Z.193.'], sample_range=[25, 25],
                                              derivative=True, pax=gax, normq=False, ltex='AIA193 LC', cur_legend=gleg,
                                              blooming_cut_off=True,cur_color='y',cmarker='-')
                    #gax.set_ylabel('Intensity Derivative')#
                    gax.set_ylabel('Normalized Intensity')#
                    #gax.set_ylim(-50, 100)
                    #ep.add_patch_to_existing_legend(legend=gleg,label = 'AIA94 LC_DER', marker='_')

                    #gax.plot(tim_plt[tidx[0]:tidx[-1] + 1], rhe_zz, c='blue', label='RHESSI 12-25 KeV')

                # goes.data['xrsb'] = 2 * (np.log10(goes.data['xrsb'])) + 26
                # print(goes.data.index)
                # xx = [str(ll) for ll in np.array(goes.data.index)]
                # print(len(tim_plt),len(xx))
                # yy = np.array(goes.data['xrsb'])
                # print(yy[0:10])
                # gax.plot(Time(xx).mjd * 24 * 3600, yy, c='yellow')
                #gax.plot(tim_plt[tidx[0]:tidx[-1] + 1], new_yy, c='k', label='GOES 1.0$-$8.0 $\AA$')
                #gax.plot(tim_plt[tidx[0]:tidx[-1] + 1], new_zz, c='', label='GOSE 0.5--4.0 $\AA$')
                #gax.set_ylim(0.0, 1.3)

                # gax.xaxis.set_major_locator(plt.MaxNLocator(5))
                # plt.colorbar()
                # rightaxis_label_time = Time(xx[-1]).mjd * 24 * 3600
                # print('here we go')
                # gax.text(rightaxis_label_time, 9.6, 'A', fontsize='15')
                # gax.text(rightaxis_label_time, 11.6, 'B', fontsize='15')
                # gax.text(rightaxis_label_time, 13.6, 'C', fontsize='15')
                # gax.text(rightaxis_label_time, 15.6, 'M', fontsize='15')
                # gax.text(rightaxis_label_time, 17.6, 'X', fontsize='15')
                font = {'weight': 'bold',
                        'size': 12}
                mpl.rc('font', **font)
                mpl.rcParams.update({'font.size': 12})
                mpl.rc('font', **font)
                def format_coord(x, y):
                    col = np.argmin(np.absolute(tim_plt - x))
                    row = np.argmin(np.absolute(freqghz - y))
                    if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                        timstr = tim_[col].isot
                        flux = spec_plt[row, col]
                        return 'time {0} = {1}, freq = {2:.3f} GHz, flux = {3:.2f} Jy'.format(col, timstr, y, flux)
                    else:
                        return 'x = {0}, y = {1:.3f}'.format(x, y)

                ax.format_coord = format_coord
                ax.set_ylabel('Frequency [GHz]')
                # gax.set_ylabel('10$^{-6}$ Watts m$^{-2}$')
                gax.set_ylabel('Normalized\n Intensity')
                gax.set_xlabel('Time [UT]')
                #ax.xaxis.set_tick_params(rotation=30)
                gax.yaxis.label.set_size(12)
                gax.tick_params(axis='y', labelsize=12)
                #gax.set_ylabel('Intensity Derivative')
                #gax.set_ylim(-50, 100)
                # labels = [item.get_text() for item in gax.get_yticklabels()]
                # labels[5] = 'M Class'
                gax.ticklabel_format(axis='y', style='sci')
                if isinstance(leg_loc, int):
                    gax.legend(loc=leg_loc,fontsize='10')
                else:
                    #gax.legend(loc='best',bbox_to_anchor=leg_loc,fontsize='10')
                    gax.legend(loc='best', fontsize='10')
                if bl:
                    ax.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + pol)
                else:
                    # ax.set_title('Medium dynamic spectrum')
                    ax.set_title('')
                if timestr:
                    # date_format = mdates.DateFormatter('%H:%M:%S.%f')
                    # ax.xaxis_date()
                    # ax.xaxis.set_major_formatter(date_format)
                    locator = AutoDateLocator()
                    ax.xaxis.set_major_locator(locator)
                    ax.xaxis.set_major_formatter(AutoDateFormatter(locator))
                ax.set_autoscale_on(False)
                formatter = mpl.dates.DateFormatter('%H:%M:%S')
                ax.xaxis.set_major_formatter(formatter)
                ax.fmt_xdata = formatter
        else:
            f = plt.figure(figsize=(8, 6), dpi=100)
            R_plot = np.absolute(spec[0, b, :, :])
            L_plot = np.absolute(spec[1, b, :, :])
            I_plot = (R_plot + L_plot) / 2.
            V_plot = (R_plot - L_plot) / 2.
            if pol == 'RRLL':
                spec_plt_1 = R_plot
                spec_plt_2 = L_plot
                polstr = ['RR', 'LL']
            if pol == 'IV':
                spec_plt_1 = I_plot
                spec_plt_2 = V_plot
                polstr = ['I', 'V']

            ax1 = f.add_subplot(211)
            freqghz = freq / 1e9
            pcol=ax1.pcolormesh(tim_plt, freqghz, spec_plt_1, cmap='jet', vmin=dmin, vmax=dmax,linewidth=0,rasterized=True,antialiased=True)
            pcol.set_edgecolor('face')
            ax1.set_xlim(tim_plt[tidx[0]], tim_plt[tidx[-1]])
            ax1.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])

            def format_coord(x, y):
                col = np.argmin(np.absolute(tim_plt - x))
                row = np.argmin(np.absolute(freqghz - y))
                if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                    timstr = tim_[col].isot
                    flux = spec_plt[row, col]
                    return 'time {0} = {1}, freq = {2:.3f} GHz, flux = {3:.2f} Jy'.format(col, timstr, y, flux)
                else:
                    return 'x = {0}, y = {1:.3f}'.format(x, y)

            ax1.format_coord = format_coord
            ax1.set_ylabel('Frequency [GHz]')
            if timestr:
                # date_format = mdates.DateFormatter('%H:%M:%S.%f')
                # ax1.xaxis_date()
                # ax1.xaxis.set_major_formatter(date_format)
                locator = AutoDateLocator()
                ax1.xaxis.set_major_locator(locator)
                ax1.xaxis.set_major_formatter(AutoDateFormatter(locator))
            ax1.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + polstr[0])
            #ax1.tick_params(axis='y',direction='in')
            ax1.set_autoscale_on(False)
            ax2 = f.add_subplot(212)
            #ax2.pcolormesh(tim_plt, freqghz, spec_plt_2, cmap='jet', vmin=dmin, vmax=dmax)
            #ax2.set_xlim(tim_plt[tidx[0]], tim_plt[tidx[-1]])
            #ax2.set_ylim(freqghz[fidx[0]], freqghz[fidx[-1]])
            if timestr:
                # date_format = mdates.DateFormatter('%H:%M:%S.%f')
                # ax2.xaxis_date()
                # ax2.xaxis.set_major_formatter(date_format)
                locator = AutoDateLocator()
                ax2.xaxis.set_major_locator(locator)
                ax2.xaxis.set_major_formatter(AutoDateFormatter(locator))

            def format_coord(x, y):
                col = np.argmin(np.absolute(tim_plt - x))
                row = np.argmin(np.absolute(freqghz - y))
                if col >= 0 and col < ntim and row >= 0 and row < nfreq:
                    timstr = tim_[col].isot
                    flux = spec_plt[row, col]
                    return 'time {0} = {1}, freq = {2:.3f} GHz, flux = {3:.2f} Jy'.format(col, timstr, y, flux)
                else:
                    return 'x = {0}, y = {1:.3f}'.format(x, y)

            ax2.format_coord = format_coord
            ax2.set_ylabel('Frequency (GHz)')
            ax2.set_title('Dynamic spectrum @ bl ' + bl.split(';')[b] + ', pol ' + polstr[1])
            #ax2.tick_params(axis='y',direction='in')
            ax2.set_autoscale_on(False)
    font = {'weight': 'bold',
        'size': 12}
    mpl.rc('font', **font)
    mpl.rcParams.update({'font.size': 12})
    mpl.rc('font', **font)
    if pax is None:
        plt.show()

'''
def wrt_dspec(specfile=None, specdat=None):
    try:
        specfile
    except NameError:
        print 'No input centfile specified for reading. Aborting...'
    if not specdat:
        print 'Output file name is not specified, use the default convention'
        specdat = specfile.replace('npz', 'dat')
    specdata = np.load(specfile)
    spec = np.copy(specdata['spec'][:, :, :, :])
    npl, nbl, nf, nt = spec.shape
    print 'Dimension of the data cube -- # of pol, # of baseline, # of frequency, # of time:'
    print npl, nbl, nf, nt
    nele = npl * nbl * nf * nt
    # need to transpose to nf, nt, npl, nbl
    spec = np.transpose(spec, (2, 3, 0, 1))
    # Turn fmx into a 1d array for writing
    spec = spec.flatten()
    # Pack the data into a byte buffer.  Note the asterisks in the next three lines!
    # If your data includes integers or other data types, use 'i' or other instead of 'f'
    buf = struct.pack(str(nf) + 'f', *specdata['freq'])
    buf += struct.pack(str(nt) + 'd', *specdata['tim'])
    buf += struct.pack(str(nele) + 'f', *spec)
    with open(specdat, 'wb') as f:
        f.write(buf)
    f.close()
'''

#def main():
#    plt_dspec(specdata='/Volumes/Data/20170820/eovsa/image_dspec/part/10s_dspec_960_94_926_136.p')
#if __name__ == '__main__':
#    main()
