def lightcurves(timerange, outdir='./', specfile=None, goes=True, hessifile=None, fermifile=None, ylog=False, hessi_smoth=0, dspec_cmap='cubehelix',
                vmax=None, vmin=None,custom_ax=None,only_dspec=True,nde=2, goes_derivative=False):
    #from sunpy.lightcurve import GOESLightCurve
    from sunpy.time import TimeRange, parse_time
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from astropy.time import Time
    import numpy as np
    import numpy.ma as ma
    from scipy.signal import medfilt
    import matplotlib.colors as colors
    from scipy.io import readsav
    from md_utils import DButil
    import os
    import copy
    #from taskinit import qa
    import ex_plot_1 as ep
    import pickle
    import sys

    timerange = Time(timerange)
    if hessifile:
        if not os.path.exists(hessifile):
            hessi_script = 'HESSI_lc.pro'
            print('Run the script {} in SSWIDL to download RHESSI summary file first!'.format(hessi_script))
            fi = open(hessi_script, 'wb')
            fi.write(b'search_network, /enable \n')
            text = "time_range = ['{}','{}'] \n".format(timerange[0].datetime.strftime('%Y-%b-%d %H:%M:%S'),
                                                          timerange[1].datetime.strftime('%Y-%b-%d %H:%M:%S'))

            fi.write(text.encode())
            fi.write(b'obs_obj = hsi_obs_summary(obs_time_interval=time_range) \n')
            fi.write(b'data = obs_obj->getdata() \n')
            fi.write(b'info = obs_obj->get( / info) \n')
            fi.write(b'obs_data = data.countrate \n')
            fi.write(b'obs_times = obs_obj->getdata(/time) \n')
            fi.write(b'obs_times_str  = anytim(obs_times,/CCSDS) \n')
            fi.write(b'obs_energies  = fltarr(2,n_elements(info.energy_edges)-1) \n')
            fi.write(b'for ll=0,n_elements(info.energy_edges)-2 do begin \n')
            fi.write(b'    obs_energies[0,ll]=info.energy_edges[ll] \n')
            fi.write(b'    obs_energies[1,ll]=info.energy_edges[ll+1] \n')
            fi.write(b'endfor \n')
            text = 'save,filename="{}",OBS_DATA,OBS_ENERGIES,obs_times_str \n'.format(hessifile)
            fi.write(text.encode())
            fi.write(b'end \n')
            fi.close()
            return -1
        hessi = readsav(hessifile)
        hessi_tim = Time(list(hessi['obs_times_str']))
        hessi_tim_plt = hessi_tim.plot_date

    if specfile[-2:] == '.p':
        with open(specfile, 'rb') as spec_open:
           specdata = pickle.load(spec_open)
        spec_open.close()
    else:
        specdata = np.load(specfile)
    #with open(specfile, 'rb') as spec_open:
    #    specdata = pickle.load(spec_open)
    #spec_open.close()
    spec = specdata['spec']
    if len(spec.shape) == 4:
        (npol, nbl, nfreq, ntim) = spec.shape
    else:
        (nfreq, ntim) = spec.shape
    spec = np.mean(np.mean(spec, axis=0), axis=0)
    freq = specdata['freq']
    freqghz = freq / 1e9
    spec_tim = Time(specdata['tim'] / 3600. / 24., format='mjd')

    #fidx_plt = np.linspace(0, nfreq - nfreq / 8.0 / 2.0, 8).astype(np.int) + nfreq / 8.0 / 2.0
    fidx_plt = np.linspace(0, nfreq - nfreq / 4.0 / 2.0, 4).astype(np.int) + nfreq / 4.0 / 2.0

    try:
        pass
        #plt.style.use('seaborn-bright')
        #params = {'font.size': 8, 'axes.grid': False, 'axes.facecolor': 'w', 'xtick.color': '#555555', 'ytick.color': '#555555',
        #          'xtick.major.size': 2.0, 'ytick.major.size': 2.0, 'xtick.minor.size': 1.0, 'ytick.minor.size': 1.0, 'axes.axisbelow': False,
        #          'axes.xmargin': 0.0, 'axes.ymargin': 0.0, 'axes.linewidth': 0.5, 'xtick.major.pad': 1.5, 'ytick.major.pad': 1.5,
        #          'lines.linewidth': 1.0}
        #mpl.rcParams.update(params)
    except:
        pass
    '''
    if goes:
        #try:
        #    pass
        #    from sunpy import lightcurve as lc
        #    from sunpy.time import TimeRange
        #    tr = TimeRange(timerange.iso)
        #    goesd = lc.GOESLightCurve.create(tr)
        #    dates = mpl.dates.date2num(parse_time(goesd.data.index))
        #except:
        goesdatafile = '/Volumes/Data/20170820/goes/goes.dat'
        if not os.path.exists(goesdatafile):
            workdir = '/Volumes/Data/20170820/goes/'
            btgoes = Time(timerange[0]).iso
            etgoes = Time(timerange[1]).iso
            goesscript = os.path.join(workdir, 'goes.py')
            goesdatafile = os.path.join(workdir, 'goes.dat')
            gt1 = mpl.dates.date2num(parse_time(btgoes))
            gt2 = mpl.dates.date2num(parse_time(etgoes))
            os.system('rm -rf {}'.format(goesscript))
            os.system('rm -rf {}'.format(goesdatafile))
            fi = open(goesscript, 'wb')
            fi.write('import os \n')
            fi.write('import shutil \n')
            fi.write('from astropy.io import fits \n')
            fi.write('import pickle \n')
            fi.write('import urllib2 \n')
            fi.write('from astropy.time import Time \n')
            fi.write('import numpy as np \n')
            fi.write('import ssl \n')
            fi.write('goesplottim = ["{0}", "{1}",{2}, {3}] \n'.format(btgoes, etgoes, gt1, gt2))
            fi.write('yr = goesplottim[0][:4] \n')
            fi.write('datstr = goesplottim[0][:4]+goesplottim[0][5:7]+goesplottim[0][8:10] \n')
            fi.write('context = ssl._create_unverified_context() \n')
            fi.write('f = urllib2.urlopen(\'https://umbra.nascom.nasa.gov/goes/fits/\'+yr, context=context) \n')
            fi.write('lines = f.readlines() \n')
            fi.write('sat_num = [] \n')
            fi.write('for line in lines: \n')
            fi.write('    idx = line.find(datstr) \n')
            fi.write('    if idx != -1: \n')
            fi.write('        sat_num.append(line[idx-2:idx]) \n')
            fi.write('if type(sat_num) is int: \n')
            fi.write('    sat_num = [str(sat_num)] \n')
            fi.write('filenames = [] \n')
            fi.write('for sat in sat_num: \n')
            fi.write('    filename = \'go\'+sat+datstr+\'.fits\' \n')
            fi.write('    url = \'https://umbra.nascom.nasa.gov/goes/fits/\'+yr+\'/\'+filename \n')
            fi.write('    f = urllib2.urlopen(url, context=context) \n')
            fi.write('    workdir=os.getcwd() \n')
            fi.write('    with open(workdir+\'/\'+filename,\'wb\') as g: \n')
            fi.write('        shutil.copyfileobj(f,g) \n')
            fi.write('    filenames.append(workdir+\'/\'+filename) \n')
            fi.write('pmerit = 0 \n')
            fi.write('for file in filenames: \n')
            fi.write('    gfits = fits.open(file) \n')
            fi.write('    tsecs = gfits[2].data[\'TIME\'][0] \n')
            fi.write('    date_elements = gfits[0].header[\'DATE-OBS\'].split(\'/\') \n')
            fi.write(
                '    temp_t = Time(date_elements[2]+\'-\'+date_elements[1]+\'-\'+date_elements[0]).plot_date + tsecs/86400. \n')
            fi.write('    i=0 \n')
            fi.write('    while temp_t[i]<goesplottim[2]: \n')
            fi.write('        i=i+1 \n')
            fi.write('    st=i-1 \n')
            fi.write('    while temp_t[i]<goesplottim[3]: \n')
            fi.write('        i=i+1 \n')
            fi.write('    et=i \n')
            fi.write('    data = gfits[2].data[\'FLUX\'][0][:,0][st:et] \n')
            fi.write('    good = np.where(data > 1.e-8) \n')
            fi.write('    merit = len(good) \n')
            fi.write('    if merit > pmerit: \n')
            fi.write('        pmerit = merit \n')
            fi.write('        goes_data = data \n')
            fi.write('        goes_t = temp_t[st:et] \n')
            fi.write('goes={\'time\':goes_t,\'xrsb\':data} \n')
            fi.write('fi2 = open("{}", "wb") \n'.format(goesdatafile))
            fi.write('pickle.dump(goes, fi2) \n')
            fi.write('fi2.close()')
            fi.close()

            try:
                os.system('python {}'.format(goesscript))
                #os.system('rm -rf {}'.format(goesscript))
            except NameError:
                print("Bad input names"
            except ValueError:
                print("Bad input values"
            except:
                print("Unexpected error:", sys.exc_info()[0]
                print("Error in generating GOES light curves. Proceed without GOES..."

        if os.path.exists(goesdatafile):
            fi1 = file(goesdatafile, 'rb')
            goesd = pickle.load(fi1)
            fi1.close()
            dates = goesd['time']
            goesdata = goesd['xrsb']
    #print(goesdata)
    '''
    if custom_ax is None:
        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(10, 6))
    else:
        axs=custom_ax
    #ax = axs[1]
    ax = axs[0]
    #gax = axs[0].twinx()
    if goes:
        goessav = '/Volumes/Data/20170820/20220511/goes/idlsave_goes.sav'
        # if goes is in sav file
        if os.path.exists(goessav):
            goes = readsav(goessav)
        #utbase = qa.convert(qa.quantity('0001/01/01/00:00:00'), 'd')['value'] + 1
        #utbase = Time('0001-01-01 00:00:00.000',format='iso').mjd * 3600.0*24.0
        anytimbase = Time('1979-01-01 00:00:00.000',format='iso').mjd*3600.0*24.0
        mjdbase = goes['utbase'] + anytimbase
        #mjdbase = goes['utbase']
        ts = goes['tarray'] + mjdbase
        #ts = goes['tarray']
        lc0 = goes['yclean'][0, :]
        lc1 = goes['yclean'][1, :]
        tr_plt = Time(timerange)
        idx = np.where((Time(ts/ 3600. / 24., format='mjd') >= tr_plt[0]) & ( Time(ts/ 3600. / 24., format='mjd')<= tr_plt[1]))
        ts_plt = ts[idx]
        lc0_plt = lc0[idx]
        lc1_plt = lc1[idx]
        ts_plt_d = Time(ts_plt / 3600. / 24., format='mjd').plot_date
        #print(type(ts_plt_d))
        #tidx_goes, = np.where((ts_plt_d >= tr_plt[0]) & (ts_plt_d <= tr_plt[1]))
        if goes_derivative:
            lc0_plt_dr = copy.deepcopy(lc0_plt)
            lc1_plt_dr = copy.deepcopy(lc1_plt)
            lc0_plt_dr = np.gradient(lc0_plt_dr)
            lc1_plt_dr = np.gradient(lc1_plt_dr)
            lc0_plt_dr = DButil.smooth(lc0_plt_dr, 20)
            lc1_plt_dr = DButil.smooth(lc1_plt_dr, 20)
            ax.plot_date(ts_plt_d, lc0_plt_dr, 'b--', label='GOES 1 - 8A deri')
            ax.plot_date(ts_plt_d, lc1_plt_dr, 'r--', label='GOES 0.5 - 4A deri')
            ax.set_ylim(-0.3e-7, 1.e-7)
        else:
            ax.plot_date(ts_plt_d, lc0_plt, 'b--',label='GOES 1 - 8A')
            ax.plot_date(ts_plt_d, lc1_plt, 'r--',label='GOES 0.5 - 4A')
        ax.legend(fontsize='8', loc='upper left')

        #print(Time(ts_plt_d[100],format='mjd').iso)
        #print(lc0_plt[100])
        #print(Time(mjdbase,format='mjd').iso)
    # tidx_hessi, = np.where((hessi_tim >= tr_plt[0]) & (hessi_tim <= tr_plt[1]))
    # for idx, eg in enumerate(hessi['obs_energies']):
    #     flux_hessi = ma.masked_array(hessi['obs_data'][tidx_hessi, idx])
    #     if hessi_smoth > 0:
    #         flux_hessi = DButil.smooth(flux_hessi, 20)
    #         flux_hessi = flux_hessi / np.nanmax(flux_hessi)
    #         ax.step(hessi_tim_plt[tidx_hessi], DButil.smooth(flux_hessi, hessi_smoth), label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))
    #     else:
    #         if idx <= 3:
    #             ax.step(hessi_tim_plt[tidx_hessi], flux_hessi, label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))
    #         else:
    #             ax.step(hessi_tim_plt[tidx_hessi], flux_hessi)
    if ylog:
        ax.set_yscale("log", nonposy='clip')
    else:
        ax.set_yscale("linear", nonposy='clip')
    ax.legend(fontsize='8')
    #ax.legend(prop={'size': 7})
    ax.set_ylabel('Count Rate\n [s$^{-1}$ detector$^{-1}$ ]')
    ax.xaxis.set_visible(False)
    #ax = axs[0]
    #with open('/Volumes/WD6T/working/20170703/fermi/1s_tres_dect{0}.p'.format(nde),'rb') as fdspec_file:
    #    fdspec = pickle.load(fdspec_file)
    #fdspec_file.close()
    #ep.plot_fermi_dynamic_spec_in_ax(cax=ax, spec_dict=fdspec, timerange=timerange,erg_max=100.0,nde=nde)
    #ax.xaxis.set_visible(False)
    #ax = axs[2]
    ax = axs[1]
    tidx_spec, = np.where((spec_tim >= tr_plt[0]) & (spec_tim <= tr_plt[1]))
    spec_tim_plt = spec_tim[tidx_spec[0]:tidx_spec[-1]].plot_date
    spec_plt = spec[:, tidx_spec[0]:tidx_spec[-1]]
    ax.pcolormesh(spec_tim_plt, freqghz, spec_plt, cmap=dspec_cmap, norm=colors.LogNorm(vmax=vmax, vmin=vmin))
    #print(freqghz)
    ax.set_ylabel('Frequency [GHz]')
    formatter = mpl.dates.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(formatter)
    ax.fmt_xdata = formatter
    axs[1].set_xlim(tr_plt.plot_date)
    # ax.set_ylim([0.0,20.0])
    ax.set_xlabel('Start time ({})'.format(tr_plt[0].datetime.strftime('%d-%b-%y %H:%M:%S')))
    #gax =ax.twinx()
    #tidx_goes, = np.where((dates >= tr_plt[0].plot_date) & (dates <= tr_plt[1].plot_date))
    #gax.plot(dates[tidx_goes], goesd.data['xrsb'][tidx_goes] / np.nanmax(goesd.data['xrsb'][tidx_goes]),
    #gax.plot(dates[tidx_goes], goesd['xrsb'][tidx_goes] / np.nanmax(goesd['xrsb'][tidx_goes]),
    #                  label='GOES 1.0$-$8.0 $\AA$')
        #ax.plot(dates[tidx_goes], goes.data['xrsa'][tidx_goes] / np.nanmax(goes.data['xrsa'][tidx_goes]), label='GOES 0.5--4.0 $\AA$')
    #tidx_spec, = np.where((spec_tim >= tr_plt[0]) & (spec_tim <= tr_plt[1]))
    '''
    if len(tidx_spec) > 1:
        spec_tim_plt = spec_tim[tidx_spec[0]:tidx_spec[-1]].plot_date
        if not only_dspec:
            flux_colors = []
            for idx, fidx in enumerate(np.round(fidx_plt).astype(np.int)):
                flux_plt = medfilt(spec[min(fidx, nfreq - 1), tidx_spec[0]:tidx_spec[-1]], 7)
                #flux_plt = medfilt(spec[min(fidx, nfreq - 1), tidx_spec[0]:tidx_spec[-1]], 5)
                p = ax.plot(spec_tim_plt, flux_plt / np.nanmax(flux_plt), label='{:.2f} GHz'.format(freqghz[min(fidx, nfreq - 1)]))
                flux_colors.append(p[0].get_color())
            ax.set_ylabel('Flux (Normalized)')
            #axs[2].set_ylim(0, 1.1)
            ax.legend(prop={'size': 7})
            #!!!!!!!!!!!!!!!!!---------------------------
            #only if horizontal
            #formatter = mpl.dates.DateFormatter('%H:%M:%S')
            #ax.xaxis.set_major_formatter(formatter)
            #ax.fmt_xdata = formatter
            #ax.set_xlim(tr_plt.plot_date)
            #ax.set_xlabel('Start time ({})'.format(tr_plt[0].datetime.strftime('%d-%b-%y %H:%M:%S')))
        ax.xaxis.set_visible(False)

        #ax = axs[3]
        ax = axs[2]
        spec_plt = spec[:, tidx_spec[0]:tidx_spec[-1]]
        ax.pcolormesh(spec_tim_plt, freqghz, spec_plt, cmap=dspec_cmap, norm=colors.LogNorm(vmax=vmax, vmin=vmin))
        print(freqghz)
        ax.set_ylabel('Frequency [GHz]')
        formatter = mpl.dates.DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(formatter)
        ax.fmt_xdata = formatter
        axs[2].set_xlim(tr_plt.plot_date)
        #ax.set_ylim([0.0,20.0])
        ax.set_xlabel('Start time ({})'.format(tr_plt[0].datetime.strftime('%d-%b-%y %H:%M:%S')))
        if not only_dspec:
            for idx, fidx in enumerate(np.round(fidx_plt).astype(np.int)):
                ax.axhline(freqghz[min(fidx, nfreq - 1)], color=flux_colors[idx], ls=':')
    else:
        print('Warning: No radio data in the timerange. Proceed without dynamic spectrum.')
    '''

    if custom_ax is None:
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.06)
        plt.show()
    else:
        return axs
    #imgdir = outdir + '/fig01_{}-{}.png'.format(tr_plt[0].datetime.strftime('%H%M%S'), tr_plt[1].datetime.strftime('%H%M%S'))
    #fig.savefig(imgdir, dpi=200)
    #print('Save image to ' + imgdir)

def lightcurves_sxr_mw(timerange, outdir='./', specfile=None, goes=True, hessifile=None, fermifile=None, ylog=False, hessi_smoth=0, dspec_cmap='cubehelix',
                vmax=None, vmin=None,custom_ax=None,only_dspec=False,nde=2):
    from sunpy.lightcurve import GOESLightCurve
    from sunpy.time import TimeRange, parse_time
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from astropy.time import Time
    import numpy as np
    import numpy.ma as ma
    from scipy.signal import medfilt
    import matplotlib.colors as colors
    from scipy.io import readsav
    from md_utils import DButil
    import os
    import ex_plot_1 as ep
    import pickle

    timerange = Time(timerange)
    specdata = np.load(specfile)
    #with open(specfile, 'rb') as spec_open:
    #    specdata = pickle.load(spec_open)
    #spec_open.close()
    spec = specdata['spec']
    if len(spec.shape) == 4:
        (npol, nbl, nfreq, ntim) = spec.shape
    else:
        (nfreq, ntim) = spec.shape
    spec = np.mean(np.mean(spec, axis=0), axis=0)
    freq = specdata['freq']
    freqghz = freq / 1e9
    spec_tim = Time(specdata['tim'] / 3600. / 24., format='mjd')

    #fidx_plt = np.linspace(0, nfreq - nfreq / 8.0 / 2.0, 8).astype(np.int) + nfreq / 8.0 / 2.0
    fidx_plt = np.linspace(0, nfreq - nfreq / 4.0 / 2.0, 4).astype(np.int) + nfreq / 4.0 / 2.0

    try:
        pass
        #plt.style.use('seaborn-bright')
        #params = {'font.size': 8, 'axes.grid': False, 'axes.facecolor': 'w', 'xtick.color': '#555555', 'ytick.color': '#555555',
        #          'xtick.major.size': 2.0, 'ytick.major.size': 2.0, 'xtick.minor.size': 1.0, 'ytick.minor.size': 1.0, 'axes.axisbelow': False,
        #          'axes.xmargin': 0.0, 'axes.ymargin': 0.0, 'axes.linewidth': 0.5, 'xtick.major.pad': 1.5, 'ytick.major.pad': 1.5,
        #          'lines.linewidth': 1.0}
        #mpl.rcParams.update(params)
    except:
        pass

    if goes:
        #tr = TimeRange(timerange.iso)
        #goes = GOESLightCurve.create(tr)
        #dates = mpl.dates.date2num(parse_time(goes.data.index))

        try:
            from sunpy import lightcurve as lc
            from sunpy.time import TimeRange
            tr = TimeRange(timerange.iso)
            goest = lc.GOESLightCurve.create(tr)
            dates = mpl.dates.date2num(parse_time(goest.data.index))
        except:
            workdir = '/Volumes/WD6T/working/20210508/goes/'
            btgoes = Time(timerange[0]).iso
            etgoes = Time(timerange[1]).iso
            goesscript = os.path.join(workdir, 'goes.py')
            goesdatafile = os.path.join(workdir, 'goes.dat')
            gt1 = mpl.dates.date2num(parse_time(btgoes))
            gt2 = mpl.dates.date2num(parse_time(etgoes))
            os.system('rm -rf {}'.format(goesscript))
            os.system('rm -rf {}'.format(goesdatafile))
            fi = open(goesscript, 'wb')
            fi.write('import os \n')
            fi.write('import shutil \n')
            fi.write('from astropy.io import fits \n')
            fi.write('import pickle \n')
            fi.write('import urllib2 \n')
            fi.write('from astropy.time import Time \n')
            fi.write('import numpy as np \n')
            fi.write('import ssl \n')
            fi.write('goesplottim = ["{0}", "{1}",{2}, {3}] \n'.format(btgoes, etgoes, gt1, gt2))
            fi.write('yr = goesplottim[0][:4] \n')
            fi.write('datstr = goesplottim[0][:4]+goesplottim[0][5:7]+goesplottim[0][8:10] \n')
            fi.write('context = ssl._create_unverified_context() \n')
            fi.write('f = urllib2.urlopen(\'https://umbra.nascom.nasa.gov/goes/fits/\'+yr, context=context) \n')
            fi.write('lines = f.readlines() \n')
            fi.write('sat_num = [] \n')
            fi.write('for line in lines: \n')
            fi.write('    idx = line.find(datstr) \n')
            fi.write('    if idx != -1: \n')
            fi.write('        sat_num.append(line[idx-2:idx]) \n')
            fi.write('if type(sat_num) is int: \n')
            fi.write('    sat_num = [str(sat_num)] \n')
            fi.write('filenames = [] \n')
            fi.write('for sat in sat_num: \n')
            fi.write('    filename = \'go\'+sat+datstr+\'.fits\' \n')
            fi.write('    url = \'https://umbra.nascom.nasa.gov/goes/fits/\'+yr+\'/\'+filename \n')
            fi.write('    f = urllib2.urlopen(url, context=context) \n')
            fi.write('    workdir=os.getcwd() \n')
            fi.write('    with open(workdir+\'/\'+filename,\'wb\') as g: \n')
            fi.write('        shutil.copyfileobj(f,g) \n')
            fi.write('    filenames.append(workdir+\'/\'+filename) \n')
            fi.write('pmerit = 0 \n')
            fi.write('for file in filenames: \n')
            fi.write('    gfits = fits.open(file) \n')
            fi.write('    tsecs = gfits[2].data[\'TIME\'][0] \n')
            fi.write('    date_elements = gfits[0].header[\'DATE-OBS\'].split(\'/\') \n')
            fi.write(
                '    temp_t = Time(date_elements[2]+\'-\'+date_elements[1]+\'-\'+date_elements[0]).plot_date + tsecs/86400. \n')
            fi.write('    i=0 \n')
            fi.write('    while temp_t[i]<goesplottim[2]: \n')
            fi.write('        i=i+1 \n')
            fi.write('    st=i-1 \n')
            fi.write('    while temp_t[i]<goesplottim[3]: \n')
            fi.write('        i=i+1 \n')
            fi.write('    et=i \n')
            fi.write('    data = gfits[2].data[\'FLUX\'][0][:,0][st:et] \n')
            fi.write('    good = np.where(data > 1.e-8) \n')
            fi.write('    merit = len(good) \n')
            fi.write('    if merit > pmerit: \n')
            fi.write('        pmerit = merit \n')
            fi.write('        goes_data = data \n')
            fi.write('        goes_t = temp_t[st:et] \n')
            fi.write('goes={\'time\':goes_t,\'xrsb\':data} \n')
            fi.write('fi2 = open("{}", "wb") \n'.format(goesdatafile))
            fi.write('pickle.dump(goes, fi2) \n')
            fi.write('fi2.close()')
            fi.close()

            try:
                os.system('python {}'.format(goesscript))
                #os.system('rm -rf {}'.format(goesscript))
            except NameError:
                print("Bad input names")
            except ValueError:
                print("Bad input values")
            except:
                print("Unexpected error:", sys.exc_info()[0])
                print("Error in generating GOES light curves. Proceed without GOES...")

            if os.path.exists(goesdatafile):
                fi1 = file(goesdatafile, 'rb')
                goest = pickle.load(fi1)
                fi1.close()
                dates = goest['time']
                goesdata = goest['xrsb']

    tr_plt = Time(timerange)
    if custom_ax is None:
        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8, 6))
    else:
        axs=custom_ax
    '''
    ax = axs[1]
    tidx_hessi, = np.where((hessi_tim >= tr_plt[0]) & (hessi_tim <= tr_plt[1]))
    for idx, eg in enumerate(hessi['obs_energies']):
        flux_hessi = ma.masked_array(hessi['obs_data'][tidx_hessi, idx])
        if hessi_smoth > 0:
            flux_hessi = DButil.smooth(flux_hessi, 20)
            flux_hessi = flux_hessi / np.nanmax(flux_hessi)
            ax.step(hessi_tim_plt[tidx_hessi], DButil.smooth(flux_hessi, hessi_smoth), label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))
        else:
            if idx <= 3:
                ax.step(hessi_tim_plt[tidx_hessi], flux_hessi, label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))
            else:
                ax.step(hessi_tim_plt[tidx_hessi], flux_hessi)
    if ylog:
        ax.set_yscale("log", nonposy='clip')
    else:
        ax.set_yscale("linear", nonposy='clip')
    ax.legend()
    #ax.legend(prop={'size': 7})
    ax.set_ylabel('Count Rate [s$^{-1}$ detector$^{-1}$ ]')
    ax.xaxis.set_visible(False)
    ax = axs[0]
    with open('/Volumes/WD6T/working/20170703/fermi/1s_tres_dect{0}.p'.format(nde),'rb') as fdspec_file:
        fdspec = pickle.load(fdspec_file)
    fdspec_file.close()
    ep.plot_fermi_dynamic_spec_in_ax(cax=ax, spec_dict=fdspec, timerange=timerange,erg_max=100.0,nde=nde)
    ax.xaxis.set_visible(False)
    '''
    ax = axs[0]
    #ax = axs[1]
    if goes:
        tidx_goes, = np.where((dates >= tr_plt[0].plot_date) & (dates <= tr_plt[1].plot_date))
        if goes and not only_dspec:
            ax.plot(dates[tidx_goes], goest.data['xrsb'][tidx_goes] / np.nanmax(goest.data['xrsb'][tidx_goes]), label='GOES 1.0$-$8.0 $\AA$')
            #ax.plot(dates[tidx_goes], goes.data['xrsa'][tidx_goes] / np.nanmax(goes.data['xrsa'][tidx_goes]), label='GOES 0.5--4.0 $\AA$')
        print(goest.data['xrsb'][tidx_goes])
    tidx_spec, = np.where((spec_tim >= tr_plt[0]) & (spec_tim <= tr_plt[1]))
    if len(tidx_spec) > 1:
        spec_tim_plt = spec_tim[tidx_spec[0]:tidx_spec[-1]].plot_date
        if not only_dspec:
            flux_colors = []
            for idx, fidx in enumerate(np.round(fidx_plt).astype(np.int)):
                flux_plt = medfilt(spec[min(fidx, nfreq - 1), tidx_spec[0]:tidx_spec[-1]], 7)
                #flux_plt = medfilt(spec[min(fidx, nfreq - 1), tidx_spec[0]:tidx_spec[-1]], 5)
                p = ax.plot(spec_tim_plt, flux_plt / np.nanmax(flux_plt), label='{:.2f} GHz'.format(freqghz[min(fidx, nfreq - 1)]))
                flux_colors.append(p[0].get_color())
            ax.set_ylabel('Normalized Intensity')
            ax.set_ylim(0, 1.1)
            ax.legend(prop={'size': 7})
            #!!!!!!!!!!!!!!!!!---------------------------
            #only if horizontal
            #formatter = mpl.dates.DateFormatter('%H:%M:%S')
            #ax.xaxis.set_major_formatter(formatter)
            #ax.fmt_xdata = formatter
            #ax.set_xlim(tr_plt.plot_date)
            #ax.set_xlabel('Start time ({})'.format(tr_plt[0].datetime.strftime('%d-%b-%y %H:%M:%S')))
        ax.xaxis.set_visible(False)

        ax = axs[1]
        #ax = axs[2]
        spec_plt = spec[:, tidx_spec[0]:tidx_spec[-1]]
        ax.pcolormesh(spec_tim_plt, freqghz, spec_plt, cmap=dspec_cmap, norm=colors.LogNorm(vmax=vmax, vmin=vmin))
        ax.set_ylabel('Frequency [GHz]')
        formatter = mpl.dates.DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(formatter)
        ax.fmt_xdata = formatter
        ax.set_xlim(tr_plt.plot_date)
        ax.set_xlabel('Start time ({})'.format(tr_plt[0].datetime.strftime('%d-%b-%y %H:%M:%S')))
        if not only_dspec:
            for idx, fidx in enumerate(np.round(fidx_plt).astype(np.int)):
                ax.axhline(freqghz[min(fidx, nfreq - 1)], color=flux_colors[idx], ls=':')
    else:
        print('Warning: No radio data in the timerange. Proceed without dynamic spectrum.')

    if custom_ax is None:
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.06)
        plt.show()
    else:
        return axs

'''

def lightcurves(timerange, outdir='./', specfile=None, goes=True, hessifile=None, fermifile=None, ylog=False, hessi_smoth=0, dspec_cmap='cubehelix',
                vmax=None, vmin=None):
    from sunpy.lightcurve import GOESLightCurve
    from sunpy.time import TimeRange, parse_time
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from astropy.time import Time
    import numpy as np
    import numpy.ma as ma
    from scipy.signal import medfilt
    import matplotlib.colors as colors
    from scipy.io import readsav
    from md_utils import DButil
    import os
    import pickle

    timerange = Time(timerange)
    if hessifile:
        if not os.path.exists(hessifile):
            hessi_script = 'HESSI_lc.pro'
            print('Run the script {} in SSWIDL to download RHESSI summary file first!'.format(hessi_script))
            fi = open(hessi_script, 'wb')
            fi.write(b'search_network, /enable \n')
            text = "time_range = ['{}','{}'] \n".format(timerange[0].datetime.strftime('%Y-%b-%d %H:%M:%S'),
                                                          timerange[1].datetime.strftime('%Y-%b-%d %H:%M:%S'))

            fi.write(text.encode())
            fi.write(b'obs_obj = hsi_obs_summary(obs_time_interval=time_range) \n')
            fi.write(b'data = obs_obj->getdata() \n')
            fi.write(b'info = obs_obj->get( / info) \n')
            fi.write(b'obs_data = data.countrate \n')
            fi.write(b'obs_times = obs_obj->getdata(/time) \n')
            fi.write(b'obs_times_str  = anytim(obs_times,/CCSDS) \n')
            fi.write(b'obs_energies  = fltarr(2,n_elements(info.energy_edges)-1) \n')
            fi.write(b'for ll=0,n_elements(info.energy_edges)-2 do begin \n')
            fi.write(b'    obs_energies[0,ll]=info.energy_edges[ll] \n')
            fi.write(b'    obs_energies[1,ll]=info.energy_edges[ll+1] \n')
            fi.write(b'endfor \n')
            text = 'save,filename="{}",OBS_DATA,OBS_ENERGIES,obs_times_str \n'.format(hessifile)
            fi.write(text.encode())
            fi.write(b'end \n')
            fi.close()
            return -1
        hessi = readsav(hessifile)
        hessi_tim = Time(list(hessi['obs_times_str']))
        hessi_tim_plt = hessi_tim.plot_date

    specdata = np.load(specfile)
    #with open(specfile, 'rb') as spec_open:
    #    specdata = pickle.load(spec_open)
    #spec_open.close()
    spec = specdata['spec']
    if len(spec.shape) == 4:
        (npol, nbl, nfreq, ntim) = spec.shape
    else:
        (nfreq, ntim) = spec.shape
    spec = np.mean(np.mean(spec, axis=0), axis=0)
    freq = specdata['freq']
    freqghz = freq / 1e9
    spec_tim = Time(specdata['tim'] / 3600. / 24., format='mjd')

    fidx_plt = np.linspace(0, nfreq - nfreq / 8.0 / 2.0, 8).astype(np.int) + nfreq / 8.0 / 2.0

    try:
        plt.style.use('seaborn-bright')
        params = {'font.size': 8, 'axes.grid': False, 'axes.facecolor': 'w', 'xtick.color': '#555555', 'ytick.color': '#555555',
                  'xtick.major.size': 2.0, 'ytick.major.size': 2.0, 'xtick.minor.size': 1.0, 'ytick.minor.size': 1.0, 'axes.axisbelow': False,
                  'axes.xmargin': 0.0, 'axes.ymargin': 0.0, 'axes.linewidth': 0.5, 'xtick.major.pad': 1.5, 'ytick.major.pad': 1.5,
                  'lines.linewidth': 1.0}
        mpl.rcParams.update(params)
    except:
        pass

    if goes:
        tr = TimeRange(timerange.iso)
        goes = GOESLightCurve.create(tr)
        dates = mpl.dates.date2num(parse_time(goes.data.index))

    tr_plt = Time(timerange)

    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(8, 6))
    ax = axs[0]
    tidx_hessi, = np.where((hessi_tim >= tr_plt[0]) & (hessi_tim <= tr_plt[1]))
    for idx, eg in enumerate(hessi['obs_energies']):
        flux_hessi = ma.masked_array(hessi['obs_data'][tidx_hessi, idx])
        if hessi_smoth > 0:
            flux_hessi = DButil.smooth(flux_hessi, 20)
            flux_hessi = flux_hessi / np.nanmax(flux_hessi)
            ax.step(hessi_tim_plt[tidx_hessi], DButil.smooth(flux_hessi, hessi_smoth), label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))
        else:
            ax.step(hessi_tim_plt[tidx_hessi], flux_hessi, label='{:.0f}-{:.0f} keV'.format(eg[0], eg[1]))

    if ylog:
        ax.set_yscale("log", nonposy='clip')
    else:
        ax.set_yscale("linear", nonposy='clip')
    ax.legend()
    ax.set_ylabel('Count Rate [s$^{-1}$ detector$^{-1}$ ]')

    ax = axs[1]
    tidx_goes, = np.where((dates >= tr_plt[0].plot_date) & (dates <= tr_plt[1].plot_date))
    if goes:
        ax.plot(dates[tidx_goes], goes.data['xrsb'][tidx_goes] / np.nanmax(goes.data['xrsb'][tidx_goes]), label='GOES 1.0--8.0 $\AA$')
        ax.plot(dates[tidx_goes], goes.data['xrsa'][tidx_goes] / np.nanmax(goes.data['xrsa'][tidx_goes]), label='GOES 0.5--4.0 $\AA$')
    tidx_spec, = np.where((spec_tim >= tr_plt[0]) & (spec_tim <= tr_plt[1]))
    if len(tidx_spec) > 1:
        spec_tim_plt = spec_tim[tidx_spec[0]:tidx_spec[-1]].plot_date
        flux_colors = []
        print(np.round(fidx_plt).astype(np.int))
        for idx, fidx in enumerate(np.round(fidx_plt).astype(np.int)):
            flux_plt = medfilt(spec[min(fidx, nfreq - 1), tidx_spec[0]:tidx_spec[-1]], 7)
            p = ax.plot(spec_tim_plt, flux_plt / np.nanmax(flux_plt), label='{:.2f} GHz'.format(freqghz[min(fidx, nfreq - 1)]))
            flux_colors.append(p[0].get_color())
        ax.set_ylabel('Flux (Normalized)')
        ax.set_ylim(0, 1.1)
        ax.legend()

        ax = axs[2]
        spec_plt = spec[:, tidx_spec[0]:tidx_spec[-1]]
        ax.pcolormesh(spec_tim_plt, freqghz, spec_plt, cmap=dspec_cmap, norm=colors.LogNorm(vmax=vmax, vmin=vmin))
        ax.set_ylabel('Frequency [GHz]')
        formatter = mpl.dates.DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(formatter)
        ax.fmt_xdata = formatter
        ax.set_xlim(tr_plt.plot_date)
        ax.set_xlabel('Start time ({})'.format(tr_plt[0].datetime.strftime('%d-%b-%y %H:%M:%S')))
        for idx, fidx in enumerate(np.round(fidx_plt).astype(np.int)):
            ax.axhline(freqghz[min(fidx, nfreq - 1)], color=flux_colors[idx], ls=':')
    else:
        print('Warning: No radio data in the timerange. Proceed without dynamic spectrum.')

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.06)
    #imgdir = outdir + '/fig01_{}-{}.png'.format(tr_plt[0].datetime.strftime('%H%M%S'), tr_plt[1].datetime.strftime('%H%M%S'))
    #fig.savefig(imgdir, dpi=200)
    #print('Save image to ' + imgdir)
'''