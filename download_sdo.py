import drms
from datetime import datetime, timedelta
import os
import shutil
import sys
import pandas as pd 


def datetime_to_str(time):
    return time.strftime("%Y-%m-%d")

#my_email = 'yw633@njit.edu'
my_email = 'nh72@njit.edu'

print( 'start downloading AIA from JSOC...')

#out_dir = '/Volumes/Data/20170715/aia/py_download/'
#out_dir = '/Volumes/WD6T/working/20170703/aia_extra/'
#out_dir = '/Volumes/WD6T/working/20220511/aia/'
#if not os.path.exists(out_dir):
#    os.makedirs(out_dir)
WAVELNTH = [131,94,171,193,211,304,335,1600,1700]
#WAVELNTH = [94,131,171,211,335]
#WAVELNTH = [335]
#WAVELNTH = [1600,1700]
series_name = 'aia.lev1_euv_12s'
#series_name = 'aia.lev1_uv_24s'
#series_name = 'aia.lev1_euv_12s'
#start_date = datetime(2017,7,3)
#start_date = datetime(2017,7,3)
#start_date = datetime(2017,8,20)
start_date = datetime(2022,5,11)
#start_date = datetime(2017,7,15)
#ending_date = datetime(2017,7,3)
#ending_date = datetime(2017,8,20)
ending_date = datetime(2022,5,11)
#ending_date = datetime(2017,7,15)
date_counter = start_date
print('start_date:', start_date)
index = 1
delta_date = ending_date - start_date
num_days = delta_date.days
 
while date_counter <= ending_date:
    try:
        for w in WAVELNTH:
            #out_dir = '/Volumes/Data/20170715/aia/py_download/'
            #out_dir = '/Volumes/WD6T/working/20170703/aia/'
            #out_dir = '/Volumes/WD6T/working/20170820/aia/'
            #out_dir = '/Volumes/Data/from_qq/DEM/'
            #out_dir = '/Volumes/WD6T/working/20220511/aia/'
            out_dir = '/Volumes/Data/20170820/20220511/aia/'
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            out_dir = out_dir + "% s" % w
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            #for h in range(19,20):
            for h in range(0,1):
                #q = series_name + '[' + datetime_to_str(date_counter) + 'T' + '%0*d' % (2,h) +':00/1h]['+ "% s" % w +']'
                #q = series_name + '[' + datetime_to_str(date_counter) + 'T15:30:00/30m]['+ "% s" % w +']'
                #q = series_name + '[' + datetime_to_str(date_counter) + 'T16:17:00/2h]['+ "% s" % w +']'
                q = series_name + '[' + datetime_to_str(date_counter) + 'T18:30:00/30m]['+ "% s" % w +']'
                print('Downloading time:', q)
            #q = q + '{' + "% s" % w + '}'
                c = drms.Client(email=my_email, verbose=True)
                r = c.export(q, method='url', protocol='fits')
                r.wait()
                r.download(out_dir)
            print('Working on index:', index, 'of', num_days)
            index = index + 1
            print('Query:' + q)
    except Exception as e:
        print('unable to download at time:', q, e)
    date_counter = date_counter + timedelta(days=1)
print('done...')