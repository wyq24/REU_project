#from suncasa.suncasatasks import subvs
#from suncasa.suncasatasks import concateovsa
import os
def makelist(tdir='', keyword1='', keyword2='', exclude=None):
    li = []
    # for root, dirs, files in os.walk(tdir):
    root = os.getcwd()
    files = os.listdir(tdir)
    for file in files:
        if exclude is None:
            if keyword1 in file and keyword2 in file:
                li.append(os.path.join(tdir, file))
        else:
            if keyword1 in file and keyword2 in file and exclude not in file:
                li.append(os.path.join(tdir, file))

    return li
# sub_tr = '2022/05/11/18:05:00~2022/05/11/18:08:00'
# tr = '2022/05/11/19:40:00~2022/05/11/19:50:00'
# opv = '/Volumes/Data/20170820/20220511/eovsa/bkg_subed/msdata/IDB20220511_1940_1950_sub_1805_1808.ms.XXYY.slfcaled'
# ipv = '/Volumes/Data/20170820/20220511/eovsa/eovsa_full/msdata/IDB20220511_1800-2000.ms.XXYY.slfcaled'
# subvs(vis = ipv, outputvis=opv,timerange=tr,subtime1=sub_tr)
vis_list = makelist(tdir='/Volumes/Data/20170820/20220511/eovsa/bkg_subed/msdata/',keyword1='IDB',keyword2='led',
                    exclude='/Volumes/Data/20170820/20220511/eovsa/bkg_subed/msdata/IDB20220511_1841_1843_sub_1805_1808.ms.XXYY.slfcaled')
print(vis_list.sort())
#concateovsa()