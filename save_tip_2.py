import numpy as np
import jdcal
import idlsave
import os.path
import matplotlib.dates as dates
import julday

#****************************************************************************

def get_time_series(misr_sav, n, date_step=2):
    date=[]
    julday = np.zeros(n)
    eff_lai = np.zeros(n)
    uncert = np.zeros(n)
    fapar = np.zeros(n)
    obscovarfluxes = np.zeros((n,6))
    stddevparams = np.zeros((n,7))
    effasym_vis = np.zeros(n)
    effssa_vis = np.zeros(n)
    truebkgdalbedo_vis = np.zeros(n)

    j=0
    for i in range(len(misr_sav.jdates)):
        if misr_sav.tipstrucarr[i].efflai>0:
            date_tmp = jdcal.jd2jcal(2400000, misr_sav.jdates[i] - 2400000)
            julday[j] = misr_sav.jdates[i]
            date.append('.'.join("%d"%d for d in date_tmp[0:3]))
            eff_lai[j] = misr_sav.tipstrucarr[i].efflai
            obscovarfluxes[j,:] = misr_sav.tipstrucarr[i].obscovarfluxes
            stddevparams[j,:] = misr_sav.tipstrucarr[i].stddevparams
            effasym_vis[j] = misr_sav.tipstrucarr[i].effasym_vis
            fapar[j] = misr_sav.tipstrucarr[i].abs_vis
            effssa_vis[j] = misr_sav.tipstrucarr[i].effssa_vis
            truebkgdalbedo_vis[j] = misr_sav.tipstrucarr[i].truebkgdalbedo_vis
            j+=1

    date_step = int(np.round(n/(n/3.)))
    return julday, eff_lai, obscovarfluxes, effasym_vis, fapar, stddevparams, effssa_vis, truebkgdalbedo_vis

#******************************************************************************


sav_toc_f = ['/media/sf_Satellite/S3VT/MISR/TIP_timeline_P027_B058_S3VT_US_V1.04-3_20150305T194214.sav',\
             '/media/sf_Satellite/S3VT/MISR/TIP_timeline_P028_B058_S3VT_US_V1.04-3_20150305T194232.sav',\
             '/media/sf_Satellite/S3VT/MISR/TIP_timeline_P029_B057_S3VT_US_V1.04-3_20150305T194251.sav']

sav_toc=[]
for i in range(len(sav_toc_f)):
        print 'Now reading ', sav_toc_f[i]
        sav_toc.append(idlsave.read(sav_toc_f[i]))

misr_tip=[]
for i in range(len(sav_toc)):
    for j in range(len(sav_toc[i].tlstruc)):
        if sav_toc[i].tlstruc[j].desc[-9:] == '+000_+000':
            misr_tip.append(sav_toc[i].tlstruc[j])
            print i, j, sav_toc[i].tlstruc[j].desc, sav_toc[i].tlstruc[j].lat, sav_toc[i].tlstruc[j].lon

n_dates=[]
for i in range(len(misr_tip)):
    k=0
    for j in range(len(misr_tip[i].jdates)):
        if misr_tip[i].tipstrucarr[j].efflai>0:
            k+=1
    n_dates.append(k)
    #If there are sites with negative spectra only show them:
    if k == 0: print misr_tip[i].desc

print len(n_dates)
print n_dates

j=1
for i in range(len(misr_tip)):
    if n_dates[i] != 0:
        (julday, eff_lai, obscovarfluxes, effasym_vis, fapar, stddevparams, effssa_vis, truebkgdalbedo_vis) = get_time_series(misr_tip[i], n_dates[i])
        np.savez('/home/max/s3vt_ng/data/misr_tip_'+misr_tip[i].desc[1:-10]+'_%d'\
                %(misr_tip[i].path[0]), julday, eff_lai, obscovarfluxes, effasym_vis, fapar, stddevparams, effssa_vis, truebkgdalbedo_vis)
        j+=1
