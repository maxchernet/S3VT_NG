import numpy as np
import pickle as pkl
import os
import semidiscrete as rt
from gp_emulator import *


def run_semidiscrete(xlai,xhc,rpl,xkab,scen,xkw,xkm,xleafn,xs1,xs2,xs3,xs4,lad,vza,vaa,sza,saa, nv=1):
        nv = 1
        wl = np.identity(2101)
        bands = np.arange(1,2102)
        params = [xlai,xhc,rpl,xkab,scen,xkw,xkm,xleafn,xs1,xs2,xs3,xs4,lad]
        rt.rt_modelpre(bands, 2101)
        return rt.rt_model(1, params, vza, vaa, sza, saa, wl,nv=nv)

gamma = [100, 1000, 1000000, 2000000, 3000000, 5000000, 8000000, 10000000, 50000000, 80000000, 1000000000,\
                10000000000, 100000000000, 1000000000000, 10000000000000, 100000000000000]
print 'gamma length:', len(gamma)

d= np.loadtxt('/home/max/s3vt/data_misr/US_Ne1_2007_ang7.brf', skiprows=1)
c=3
doy_misr = d[::7,0].astype(int)
vza = np.round(d[c::7, 2])
vaa = np.round(d[c::7, 3])
sza = np.round(d[c::7, 4])
saa = np.round(d[c::7, 5])
raa = np.round(abs(vaa - saa))
brf_misr = d[c::7, 6:]
print vza.shape
print doy_misr
params = np.zeros((doy_misr.shape[0], len(gamma), 10))
print params.shape

#xlai,xhc,rpl,xkab,scen,xkw,xkm,xleafn,xs1,xs2
doy = np.arange(1,366)
lai = np.zeros((len(gamma), 365))
for i in range(len(gamma)):
    fname = '/home/max/s3vt_ng/output_ng/misr_all_2007_ang7_gamma%d.pkl'%gamma[i]
    if os.path.isfile(fname):
        retval = pkl.load( open(fname, 'rb') )
        
        params[:,i,0] = retval['transformed_map']['xlai'][doy_misr-1]
        params[:,i,1] = retval['real_map']['xhc'][doy_misr-1]
        params[:,i,2] = retval['real_map']['rpl'][doy_misr-1]
        params[:,i,3] = retval['transformed_map']['xkab'][doy_misr-1]
        params[:,i,4] = retval['real_map']['scen'][doy_misr-1]
        params[:,i,5] = retval['transformed_map']['xkw'][doy_misr-1]
        params[:,i,6] = retval['transformed_map']['xkm'][doy_misr-1]
        params[:,i,7] = retval['real_map']['xleafn'][doy_misr-1]
        params[:,i,8] = retval['real_map']['xs1'][doy_misr-1]
        params[:,i,9] = retval['real_map']['xs2'][doy_misr-1]
        
        lai[i, :] = retval['real_map']['xlai']
#Load emulators
emulators = {}
emu_file = '/media/sf_MISR_EOLDAS/emul/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz'
for j in range(len(doy_misr)):
    #print file_emu
    file_emu = emu_file%(sza[i], vza[i], raa[i])
    if os.path.isfile(file_emu):
        #print tag
        emulators[j] = MultivariateEmulator (dump=file_emu)
#run model for all gammas
brf_ldas = np.zeros((doy_misr.shape[0], 2101))
gp_brf_ldas = np.zeros((doy_misr.shape[0], 2101))
for i in range(len(gamma)):
    print 'gamma=%d'%gamma[i]
    for j in range(len(doy_misr)):
        #print params[i,:]
        p = np.concatenate([params[j,i,:], [0.,0.,2.,vza[j],vaa[j], sza[j],saa[j]]])
        #print p.shape
        brf_ldas[j,:] = run_semidiscrete(*np.concatenate([params[j,i,:], [0.,0.,2.,vza[j],vaa[j], sza[j],saa[j]]]))[0]
        gp_brf_ldas[j,:] = emulators[j].predict(params[j,i,:])[0]
    np.savez('/home/max/s3vt_ng/output_ng/cross_val_gamma_%d'%gamma[i], gamma, brf_ldas, gp_brf_ldas)
