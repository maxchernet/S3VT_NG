import numpy as np
from gp_emulator import *
import semidiscrete as rt
import os
import scipy.stats as ss

def run_semidiscrete(xlai,xhc,rpl,xkab,scen,xkw,xkm,xleafn,xs1,xs2,xs3,xs4,lad,vza,vaa,sza,saa, nv=1):
        nv = 1
        wl = np.identity(2101)
        bands = np.arange(1,2102)
        params = [xlai,xhc,rpl,xkab,scen,xkw,xkm,xleafn,xs1,xs2,xs3,xs4,lad]
        rt.rt_modelpre(bands, 2101)
        return rt.rt_model(1, params, vza, vaa, sza, saa, wl,nv=nv)

#**************************************************
#Read bounds of the parameters from .conf file
def read_bounds(conf_file):
        bounds = np.zeros((10,2))
        f = open(conf_file)
        list_opt = f.read().split('\n')
        for i in range(0,len(list_opt)):
                if '[parameter.x.assoc_bounds]' in list_opt[i]:
                        tmp_lst =  [ss.split('=') for ss in list_opt[i+2:i+12]]
                        j=0
                        for line in tmp_lst:
                                #print line
                                bounds[j,:] =  np.array(line[1].split(',')).astype(float)
                                j+=1
        return bounds

#****************************************************

def do_emul(emul_dir, is_gp=False, site=1):
        #time_total_1 = time.clock()
        year = np.arange(2000, 2015)
        cam = ['Cf','Bf','Af','An','Aa','Ba','Ca']
        n_ang = 7
        train_file = 'nad_%03d_sza_%03d_vza_%03d_raa_train'
        angels=[0,0,0,0,0]
        for y in year:
                d = np.loadtxt('/home/max/s3vt/data_misr/US_Ne%d_%d_ang7.brf'%(site,y), skiprows=1)
                #print 'number of days in a year %d:'%y, d.shape[0]/n_ang
                vza = np.zeros((int(d.shape[0]/n_ang), n_ang))
                vaa = np.zeros((int(d.shape[0]/n_ang), n_ang))
                sza = np.zeros(int(d.shape[0]/n_ang))
                saa = np.zeros(int(d.shape[0]/n_ang))
                doy = np.zeros(int(d.shape[0]/n_ang))
                j=0
                for i in range(0,d.shape[0],n_ang):
                        vza[j,:] = d[i:i+n_ang,2]
                        vaa[j,:] = d[i:i+n_ang,3]
                        sza[j] = d[i,4]
                        saa[j] = d[i,5]
                        doy[j] = d[i,0]
                        j=j+1

                for j in range(sza.shape[0]):
                        #print 'DOY: ',doy[j]
                        #timeBefore = time.clock()
                        #Only train a model
                        if is_gp == False:
                                #Check if we already run model for this combinations of angels
                                vza_n=0 #A number of already calculated view angles for one day
                                for v in range(vza.shape[1]):
                                        file_train = emul_dir + train_file%\
                                                (np.round(sza[j]), np.round(vza[j,v]), np.round(abs(vaa[j,v]-saa[j])))
                                        if os.path.isfile(file_train+'.npz'): vza_n = vza_n + 1
                                #If at least one combination of angles for this day is not exists
                                if vza_n < n_ang:
                                        #print 'Year %d, DOY %d, Site %d. Do RT run for angles:'%(y,doy[j],site)
                                        #print vza[j,:], vaa[j,:], sza[j], saa[j]
                                        #Create samples by running a model
                                        #(y_obs, x_train) = create_emulator(vza[j,:], vaa[j,:], sza[j], saa[j], nv=len(vza))
                                        #print vza[j,:]
                                        #Save these samples
                                        for v in range(vza.shape[1]):
                                                file_train = emul_dir + train_file%\
                                                        (np.round(sza[j]), np.round(vza[j,v]), np.round(abs(vaa[j,v]-saa[j])))
                                                #np.savez(file_train, y_obs[v,:,:], x_train)
                                #else:
                                        #print '*****'
                                        #print 'These combinations of angles already exists (year %d, DOY %d, site %d):'%(y,doy[j],site)
                                        #print vza[j,:], vaa[j,:], sza[j], saa[j]
                                        #print '*****'
                        #Train emulators
                        else:
                                for v in range(vza.shape[1]):
                                        file_train = emul_dir+'nad_%03d_sza_%03d_vza_%03d_raa_train'%\
                                                (np.round(sza[j]), np.round(vza[j,v]), np.round(abs(vaa[j,v]-saa[j])))
                                        angels = np.vstack(( angels, [np.round(abs(vaa[j,v]-saa[j])), np.round(vza[j,v]),\
                                            np.round(vaa[j,v]), np.round(sza[j]), np.round(saa[j])] ))
        return angels

#***********************************************************

emul_dir = '/media/sf_MISR_EOLDAS/emul/'
#for site in [1]:
angels = do_emul(emul_dir=emul_dir, is_gp=True, site=1)
#do_emul(emul_dir=emul_dir, is_gp=True, site=site)
print angels.shape

#*********************************************************

bounds = read_bounds('/home/max/misr/conf/misr_obs1.conf')
min_vals = np.append(bounds[:,0], [15,10,10])
max_vals = np.append(bounds[:,1], [70,70,400])
n_params=min_vals.shape[0]
print n_params
n_train = 1000000
dist=[]
#Get a distrubution for sampling
for k in xrange(n_params):
    dist.append(ss.uniform(loc=min_vals[k], scale=max_vals[k]-min_vals[k]))
# The training dataset is obtaiend by a LatinHypercube Design
x_train = lhd(dist=dist, size=n_train )
#print np.round(x_train[:,10:13])

#***************************************************************

j=0
emu={}
brf_nadim = np.zeros(2101)
brf_emul = np.zeros(2101)
for i in range(n_train):
    s = np.round(x_train[i,10]).astype(int)
    v = np.round(x_train[i,11]).astype(int)
    r = np.round(x_train[i,12]).astype(int)
    f_name = '/media/sf_MISR_EOLDAS/emul/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz'%(s,v,r)
    #print f_name
    if os.path.isfile(f_name):
        j=j+1
        print f_name
        tag = tuple([s,v])
        #vza_t.append(v)
        emu[tag] = MultivariateEmulator (dump=f_name)
        brf_pred = emu[tag].predict(x_train[i,0:10])
        
        print x_train[i,10], x_train[i,11], x_train[i,12]
        
        ind = np.where( np.logical_and(np.round(x_train[i,12]) == angels[:,0], np.round(x_train[i,10]) == angels[:,3]) )
        print angels[ind[0][0],:]
        
        x_params = np.append(x_train[i,0:10], [0,0,2])
        x_params = np.append(x_params, [x_train[i,11], angels[ind[0][0],2], x_train[i,10], angels[ind[0][0],4]])
        
        brf_full = run_semidiscrete(*x_params)
        
        brf_nadim = np.vstack(( brf_nadim,  brf_full[0]))
        brf_emul = np.vstack(( brf_emul,  brf_pred[0]))
print j
print brf_nadim.shape
print brf_emul.shape
np.savez('/home/max/s3vt_ng/data/gp-validation', brf_nadim, brf_emul)
