import semidiscrete as rt
from gp_emulator import *
import scipy.stats as ss
from collections import OrderedDict
from eoldas_ng import *
import os
from read_bounds import *
#from first_guess import *
import time
import pickle
from fapar_3 import *


timeBefore1 = time.clock()

def get_state():
        #Read bounds
        bounds = read_bounds('/home/max/s3vt/conf/misr_02')
        bounds[0,:] = -2.*np.log(bounds[0,:])[::-1]
        bounds[3,:] = -100.*np.log(bounds[3,:])[::-1]
        bounds[6,:] = -(1./100.)*np.log(bounds[6,:])[::-1]
        print bounds[0,:]
        print bounds[3,:]
        print '%.5f %.5f'%(bounds[5,0], bounds[5,1])
        print '%.5f %.5f'%(bounds[6,0], bounds[6,1])
        #**************************************************************
        #Define state and its parameters
        #***************************************************************
        #xlai,xhc,rpl,xkab,scen,xkw,xkm,xleafn,xs1,xs2
        state_config = OrderedDict ()

        state_config['xlai'] = VARIABLE
        state_config['xhc'] = VARIABLE
        state_config['rpl'] = VARIABLE
        state_config['xkab'] = VARIABLE
        state_config['scen'] = VARIABLE
        state_config['xkw'] = VARIABLE
        state_config['xkm'] = VARIABLE
        state_config['xleafn'] = VARIABLE
        state_config['xs1'] = VARIABLE
        state_config['xs2'] = VARIABLE

        default_par = OrderedDict ()
        default_par['xlai'] = 2.
        default_par['xhc'] = 1.
        default_par['rpl'] = 0.01
        default_par['xkab'] = 40.
        default_par['scen'] = 0.0
        default_par['xkw'] = (-1/50.)*np.log ( 0.995 )
        default_par['xkm'] = (-1/100.)*np.log ( 0.995 )
        default_par['xleafn'] = 1.5
        default_par['xs1'] = 1.0
        default_par['xs2'] = 0.0

        parameter_min = OrderedDict()
        parameter_max = OrderedDict()

        #bounds = read_bounds('/home/max/s3vt/conf/misr')
        min_vals = bounds[:,0]
        max_vals = bounds[:,1]

        for i, param in enumerate ( state_config.keys() ):
            parameter_min[param] = min_vals[i]
            parameter_max[param] = max_vals[i]

        # Define parameter transformations
        transformations = {
                    'xlai': lambda x: np.exp ( -x/2. ), \
                    'xkab': lambda x: np.exp ( -x/100. ), \
                    'xkw': lambda x: np.exp ( -50.*x ), \
                    'xkm': lambda x: np.exp ( -100.*x ) }
        inv_transformations = {
                    'xlai': lambda x: -2.*np.log ( x ), \
                    'xkab': lambda x: -100.*np.log ( x ), \
                    'xkw': lambda x: (-1/50.)*np.log ( x ), \
                    'xkm': lambda x: (-1/100.)*np.log ( x ) }

        # Define the state grid. In time in this case
        state_grid = np.arange ( 1, 366)

        state = State ( state_config, state_grid, default_par, parameter_min, parameter_max, verbose=False)

        # Set the transformations
        state.set_transformations ( transformations, inv_transformations )
        # Define the state
        return state

#*********************************************
#Cross validation
#********************************************
def cross_valid(state_file, cross_val=0):
        #Load data a file
        d = np.loadtxt(state_file, skiprows=1)
        c=3
        rho = d[c::7, 6:]
        #doys = d[::7,0].astype(int)
        print 'rho.shape:', rho.shape
        #Cross validation if cross_val > 0
        cross = 100/float(100-cross_val)
        #generate an array of integers of size rho
        a = np.arange(rho.shape[0])
        #randomly mix it up
        np.random.shuffle(a)
        #take % of initial size
        ind_cross = a[:rho.shape[0]/cross]
        #sort index incrementaly
        ind_cross = np.sort(ind_cross)
        return ind_cross


#***************************************************************************
#Define observational operator
#***************************************************************************

def get_obs_operator(emu_file, state, year, ind_cross, state_file, cam='An'):
    
    wl_misr = [443, 555, 670, 865]
    wl_width=[21,15,11,20]
    wl_full = np.arange(400,2501)
    #Load data a file
    d = np.loadtxt(state_file, skiprows=1)
    if cam=='An': c=3
    if cam=='Af': c=2
    if cam=='Aa': c=4
    if cam=='Bf': c=1
    if cam=='Ba': c=5
    if cam=='Cf': c=0
    if cam=='Ca': c=5
    doys = d[::7,0].astype(int)
    vza = np.round(d[c::7, 2])
    sza = np.round(d[c::7, 4])
    raa = np.round(abs(d[c::7, 3] - d[c::7, 5]))
    rho = d[c::7, 6:]
    print 'rho.shape:', rho.shape
    #446nm +-21, 558nm +-15, 672nm +-11, 866nm +-20
    b_min = np.array(wl_misr) - np.array(wl_width)
    b_max = np.array(wl_misr) + np.array(wl_width)
    
    #if cros_vall==0 they will be the same size
    rho = rho[ind_cross,:]
    vza = vza[ind_cross]
    sza = sza[ind_cross]
    raa = raa[ind_cross]
    doys = doys[ind_cross]
    print 'rho.shape:', rho.shape

    n_bands = b_min.shape[0]
    band_pass = np.zeros(( n_bands,2101), dtype=np.bool)
    bw = np.zeros( n_bands )
    bh = np.zeros( n_bands )
    for i in xrange( n_bands ):
        band_pass[i,:] = np.logical_and ( wl_full >= b_min[i], wl_full <= b_max[i] )
        bw[i] = b_max[i] - b_min[i]
        bh[i] = ( b_max[i] + b_min[i] )/2.
    
    sigma_min = 0.001
    sigma_max = 0.04
    sigma_obs = (sigma_max - sigma_min)*(bh-bh.min())/(bh.max()-bh.min())
    sigma_obs += sigma_min
    bu = sigma_obs
    print bu
    bu = [0.0040, 0.0044, 0.0047, 0.0054]
    print bu
    rho_big = np.zeros(( 365,n_bands ))
    mask = np.zeros(( 365, 4))

    for i in state.state_grid:
        if i in doys:
            r = rho[doys==i, :].squeeze()
            rho_big[i-1,:] = r
            #mask[ i, :] = [ i, raa[doys==i], np.round(sza[doys==i]/5.)*5, np.round(vza[doys==i]/5.)*5 ]
            mask[ i-1, :] = [ i, raa[doys==i], np.round(sza[doys==i]), np.round(vza[doys==i]) ]
        
    ind = np.where(mask!=0)
    #print mask[ind]

    emulators = {}
    for i in range(doys.shape[0]):
        
        #tag = tuple([ np.round(sza[i]/5.)*5, np.round(vza[i]/5.)*5 ])
        tag = tuple([ np.round(sza[i]), np.round(vza[i]) ])
    
        file_emu = emu_file%(sza[i], vza[i], raa[i])
        print sza[i], vza[i], raa[i]
        if os.path.isfile(file_emu):
            print file_emu
            print tag
            emulators[tag] = MultivariateEmulator ( dump=file_emu)
        else: print 'emulator file '+file_emu+' does not exist'
    obs = ObservationOperatorTimeSeriesGP(state.state_grid, state, rho_big, mask, emulators, bu, band_pass, bw)
    return obs, doys
    #return rho_big[:,2]

#******************************************************************************************************************
#Define prior
#******************************************************************************************************************
def get_prior(state):
        mu_prior = OrderedDict ()
        prior_inv_cov = OrderedDict ()
        prior_inv_cov['xlai'] = np.array([1.0])
        prior_inv_cov['xhc'] = np.array([1.0])
        prior_inv_cov['rpl'] = np.array([1.0])
        prior_inv_cov['xkab'] = np.array([1.0])
        prior_inv_cov['scen'] = np.array([1.0])
        prior_inv_cov['xkw'] = np.array([1.0])
        prior_inv_cov['xkm'] = np.array([1.0])
        prior_inv_cov['xleafn'] = np.array([1.0])
        prior_inv_cov['xs1'] = np.array([2.0])
        prior_inv_cov['xs2'] = np.array([2.0])
            
        for param in state.parameter_min.iterkeys():
            if state.transformation_dict.has_key ( param ):
                mu_prior[param] = state.transformation_dict[param](np.array([state.default_values[param]]) )
            else:
                mu_prior[param] = np.array([state.default_values[param]])
            prior_inv_cov[param] = 1./prior_inv_cov[param]**2

        #load lai "phenological model"
        #prior_lai = np.load('/home/max/s3vt_ng/data/fit_ts_mean.npz')['arr_0']
        #mu_prior['xlai'] = state.transformation_dict['xlai'](prior_lai)

        return Prior ( mu_prior, prior_inv_cov )

#*******************************************************************
#Do a first guess
#*******************************************************************
def get_first_guess(state, obs, year, save_dir, n_site=1):
        xlai1=[] 
        xhc1=[] 
        rpl1=[] 
        xkab1=[]
        scen1=[]
        xkw1=[]
        xkm1=[]
        xleafn1=[] 
        xs11=[]
        xs21=[]
        wl = np.array([443, 555, 670, 865])
        #Make a distribution
        distr_file = '/home/max/s3vt/data/US_Ne%d_misr_distr'%1
        #read train dataset from a file
        npzfile = np.load(distr_file+'.npz')
        #load train state parameters
        x_train = npzfile['arr_0']
        #train reflectance
        brf_train = npzfile['arr_1']
        train_misr=[]
        #Reduce number of bands of trained dataset to 
        #a nummber of bands of satellite data
        for i in xrange(wl.size):
            tmp = brf_train[:,wl[i]-400]
            train_misr.append(tmp)
        train_misr = np.array(train_misr)
        n_train = brf_train.shape[0]

        for i in state.state_grid:
                if obs.mask[i-1, 0]!=0:
                        #print 'obs.observations[i-1,:]:', obs.observations[i-1,:]
                        #Find sum of difference
                        sum_misr=[]
                        for j in xrange(n_train):
                                sum_misr.append(np.sum(abs(obs.observations[i-1,:] - train_misr[:,j])))
                                #print 'obs.observations[i-1,:]:',obs.observations[i-1,:]
                                #print '**'
                                #print 'train_misr[:,i]:', train_misr[:,j]
                                #print'***'
                        #Find a minimum and this is an end of first guess
                        min_i = np.argmin(sum_misr)

                        #Initial values = first guess
                        
                        xlai1.append(-2.*np.log(x_train[min_i,0]))
                        xhc1.append(x_train[min_i,1])
                        rpl1.append( x_train[min_i,2] )
                        xkab1.append(-100*np.log( x_train[min_i,3] ))
                        scen1.append( x_train[min_i,4] )
                        xkw1.append(-(1./50.)*np.log( x_train[min_i,5] ))
                        xkm1.append(-(1./100.)*np.log( x_train[min_i,6] ))
                        xleafn1.append( x_train[min_i,7] )
                        xs11.append( x_train[min_i,8] )
                        xs21.append( x_train[min_i,9] )
                
        x_dict = {}
        x_dict['xlai'] = np.interp(state.state_grid, doys, xlai1)
        x_dict['xhc'] = np.interp(state.state_grid, doys, xhc1)
        x_dict['rpl'] = np.interp(state.state_grid, doys, rpl1)
        x_dict['xkab'] = np.interp(state.state_grid, doys, xkab1)
        x_dict['scen'] = np.interp(state.state_grid, doys, scen1)
        x_dict['xkw'] = np.interp(state.state_grid, doys, xkw1)
        x_dict['xkm'] = np.interp(state.state_grid, doys, xkm1)
        x_dict['xleafn'] = np.interp(state.state_grid, doys, xleafn1)
        x_dict['xs1'] = np.interp(state.state_grid, doys, xs11)
        x_dict['xs2'] = np.interp(state.state_grid, doys, xs21)
        pickle.dump( x_dict, open( save_dir+'first_guess_Ne%d_%d.pkl'%(n_site, year), "wb" ) )
        return x_dict

#.............................................................................................................................
#*MAIN************************************************************************************************************************
#.............................................................................................................................
#US_Ne1_%d_ang7.brf
#data_dir = '/home/max/s3vt/data_misr/'

lad=5
save_dir = '/home/max/s3vt_ng/output_ng/no_prior_lad%d/'%lad

for n_site in [2,3]:
        for year in [2002, 2004, 2006, 2008]: #range(2001,2009):
                state_file = '/home/max/s3vt/data_misr/US_Ne%d_%d_ang7.brf'%(n_site, year)
                state = get_state()
                #% of taken off data for cross validation
                cross_val = 0
                ind_cross = cross_valid(state_file, cross_val=cross_val)
                #print 'An'
                obs_misr_an, doys = get_obs_operator('/media/sf_MISR_EOLDAS/emul_lad'+str(lad)+'/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz', state, year, ind_cross, state_file, cam='An')
                #print 'Af'
                obs_misr_af, doys = get_obs_operator('/media/sf_MISR_EOLDAS/emul_lad'+str(lad)+'/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz', state, year, ind_cross, state_file, cam='Af')
                #print 'Aa'
                obs_misr_aa, doys = get_obs_operator('/media/sf_MISR_EOLDAS/emul_lad'+str(lad)+'/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz', state, year, ind_cross, state_file, cam='Aa')
                #print 'Bf'
                obs_misr_bf, doys = get_obs_operator('/media/sf_MISR_EOLDAS/emul_lad'+str(lad)+'/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz', state, year, ind_cross, state_file, cam='Bf')
                #print 'Ba'
                obs_misr_ba, doys = get_obs_operator('/media/sf_MISR_EOLDAS/emul_lad'+str(lad)+'/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz', state, year, ind_cross, state_file, cam='Ba')
                #print 'Cf'
                obs_misr_cf, doys = get_obs_operator('/media/sf_MISR_EOLDAS/emul_lad'+str(lad)+'/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz', state, year, ind_cross, state_file, cam='Cf')
                #print 'Ca'
                obs_misr_ca, doys = get_obs_operator('/media/sf_MISR_EOLDAS/emul_lad'+str(lad)+'/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz', state, year, ind_cross, state_file, cam='Ca')

                prior = get_prior(state)

                for gamma in [10000000]:
                #[50000000, 80000000, 100000000000000]:#[100000000000, 1000000000000, 10000000000000, 100000000000000]:#[1000000000, 10000000000]:#[100, 1000, 1000000, 2000000, 3000000, 5000000, 8000000, 10000000]:
                        #***************************************************************************
                        #Define temporal constraint
                        #***************************************************************************
                        temporal = TemporalSmoother ( state.state_grid, gamma=gamma, required_params=["xlai", "xhc", "rpl", "xkab", "scen", "xkw", "xkm", "xleafn", "xs1", "xs2"] )

                        #********************************************************************
                        #Add all operators to state
                        #********************************************************************
                        state.add_operator("Prior", prior )
                        state.add_operator("Temporal", temporal)
                        state.add_operator("MISR_An", obs_misr_an )
                        state.add_operator("MISR_Af", obs_misr_af )
                        state.add_operator("MISR_Aa", obs_misr_aa )
                        state.add_operator("MISR_Bf", obs_misr_bf )
                        state.add_operator("MISR_Ba", obs_misr_ba )
                        state.add_operator("MISR_Cf", obs_misr_cf )
                        state.add_operator("MISR_Ca", obs_misr_ca )

                        
                        x_dict = get_first_guess(state, obs_misr_bf, year, save_dir, n_site=n_site)

                        #******************************************************************
                        #Do Optimization
                        #******************************************************************
                        timeBefore2 = time.clock()
                        retval = state.optimize ( x_dict, do_unc=True )
                        timeAfter2 = time.clock()
                        elapsed_time2 = timeAfter2 - timeBefore2
                        print 'optimization elapsed time (s): ', elapsed_time2
                        f_retval = save_dir+'misr_all_Ne%d_%d_ang7.pkl'%(n_site, year)

                        # Save a dictionary into a pickle file.
                        pickle.dump( retval, open( f_retval, "wb" ) )

                        # Write files which will be used by Fortran for calculation of fAPAR
                        dir_params = '/home/max/s3vt_ng/output_ng/no_prior_lad%d/'%lad
                        dir_fapar_mean = '/home/max/s3vt_ng/output_fapar_mean/no_prior_lad%d/'%lad
                        dir_fapar_sd = '/home/max/s3vt_ng/output_fapar_sd/no_prior_lad%d/'%lad
                        try:
                                write_fapar_in(f_retval, year, dir_params, dir_fapar_mean, dir_fapar_sd, n_train=50, n_site=n_site, lad=lad)
                        except ValueError:
                                print 'Something wrong with write_fapar_in()...'

timeAfter1 = time.clock()
elapsed_time1 = timeAfter1 - timeBefore1
print 'Total elapsed time (s): ', elapsed_time1
print 'Done!!!'
