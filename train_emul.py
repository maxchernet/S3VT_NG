import numpy as np
from gp_emulator import *

#load train samples
file_train = '/media/sf_MISR_EOLDAS/emul/nad_025_sza_012_vza_212_raa_train.npz'
y_obs_gp = np.load(file_train)['arr_0']
x_train =  np.load(file_train)['arr_1']
file_gp = file_train.replace('train', 'gp')
#Create emulators
em = MultivariateEmulator (X=y_obs_gp, y=x_train)
em.dump_emulator(file_gp)
