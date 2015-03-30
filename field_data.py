import numpy as np
from julday import *

def get_field(f_name):
    #Load field data
    ground = np.loadtxt(f_name, skiprows=1, usecols=(0,5,6,7))
    #Find julian days and actual dates of field measurements
    date_field=[]
    date_field_jul=[]
    for i in range(len(ground[:,0])):
        date_tmp = doy2date(int(ground[i,0]), int(ground[i,1]))
        date_field_jul.append(date2jul(date_tmp))
        date_field.append('%d.%d.%d'%(date_tmp.year, date_tmp.month, date_tmp.day))
    #Transform fAPAR LAI by interception
    lai_field = -np.log(1.-ground[:,3])*np.cos(30.*np.pi/180.)/0.5
    fapar = ground[:,3]
    return date_field, date_field_jul, lai_field, fapar

#Match field and satellite dates for comparison
def field_date(date1, date2):
    ind_1=[]
    ind_2=[]
    for d in range(len(date2)):
        dt = np.min(abs(date1 - date2[d]))
        ind = np.argmin(abs(date1 - date2[d]))
        if dt <= 1.: 
            ind_1 = np.append(ind_1, ind).astype(int)
            ind_2 = np.append(ind_2, d).astype(int)
    
    tmp_arr = np.array(date1)
    date1_2 = tmp_arr[ind_1]
    tmp_arr = np.array(date2)
    date2_2 = tmp_arr[ind_2]
    return date1_2, date2_2, ind_1, ind_2

if __name__=="__main__":
        (date_field_us2, date_field_jul_us2, lai_field_us2, fapar_field_us2) = get_field('/media/sf_MISR_EOLDAS/fAPARgreen_NE2.txt')
        print date_field_us2
