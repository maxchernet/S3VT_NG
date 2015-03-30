import numpy as np

def get_tip(n_site=1):
        julday_tip_us1=[]
        julday_tip_us2=[]
        julday_tip_us3=[]

        lai_tip_us1=[]
        lai_tip_us2=[]
        lai_tip_us3=[]

        fapar_tip_us1=[]
        fapar_tip_us2=[]
        fapar_tip_us3=[]

        lai_std_tip_us1=[]
        lai_std_tip_us2=[]
        lai_std_tip_us3=[]

        fapar_std_tip_us1=[]
        fapar_std_tip_us2=[]
        fapar_std_tip_us3=[]

        effssa_vis_us1=[]
        effssa_vis_us2=[]
        effssa_vis_us3=[]

        effssa_std_us1=[]
        effssa_std_us2=[]
        effssa_std_us3=[]

        truebkgdalbedo_vis_us1=[]
        truebkgdalbedo_vis_us2=[]
        truebkgdalbedo_vis_us3=[]

        truebkgdalbedo_std_us1=[]
        truebkgdalbedo_std_us2=[]
        truebkgdalbedo_std_us3=[]

        for i in [27,28,29]:
                julday_tip_us1 = np.concatenate((julday_tip_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne1_%d.npz'%i)['arr_0']))
                julday_tip_us2 = np.concatenate((julday_tip_us2, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne2_%d.npz'%i)['arr_0']))
                julday_tip_us3 = np.concatenate((julday_tip_us3, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne3_%d.npz'%i)['arr_0']))

                lai_tip_us1 = np.concatenate((lai_tip_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne1_%d.npz'%i)['arr_1']))
                lai_tip_us2 = np.concatenate((lai_tip_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne2_%d.npz'%i)['arr_1']))
                lai_tip_us3 = np.concatenate((lai_tip_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne3_%d.npz'%i)['arr_1']))

                lai_std_tip_us1 = np.concatenate((lai_std_tip_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne1_%d.npz'%i)['arr_5'][:,0]))
                lai_std_tip_us2 = np.concatenate((lai_std_tip_us2, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne2_%d.npz'%i)['arr_5'][:,0]))
                lai_std_tip_us3 = np.concatenate((lai_std_tip_us3, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne3_%d.npz'%i)['arr_5'][:,0]))

                effssa_std_us1 = np.concatenate((effssa_std_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne1_%d.npz'%i)['arr_5'][:,1]))
                effssa_std_us2 = np.concatenate((effssa_std_us2, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne2_%d.npz'%i)['arr_5'][:,1]))
                effssa_std_us3 = np.concatenate((effssa_std_us3, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne3_%d.npz'%i)['arr_5'][:,1]))

                truebkgdalbedo_std_us1 = np.concatenate((truebkgdalbedo_std_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne1_%d.npz'%i)['arr_5'][:,3]))
                truebkgdalbedo_std_us2 = np.concatenate((truebkgdalbedo_std_us2, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne2_%d.npz'%i)['arr_5'][:,3]))
                truebkgdalbedo_std_us3 = np.concatenate((truebkgdalbedo_std_us3, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne3_%d.npz'%i)['arr_5'][:,3]))

                fapar_std_tip_us1 = np.concatenate((fapar_std_tip_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne1_%d.npz'%i)['arr_2'][:,1]))
                fapar_std_tip_us2 = np.concatenate((fapar_std_tip_us2, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne2_%d.npz'%i)['arr_2'][:,1]))
                fapar_std_tip_us3 = np.concatenate((fapar_std_tip_us3, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne3_%d.npz'%i)['arr_2'][:,1]))

                fapar_tip_us1 = np.concatenate((fapar_tip_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne1_%d.npz'%i)['arr_4']))
                fapar_tip_us2 = np.concatenate((fapar_tip_us2, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne2_%d.npz'%i)['arr_4']))
                fapar_tip_us3 = np.concatenate((fapar_tip_us3, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne3_%d.npz'%i)['arr_4']))
                
                effssa_vis_us1 = np.concatenate((effssa_vis_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne1_%d.npz'%i)['arr_6']))
                effssa_vis_us2 = np.concatenate((effssa_vis_us2, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne2_%d.npz'%i)['arr_6']))
                effssa_vis_us3 = np.concatenate((effssa_vis_us3, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne3_%d.npz'%i)['arr_6']))
                
                truebkgdalbedo_vis_us1 = np.concatenate((truebkgdalbedo_vis_us1, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne1_%d.npz'%i)['arr_7']))
                truebkgdalbedo_vis_us2 = np.concatenate((truebkgdalbedo_vis_us2, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne2_%d.npz'%i)['arr_7']))
                truebkgdalbedo_vis_us3 = np.concatenate((truebkgdalbedo_vis_us3, np.load('/home/max/s3vt_ng/data/misr_tip_US_Ne3_%d.npz'%i)['arr_7']))

        ind1 = np.argsort(julday_tip_us1)
        julday_tip_us1 = julday_tip_us1[ind1]
        lai_tip_us1 = lai_tip_us1[ind1]
        fapar_tip_us1 = fapar_tip_us1[ind1]
        lai_std_tip_us1 = lai_std_tip_us1[ind1]
        effssa_vis_us1 = effssa_vis_us1[ind1]
        truebkgdalbedo_vis_us1 = truebkgdalbedo_vis_us1[ind1]
        effssa_std_us1 = effssa_std_us1[ind1]
        truebkgdalbedo_std_us1 = truebkgdalbedo_std_us1[ind1]

        ind1 = np.argsort(julday_tip_us2)
        julday_tip_us2 = julday_tip_us2[ind1]
        lai_tip_us2 = lai_tip_us2[ind1]
        fapar_tip_us2 = fapar_tip_us2[ind1]
        lai_std_tip_us2 = lai_std_tip_us2[ind1]
        effssa_vis_us2 = effssa_vis_us2[ind1]
        truebkgdalbedo_vis_us2 = truebkgdalbedo_vis_us2[ind1]
        effssa_std_us2 = effssa_std_us2[ind1]
        truebkgdalbedo_std_us2 = truebkgdalbedo_std_us2[ind1]

        ind1 = np.argsort(julday_tip_us3)
        julday_tip_us3 = julday_tip_us3[ind1]
        lai_tip_us3 = lai_tip_us3[ind1]
        fapar_tip_us3 = fapar_tip_us3[ind1]
        lai_std_tip_us3 = lai_std_tip_us3[ind1]
        effssa_vis_us3 = effssa_vis_us3[ind1]
        truebkgdalbedo_vis_us3 = truebkgdalbedo_vis_us3[ind1]
        effssa_std_us3 = effssa_std_us3[ind1]
        truebkgdalbedo_std_us3 = truebkgdalbedo_std_us3[ind1]

        #julday_tip = np.vstack((julday_tip_us1, julday_tip_us2))
        #print julday_tip_us2.shape, lai_tip_us2.shape
        if n_site==1:
                return julday_tip_us1, lai_tip_us1, fapar_tip_us1, lai_std_tip_us1, fapar_std_tip_us1,\
                        effssa_vis_us1, truebkgdalbedo_vis_us1, effssa_std_us1, truebkgdalbedo_std_us1
        if n_site==2:
                return julday_tip_us2, lai_tip_us2, fapar_tip_us2, lai_std_tip_us2, fapar_std_tip_us2,\
                        effssa_vis_us2, truebkgdalbedo_vis_us2, effssa_std_us2, truebkgdalbedo_std_us2
        if n_site==3:
                return julday_tip_us3, lai_tip_us3, fapar_tip_us3, lai_std_tip_us3, fapar_std_tip_us3,\
                        effssa_vis_us3, truebkgdalbedo_vis_us3, effssa_std_us3, truebkgdalbedo_std_us3

if __name__=="__main__":
        get_tip()
