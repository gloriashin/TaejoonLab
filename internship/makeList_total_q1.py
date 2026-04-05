import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import os

DGRP_dir = "/work/kkwon/20220215_Drosophila_cluster_total/"  #directory with file
DGRPlist = os.listdir(DGRP_dir)  #make a list with DGRP name

for so in DGRPlist:
    if os.path.isdir(DGRP_dir+so) :
        arenadir = DGRP_dir + "%s/" %so  #directory with arena folder
        
        arenalist = os.listdir(arenadir)
        arenalist = sorted(arenalist)
        df = pd.DataFrame()

        for arena in arenalist:
            
            if os.path.isdir(arenadir+arena):
                dir = arenadir + arena + "/distances/"
                flylist = os.listdir(dir)
                flylist = sorted(flylist)
                num = len(flylist)
                
                for i in range(num):
                    dist1 = pd.read_csv(dir+flylist[i])

                    Q1 = len(dist1)//4
                    dist1 = dist1.iloc[1:Q1,6+i:5+num]
                    #dist1 = dist1.drop(dist1.columns[i], axis=1) 
                    df =  pd.concat([df,dist1],axis=1)  



            
            else:
                continue
        df.to_csv('reads_Q1_merged_%s.csv' %so, sep=',', na_rep='NaN', index = None, header = False)
    else:
        continue


