import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import os
social = [73,563,370]
nonsocial = [317,360,707]
DGRP_dir = "/work/kkwon/20220215_Drosophila_cluster_total/"
val_data = []
for so in social:
    arenadir = DGRP_dir + "DGRP%s/" %so
    
    arenalist = os.listdir(arenadir)
    

    for arena in arenalist:
        if os.path.isdir(arenadir+arena):
            dir = arenadir + arena + "/distances/"
            flylist = os.listdir(dir)
            num = len(flylist)
            
            for i in range(num):
                fly1 = pd.read_csv(dir+flylist[i])
                Q1 = len(fly1)//4
                dist1 = fly1.iloc[1:Q1,5:5+num]


                for j in range(num):
                    if dist1.iloc[0,j]!= 0:
           
        
                        val_data_j = dist1.iloc[:,j].tolist()
                        val_data.extend(val_data_j)
                    else: 
                        continue
        else:
            continue


plt.hist(val_data, bins=100)
plt.xlabel('distance')
plt.ylabel('frequency')
plt.title('Distance Histogram: DGRP social')
plt.savefig('DGRP_social_total.png')
