import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import os
import csv

DGRP_dir = "."  #directory path with file
DGRPlist = os.listdir(DGRP_dir) # make a list with the file name
del DGRP_dir    #delete variable for memory
DGRPlist=sorted(DGRPlist)
#arr = np.array([]) #make empty array
arr1 = np.array([])
for so in DGRPlist:
    if 'csv' in so:  #only csv file
        DGRP = pd.read_csv(so)
        arr = np.array([])
        for i in range(len(DGRP)):
            array = DGRP.iloc[[i]][DGRP < 100].count().sum()
            arr = np.append(arr, array)
        if len(arr) < len(arr1):
            arr1[:len(arr)] += arr
        else:
            c = arr.copy()
            c[:len(arr1)] += arr1
            arr1=c
#plt.figure(figsize=(5,4))
plt.plot(arr1)
plt.xlim(0,250)
plt.xlabel('Frame')
plt.ylabel('Frequency(Counts)')
plt.title('Total Distance Under 100 change for Time')
plt.savefig('DGRP_time_total.png')
        #array = np.genfromtxt(so, delimiter=',')
        #arr = np.append(arr, array)
##filter under 100
#filter = np.where(arr<=100)
#arr2 = arr[filter]


##print max and min
#print(arr.shape) #to check the size of the array
#print(np.min(arr))
#print(np.min(arr[arr>0]))
#print(np.max(arr[arr<800]))


#plt.hist(arr, bins=400)
#plt.xlabel('distance')
#plt.ylabel('frequency')
#plt.title('Distance Histogram: DGRP last 1min filtered')
#plt.savefig('DGRP_last1min_total_distance.png')
