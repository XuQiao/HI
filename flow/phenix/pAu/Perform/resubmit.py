#!/bin/python3

import os
import re

dir = "/scratch/xuq7/phenix/pAu/Perform/"
arr = os.listdir(dir)
arr = [x for x in arr if re.match(r'AnapAumbst\_\d{1,3}.root',x)]
#arr = [x for x in arr]
arr = list(map(lambda x:dir+x,arr))
bad = []
for i in range(len(arr)):
#    print(os.path.getsize(arr[i]))
     size = os.path.getsize(arr[i])
     if size < 1000000:
        badnumber = int(arr[i].split('_')[1].split('.')[0])
        print('./submit '+str(badnumber))
        os.system('./submit.sh '+str(badnumber))
        bad.append(badnumber)
#print(bad)
print("bad=%d" % len(bad))

