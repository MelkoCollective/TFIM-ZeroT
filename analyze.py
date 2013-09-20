#!/usr/bin/python

import glob
import numpy as np

def sortNames(fold_list,num_skip=1):
    return sorted(fold_list,key=lambda k: int(k[num_skip:]))

def extractData(filelist):
    y = 0
    y2 = 0
    N = 0
    count = 0
    for f in filelist:
        data = np.loadtxt(f)
        y += np.mean(data)
        y2 += np.mean(data**2)
        N += 1
        count += len(data)
    y /= N
    y2 /= N
    err = (y2 - y**2)/(count**0.5)
    return y,err

# Assuming max_add = 8, define what counting number refers to complete circles
# as a function of their radius.
# Using a crude assumption that radius = max width/2
circle_map = { 0:0.5,
1:1.5,
3:2.5,
5:3.5,
8:4.5
}

def main():
    folder_fmt = 'L%03d/r%03d/R%03d'
    filename = '01.data'
    L_folders = glob.glob('L*')
    L_folders = sortNames(L_folders)
    all_data = []
    for iL in L_folders:
        nL = int(iL[1:])
        r_folders = glob.glob(iL + '/r*')
        r_folders = sortNames(r_folders,6)
        y_total = 0
        err2_total = 0
        for ir in r_folders:
            nr = int(ir[6:])
            samples = glob.glob(ir + '/R*/01.data')
            y,err = extractData(samples)
            y_total += -1.*np.log(y)
            err2_total += (err/y)**2
            try:
                radius = circle_map[nr]
                all_data.append([nL,radius,y_total,err2_total**(0.5)])
            except KeyError:
                pass

if __name__ == '__main__':
    main()
