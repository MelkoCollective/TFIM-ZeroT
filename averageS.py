#!/usr/bin/env python

#import matplotlib.pyplot as plt
#import os
import numpy as np

data= '01.data' 
fin = open(data, 'r')
x = 0
x2 = 0
count = 0
for lines in fin:
    line = lines.split()
    temp = float(line[0])
    x +=  temp
    x2 +=  temp*temp
    count += 1

x = x/count
x2 = x2/count

#variance of the raw swap
Var = x2 - x*x
StDev = np.sqrt(Var/(1.0 * count - 1))

#error propagation: see e.g.
#http://en.wikipedia.org/wiki/Propagation_of_uncertainty#cite_note-harris2003-12

Entropy = -np.log(x)
E_stdev = StDev/x

print Entropy,E_stdev



