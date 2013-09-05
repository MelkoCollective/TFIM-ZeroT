#!/usr/bin/python
'''
Python code to analyze the results in 01.data and return:
-Expectation <swap>
-Error in <swap> (std dev of the mean)
-Expectation -ln(<swap>)
-Error in -ln(<swap>) (propogation of the above error)

Useage:
./analyze.py [filename]
Analyzes the data file [filename] (01.data if none given) and produces the
above estimators and errors

author: Stephen Inglis
2012-03-06
'''

from math import log,fabs
from sys import argv

def loadfile(filename):
	file = open(filename)
	data = []
	for lines in file:
		temp = lines.split()
		data.append([float(i) for i in temp])
	return data

def main1():
	if len(argv) > 1:
		filename = argv[1]
	else:
		filename = '01.data'
	try:
		data = loadfile(filename)
	except IOError:
		print "File %s not found" % filename
		return
	numBin = len(data)
	numSize = len(data[0])
	y = [0 for i in range(len(data[0]))]
	y2 = [0 for i in range(len(data[0]))]
	err = [0 for i in range(len(data[0]))]
	f = [0 for i in range(len(data[0]))]
	fErr = [0 for i in range(len(data[0]))]
	for i in data:
		for j in range(len(i)):
			y[j] += i[j]
			y2[j] += i[j]*i[j]
	for i in range(numSize):
		err[i] = ((numBin*y2[i] - y[i]*y[i])**0.5)/(numBin**1.5)
		y[i] = y[i]/numBin
	for i in range(numSize):
		f[i] = -1.0*log(y[i])
		fErr[i] = fabs(-1.0/y[i] * err[i])
	print "#%13s %14s %14s %14s" % ('<Swap>','d<Swap>','-ln(<Swap>)','d[-ln(<Swap>)]')
	for i in range(numSize):
		print "%0.12f %0.12f %0.12f %0.12f" % (y[i],err[i],f[i],fErr[i])

if __name__ == "__main__":
	main1()
