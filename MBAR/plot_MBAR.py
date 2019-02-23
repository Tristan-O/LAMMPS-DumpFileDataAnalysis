#!/usr/bin/env python
import sys
import numpy as np
from pymbar import timeseries
from pymbar import mbar
import matplotlib.pyplot as plt
import argparse

#Get data file (in .npz format)
parser = argparse.ArgumentParser()
parser.add_argument('-f','--infile', type=str)
parser.add_argument('-b','--binsize', type=float,default = 0.5)
parser.add_argument('-t','--title', type=str, default='noTitleProvided')
args = parser.parse_args()


with open(args.infile, 'rb') as f:
  npzfile = np.load(f)

  #Read data file
  Umat = npzfile['Umat']
  numSamples = npzfile['numSamples']
  allDists1D = npzfile['allDists1D']

#Now create and plot PMF
mbarobj = mbar.MBAR(Umat, numSamples,verbose=True)

deltaGs, deltaGerr, thetaStuff = mbarobj.getFreeEnergyDifferences()

print("\nFree energies between states:")
print(deltaGs[0])
print("\nwith uncertainties:")
print(deltaGerr[0])

#Now also calculate pmf
pmfBins = np.arange(np.min(allDists1D)-(1E-20), np.max(allDists1D), args.binsize)
pmfBinInds = np.digitize(allDists1D, pmfBins) - 1
pmfBinCents = 0.5*(pmfBins[:-1] + pmfBins[1:])
nBins = len(pmfBinCents)

allU = np.zeros( np.shape(Umat) ) 

pmfVals, pmfErr = mbarobj.computePMF(allU, pmfBinInds, nBins)

print('')
print('bin centers')
print(pmfBinCents)
print('pmf Values')
print(pmfVals)
print('Uncertainties')
print(pmfErr)

#Quickly plot and save pmf
pmfVals += 2*np.log(pmfBinCents)
plt.errorbar([ [list(pmfBinCents), list(pmfVals), list(pmfErr)] ], args.title)
plt.show()