
import io
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def printname(name):
	print(name)
###############################################
print("-----------")
print("read")

f = h5py.File("test.hdf5", "r")
print(list(f.keys()))
f.visit(printname)



print("-----------")
dset = f['daqdata/daqevtno']
dset2 = f['daqdata/ch']
dset3 = f['daqdata/decimal_timestamp']
dset4 = f['daqdata/psa_energy']
for (i, j,k) in zip(dset[0:10],dset2,dset4):
	print(i," ",j," ",k)
print(dset.dtype)

plt.hist(dset4, bins = np.arange(0,3,0.01)) 
plt.show()
