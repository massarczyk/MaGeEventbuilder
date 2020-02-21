
import io
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

###############################################
print("-----------")
print("read")

f = h5py.File("test.hdf5", "r")
#f = h5py.File("daq_testdata.h5", "r")

members = []
f.visit(members.append)

df = {}
for i in range(len(members)):
	dset =f[members[i]]
	print(members[i]) #, " ",dset," ",isinstance(dset,h5py.Dataset))
	#	if isinstance(dset,h5py.Dataset) :
	#		df[members[i]] = dset

#dataframe = pd.DataFrame(df)
#print(dataframe)


print("-----------")
dset = f['daqdata/daqevtno']
dset2 = f['daqdata/waveform_lf/values/flattened_data']
dset3 = f['daqdata/waveform_lf/values/cumulative_length']
for (i, j,k) in zip(dset[0:10],dset2,dset3):
	print(i," ",j," ",k)
print(len(dset), " ", len(dset2))

"""
plt.hist(dset2, bins = np.arange(0,3,0.010)) 
plt.show()
"""

for i in range(len(dset3)):
	print(i," ",dset3[i])
	upperlimit = int(dset3[i])
	if i:
		lowerlimit = int(dset3[i-1])
	else:
		lowerlimit = 0
	wf = dset2[lowerlimit:upperlimit]
	xaxis = np.arange(0,len(wf))
	plt.scatter(xaxis,wf)
	plt.show()
	input("Press Enter to continue...")



