#####################################################
#use as 
#python converter.py [inDir] [inFile] [outDir]
#####################################################

#!/usr/bin/env python
import sys, time, os, pywt
import numpy as np
import pandas as pd
import warnings
import h5py
import uproot
import argparse
import matplotlib.pyplot as plt
import itertools


#####################################################
def main(argv):
  start_time = time.time()
  par = argparse.ArgumentParser(description="MaGe ROOT->HDF5 converter")
  par.add_argument("input_file",action="store",type=str, help="input file (required)")
  par.add_argument("out_dir", action="store", type=str, help="output directory (required)")
  par.add_argument("-t", "--test", action="store_true", help="test mode")
  args = par.parse_args()

  inFileName = args.input_file
  outDir = args.out_dir
  outFileName = outDir + '/' + inFileName.split('.')[0] + '.hdf5'
  print ("Converting file  : {}".format(inFileName))
  print ("into file        : {}".format(outFileName))

  # default (no options) -- run conversion
  pandas_frame = uproot_file(inFileName)
  print("# start Event building --- %s seconds ---" % (time.time() - start_time))
  

  
  
  
#####################################################  
def uproot_file(in_file):
  tt = uproot.open(in_file)["fTree/eventSteps"]
  print("... Tree has ", tt.numentries, " entries")
  tt_start = 2 #2
  tt_stop =  10 #tt.numentries-2

  #empty data frame
  df={}
  ###go through the keys
  for key1 in tt.keys():
    if key1.decode("utf-8").find("MGTDataObject") >= 0:
      continue #ignore the data object
    #for steps we wanna have the details added to the pandas frame
    if key1.decode("utf-8").find("fSteps") >= 0:
      for key2 in tt[key1].keys():
        ignoreList = ['MGTDataObject','fProcessName','fTimeOffset','fIsPreStep','fParticleID','fTrackID','fParentTrackID',
        'fCopyNo','fSensVolID','fKineticEnergy','fStepLength','fTotalTrackLength',
        'fX','fY','fZ','fPx','fPy','fPz','fStepNumber','fTrackWeight'] #ignore these keys
        if any(c in key2.decode("utf-8") for c in ignoreList):
          continue  #ignore the data object
        df[key2.decode("utf-8").replace('.','')] = tt[key2].array(entrystart=tt_start, entrystop=tt_stop) #first two and last two entries cause volname to crash 
      continue     #dont add the frame fsteps itself
    #add add NSteps and EventID to df
    df[key1.decode("utf-8")] = tt[key1].array(entrystart=tt_start, entrystop=tt_stop)
  
  
  print("... read data")
  #make it a DataFrame
  dataframe = pd.DataFrame(df)
    
    
  print("... converted to pandas")
  #split the one long byte entry into the individual Processes for each event
  if 'fStepsfProcessName' in dataframe.columns:
    test = dataframe['fStepsfProcessName'].values
    dataframe['fStepsfProcessName'] = splitProcesses(test)
    print("... modified process names")

  #split the one long byte entry into the individual DetectorNumbers for each event
  test = dataframe['fStepsfPhysVolName'].values
  dataframe['fStepsfPhysVolName'] = splitVolumes(test)
  print("... modified volume IDs")
  

  #flatten the panda frame
  dataframe2 = dataframe.set_index(['fEventID','fNSteps']) \
                        .apply(lambda x: x.explode()) \
                        .reset_index() 
  #remove empty steps
  dataframe2 = dataframe2[dataframe2.fStepsfEdep != 0]
  print("... removed Edep = 0 entries and flattened dataframe")
  
  #modify entires into floats and integers
  for column in dataframe2:
    if column.find("fSteps") < 0:
      dataframe2[column] = dataframe2[column].astype('int')
    else:
      dataframe2[column] = dataframe2[column].astype('float64')
  print("... change data types")

   #set the the time to zero (avoids big numbers, not really necessary)
  dataframe2['fStepsfT'] = dataframe2['fStepsfT'] - min(dataframe2['fStepsfT'].values)


  #plt.figure(figsize=(9, 3))
  #hist = dataframe2.hist(bins=100)
  #plt.show()


  print(dataframe2)

  
  return dataframe

  


#####################################################
def splitProcesses(name):  
  #empty array
  newname = []
  n = 0
  for i in name:
    newname.append([])
    first = 0
    last = 0
    for j in i:
      last = last+1
      if j < 31:   #if a number, then its time for the next process
        newname[n].append(i[first:last-1])
        first = last        
    newname[n].append(i[first:last]) #dont forget the last one
    n = n +1 
  return newname

#####################################################
def splitVolumes(name):  
  #empty array
  detectorID = 0
  newname = []
  n = 0
  for i in name:
    newname.append([])
    j = i.split(b'\x3E')
    for detector in j:
      DetectorID = detector.split(b'\x5F')
      id = int(DetectorID[1])*10000 + int(DetectorID[3])*100 + int(DetectorID[5])   #== cryostat*10000 + string*100 + detector, need to add LAr
      newname[n].append(id)
    n = n +1 
  return newname
  




#####################################################

if __name__ == "__main__":
  main(sys.argv[1:])

