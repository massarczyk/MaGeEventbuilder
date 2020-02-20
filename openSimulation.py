##########################################################################################################
#use as 
#python converter.py [inDir] [inFile] [outDir]
##########################################################################################################

#!/usr/bin/env python
import sys, time, os, pywt, tqdm, datetime
import numpy as np
import pandas as pd
import warnings
import h5py
import uproot
import argparse
import matplotlib.pyplot as plt
import itertools

##########################################################################################################
# main function

##########################################################################################################
def main(argv):
  start_time = time.time()
  par = argparse.ArgumentParser(description="MaGe ROOT->HDF5 converter")
  par.add_argument("input_file",action="store",type=str, help="input file (required)")
  par.add_argument("out_dir", action="store", type=str, help="output directory (required)")
  par.add_argument("-t", "--test", action="store_true", help="test mode")
  par.add_argument("-type", action="store", help="cal or bkg", dest='data_type')
  args = par.parse_args()



  inFileName = args.input_file
  outDir = args.out_dir
  outFileName = outDir + '/' + inFileName.split('.')[0] + '.hdf5'
  print ("# Converting file  : {}".format(inFileName))
  print ("# into file        : {}".format(outFileName))

  # read the file and get needed data
  print("# %s sec : start Reading file" % (time.time() - start_time))
  pandas_inputframe = uproot_file(inFileName,args.test) #return a step-wise panda data frame from the file

  # build individual events according to data structure
  print("# %s sec : start Event building" % (time.time() - start_time))  
  pandas_outputframe = event_building(pandas_inputframe,args.test,args.data_type)


   # write pandas frame in a hdf 5 file (! outsource to a new function)
  f = h5py.File(outFileName, "w")
  for column in pandas_outputframe:

    if column.find('waveform_lf') > 0:
      outputname = 'daqdata/waveform_lf/' + column
    elif column.find('waveform_hf') > 0:
      outputname = 'daqdata/waveform_hf/' + column
    else:
      outputname = 'daqdata/' + column
    print(column, " ",outputname)
    f.create_dataset(outputname,data=pandas_outputframe[column].values)

  f.create_dataset('daqdata/index',data=pandas_outputframe.index)
  f.close() 
    
  print("# %s sec :  end" % (time.time() - start_time))  

  
########################################################################################################## 
# this function opens the root file with uproot
# it transfers the data in a pandas dataframe
# processes and volumenames are one long string even the events have a lot of steps
# I split up the string array outside the panda frame as individual array (was the fastest)
# if additional MaGe values flags are needed/wished remove the items from the ignoreList 
# (keep the frame as small as possible to speed up)
# Event data frame is flattened at the end to optimize use of pandas 
# (fixed structure over flexible vector with length = number of steps)
##########################################################################################################  
def uproot_file(in_file,testmode):
  tt = uproot.open(in_file)["fTree/eventSteps"]
  print("... Tree has ", tt.numentries, " Events")
  tt_start = 2 #2
  tt_stop =  tt.numentries-2
  if testmode:
    tt_stop = min(1000,tt.numentries-2)

    
  #empty data frame
  df={}
  ###go through the keys
  for key1 in tt.keys():
    ignoreList = ['MGTDataObject','fNSteps']
    if any(c in key1.decode("utf-8") for c in ignoreList):
      continue #ignore the data object and dont need steps
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
    name_string = dataframe['fStepsfProcessName'].values
    dataframe['fStepsfProcessName'] = splitProcesses(name_string)
    print("... modified process names")

  #split the one long byte entry into the individual DetectorNumbers for each event
  name_string = dataframe['fStepsfPhysVolName'].values
  dataframe['fStepsfPhysVolName'] = splitVolumes(name_string)
  print("... modified volume IDs")
  

  #flatten the panda frame
  dataframe2 = dataframe.set_index(['fEventID']) \
                        .apply(lambda x: x.explode()) \
                        .reset_index()                       
                        
  #remove empty steps
  dataframe2 = dataframe2[dataframe2.fStepsfEdep != 0].reset_index(drop=True)
  print("... removed Edep = 0 entries and flattened dataframe")
  
  
  
  #modify entires into floats and integers
  for column in dataframe2:
    if column.find("fSteps") < 0:
      dataframe2[column] = dataframe2[column].astype('int64')
    elif column.find("fStepsfPhysVolName") == 0:
      dataframe2[column] = dataframe2[column].astype('int64')
    else:
      dataframe2[column] = dataframe2[column].astype('float64')
    #dataframe2[column].reset_index()[column].unique()
  print("... change data types")

   #set the the time to zero (avoids big numbers, not really necessary)
  dataframe2['fStepsfT'] = dataframe2['fStepsfT'] - min(dataframe2['fStepsfT'].values)
  


  if testmode:
    pd.set_option('display.max_rows', dataframe2.shape[0]+1)
    print(dataframe2)
  
  return dataframe2

  


##########################################################################################################
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

##########################################################################################################
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

########################################################################################################## 
# this function first breaks the individual (known) MaGe Events in to individual steps
# Steps within one Geant4-Event are time sorted
# the timesorted array is used to take all the steps in a coincidence window (makeEvents function) 
# and push them in a smaller panda frame
# all steps within this smaller frame are sorted by detector
# then all a sub-sub frame is created only consisting events in a timewindow of one detector 
# (inner while structure)
# this sub-sub frame is pushed towards the makeWafevorm function in order to return a waveform
##########################################################################################################  
def event_building(data_frame,testmode,CalOrBkg):
  # as in https://github.com/orgs/legend-exp/teams/simulations-and-analysis/discussions/12?from_comment=12
  
  output_frame = pd.DataFrame(columns=['ch','daqclk','daqevtno','decimal_timestamp','evtno','evttype','inverted','muveto','muveto_sample','psa_energy', 'trigno','unixtime']) 
  if CalOrBkg.find("cal"):
    evttype = 3
  else:
    evttype = 3
    
  
  print("... processing events")
  pbar = tqdm.tqdm(total=data_frame.shape[0]) #,miniters=int(1+0.005*data_frame.shape[0]))
  
  index = 0
  while True:
    #pick only the same event number and sort it by time
    subset_frame = data_frame.loc[(data_frame.fEventID == data_frame.fEventID[index])].sort_values(by=['fStepsfT']).reset_index(drop=True)
    #set time to 0 , not so necesary
    subset_frame['fStepsfT'] = subset_frame['fStepsfT'] - min(subset_frame['fStepsfT'].values)
    
    output_frame = output_frame.append(makeEvents(subset_frame,output_frame,evttype), ignore_index=True) #attach the individual Events to the total output list
    pbar.update(subset_frame.shape[0])   
    
    index += subset_frame.shape[0]    

    
    #return when we went throuh all Events
    if(index >= data_frame.shape[0]):
      break

   
  return output_frame

##########################################################################################################
def makeEvents(stepList,outFrame,COrB): 
  Event_frame = pd.DataFrame().reindex_like(outFrame).dropna() #copy the structure of the outputframe in an empty EVENT frame

  ######paramter of the wf / Event building
  dtime = 150000 # 160 us as in GERDA - 10us offset
  daqclock = 10  #10 daqticks

  
  time = 0
  index = 0

  while True:  
    #split into subsets that are in coincidence  
    subset_frame = stepList.loc[(stepList.fStepsfT >= time) & (stepList.fStepsfT < time+dtime)].sort_values(by=['fStepsfPhysVolName'])
    subset_frame = subset_frame.reset_index(drop=True)
    #now handle the indvidual detectors (det sorted!)
    index2 = 0
    #print(subset_frame)
    n=0
    while True:
      detector = subset_frame['fStepsfPhysVolName'][index2]
      eventno = subset_frame['fEventID'][index2]
      eventtime = subset_frame['fStepsfT'][index2]
      timer = datetime.datetime.utcnow().timestamp()
      #print("----",detector,"-",eventno," ", timer )  
      #this frame is all the steps of one detector within the timewindow
      one_det_frame = subset_frame.loc[subset_frame.fStepsfPhysVolName == detector].reset_index(drop=True)

      eventenergy = one_det_frame['fStepsfEdep'].sum()
      index2 += one_det_frame.shape[0]   
      
      #here we add the data for each event to the panda ferame of this detector at this time in this event
      Event_frame = Event_frame.append({'ch':detector,
                                        'daqclk':daqclock,
                                        'daqevtno':eventno,
                                        'decimal_timestamp':eventtime,
                                        'evtno':n,
                                        'evttype':COrB,
                                        'inverted':0,
                                        'muveto':0,
                                        'muveto_sample':0,
                                        'psa_energy':eventenergy,
                                        'trigno':n,
                                        'unixtime':timer
                                        }
                                        ,ignore_index=True)

      #return when we went throuh all steps of this detecotr at this time
      n = n + 1
      if(index2 >= subset_frame.shape[0]): 
        break
    
    #update index
    index += subset_frame.shape[0]  
    #break if we are at the end of the list
    if(index >= stepList.shape[0]):
      break 
    #update time if not at the end of list   
    time = stepList['fStepsfT'][index] #set time to the next event outside the coincidence window (time sorted array)
  
  return Event_frame
##########################################################################################################
def makeWaveform(stepFrame, daqclock, ):
  
  return 


##########################################################################################################
##########################################################################################################

if __name__ == "__main__":
  main(sys.argv[1:])

