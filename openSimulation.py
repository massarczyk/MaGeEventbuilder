##########################################################################################################
#use as 
#python converter.py [inFile] [outDir] -type [type]
#e.g.  python openSimulation.py test.root . -type cal -t
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

  # write to output file
  print("# %s sec : write to output file" % (time.time() - start_time))  
  writeDataToFile(outFileName,pandas_outputframe,args.test)
    
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
    tt_stop = min(18,tt.numentries-2)

    
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
    name_string = dataframe['fStepsfProcessName'].to_numpy()
    dataframe['fStepsfProcessName'] = splitProcesses(name_string)
    print("... modified process names")

  #split the one long byte entry into the individual DetectorNumbers for each event
  name_string = dataframe['fStepsfPhysVolName'].to_numpy()
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
  dataframe2['fStepsfT'] = dataframe2['fStepsfT'] - min(dataframe2['fStepsfT'].to_numpy())
  


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
   # cant uses slashes, naming != naming in hdf5 file
  output_frame = pd.DataFrame(columns=['ch','daqclk','daqevtno','decimal_timestamp','evtno','evttype','inverted',
                                       'muveto','muveto_sample','psa_energy', 'trigno','unixtime',
                                       'waveform_hf_t0','waveform_hf_dt','waveform_hf_length','waveform_hf_data',
                                       'waveform_lf_t0','waveform_lf_dt','waveform_lf_length','waveform_lf_data']) 
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
    subset_frame['fStepsfT'] = subset_frame['fStepsfT'] - min(subset_frame['fStepsfT'].to_numpy())
    
    output_frame = output_frame.append(makeEvents(subset_frame,output_frame,evttype), ignore_index=True) #attach the individual Events to the total output list
    pbar.update(subset_frame.shape[0])   
    
    index += subset_frame.shape[0]    

    
    #return when we went throuh all Events
    if(index >= data_frame.shape[0]):
      break

  #get the cumulatative lengths as needed for the output format
  output_frame['waveform_hf_length'] = output_frame.waveform_hf_length.cumsum()
  output_frame['waveform_lf_length'] = output_frame.waveform_lf_length.cumsum()


  return output_frame

##########################################################################################################
def makeEvents(stepList,outFrame,COrB): 
  Event_frame = pd.DataFrame().reindex_like(outFrame).dropna() #copy the structure of the outputframe in an empty EVENT frame

  ######paramter of the wf / Event building
  
  daqparameter = {
    'daq_hf_clock'   : 10,    #10 ns daqticks
    'daq_hf_t0'      : 76420, #where wf starts
    'daq_hf_length'  : 1000,  #10 daqticks
    'daq_lf_clock'   : 40,    #10 ns daqticks
    'daq_lf_t0'      : 0,     #where wf starts
    'daq_lf_length'  : 4000}  #10 daqticks
  dtime = daqparameter['daq_lf_length']*daqparameter['daq_lf_clock']/2 # 160 us as in GERDA, signal in the middle
  
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
      unixtime = subset_frame['fStepsfT'][index2]
      eventtime = int(unixtime / 1E9) #ns to s
      eventtimedec = unixtime - eventtime*1E9
      unixtime = unixtime / 1E9 #bring it to s for same format

      #print("----",detector,"-",index," ",index2, " ",n)  
      #this frame is all the steps of one detector within the timewindow
      one_det_frame = subset_frame.loc[subset_frame.fStepsfPhysVolName == detector].reset_index(drop=True)
      wf_hf, wf_lf = makeWaveform(one_det_frame,daqparameter)
      eventenergy = one_det_frame['fStepsfEdep'].sum() * 1000 #MeV to keV
      index2 += one_det_frame.shape[0]   
      
      #here we add the data for each event to the panda ferame of this detector at this time in this event
      Event_frame = Event_frame.append({'ch':detector,
                                        'daqclk':eventtime,
                                        'daqevtno':eventno,
                                        'decimal_timestamp':eventtimedec,
                                        'evtno':n,
                                        'evttype':COrB,
                                        'inverted':0,
                                        'muveto':0,
                                        'muveto_sample':0,
                                        'psa_energy':eventenergy,
                                        'trigno':n,
                                        'unixtime':unixtime,
                                        'waveform_hf_t0':daqparameter['daq_hf_t0'],
                                        'waveform_hf_dt':daqparameter['daq_hf_clock'],
                                        'waveform_hf_length':daqparameter['daq_hf_length'],
                                        'waveform_hf_data':wf_hf,
                                        'waveform_lf_t0':daqparameter['daq_lf_t0'],
                                        'waveform_lf_dt':daqparameter['daq_lf_clock'],
                                        'waveform_lf_length':daqparameter['daq_lf_length'],
                                        'waveform_lf_data':wf_lf
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
def makeWaveform(stepFrame, daqpar):
  #replace the global time with the travel time in the detector
  stepFrame['fStepsfT'] = stepFrame.apply(distance, axis=1)
  signalstarttime = daqpar['daq_lf_length']*daqpar['daq_lf_clock']/2
  
  print(stepFrame['fStepsfT'])

  
  wf1 = np.zeros(daqpar['daq_hf_length'])
  
  wf2 = np.zeros(daqpar['daq_lf_length'])
  
  
  return wf1,wf2
##########################################################################################################
def distance(x):
  #distance in mm , speed ~0.05 mm/ns
  return (np.sqrt( x['fStepsfLocalX']**2 + x['fStepsfLocalY']**2 + (x['fStepsfLocalZ']+40)**2 ) / 0.05 )

##########################################################################################################
##########################################################################################################
def writeDataToFile(outFile,outputframe,testmode):
   # write pandas frame in a hdf 5 file (! outsource to a new function)
  f = h5py.File(outFile, "w")

  for column in outputframe:
    datatype = 'int64'
    if column.find('waveform_lf') >= 0:
      if column.find('data') > 0:  
        outputname = 'daqdata/waveform_lf/values/flattened_data'
      elif column.find('length') > 0:
        outputname = 'daqdata/waveform_lf/values/cumulative_length'
      else:
        outputname = 'daqdata/waveform_lf/' + column[12:]

    elif column.find('waveform_hf') >= 0:
      if column.find('data') > 0:  
        outputname = 'daqdata/waveform_hf/values/flattened_data'
      elif column.find('length') > 0:
        outputname = 'daqdata/waveform_hf/values/cumulative_length'
      else:
        outputname = 'daqdata/waveform_hf/' + column[12:]

    else:
      outputname = 'daqdata/' + column
      if column.find('psa_energy') >= 0:
        datatype='float64'
      elif column.find('unixtime') >= 0:
        datatype='float64'
      
    
    
    #flatten the wf to data, doesnt do anything to the other columns since they are already
    dataset = outputframe.explode(column).reset_index(drop=True)[column].astype(datatype).to_numpy()
    #dataset = outputframe.explode(column)[column].to_numpy()  
    #dataset = np.arange(len(dataset))
    if testmode:
      print("----------")
      print(column, " ",outputname," ",len(dataset)," ", dataset.shape, " ",datatype)
      print(dataset)


    f.create_dataset(outputname,data=dataset)

  f.create_dataset('daqdata/index',data=outputframe.index)
  f.close() 
  return

##########################################################################################################
##########################################################################################################

if __name__ == "__main__":
  main(sys.argv[1:])

