 ### Version history
 # 02/2020 created R. Massarczyk (massarczyk@lanl.gov)

call:
python openSimulation.py test.root . -type cal 

  # for monte Carlo
  # - it the timestamp the time of the trigger or minus the baselineoffset (10us) ?
  # - do we wanna split the event type 3 (evt_mc) in calibration and background ?
  # - psa_energy == online energy ?
  # - unixtime in simulation == 0 for all immediate events and the tiem for delayed events
  # - online energy == total energy of detector energy (no wf/ dead layer correction)
  # - why cummulative length ? is the length necessary at all ?
  # in Olis example hf t0 is 76420, this is 76us baseline ???
  
  
  # in general
  # why unixtime, daqclk and decimal_timestamp, for high precision ??? why not only use the sum of the latter two, or unixtime in ns
  # why do we save t0 dt for each event
  # why do we save cumlative length at all just split array in equal parts (or is this not equal)
  # what is the windowing




ch                  : Cryostat*10000 + String*1000 + Detector  : e.g 11207 = in Cryostat 1 the 7th detector in string 12
daqclk              : unixtime[s] without decimals
daqevtno            : eventnumber (the MaGe evntno)
decimal_timestamp   : decimals of the unixtimestamp[s]
evtno               : always 0 ??? (I use the event in this coincidnce window)
evttype             : evt_undef=0,evt_real=1,evt_pulser=2,evt_mc=3,evt_baseline=4
index               : always 0 ??? (I use now overal index of these event (!= evntno)
inverted            : wf inverted (1) or not (0)
muveto              : this channel is in the mu veto ???
muveto_sample       : idk what this is
psa_energy          : total energy of detector energy (in sim pure sum of all steps in this detector and timeframe, dead layer correction)
trigno              : always 1 ???
unixtime            : unixtime
waveform:
  dt                : in ns, sampling rate (10 for hf, 40 for lf)
  t0                : in ns, wf start (7640 for hf, 10 for lf)
    values          
  cumulative_length : in samples, in 1000 (4096) steps for hf(lf)
  flatten_data      : why ???? the website says array<1>{array<1>{real}} but now its an array<1>{real}, thats very unpractical to load into pandas.
                      array<1>{array<1>{real}} would also save the cumulative length


thing to do:
- chuncksize and max files size 
- compression
- implement deadlayer 
- implement siggen library readout
