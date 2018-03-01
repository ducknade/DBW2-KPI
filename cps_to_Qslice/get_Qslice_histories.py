#!/usr/bin/python
import analyze
import os
import sys
import time
import math

# argument format: --open/periodic #ens_name #therm_trajs #wflow_dir
# e.x., ./get_Qslices_history --open open-10x20 300 dir/to/wflow

if len(sys.argv) != 5: raise Exception('need 4 arguments but %d are given.' % len(sys.argv))
if sys.argv[1] == '--periodic': bc_open = False
elif sys.argv[1] == '--open': bc_open = True
else: raise Exception('unknown command line option %s' % sys.argv[1])

ens_name = sys.argv[2]
print "ens_name =", ens_name

#TODO: consider this
therm_trajs = int(sys.argv[3])

print "Skipping first %d trajectories" % therm_trajs

folder = sys.argv[4]
print "reading from %s" % folder

# Jiqun Tu
if not os.path.exists('%s' % ens_name):
    os.makedirs('%s' % ens_name)

def read_results(filename):
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()

  wft = 0

  wfts = []
  Es = []
  sliceEs = []
  Qs = []
  sliceQs = []

  c5 = 1./20. 

  i = 0
  while i < len(lines):
    if "AlgWilsonFlow" in lines[i]:
      dwft = float(lines[i].split()[3])
      wft += dwft
      i += 1
      continue

    wfts.append(wft)

    E = float(lines[i].strip())

    Es.append(E)
    sliceE = [float(x) for x in lines[i+1].split()]
    sliceEs.append(sliceE)

    Q_11 = float(lines[i+9].split()[3])
    Q_12 = float(lines[i+10].split()[3])
    Q_22 = float(lines[i+11].split()[3])
    Q_33 = float(lines[i+12].split()[3])
    Q_13 = float(lines[i+13].split()[3])
    Q = analyze.improved_Q(Q_11, Q_12, Q_22, Q_33, Q_13, c5)
    Qs.append(Q)

    sQ_11 = [float(x) for x in lines[i+14].split()[2:]]
    sQ_12 = [float(x) for x in lines[i+15].split()[2:]]
    sQ_22 = [float(x) for x in lines[i+16].split()[2:]]
    sQ_33 = [float(x) for x in lines[i+17].split()[2:]]
    sQ_13 = [float(x) for x in lines[i+18].split()[2:]]
    sliceQ = [analyze.improved_Q(sQ_11[a], sQ_12[a], sQ_22[a], sQ_33[a], sQ_13[a], c5) for a in range(len(sQ_11))]
    sliceQs.append(sliceQ)

    i += 19

  return (wfts, Es, sliceEs, Qs, sliceQs)

now = time.time()
def file_age(filename):
  return now - os.path.getmtime(filename)

print "Counting results files..."
# folder = '/home/gregm/DBW2/%s/results/alg_wflow' % ens_name
confs = [int(f[6:]) for f in os.listdir(folder) if f[:6] == 'wflow.' and file_age(folder + '/' + f) > 60]
confs = [x for x in confs if x >= therm_trajs] #skip thermalization
confs.sort()

if len(confs) == 0:
  print "No configurations yet; skipping\n"
  exit()

dconf = confs[1] - confs[0] #separation in configurations between adjacent measurements
assert all([confs[i+1] - confs[i] == dconf for i in range(len(confs)-1)])

print "%d results files from %d configurations" % (len(confs), max(confs))

print "Reading results files..."
all_results = [read_results(folder + '/wflow.%d' % conf) for conf in confs]
print "Read results files"

all_wfts, all_Es, all_sliceEs, all_Qs, all_sliceQs = zip(*all_results)

wfts = all_wfts[0] #The wilson flow times at which measurements have been made
Nt = len(all_sliceEs[0][0]) #number of time slices
T = Nt
L = Nt/2
ts = [i for i in range(Nt)]

Ehists = dict(zip(wfts, zip(*all_Es))) #Ehists[t] is the history of E measurements at flow time t
Qhists = dict(zip(wfts, zip(*all_Qs))) #Qhists[t] is the history of Q measurements of flow time Q

sliceEhists = {} #sliceEhists[wft][t][c] is the E measured at wilson flow time wft on slice t on conf c
sliceQhists = {} #sliceQhists[wft][t][c] is the Q measured at wilson flow time wft on slice t on conf c
for w in range(len(wfts)):
  sliceEhists[wfts[w]] = [[all_sliceEs[c][w][t] for c in range(len(confs))] for t in range(Nt)]
  sliceQhists[wfts[w]] = [[all_sliceQs[c][w][t] for c in range(len(confs))] for t in range(Nt)]

def get_central_Qhist(wft, nslices): 
  tmin = Nt/2 - nslices/2
  tmax = nslices + tmin
  return analyze.list_sum(sliceQhists[wft][tmin:tmax])

def get_central_Ehist(wft, nslices): 
  tmin = Nt/2 - nslices/2
  tmax = nslices + tmin
  return analyze.list_mean(sliceEhists[wft][tmin:tmax])

#Compute the reference flow times t0 and w0
W = Nt/4 if bc_open else 0
mean_Es = [analyze.mean(get_central_Ehist(wft, Nt-2*W)) for wft in wfts]
t0 = analyze.flow_scale_luscher(wfts, mean_Es)
w0 = analyze.flow_scale_BMW(wfts, mean_Es)

#Pick a reference Wilson flow time, and interpolate Q and E measurements to that
#flow time
wft_ref = t0**2 #TODO: consider switching to w0**2
(wft_1, wft_2) = next((wfts[i], wfts[i+1]) for i in range(len(wfts)-1) if wfts[i+1] > wft_ref)
interp = (wft_ref - wft_1)/(wft_2 - wft_1)
def interpolate_lists(interp, data1, data2):
  return [(1-interp)*x1 + interp*x2 for (x1, x2) in zip(data1, data2)]
Ehist_ref = interpolate_lists(interp, Ehists[wft_1], Ehists[wft_2])
Qhist_ref = interpolate_lists(interp, Qhists[wft_1], Qhists[wft_2])
sliceEhists_ref = [interpolate_lists(interp, sliceEhists[wft_1][t], sliceEhists[wft_2][t]) for t in range(Nt)]
sliceQhists_ref = [interpolate_lists(interp, sliceQhists[wft_1][t], sliceQhists[wft_2][t]) for t in range(Nt)]

for t in ts:
  print "Writing history file for slice %d" % t
  analyze.write_gnuplot_file('%s/Qslice_%d.dat' % (ens_name, t), [sliceQhists_ref[t]])
