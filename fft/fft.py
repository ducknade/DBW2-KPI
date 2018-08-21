#!/usr/bin/python
import numpy as np
import math

#ens_name="periodic-24x64x24ID"
#Njack = 76
#T = 64
#tau_limit = 10
#dtau = 1
#fit_min = 6
#fit_max = 32

#ens_name="24x64x12ID-N"
#Njack = 41
#T = 64
#tau_limit = 15
#dtau = 1
#fit_min = 3
#fit_max = 32

#ens_name="periodic-12x32x12ID-B"
#Njack = 86
#T = 32
#tau_limit = 10
#dtau = 1
#fit_min = 3
#fit_max = 16

#ens_name="periodic-16x32"
#Njack = 20
#T = 32
#tau_limit = 15

#ens_name="l2464f211b600m0102m0509m635a"
#Njack = 106
#T = 64
#tau_limit = 6
#dtau = 5
#fit_min = 10
#fit_max = 26

#ens_name="l3296f211b630m0074m037m440a"
#Njack = 51
#T = 96
#tau_limit = 6
#dtau = 6
#fit_min = 0
#fit_max = 32

#ens_name="periodic-8x16"
#Njack = 97
#T = 16
#tau_limit = 10
#dtau = 10
#fit_min = 0
#fit_max = 8

#ens_name="periodic-10x20"
#Njack = 70
#T = 20
#tau_limit = 10
#dtau = 15
#fit_min = 0
#fit_max = 6

#ens_name="periodic-12x24"
#Njack = 54
#T = 24
#tau_limit = 10
#dtau = 21
#fit_min = 0
#fit_max = 8

ens_name="periodic-14x28"
Njack = 37
T = 28
tau_limit = 10
dtau = 28
fit_min = 0
fit_max = 6


filename_stem = "../Qslice_to_correlation/%s/correlations_jackknife" % (ens_name)
#filename_stem = "../Qslice_to_correlation/periodic-24x64x24ID/correlations_jackknife"
#filename_stem = "../Qslice_to_correlation/periodic-16x32/correlations_jackknife" 

#Njack = 106
#Njack = 20

def mean(data):
    return sum(data)/len(data)

def std(data):
    #  print data
  N = len(data)
  return math.sqrt(abs(N/(N-1.0) * (mean([x*x for x in data]) - mean(data)**2)))

def compute(jack):

    #filename_stem = "../Qslice_to_correlation/24x64x12ID-N/correlations_jackknife0_dt"
  #filename_stem = "../Qslice_to_correlation/periodic-24x64x24ID/correlations_jackknife0_dt"

  data = []
  for tau in range(tau_limit):
    data.append( [0. for dt in range(T) ] )

  for dt in range(T): 
    filename = '%s%d_dt%d.dat' %(filename_stem, jack, dt)
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    if (dt-T/2)**2 >= 0:
      for tau in range(tau_limit):
        data[tau][dt] = float( lines[tau].split()[1] )
#  print data  
  fft_data = [ np.fft.fft( data[tau] ) for tau in range(tau_limit) ]

#  print "data"
#  print data[0]
#  print "fft_data"
#  print fft_data
#  m = []
#  for j in range(len(fft_data[0])):
#    m.append( max( [ fft_data[k][j].real for k in range(len(fft_data)) ] ) )
#  for i in range(tau_limit):
#    for j in range(len(fft_data[0])):
#      print "%+12.8f" %(fft_data[i][j].real/m[j]),
#    print "\n" 
#  for n in range(T):
  #    print sum([ (fft_data[n+1][dt+1]-fft_data[n+1][dt])/fft_data[n+1][dt] - (fft_data[n][dt+1]-fft_data[n][dt])/fft_data[n][dt] for dt in range(3) ])
  #    print np.fft.fftfreq(64)[n]+1/128., sum([ ( (fft_data[tau+1][n+1]-fft_data[tau][n+1])/(fft_data[tau+1][n+1]+fft_data[tau][n+1]) - (fft_data[tau+1][n]-fft_data[tau][n])/(fft_data[tau+1][n]+fft_data[tau][n]) ).real for tau in range(3) ])
  #    print np.fft.fftfreq(64)[n], ([ ( (fft_data[tau+1][n]-fft_data[tau][n])*2/(fft_data[tau+1][n]+fft_data[tau][n]) ).real for tau in range(3) ])
#      eff = [ np.log(abs(( (fft_data[tau][n])/(fft_data[tau+1][n]) ).real)) for tau in range(1) ]
#  return [np.log((( (fft_data[1][n])/(fft_data[1+1][n]) ).real)) for n in range(T)]
  return [np.log((( (fft_data[0][n].real)/(fft_data[0+1][n]).real ).real)) for n in range(T)]
#  return [ np.arccosh(abs(( (fft_data[2][n]+fft_data[0][n])/(fft_data[1][n]) ).real)) for n in range(T)]

#compute(0)

lst = []
for jack in range(Njack):
  lst.append( compute(jack) )

for n in range(T):
#for n in range(T/2):
#  print "%8.4e %12.8e %12.8e" % ( np.fft.fftfreq(T)[n], lst[0][n+2]+lst[0][n]-2*lst[0][n+1], np.sqrt(Njack+1)*std( [lst[j][n+2]+lst[j][n]-2*lst[j][n+1] for j in range(1,Njack,1)]) )
  print "%8.4e %12.8e %12.8e" % ( np.fft.fftfreq(T)[n], lst[0][n], np.sqrt(Njack+1)*std([lst[j][n] for j in range(1,Njack,1)]) )
#      print "%+8.4e %12.8e %12.8e %12.8e" % ( np.fft.fftfreq(64)[n], eff[0], eff[1], eff[2] ) 

print "eff mass"
for n in range(T-1):
#for n in range(T/2):
#  print "%8.4e %12.8e %12.8e" % ( np.fft.fftfreq(T)[n], lst[0][n+2]+lst[0][n]-2*lst[0][n+1], np.sqrt(Njack+1)*std( [lst[j][n+2]+lst[j][n]-2*lst[j][n+1] for j in range(1,Njack,1)]) )
  print "%8.4e %12.8e %12.8e" % ( np.fft.fftfreq(T)[n], lst[0][n+1]-lst[0][n], np.sqrt(Njack+1)*std([ lst[j][n+1]-lst[j][n] for j in range(1,Njack,1)]) )
#      print "%+8.4e %12.8e %12.8e %12.8e" % ( np.fft.fftfreq(64)[n], eff[0], eff[1], eff[2] ) 


x = [ 2.-2.*np.cos(k*2.*math.pi) for k in np.fft.fftfreq(T)]
y = [ B/float(dtau) for B in lst[0] ]
w = [ 1./(np.sqrt(Njack+1)*std([lst[j][n] for j in range(1,Njack,1)]))**2 for n in range(T) ]

print x[fit_min:fit_max]
print y[fit_min:fit_max]
print w[fit_min:fit_max]

print np.polyfit(x[fit_min:fit_max],y[fit_min:fit_max],1,w=w[fit_min:fit_max])

# Train the model using the training sets
