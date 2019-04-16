#!/usr/bin/python
from __future__ import print_function
import numpy as np
import math

import os
import sys
import fnmatch


#ens_name="periodic-24x64x24ID"
#Njack = 76
#T = 64
#tau_limit = 10
#dtau = 1
#fit_min = 5
#fit_max = 20

#ens_name="24x64x24ID-GL"
#Njack = 19
#T = 64
#tau_limit = 10
#dtau = 1
#fit_min = 5
#fit_max = 20

#ens_name="24x64x12ID-N"
#Njack = 41
#T = 64
#tau_limit = 15
#dtau = 1
#fit_min = 2
#fit_max = 19

#ens_name="periodic-12x32x12ID-B"
#Njack = 87
#T = 32
#tau_limit = 10
#dtau = 1
#fit_min = 2
#fit_max = 9

#ens_name="periodic-16x32"
#Njack = 20
#T = 32
#tau_limit = 15
#dtau = 1
#fit_min = 5
#fit_max = 20

#ens_name="l2464f211b600m0102m0509m635a"
#Njack = 106
#T = 64
#tau_limit = 6
#dtau = 5
#fit_min = 3
#fit_max = 15

#ens_name="l3296f211b630m0074m037m440a"
#Njack = 51
#T = 96
#tau_limit = 6
#dtau = 6
#fit_min = 0
#fit_max = 23

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

#ens_name="periodic-14x28"
#Njack = 37
#T = 28
#tau_limit = 10
#dtau = 28
#fit_min = 0
#fit_max = 6

ens_name="Beta6.4_t1.0_K0.9"
Njack = 14
T = 32
tau_limit = 10
dtau = 5
fit_min = 0
fit_max = 6

def generate_gnuplot_script(ens, dt, num_col, diff=False):
  script_path="plots/%s_adaptive.gnuplot" %ens
  if diff:
    script_path="plots/%s_difference.gnuplot" %ens
  with open(script_path, 'w') as f:
    print("set terminal epslatex size 20.00cm, 10.00cm color lw 2 standalone colortext font 10 \\", file=f)
    print("	header '\usepackage{grffile} \usepackage{amssymb}'", file=f)
    print("", file=f)
    print("set xtics", file=f)
    print("set ytics", file=f)
    if diff:
      print("set ylabel '$\log[C(k,0)/C(k,%d)]/%d-\log[C(k,%d)/C(k,%d\cdot2)]/%d$'" %(dt, dt, dt, dt, dt), file=f)
    else:
      print("set ylabel '$\log[C(k,0)/C(k,%d)]/%d$'" %(dt, dt), file=f)
    print("set xlabel '$2-2\cos k$'", file=f)
    print("", file=f)
    print("set xrange [-0.1:5.8]", file=f)
    if diff:
      print("set yrange [-0.25:0.65]", file=f)
    else:
      print("set yrange [-0.05:0.85]", file=f)
    print("", file=f)
    print("set bars small", file=f)
    print("", file=f)
    if diff:
      print("outfile = 'plots/%s_difference'    # Also used for filename generation so set it" %ens, file=f)
    else:
      print("outfile = 'plots/%s_adaptive'    # Also used for filename generation so set it" %ens, file=f)
    print("", file=f)
    print("set key spacing 1.4", file=f)
    print("set key enhanced", file=f)
    print("set key top right", file=f)
    print("", file=f)
    print("set output outfile.'.tex'", file=f)
    print("", file=f)
    print("two_pi = 6.28318530718", file=f)
    print("two_pi = 6.28318530718", file=f)
    print("", file=f)
    if diff:
      print("set label 1 at graph 0.1, graph 0.9 '\\texttt{%s diff.}' " %ens.replace("_","\\_"), file=f)
      print("infile = 'data/%s_difference.dat'" %ens, file=f)
    else:
      print("set label 1 at graph 0.1, graph 0.9 '\\texttt{%s value}' " %ens.replace("_","\\_"), file=f)
      print("infile = 'data/%s_adaptive.dat'" %ens, file=f)
    print("", file=f)
    print("plot for [i=1:%d:1] infile using (2-2*cos($1*two_pi)+0.04*i):(column(2*i)):(column(2*i+1)) with ye lw 1.5 pt 7 title columnheader(i*2), \\" %num_col, file=f)
    print("", file=f)
    print("set output # Need to close the output file before latex can work with it", file=f)
    print("", file=f)
    print("system 'pdflatex -shell-escape -interaction=batchmode -output-directory=plots '.outfile.\\", file=f)
    print("	'.tex > /dev/null; rm '.outfile.'-inc.eps '.outfile.\\", file=f)
    print("		'-inc-eps-converted-to.pdf '.outfile.'.tex '.\\", file=f)
    print("			outfile.'.log '.outfile.'.aux'\\", file=f)
    print("", file=f)
    print("system 'open '.outfile.'.pdf'", file=f)

def round_sig(x, nsig):
  return round(x, nsig-int(np.floor(np.log10(abs(x))))-1)

def format_number_with_err(cv, err):
  print("%.6f %.6f" %(cv, err))
  nord = int(np.floor(np.log10(abs(cv))))
  nsig = 1-int(np.floor(np.log10(abs(err))))
  tmp_err = err*10.0**nsig
  nsig += nord + 1
  cv_str = str(round_sig(cv,nsig))
  cv_str_len = len(cv_str)
  if nord > -1:
    cv_str_len -= 1
  else:
    cv_str_len -= abs(nord) + 1
  for i in range(0,nsig-cv_str_len):
    cv_str += "0"
  return str("${0:s}({1:d})$".format(cv_str,int(round_sig(tmp_err,2))))

def mean(data):
    return sum(data)/len(data)

def std(data):
  N = len(data)
  return math.sqrt(abs(N/(N-1.0) * (mean([x*x for x in data]) - mean(data)**2)))

def cv_err_jackknife(data):
  n = len(data) - 1
  jkn_avg = mean(data[1:n+1])
  cv = n*data[0]-(n-1)*jkn_avg
  err = np.sqrt(n+1)*std(data[1:n+1])
  return [cv, err]

def compute(filename_stem, jack, T, tau_limit, diff):

  data     = []
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
  fft_data = [ np.fft.fft( data[tau] ) for tau in range(tau_limit) ]

  if diff:
    return [ 0.5*np.log(fft_data[0+0][n].real*fft_data[0+2][n].real/(fft_data[0+1][n].real)**2) if fft_data[0+0][n].real*fft_data[0+2][n].real>0. else np.nan for n in range(T)]
  else:
    return [ np.log(fft_data[0+0][n].real/fft_data[0+1][n].real) if fft_data[0+0][n].real/fft_data[0+1][n].real>0. else np.nan for n in range(T)]

def compute_ens(ens_name, Njack, T, tau_limit, dtau, fit_min, fit_max, label, diff):
  filename_stem = "../Qslice_to_correlation/%s/correlations_jackknife" % (ens_name)
  x = [ 2.-2.*np.cos(k*2.*math.pi) for k in np.fft.fftfreq(T)]

#  D = []
#  t = []
  lst = []
  for jack in range(Njack):
    jlst = compute(filename_stem, jack, T, tau_limit, diff)
    lst.append( [y/float(dtau) for y in jlst] )
#    y = [ B for B in lst[jack] ]
#    pol = np.polyfit(x[fit_min:fit_max],y[fit_min:fit_max],1)
#    D.append(pol[0])
#    t.append(1./pol[1])
#  
#  f = open(label, 'w')
#  
#  f.write( "D: %12.8e %12.8e\n" % ( D[0], np.sqrt(Njack+1)*std( D[1:Njack] ) ) )
#  f.write( "t: %12.8e %12.8e\n" % ( t[0], np.sqrt(Njack+1)*std( t[1:Njack] ) ) )
#  f.write( "[]: %.4f %.4f\n" % (x[fit_min], x[fit_max-1]) )
#  print("%s & %s & %s \\\\" %( label, format_number_with_err(D[0], np.sqrt(Njack+1)*std( D[1:Njack])), format_number_with_err(t[0], np.sqrt(Njack)*std( t[1:Njack])) ) ) 
#  print "D: %12.8e %12.8e" % ( D[0], np.sqrt(Njack+1)*std( D[1:Njack] ) ) 
#  print "t: %12.8e %12.8e" % ( t[0], np.sqrt(Njack+1)*std( t[1:Njack] ) ) 

#  for n in range(T):
#    f.write("%+8.4e %12.8e %12.8e\n" % ( np.fft.fftfreq(T)[n], lst[0][n], np.sqrt(Njack)*std([lst[j][n] for j in range(1,Njack,1)]) ) )
  return [ cv_err_jackknife([lst[j][n] for j in range(Njack)]) for n in range(T)]


def compute_ens_wfs(ens_name, Njack, T, tau_limit, dtau, fit_min, fit_max, label, diff=False):
  fft = []
  num_col = 0
  for dir_path in os.listdir('../Qslice_to_correlation'):
    if fnmatch.fnmatch(dir_path, "%s_wfs*"%(ens_name)):
      # print dir_path 
      wfs = float(dir_path.split("wfs")[1])
      # print wfs
      fft.append([wfs, compute_ens(dir_path, Njack+1, T, tau_limit, dtau, fit_min, fit_max, label, diff)])
      num_col = num_col + 1
  f = open(label, 'w')
  # print fft
  f.write("wfs     ")
  for w in range(len(fft)):
      # f.write( "%6.4f %12.8 %12.8\n" %(np.fft.fftfreq(T)[n], fft[w][1][n][0], fft[w][1][n][1]) )  
      f.write( "$t^2{\\\\langle}E{\\\\rangle}=%.2f$     $wfs%+.2f$     " % (fft[w][0], fft[w][0]) )  
  f.write("\n")
  for n in range(T):
    f.write("%+6.4f " % (np.fft.fftfreq(T)[n]))
    for w in range(len(fft)):
      # f.write( "%6.4f %12.8 %12.8\n" %(np.fft.fftfreq(T)[n], fft[w][1][n][0], fft[w][1][n][1]) )  
      f.write( "%+10.6e %+10.6e " % (fft[w][1][n][0], fft[w][1][n][1]) )  
    f.write("\n")
  return num_col

# ens_name, Njack, T, tau_limit, dtau, fit_min, fit_max
#compute_ens("periodic-12x32x12ID-B", 87, 32, 10, 1, 2, 10, "mdwf1.dat")
#compute_ens("periodic-24x64x24ID", 76, 64, 10, 1, 5, 20, "mdwf2.dat")
#compute_ens_wfs("24x64x12ID-N-adaptive", 40, 64, 15, 1, 2, 22, "mdwf3_adaptive.dat")
#compute_ens("l2464f211b600m0102m0509m635a", 106, 64, 6, 5, 4, 22, "milc1.dat")
#compute_ens("l3296f211b630m0074m037m440a", 51, 96, 6, 6, 0, 24, "milc2.dat")
#
#compute_ens("periodic-8x16", 97, 16, 10, 10, 0, 5, "gauge08.dat")
#compute_ens("periodic-10x20", 70, 20, 10, 15, 0, 6, "gauge10.dat")
#compute_ens("periodic-12x24", 54, 24, 10, 21, 0, 6, "gauge12.dat")
#compute_ens("periodic-14x28", 37, 28, 10, 28, 0, 6, "gauge14.dat")
#compute_ens_wfs("periodic-16x32", 518, 32, 15, 40, 0, 6, "gauge16_adaptive.dat")
#compute_ens("Beta6.4_t1.0_K0.9", 14, 32, 6, 5, 0, 6, "Beta6.4_t1.0_K0.9.dat")
#compute_ens("L32_beta6.4_K0.9_t1.0", 16, 32, 10, 5, 0, 5, "L32_beta6.4_K0.9_t1.0")

if len(sys.argv) != 5: raise Exception("Need at least 5 arguments but %d are given." % len(sys.argv))

ens_stem = sys.argv[1]
print("ensemble stem: %s" % ens_stem)

njack = int(sys.argv[2])
print("number of jackknife samples: %d" % njack)

ens_t = int(sys.argv[3])
print("T = %d" % ens_t)

ens_dtau = int(sys.argv[4])
print("dtau = %d" % ens_dtau)

num_col = compute_ens_wfs(ens_stem, njack, ens_t, 5, ens_dtau, 0, 6, "data/%s_adaptive.dat"%ens_stem, False)
generate_gnuplot_script(ens_stem, ens_dtau, num_col, False)
os.system("gnuplot -c plots/"+ ens_stem +"_adaptive.gnuplot")

compute_ens_wfs(ens_stem, njack, ens_t, 5, ens_dtau, 0, 6, "data/%s_difference.dat"%ens_stem, True)
generate_gnuplot_script(ens_stem, ens_dtau, num_col, True)
os.system("gnuplot -c plots/"+ ens_stem +"_difference.gnuplot")

#ens_name="periodic-16x32"
#Njack = 20
#T = 32
#tau_limit = 15
#dtau = 1
#fit_min = 5
#fit_max = 20

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

#ens_name="periodic-14x28"
#Njack = 37
#T = 28
#tau_limit = 10
#dtau = 28
#fit_min = 0
#fit_max = 6

