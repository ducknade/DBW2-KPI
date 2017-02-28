#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <cassert>

struct Correlator {
  std::vector<double> MD_time;
  std::vector<double> data;
  std::vector<double> error;
};

std::vector<Correlator> true_corrs;
char last_ens_name[250] = "";
int last_jack_block = -1;

bool write_shape = false;

std::vector<Correlator> read_true_correlations(const char* file_format, int T) {
  std::vector<Correlator> ret(T);

  for(int t = 0; t < T; t++) {
    char filename[256];
    sprintf(filename, file_format, t);

    std::ifstream f(filename);
    if(!f.is_open()) {
      mexErrMsgIdAndTxt("Fitter:file", "couldn't open %s\n", filename);
    } else {
      //mexPrintf("reading from %s\n", filename);
    }

    while(f) {
      double MD_time, corr, err;
      f >> MD_time;
      f >> corr;
      f >> err;
      ret[t].MD_time.push_back(MD_time);
      ret[t].data.push_back(corr);
      ret[t].error.push_back(err);
    }
    f.close();
  }

  return ret;
}

std::vector<double> compute_residuals_periodic(const char* ens_name, int T, 
                                          double delta_MD, int compare_interval,
                                          const double *params) {
  std::vector<double> residuals;

  double tauexp = params[0];
  double D = params[1];

  // prepare initial state:
  double x[T];
  for(int t = 0; t < T; t++) {
    x[t] = true_corrs[t].data[0];
  }

  int Ncomps = true_corrs[0].MD_time.size() - 1;
  //mexPrintf("max MD = %f, Ncomps = %f\n", true_corrs[0].MD_time[Ncomps], Ncomps);
  
  std::vector<FILE*> acfs;
  std::vector<FILE*> acfs_model;
  for(int t = 0; t < T; t++) {
    char filename[256];
    sprintf(filename, "ACFs/%s_dt%d.dat", ens_name, t);
    acfs.push_back(fopen(filename, "w"));
    sprintf(filename, "ACFs/%s_dt%d_model.dat", ens_name, t);
    acfs_model.push_back(fopen(filename, "w"));
  }


  double MD_time = 0.0;
  for(int c = 1; c <= Ncomps; c++) {
    if(write_shape) {
      char filename[256];
      sprintf(filename, "shape/%s_md%d.dat", ens_name, (int)(0.5 + MD_time));
      FILE* f = fopen(filename, "w");
      for(int t = 0; t < T; t++) {
        fprintf(f, "%d %0.10f %0.10f %0.10f\n", t, true_corrs[t].data[c-1], true_corrs[t].error[c-1], x[t]);
      }
      fclose(f);

      for(int t = 0; t < T; t++) {
        fprintf(acfs[t], "%f %0.10f %0.10f\n", MD_time, true_corrs[t].data[c-1], true_corrs[t].error[c-1]);
      }
    }

    // Evolve in MD time
    for(int i = 0; i < compare_interval; i++) {
      if(write_shape) {
        for(int t = 0; t < T; t++) {
          fprintf(acfs_model[t], "%f %0.10f\n", MD_time, x[t]);
        }
      }

      double newX[T];
      for(int t = 0; t < T; t++) {
        newX[t] = x[t] + delta_MD * (-1.0/tauexp * x[t])
                       + delta_MD * D * (x[(t+1)%T] - 2*x[t] + x[(t-1+T)%T]);
      }
      for(int t = 0; t < T; t++) x[t] = newX[t];
      MD_time += delta_MD;

    }
    
    // Do a comparison
    for(int t = 0; t < T; t++) {
      double residual = (x[t] - true_corrs[t].data[c]) / true_corrs[t].error[c];
      residuals.push_back(residual);
    }
  }

  for(int t = 0; t < T; t++) {
    fclose(acfs[t]);
    fclose(acfs_model[t]);
  }

  return residuals;
}

std::vector<double> driver(const char* ens_name, int T, int jackknife_num, double *params) {
  // Load measured correlations
  if(strcmp(ens_name, last_ens_name) != 0 || last_jack_block != jackknife_num) {
    //mexPrintf("Need to load correlations for ensemble %s jackknife block %d\n", ens_name, jackknife_num);
    char true_corr_file_format[256];
    sprintf(true_corr_file_format, "/home/gregm/DBW2/analysis/heat/fitting/correlations/%s/correlations_jackknife%d_dt%%d.dat", ens_name, jackknife_num);
    true_corrs = read_true_correlations(true_corr_file_format, T);
  } else {
    //mexPrintf("Correlations are already loaded for ensemble %s\n", ens_name);
  }
  if(T != true_corrs.size()) {
    mexErrMsgIdAndTxt("Fitter:cache", "Loaded correlators have wrong T: size = %d, expected %d\n", true_corrs.size(), T);
  }
  strcpy(last_ens_name, ens_name);
  last_jack_block = jackknife_num;

  // Compute residuals
  const double MD_meas_interval = true_corrs[0].MD_time[1] - true_corrs[0].MD_time[0];
  const int compare_interval = 40;
  const double delta_MD = MD_meas_interval / compare_interval;
  std::vector<double> residuals = compute_residuals_periodic(ens_name, T, delta_MD, compare_interval, params);

  return residuals;
}


/* The gateway function */
void mexFunction( int num_outputs, mxArray *outputs[],
    int num_inputs, const mxArray *inputs[])
{
  write_shape = false;

  /* check for proper number of arguments */
  if(num_inputs != 4) {
    if(num_inputs == 5) {
      write_shape = true;
      mexPrintf("Writing shape.\n");
    } else {
      mexErrMsgIdAndTxt("Fitter:input", "Four inputs required.");
    }
  }
  if(num_outputs != 1) {
    mexErrMsgIdAndTxt("Fitter:input", "One output required.");
  }

  if(!mxIsChar(inputs[0])) {
    mexErrMsgIdAndTxt("Fitter:input", "First input must be a string.");
  }
  const char* ens_name = mxArrayToString(inputs[0]);
  //mexPrintf("ens_name = %s\n", ens_name);

  if( !mxIsDouble(inputs[1]) || 
      mxIsComplex(inputs[1]) ||
      mxGetNumberOfElements(inputs[1])!=1 ) {
    mexErrMsgIdAndTxt("Fitter:input", "Second input must be a scalar.");
  }
  int T = (int) mxGetScalar(inputs[1]);

  if( !mxIsDouble(inputs[2]) || 
      mxIsComplex(inputs[2]) ||
      mxGetNumberOfElements(inputs[2])!=1 ) {
    mexErrMsgIdAndTxt("Fitter:input", "Third input must be a scalar.");
  }
  int jackknife_num = (int) mxGetScalar(inputs[2]);
  //mexPrintf("jackknife_num = %d\n", jackknife_num);

  if(mxGetM(inputs[3]) != 1) {
    mexErrMsgIdAndTxt("Fitter:input", "Fourth input must be a row vector.");
  }
  /* create a pointer to the real data in the input matrix  */
  double* params = mxGetPr(inputs[3]);

  // Compute the residuals
  std::vector<double> residuals = driver(ens_name, T, jackknife_num, params);

  // Create the Matlab output matrix
  outputs[0] = mxCreateDoubleMatrix(1, residuals.size(), mxREAL);

  // Get a pointer to the underlying array for the matlab output matrix
  double* output_residuals = mxGetPr(outputs[0]);

  // Copy the residuals to the Matlab matrix
  memcpy(output_residuals, residuals.data(), residuals.size() * sizeof(residuals[0]));
}

