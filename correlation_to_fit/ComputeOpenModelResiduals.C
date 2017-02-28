#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <cassert>

std::vector<double> meas_MD_times;
std::vector<std::vector<std::vector<double> > > meas_corrs;
std::vector<std::vector<std::vector<double> > > meas_corr_errors;

char last_ens_name[250] = "";
int last_jack_block = -1;

bool write_shape = false;

void read_true_correlations(const char* file_format, int T) {
  meas_corrs = std::vector< std::vector<std::vector<double> > >(T, std::vector< std::vector<double> >(T, std::vector<double>()));
  meas_corr_errors = std::vector< std::vector<std::vector<double> > >(T, std::vector< std::vector<double> >(T, std::vector<double>()));

  for(int t1 = 0; t1 < T; t1++) {
    for(int t2 = 0; t2 < T; t2++) {
      char filename[256];
      sprintf(filename, file_format, t1, t2);
      std::ifstream f(filename);
      if(f.bad()) {
        mexErrMsgIdAndTxt("Fitter:file", "couldn't open %s\n", filename);
      } else {
        //printf("reading from %s\n", filename);
      }

      meas_MD_times = std::vector<double>();

      while(f) {
        double MD_time, corr, err;
        f >> MD_time;
        f >> corr;
        f >> err;
        meas_MD_times.push_back(MD_time);
        meas_corrs[t1][t2].push_back(corr);
        meas_corr_errors[t1][t2].push_back(err);
      }
      f.close();
    }
  }
}

std::vector<double> compute_residuals_open(const char* ens_name, int T, 
                                           double delta_MD, int compare_interval, const double *params) {
  std::vector<double> residuals;

  double tauexp = params[0];
  double D[T-1];
  for(int i = 0; i < T/2; i++) {
    D[i] = D[T-2-i] = params[i+1]; // time reversal symmetry
  }
//  double D[T];
//  for(int i = 0; i < T/2; i++) {
//    D[i] = D[T-1-i] = params[i+1]; // time reversal symmetry
//  }

//  double tauexp[T];
//  for(int i = 0; i < T; i++) {
//    tauexp[i] = params[i];
//  }
//  double D[T-1];
//  for(int i = 0; i < T/2; i++) {
//    D[i] = D[T-2-i] = params[T+i];
//  }

  int Ncomps = meas_MD_times.size() - 1;

  for(int t0 = 0; t0 < T/2; t0++) { // only go to T/2 b/c of time reversal symmetry
    // prepare initial state:
    double x[T];
    for(int t = 0; t < T; t++) {
      x[t] = 0.5 * (meas_corrs[t0][t][0] + meas_corrs[T-t0-1][T-t-1][0]); // time reversal symmetry
    }
    x[0] = x[T-1] = 0; // open BC's

    std::vector<FILE*> acfs;
    std::vector<FILE*> acfs_model;
    for(int t = 0; t < T; t++) {
      char filename[256];
      sprintf(filename, "ACFs/%s_t%d_s%d.dat", ens_name, t, t0);
      acfs.push_back(fopen(filename, "w"));
      sprintf(filename, "ACFs/%s_t%d_s%d_model.dat", ens_name, t, t0);
      acfs_model.push_back(fopen(filename, "w"));
    }

    double MD_time = 0.0;
    for(int c = 1; c <= Ncomps; c++) {
      if(write_shape) {
        char filename[256];
        sprintf(filename, "shape/%s_md%d_s%d.dat", ens_name, (int)(0.5 + MD_time), t0);
        FILE* f = fopen(filename, "w");
        for(int t = 0; t < T; t++) {
          fprintf(f, "%d %0.10f %0.10f %0.10f\n", t, meas_corrs[t0][t][c-1], meas_corr_errors[t0][t][c-1], x[t]);
        }
        fclose(f);

        for(int t = 0; t < T; t++) {
          fprintf(acfs[t], "%f %0.10f %0.10f\n", MD_time, meas_corrs[t0][t][c-1], meas_corr_errors[t0][t][c-1]);
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
        for(int t = 1; t < T-1; t++) {
          newX[t] = x[t] + delta_MD * (-1.0/tauexp * x[t])
                         + delta_MD * (D[t]*(x[t+1] - x[t]) - D[t-1]*(x[t] - x[t-1]));
                         //+ delta_MD * (D[t+1]*x[t+1] + D[t-1]*x[t-1] - 2*D[t]*x[t]);
        }
        newX[0] = newX[T-1] = 0; // open BC's
        for(int t = 0; t < T; t++) x[t] = newX[t];
        MD_time += delta_MD;
      }
    
      // Do a comparison
      for(int t = 1; t < T-1; t++) {
        double residual = (x[t] - meas_corrs[t0][t][c]) / meas_corr_errors[t0][t][c];
        //mexPrintf("x = %f, meas = %f, err = %f, residual = %f\n", x[t], meas_corrs[t0][t][c], meas_corr_errors[t0][t][c], residual);
        residuals.push_back(residual);
      }
    }

    for(int t = 0; t < T; t++) {
      fclose(acfs[t]);
      fclose(acfs_model[t]);
    }
  } // t0

  return residuals;
}

std::vector<double> driver(const char* ens_name, int T, int jackknife_num, double *params) {
  // Load measured correlations
  if(strcmp(ens_name, last_ens_name) != 0 || last_jack_block != jackknife_num) {
    //mexPrintf("Need to load correlations for ensemble %s jackknife block %d\n", ens_name, jackknife_num);
    char true_corr_file_format[256];
    sprintf(true_corr_file_format, "../Qslice_to_correlation/%s/correlations_jackknife%d_t%%d_t%%d.dat", ens_name, jackknife_num);
    read_true_correlations(true_corr_file_format, T);
  } else {
    //mexPrintf("Correlations are already loaded for ensemble %s\n", ens_name);
  }
  if(T != meas_corrs.size()) {
    mexErrMsgIdAndTxt("Fitter:cache", "Loaded correlators have wrong T: size = %d, expected %d\n", meas_corrs.size(), T);
  }
  strcpy(last_ens_name, ens_name);
  last_jack_block = jackknife_num;

  const double MD_meas_interval = meas_MD_times[1] - meas_MD_times[0];
  const int compare_interval = 40;
  const double delta_MD = MD_meas_interval / compare_interval;

  //mexPrintf("computing residuals...\n");
  std::vector<double> residuals = compute_residuals_open(ens_name, T, delta_MD, compare_interval, params);
  //mexPrintf("got residuals\n");

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

