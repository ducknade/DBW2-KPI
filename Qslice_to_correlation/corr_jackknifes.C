#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <sys/stat.h>

char* ens_name;
int T;
int N;
int meas_spacing_MD;
int max_md_sep;
int err_max_md_sep;
int num_replicas;
std::vector< std::vector< std::vector<double> > > Qslice;

int jackknife_num;
int jack_rep_exclude;
int exclude_start;
int exclude_end;

bool included(int rep, int md) {
  return !((rep == jack_rep_exclude) && (md >= exclude_start) && (md < exclude_end));
}

void get_Qslice(std::vector<std::string> replica_names) {
  Qslice = std::vector<std::vector<std::vector<double> > >(num_replicas);

  for(int rep = 0; rep < num_replicas; rep++) {
    for(int t = 0; t < T; t++) {
      std::vector<double> hist;
      int rep_N = 0;

      char filename[256];
      sprintf(filename, "../cps_to_Qslice/%s/Qslice_%d.dat", replica_names[rep].c_str(), t);
      std::ifstream f(filename);

      if(!f.is_open()) {
        printf("Couldn't read from file %s\n", filename);
        exit(-1);
      }

      while(f) {
        double Q;
        f >> Q;
        hist.push_back(Q);
        rep_N++;
      }
      f.close();

      Qslice[rep].push_back(hist);

      //printf("Read time slice %d of replica %s: N = %d\n", t, replica_names[rep].c_str(), rep_N);

      // Find the length of the smallest run
      if(rep == 0 || rep_N < N) {
        N = rep_N;
      }
    }
  }

  printf("Global N = %d\n", N);
}

std::vector<std::vector<double> > compute_correlations_periodic() {
  std::vector<std::vector<double> > correlations;

  for(int dt = 0; dt < T; dt++) {
    //printf("dt = %d\n", dt);
    std::vector<double> correlations_dt;

    if(dt > T/2) {
      // use symmetry
      correlations_dt = correlations[T-dt];
    } else {
      for(int md_sep = 0; md_sep <= max_md_sep; md_sep++) {
        double corr = 0;
        int n = 0;
        for(int rep = 0; rep < num_replicas; rep++) {
          for(int t1 = 0; t1 < T; t1++) {
            int t2 = (t1 + dt) % T;
            for(int md = 0; md < N - md_sep; md++) {
              if(included(rep, md) && included(rep, md+md_sep)) {
                corr += Qslice[rep][t1][md] * Qslice[rep][t2][md+md_sep];
                n++;
              }
            }
          }
        }
        corr /= n;
        correlations_dt.push_back(corr);
      }
    }

    correlations.push_back(correlations_dt);
  }

  return correlations;
}

std::vector<std::vector<double> > compute_errors_periodic(const std::vector<std::vector<double> > &correlations) {
  std::vector<std::vector<double> > correlation_errors;

  for(int dt = 0; dt < T; dt++) {
    std::vector<double> err_dt;
    if(dt > T/2) {
      err_dt = correlation_errors[T-dt];
    } else {
      for(int md_sep = 0; md_sep <= max_md_sep; md_sep++) {
        double err = 0;
        for(int v = 0; v < T; v++) {
          for(int k = -err_max_md_sep; k <= err_max_md_sep; k++) {
            err += correlations[v].at(abs(k)) * correlations[v].at(abs(k)); 
            if(abs(k+md_sep) <= err_max_md_sep && abs(k-md_sep) <= err_max_md_sep) {
              err += correlations[(v+dt)%T].at(abs(k+md_sep)) * correlations[(v-dt+T)%T].at(abs(k-md_sep));
            }
          }
        }
        err /= 2*N*T*num_replicas; // there is a 2 because k is double counted?
        err = std::sqrt(err);
        err_dt.push_back(err);
      }
    }
    correlation_errors.push_back(err_dt);
  }

  return correlation_errors;
}

void write_correlations_periodic(std::vector<std::vector<double> > correlations, 
                                 std::vector<std::vector<double> > correlation_errors) {
// Jiqun Tu  
  char filestem[256];
  sprintf(filestem, "./%s", ens_name);
  mkdir(filestem, 0777);

  for(int dt = 0; dt < T; dt++) {
    char filename[256];
    sprintf(filename, "./%s/correlations_jackknife%d_dt%d.dat", ens_name, jackknife_num, dt);
    FILE* f = fopen(filename, "w");
    if(!f) {
      printf("Couldn't open %s to write!\n", filename);
      exit(-1);
    }
    for(int md_sep = 0; md_sep <= max_md_sep; md_sep++) {
      int MD = md_sep * meas_spacing_MD;
      fprintf(f, "%d\t%0.10f\t%0.10f\n", MD, correlations[dt][md_sep], correlation_errors[dt][md_sep]);
    }
    fclose(f);
  }
}


std::vector<std::vector<std::vector<double> > > compute_correlations_open() {
  std::vector<std::vector<std::vector<double> > > correlations;

  for(int t1 = 0; t1 < T; t1++) {
    //printf("t1 = %d\n", t1);
    correlations.push_back(std::vector< std::vector<double> >());
    for(int t2 = 0; t2 < T; t2++) {
      std::vector<double> correlations_t1t2;

      for(int md_sep = 0; md_sep <= max_md_sep; md_sep++) {
        double corr = 0;
        int n = 0;
        for(int rep = 0; rep < num_replicas; rep++) {
          for(int md = 0; md < N - md_sep; md++) {
            if(included(rep, md) && included(rep, md+md_sep)) {
              // fancy estimator: average over all symmetries
              corr += Qslice[rep][t1][md] * Qslice[rep][t2][md+md_sep];
              corr += Qslice[rep][t2][md] * Qslice[rep][t1][md+md_sep];
              corr += Qslice[rep][T-t1-1][md] * Qslice[rep][T-t2-1][md+md_sep];
              corr += Qslice[rep][T-t2-1][md] * Qslice[rep][T-t1-1][md+md_sep];
              n += 4;
            }
          }
        }
        corr /= n;
        correlations_t1t2.push_back(corr);
      }
      correlations[t1].push_back(correlations_t1t2);
    }
  }

  return correlations;
}

// covar(corr(s, t), corr(u, v))
double covariance(const std::vector<std::vector<std::vector<double> > > &correlations,
                  int s, int t, int u, int v, int md_sep) {
  double ret = 0;
  for(int k = -err_max_md_sep; k <= err_max_md_sep; k++) {
    //printf("k = %d, correlations[s][u].size() = %d\n", k, correlations[s][u].size());
    ret += correlations[s][u].at(abs(k)) * correlations[t][v].at(abs(k));
    if(abs(k+md_sep) <= err_max_md_sep && abs(k-md_sep) <= err_max_md_sep) {
      //printf("abs(k+md_sep) = %d, abs(k-md_sep) = %d, correlations[s][v].size() = %d, correlations[t][u].size() = %d\n", 
      //       abs(k+md_sep), abs(k-md_sep), correlations[s][v].size(), correlations[t][u].size());
      ret += correlations[s][v].at(abs(k+md_sep)) * correlations[t][u].at(abs(k-md_sep));
    }
  }
  ret /= N*num_replicas;
  return ret;
}

// var( (1/4) (corr(s, t) + corr(t, s) + corr(T-s-1, T-t-1) + corr(T-t-1, T-s-1)) )
double var_fancy_estimator(const std::vector<std::vector<std::vector<double> > > &correlations,
                           int s, int t, int md_sep) {
  double ret = 0;
  int ss[] = {s, t, T-s-1, T-t-1};
  int ts[] = {t, s, T-t-1, T-s-1};
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      int s1 = ss[i];
      int t1 = ts[i];
      int s2 = ss[j];
      int t2 = ts[j];
      ret += covariance(correlations, s1, t1, s2, t2, md_sep);
    }
  }
  ret /= 16;
  return ret;
}


std::vector<std::vector<std::vector<double> > > compute_errors_open(const std::vector<std::vector<std::vector<double> > > &correlations) {
  std::vector<std::vector<std::vector<double> > > correlation_errors;

  for(int t1 = 0; t1 < T; t1++) {
    //printf("t1 = %d\n", t1);
    correlation_errors.push_back(std::vector< std::vector<double> >());
    for(int t2 = 0; t2 < T; t2++) {
      std::vector<double> err_t1t2;
      for(int md_sep = 0; md_sep <= max_md_sep; md_sep++) {
        //double err = 0;
        //for(int k = -err_max_md_sep; k <= err_max_md_sep; k++) {
        //  err += correlations[t1][t1].at(abs(k)) * correlations[t2][t2].at(abs(k));
        //  if(abs(k+md_sep) <= err_max_md_sep && abs(k-md_sep) <= err_max_md_sep) {
        //    err += correlations[t1][t2].at(abs(k+md_sep)) * correlations[t1][t2].at(abs(k-md_sep));
        //  }
        //}
        //err /= N;
        //err = std::sqrt(err);
        double err = std::sqrt(var_fancy_estimator(correlations, t1, t2, md_sep));
        err_t1t2.push_back(err);
      }
      correlation_errors[t1].push_back(err_t1t2);
    }
  }

  return correlation_errors;
}

void write_correlations_open(std::vector<std::vector<std::vector<double> > > correlations,
                             std::vector<std::vector<std::vector<double> > > correlation_errors) {
  
// Jiqun Tu  
  char filestem[256];
  sprintf(filestem, "./%s", ens_name);
  mkdir(filestem, 0777);

  for(int t1 = 0; t1 < T; t1++) {
    //printf("t1 = %d\n", t1);
    for(int t2 = 0; t2 < T; t2++) {
      char filename[256];
      sprintf(filename, "./%s/correlations_jackknife%d_t%d_t%d.dat", ens_name, jackknife_num, t1, t2);
      FILE* f = fopen(filename, "w");
      if(!f) {
        printf("Couldn't open %s to write!\n", filename);
        exit(-1);
      }
      for(int md_sep = 0; md_sep <= max_md_sep; md_sep++) {
        int MD = md_sep * meas_spacing_MD;
        fprintf(f, "%d\t%0.10f\t%0.10f\n", MD, correlations[t1][t2][md_sep], correlation_errors[t1][t2][md_sep]);
      }
      fclose(f);
    }
  }
}

void run_jackknife_block(bool bc_open) {
  if(bc_open) {
    //printf("Computing correlations...\n");
    std::vector<std::vector<std::vector<double> > > correlations = compute_correlations_open();
    //printf("Computing error bars...\n");
    std::vector<std::vector<std::vector<double> > > errors = compute_errors_open(correlations);
    //printf("Writing correlations...\n");
    write_correlations_open(correlations, errors);
  } else {
    //printf("Computing correlations...\n");
    std::vector<std::vector<double> > correlations = compute_correlations_periodic();
    //printf("Computing error bars...\n");
    std::vector<std::vector<double> > errors = compute_errors_periodic(correlations);
    //printf("Writing correlations...\n");
    write_correlations_periodic(correlations, errors);
  }
}

int main(int argc, char **argv) {
  if(argc < 10) {
    printf("Usage: ./corr_jackknifes bc ens_name T meas_spacing_MD block_size_MD max_md_sep err_max_md_sep num_reps <rep names>\n");
    exit(-1);
  }
  bool bc_open;
  if(strcmp(argv[1], "open") == 0) {
    bc_open = true;
  } else if(strcmp(argv[1], "periodic") == 0) {
    bc_open = false;
  } else {
    printf("Unrecognized bc %s\n", argv[1]);
  }
  ens_name = argv[2];
  T = atoi(argv[3]);
  meas_spacing_MD = atoi(argv[4]);
  int block_size = atoi(argv[5]) / meas_spacing_MD;
  max_md_sep = atoi(argv[6]) / meas_spacing_MD;
  err_max_md_sep = atoi(argv[7]) / meas_spacing_MD;
  num_replicas = atoi(argv[8]);
  if(argc != 9 + num_replicas) {
    printf("Got %d replica names; expected %d\n", argc - 9, num_replicas);
    exit(-1);
  }

  printf("Analyzing correlations on %s (T = %d, boundary conditions = %s).\n", ens_name, T, bc_open ? "open" : "periodic");
  printf("Measurment spacing is %d MD units\n", meas_spacing_MD);
  printf("Block size is %d MD units = %d measurements\n", atoi(argv[5]), block_size);
  printf("Max separation to compute correlations is %d MD units = %d measurements\n", atoi(argv[6]), max_md_sep);
  printf("Max separation in error sum is %d MD units = %d measurements\n", atoi(argv[7]), err_max_md_sep);

  printf("Replicas:\n");
  std::vector<std::string> replica_names;
  for(int i = 9; i < argc; i++) {
    replica_names.push_back(std::string(argv[i]));
    printf("%s\n", argv[i]);
  }

  printf("Reading Qslice data...\n");
  get_Qslice(replica_names);

  printf("Num jackknife blocks per replica = %d\n", N / block_size);

  jackknife_num = 0;
  jack_rep_exclude = -1;
  printf("Doing jackknife block #%d - no exclusion\n", jackknife_num);
  run_jackknife_block(bc_open);

  for(jack_rep_exclude = 0; jack_rep_exclude < num_replicas; jack_rep_exclude++) {
    for(int jack_block_exclude = 0; jack_block_exclude < N / block_size; jack_block_exclude++) {
      jackknife_num++;
      exclude_start = block_size * jack_block_exclude;
      exclude_end = block_size * (jack_block_exclude + 1);
      printf("Doing jackknife block #%d - excluding [%d, %d) from replica %d\n", jackknife_num, exclude_start, exclude_end, jack_rep_exclude);
      run_jackknife_block(bc_open);
    }
  }
}



  
