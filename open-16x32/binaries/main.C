#include <config.h>
#include <util/lattice.h>
#include <alg/do_arg.h>
#include <alg/no_arg.h>
#include <alg/alg_plaq.h>
#include <comms/sysfunc_cps.h>
#include <alg/alg_hmc.h>
#include <alg/alg_int.h>
#include <alg/hmc_arg.h>
#include <alg/int_arg.h>
#include <alg/hmd_arg.h>
#include <util/time_cps.h>
#include <alg/alg_smear.h>
#include <alg/alg_tcharge.h>
#include <util/lat_cont.h>
#include <alg/alg_actiondensity.h>
#include <alg/alg_wilsonflow.h>
#include <fenv.h>

#include <unistd.h>

USING_NAMESPACE_CPS

static const char* cname = "";

//a macro to load an argument struct variable from a VML file of the same name
#define decode_vml(arg_name)  do{                                      \
    if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )               \
      ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");           \
  } while(0)

//a macro to save an argument struct variable to a VML file of the same name
#undef encode_vml
#define encode_vml(arg_name, traj) do{                                  \
    char vml_file[256];                                                 \
    sprintf(vml_file, #arg_name".%d", traj);                            \
    if( !arg_name.Encode(vml_file, #arg_name) ){                        \
      ERR.General(cname, fname, #arg_name " encoding failed.\n");       \
    }                                                                   \
  }while(0)


// global argument structures
DoArg do_arg;
EvoArg evo_arg;
HmcArg hmc_arg;
ActionGaugeArg gauge_arg;
IntABArg ab1_arg;

// where to put measurement output files
static const char* hmc_stem = "../results/alg_hmc/hmc";
static const char* plaq_stem = "../results/alg_plaq/plaq";
static const char* wflow_stem = "../results/alg_wflow/wflow";

const int trajs_per_run = 100;

// controls the Wilson flow measurements
const int wflow_trajs_per_measurement = 20;
const Float wflow_dt = 0.05;
const int wflow_nsteps = 80;
const int wflow_measure_start_step = 0;
const int wflow_measure_interval = 1;


void setup(int *argc, char **argv[]);
void do_all_measurements(int traj);
void measure_plaq(int traj);
void do_wflow_measurements(int traj,
                          Float dt,
                          int num_time_steps,
                          int measure_start_step,
                          int measure_interval);
void run_hmc_trajectory(int traj, AlgIntAB &top_level_integrator);
void checkpoint(int traj);
void erase_file(const char* filename);


int main(int argc, char *argv[])
{
  const char* fname = "main()";

  setup(&argc, &argv);

  GJP.TopenBc(true); //Use open boundary conditions
  if(GJP.TopenBc()) {
    GnoneFnone lat;
    lat.ZeroTboundary();
  }

  VRB.Result(cname, fname, "\n\n\n ************ Starting to create actions and integrators **********\n\n\n");

  AlgMomentum mom;
  AlgActionGauge gauge(mom, gauge_arg);
  AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge, ab1_arg);

  VRB.Result(cname, fname, "\n\n\n ************ Finished creating actions and integrators ***********\n\n\n");

  //main configuration-generating loop
  int traj = evo_arg.traj_start;
  const int start_traj = traj;
  while(1) {
    do_all_measurements(traj);
    
    run_hmc_trajectory(traj, ab1);
    traj++;
    
    if((traj % wflow_trajs_per_measurement) == 0) {
      checkpoint(traj);
      if(traj - start_traj >= trajs_per_run) break;
    }
  }

  AlgIntAB::Destroy(ab1);

  VRB.Result(cname, fname, "Program exiting normally\n");

  End();

  return 0;
}


void setup(int *argc, char **argv[])
{
  const char* fname = "setup()";

  Start(argc, argv);

  if(chdir("../vmls") != 0) {
    ERR.General(cname, fname, "Changing director to ../vmls failed!\n");
  }

  decode_vml(do_arg);
  decode_vml(evo_arg);
  decode_vml(hmc_arg);
  decode_vml(gauge_arg);
  decode_vml(ab1_arg);

  if(chdir(evo_arg.work_directory) != 0) {
    ERR.General(cname, fname, "Changing directory to %s failed.\n", evo_arg.work_directory);
  }
  
  VRB.Result(cname, fname, "Successfully read VML files.\n");

  GJP.Initialize(do_arg);
  LRG.Initialize();
  LRG.setSerial();

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
}

void erase_file(const char* filename)
{
  FILE* f = Fopen(filename, "w");
  Fclose(f);
}
    
void run_hmc_trajectory(int traj, AlgIntAB &top_level_integrator)
{
  const char* fname = "run_hmc_trajectory()";
  
  VRB.Result(cname, fname, "\n\n\n ************** Starting HMC trajectory ***************\n\n\n");
  Float timer = dclock();

  //Decide whether to do a reproducibility test on this trajectory
  if( (evo_arg.reproduce_interval > 0) && 
      ((traj % evo_arg.reproduce_interval) == 0)) {
    VRB.Result(cname, fname, "Running trajectory %d with reproduction\n", traj);
    hmc_arg.reproduce = REPRODUCE_YES;
  } else {
    VRB.Result(cname, fname, "Running trajectory %d without reproduction\n", traj);
    hmc_arg.reproduce = REPRODUCE_NO;
  }

  //for the first few trajectories, automatically turn off the metropolis step
  MetropolisType saved_metropolis = hmc_arg.metropolis;
  if(traj < 30) {
    VRB.Result(cname, fname, "Disabled Metropolis step for early trajectory (traj = %d)\n", traj);
    hmc_arg.metropolis = METROPOLIS_NO;
  }
  
  CommonArg common_arg;
  char filename[512];
  sprintf(filename, "%s.%d", hmc_stem, traj);
  erase_file(filename);
  common_arg.set_filename(filename);

  //run the HMC algorithm for one trajectory
  AlgHmc hmc(top_level_integrator, common_arg, hmc_arg);
  hmc.run();

  hmc_arg.metropolis = saved_metropolis;

  VRB.Result(cname, fname, "\n\n\n ************** Finished HMC trajectory in %e seconds ***************\n\n\n", dclock() - timer);
}


void do_all_measurements(int traj)
{
  measure_plaq(traj);

  if((traj % wflow_trajs_per_measurement) == 0) {
    do_wflow_measurements(traj, wflow_dt, wflow_nsteps, wflow_measure_start_step, wflow_measure_interval);
  }
}



void measure_plaq(int traj)
{
  const char* fname = "measure_plaq()";

  VRB.Result(cname, fname, "\n\n\n ************ Starting plaquette measurement *************\n\n\n");
  Float timer = dclock();

  GwilsonFnone lat;

  CommonArg common_arg;
  char filename[512];
  sprintf(filename, "%s.%d", plaq_stem, traj);
  erase_file(filename);
  common_arg.set_filename(filename);
  
  NoArg no_arg;

  AlgPlaq plaq(lat, &common_arg, &no_arg);
  plaq.run();

  VRB.Result(cname, fname, "\n\n\n ************ Finished plaquette measurement in %e seconds *************\n\n\n", dclock() - timer);
}
  

void do_wflow_measurements(int traj,
                           Float dt,
                           int num_time_steps,
                           int measure_start_step,
                           int measure_interval)
{
  const char* fname = "do_wflow_measurements()";
  VRB.Result(cname, fname, "\n\n\n ************ Starting Wilson flow measurement *************\n\n\n");
  Float timer = dclock();

  GwilsonFnone lat;

  // save the unsmeared lattice
  LatticeContainer lat_cont;
  lat_cont.Get(lat);
  
  CommonArg common_arg;
  char filename[512];
  sprintf(filename, "%s.%d", wflow_stem, traj);
  erase_file(filename);
  common_arg.set_filename(filename);

  AlgWilsonFlow wilson_flow(lat, &common_arg, dt);
  AlgActionDensity action_density(lat, &common_arg);
  AlgTcharge tcharge (lat, &common_arg);

  int step;
  for(step = 0; step < num_time_steps; step++) {
    if(step >= measure_start_step && ((step % measure_interval) == 0)) {
      VRB.Result(cname, fname, "Smart-Running AlgActionDensity %d...\n", step);
      action_density.run();
      tcharge.run();
    }

    wilson_flow.run(); //integrate the wilson flow for one time step
  }

  if(step >= measure_start_step && ((step % measure_interval) == 0)) {
    VRB.Result(cname, fname, "Smart-Running AlgActionDensity %d...\n", step);
    action_density.run();
    tcharge.run();
  }

  // restore the unsmeared lattice
  lat_cont.Set(lat);

  VRB.Result(cname, fname, "\n\n\n ************ Finished Wilson flow measurement in %e seconds *************\n\n\n", dclock() - timer);
}



void checkpoint(int traj)
{
  const char *fname="checkpoint()";
  
  VRB.Result(cname, fname, "\n\n\n ***************** Starting checkpoint **************** \n\n\n");

  char lat_file[256];
  char rng_file[256];

  Float time = -dclock();

  // Save this config to disk
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

  sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
  QioArg wt_arg(lat_file,0.001);

  wt_arg.ConcurIONumber=evo_arg.io_concurrency;
  WriteLatticeParallel wl;
  wl.setSerial();
  wl.setHeader(evo_arg.ensemble_id,evo_arg.ensemble_label,traj);
  wl.write(lat,wt_arg);

  if(!wl.good())
    ERR.General(cname,fname,"Failed write lattice %s",lat_file);

  LatticeFactory::Destroy();

  // Save the RNG's
  sprintf(rng_file,"%s.%d",evo_arg.rng_file_stem,traj);
  LRG.setSerial();
  if ( !LRG.Write(rng_file) )
    ERR.General(cname,fname,"Failed write RNG file %s",rng_file);

  // Update the parameter files for restart
  do_arg.start_seed_filename = rng_file;
  do_arg.start_seed_kind = START_SEED_FILE;
  do_arg.start_conf_filename = lat_file;
  do_arg.start_conf_kind = START_CONF_FILE;
  evo_arg.traj_start = traj;

  encode_vml(do_arg, traj);
  encode_vml(evo_arg, traj);

  time += dclock();
  print_flops("","checkpoint()",0,time);
    
  VRB.Result(cname, fname, "\n\n\n ***************** Finished checkpoint **************** \n\n\n");
}


