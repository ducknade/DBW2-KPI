#include <fstream>

#include <sys/stat.h>
#include <unistd.h>
#include <cstdio>

#include <qmp.h>
#include <config.h>
#include <util/lattice.h>
#include <util/lat_cont.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/random.h>
#include <alg/alg_pbp.h>
#include <alg/do_arg.h>
#include <alg/common_arg.h>
#include <alg/pbp_arg.h>
#include <alg/alg_actiondensity.h>
#include <alg/alg_wilsonflow.h>
#include <alg/alg_tcharge.h>
#include <alg/alg_smear.h>

#include <alg/alg_meas.h>
#include <util/ReadLatticePar.h>
#include <util/qioarg.h>

USING_NAMESPACE_CPS

static const char* cname = "";

bool mkdir_recursively(const char *dir){
        char tmp[256];
        char *p = NULL;
        size_t len;

        snprintf(tmp, sizeof(tmp),"%s",dir);
        len = strlen(tmp);
        if(tmp[len - 1] == '/') tmp[len - 1] = 0;
        for(p = tmp + 1; *p; p++){
            if(*p == '/'){
                *p = 0;
                mkdir(tmp, 0777);
                *p = '/';
            }

        }
        return(!mkdir(tmp, 0777));
}

void set_do_arg(DoArg& do_arg, const int total_size[4])
{
	do_arg.x_sites = total_size[0];
	do_arg.y_sites = total_size[1];
	do_arg.z_sites = total_size[2];
	do_arg.t_sites = total_size[3];
	do_arg.s_sites = 2;
	do_arg.dwf_height = 1.0;
	do_arg.x_bc = BND_CND_PRD;
	do_arg.y_bc = BND_CND_PRD;
	do_arg.z_bc = BND_CND_PRD;
	do_arg.t_bc = BND_CND_PRD;
	do_arg.start_conf_kind = START_CONF_ORD;
	do_arg.start_seed_kind = START_SEED_INPUT;
	do_arg.start_seed_value = 123121;
	do_arg.x_nodes = 0;
	do_arg.y_nodes = 0;
	do_arg.z_nodes = 0;
	do_arg.t_nodes = 0;
	do_arg.s_nodes = 0;
	do_arg.x_node_sites = 0;
	do_arg.y_node_sites = 0;
	do_arg.z_node_sites = 0;
	do_arg.t_node_sites = 0;
	do_arg.s_node_sites = 0;
	do_arg.gfix_chkb = 1;
}

void load_config(char* filename)
{
    const char *cname = "cps";
    const char *fname = "load_config";

	VRB.Result(cname, fname, "-------- starting to load configuration: %s\n", filename);	
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    QioArg rd_arg(filename, 0.001);
    rd_arg.ConcurIONumber = 1;
    ReadLatticeParallel rl;
    rl.read(lat, rd_arg);
    if(!rl.good()) ERR.General(cname, fname, "failed to read configuration %s\n", filename);
    LatticeFactory::Destroy();
	
	VRB.Result(cname, fname, "-------- finished to load configuration: %s\n", filename);	
}

int main(int argc, char* argv[]){
	
	const char* fname = "main()";

// -------- INPUT DECK!!!

// input lattice configuration directory
	const char config_dir_stem[] = "/home/gregm/DBW2/open-10x20/configurations/ckpoint_lat";
// 4d lattice size
	const int total_size[] = {10, 10, 10, 20}; // {x, y, z, t}
// start configuration number
	const int config_num_start = 10008;
// total number of configurations
	const int config_num_count = 9619;
// interval between two consecutive configurations
	const int config_num_incre = 12;
// dt for wilson flow 
	const double wflow_dt =      0.05;
// number of wilson flow steps for each loaded configuration
	const int wflow_steps =      30;
// ensemble name
	const char* ens_name =       "open-10x20";

// -------- END INPUT DECK!!!

// setup
	Start(&argc, &argv);
	
	DoArg do_arg;
	set_do_arg(do_arg, total_size);	
	GJP.Initialize(do_arg);
 	LRG.Initialize();

	mkdir("../results", 0777);
	char ens_dir[512]; sprintf(ens_dir, "../results/%s", ens_name);
	mkdir(ens_dir, 0777);
	char wflow_dir_stem[512]; sprintf(wflow_dir_stem, "../results/%s/wflow", ens_name);
	char config_dir[512];

	int config_num;
	for(int i = 0; i < config_num_count; i++){
		config_num = i * config_num_incre + config_num_start;
		sprintf(config_dir, "%s.%d", config_dir_stem, config_num);
		load_config(config_dir);
		
		Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

		char wflow_dir[512]; sprintf(wflow_dir, "%s.%d", wflow_dir_stem, config_num);
		CommonArg common_arg("wilson_flow", wflow_dir);

		AlgWilsonFlow wflow(lat, &common_arg, wflow_dt);
		AlgActionDensity density(lat, &common_arg);
		AlgTcharge tcharge(lat, &common_arg);

		for(int j = 0; j < wflow_steps; j++){
			density.run();
			tcharge.run();
			wflow.run();
		}

		density.run();
		tcharge.run();

		LatticeFactory::Destroy(); 
	}

	End();
	return 0;
}

