#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>
#include <string>
#include <OpenMesh/Apps/Assembly/energy.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;


struct RunParameters{
	//number of timesteps to run
	long timesteps;

	//frequency of dumping configurations
	long dtsave;

	//don't dump configurations at the beginning, it can be slow
	long dtskip;

	//full path to where to save data
	std::string data_folder;

	//ensemble ID used to seed the random number generator
	long ensemble;

	//spring constant for umbrella sampling
	double spring_const;

	//number of steps in an umbrella window
	long nsteps_umbrella;

	//shift width of umbrella window; windows are centered at 0, dN_umbrella, 2*dN_umbrella, ...
	int dN_umbrella;

	// .om file for initial configuration
	std::string init_file;

	// byte position in the .om file
	long long init_file_pos;

	//lower k_fusion, k_fusion2, k_fusion_edge, k_insertion if p_propose>1
	bool adaptive_rates;

	//convet to lammps files at the end?
	bool convert_to_lammps;

	//convert to vtk files at the end?
	bool convert_to_vtk;

	//convert to lammps traj files at the end?
	bool convert_to_lammps_trajectory;

	std::string path_to_source;

	//convert to Carlos' vertices.dat, faces.dat, edges.dat at the end?
	bool convert_to_dat;

	//burn-in time for umbrella windows; statistics is not taken during burnin time
	long dt_burnin;

	//stop the simulation if the structure reached this size
	int Nmax;

	//checkpoint creation frequency
	long checkpoint_freq;

	//continue simulation if possible (if checkpoint file is present)
	bool continue_if_possible;

	//the maximum amount of physical time to run; set it to less than 2 days for Stampede
	//so that everything gets saved and copied properly before the scheduler cancels 
	long max_seconds;

	//wall potential is amplitude*exp(-hardness*dz)
	double wall_hardness;

	//wall potential is amplitude*exp(-hardness*dz)
	double wall_amplitude;

	//dz / timesteps
	double wall_speed;

	//number of timesteps in a sawtooth period
	long wall_period;

	//frequency of dumping snapshots; if not set, the value of dt_save will be used
	//the point here is to be able to log more frequently and save configurations less frequently to save space
	//it is assume that dtsave_snapshot is a multiple of dtsave
	long dtsave_snapshot;

	//whether to stop the simulation when the structure is fully closed
	bool stop_at_closure;	


	std::string externalpotential;
	double extparam1;
	double extparam2;
	double extparam3;
	double extparam4;
	double extparam5;	
	double extparam6;


	//einstein spring constant for thermodynamic integration
	double k_einstein;

	//find optimal spring constant from rms of displacement fluctuations
	//bool find_optimal_k_einstein;

};

std::string getenv_str(std::string env_name);

void read_set_mesh_params(MyMesh & mesh, std::string input_filename);

void read_set_run_params(RunParameters & rp, std::string input_filename);

void init_mesh_from_om(MyMesh & mesh, std::string om_file, long long om_byte_position, std::string input_filename);

long long dump_om(MyMesh mesh, long t, std::ostream& os);

int dump_om(MyMesh mesh, std::string om_filename);

int dump_data(MyMesh mesh, long timestep, std::ostream& os, long long key_to_snapshot=0);

int dump_data(MyMesh mesh, long timestep, std::ostream& os, long long key_to_snapshot, std::vector<double> r_acc);

int dump_lammps_snapshot(MyMesh mesh, std::string outfile_name);

int convert_om_to_lammps_snapshots(std::string om_filename, std::string lammps_base);

//void dump_mesh(MyMesh mesh, std::string om_filename);

int convert_om_to_lammps_trajectory(std::string om_filename, std::string lammps_base="");

int convert_om_to_VTK_snapshots(std::string om_filename, std::string vtk_base, std::string vtk_conf_file);

int dump_vtk_snapshot(MyMesh mesh, std::string outfile_name, std::string vtk_conf_file);

void print_mesh(MyMesh & mesh);

//void print_mesh_stripped(MyMesh & mesh);

void create_single_subunit_om(std::string om_filename, int face_type, int edge_type1, int edge_type2, int edge_type3);

void create_hexamer_subunit_om(std::string om_filename, int face_type, int edge_type1, int edge_type2, int edge_type3);

void convert_dat_to_om(std::string vertices_dat, std::string faces_dat, std::string edges_dat, std::string omfile);

void dump_dat_snapshot(MyMesh mesh, std::string vertices_dat, std::string faces_dat, std::string edges_dat);

int convert_om_to_dat_snapshots(std::string om_filename, std::string dat_base);

int convert_om_to_dat_snapshots(std::string om_filename, std::string dat_base, long long om_byte_position);

int create_checkpoint(MyMesh & mesh, UmbrellaWindow & uw, long timestep, std::string ckp_filename);

int load_from_checkpoint(MyMesh & mesh, UmbrellaWindow & uw, long & timestep, std::string ckp_filename, std::string input_filename);

int create_checkpoint_wall(MyMesh & mesh, Wall & wall, long timestep, std::string ckp_filename);

int load_from_checkpoint_wall(MyMesh & mesh, Wall & wall, long & timestep, std::string ckp_filename, std::string input_filename);

bool file_exists (std::string name);

MyMesh get_stripped_mesh(MyMesh mesh);

int dump_lammps_snapshot_double_edge(MyMesh mesh, std::string outfile_name);

int convert_om_to_lammps_snapshots_double_edge(std::string om_filename, std::string lammps_base);

void dump_error(std::string error_msg, std::string filename);

long long append_snapshot(std::string om_filename, long long om_byte_position, std::string omfile_to_append_to);

int dump_lammps_typed_edges_snapshot(MyMesh mesh, std::string outfile_name);