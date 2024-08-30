#include <iostream>
#include <stdlib.h>
#include <getopt.h>
#include <ctime>


#include <OpenMesh/Apps/Assembly/IO.hh>
#include <OpenMesh/Apps/Assembly/run.hh>
#include <OpenMesh/Apps/Assembly/analyze.hh>


// ./assembly -d -i input_carlos2.json --> dynamical run
// ./assembly -u -i input_carlos2.json  --> umbrella sampling
// ./assembly -p -i pt_input.json  --> parallel tempering
// ./assembly -w -i input_carlos2.json  --> wham only on existing data
// ./assembly -W -i PT_input.json  --> wham only on existing PT data
// ./assembly -c snapshots.om -f conversions/ --> convert .om 
// ./assembly -c snapshots.om  -b 803215 -f conversions/ --> convert .om, only the specified byte position 
// ./assembly -g datfiles/ -o structure.om --> generate .om from datfiles/vertices.dat, faces.dat, edges.dat
// ./assembly -a -i input_carlos2.json --> analyze; snapshots.om already written in the .json file
// ./assembly -t -i snapshots.om --> analyze topology only based on what's in snapshots.om alone; only simplified mesh can be accessed
// ./assembly -q -i input_carlos2.json --> quench structure to kT=0
// ./assembly -S -i input_carlos2.json --> dynamical run with wall squishing
// ./assembly -T -i input_carlos2.json --> thermodynamic integration
// ./assembly -i omfile.om -b 803215 -A omfile_to_append_to.om --> append snapshot from position b in omfile.om to the file in omfile_to_append_to.om
int main(int argc, char **argv){
	extern char *optarg;
	extern int optind;
	int c, dynamic=0, umbrella=0, convert=0, analyze=0, wham=0, pt=0, debug=0, generate=0, WHAM=0, topology_analyze=0, quench=0, wall_squish=0, extpotential=0, thermodynamic_integration=0, append_snapshot_flag=0;
	int const_load_wall=0;
	int err;
	long long b=-1;
	std::string input_fname, om_snapshots_fname, conversion_foldername, datfiles, structurename;
	while ((c = getopt(argc, argv, "duai:c:s:f:wpDg:o:WtqSelb:TA:")) != -1)
		switch (c) {
		case 'd':
			dynamic = 1;
			break;
		case 'u':
			umbrella = 1;
			break;
		case 'w':
			wham = 1;
			break;	
		case 'W':
			WHAM=1;
			break;			
		case 'a':
			analyze = 1;
			break;
		case 'i':
			input_fname = optarg;
			break;
		case 'c':
			convert = 1;
			om_snapshots_fname = optarg;
		
			break;

		case 'f':
			conversion_foldername = optarg;	
			break;
		case 's':
			om_snapshots_fname = optarg;
			break;
		case 'p':
			pt=1;
			break;
		case 'D':
			debug=1;
			break;
		case 'g':
			generate=1;
			datfiles = optarg;
			break;
		case 'o':
			structurename = optarg;
			break;
		case 't':
			topology_analyze=1;
			break;
		case 'q':
			quench=1;
			break;
		case 'S':
			wall_squish=1;
			break;
		case 'e':
			extpotential=1;
			break;
		case 'l':
			const_load_wall=1;
			break;
		case 'T':
			thermodynamic_integration=1;
			break;
		case 'b':
			b=atoll(optarg);
			break;
		case 'A':
			append_snapshot_flag=1;
			om_snapshots_fname = optarg;
			break;
		case '?':
			err = 1;
			break;
		}

	//std::time_t now = std::time(nullptr);
	//char * t_start = std::asctime(std::localtime(&now));
	std::time_t t_start = std::time(nullptr);
	std::time_t t_end;

	if (dynamic){
		std::cout<<"Starting dynamical run at "<<std::asctime(std::localtime(&t_start))<<std::endl;
		run_dynamical(input_fname);
		t_end = std::time(nullptr);		
		std::cout<<"Dynamical run done at "<<std::asctime(std::localtime(&t_end))<<std::endl;		
	}


	if (extpotential){
		std::cout<<"Starting external potential dynamical run at "<<std::asctime(std::localtime(&t_start))<<std::endl;
		run_dynamical_external_potential(input_fname);
		t_end = std::time(nullptr);		
		std::cout<<"External potential run done at "<<std::asctime(std::localtime(&t_end))<<std::endl;		
	}

	if (const_load_wall){
		std::cout<<"Starting wall potential const load dynamical run at "<<std::asctime(std::localtime(&t_start))<<std::endl;
		run_wall_potential_constant_load(input_fname);
		t_end = std::time(nullptr);		
		std::cout<<"Wall potential const load run done at "<<std::asctime(std::localtime(&t_end))<<std::endl;		
	}	

	if (umbrella){
		std::cout<<"Starting umbrella run at "<<std::asctime(std::localtime(&t_start))<<std::endl;
		run_umbrella(input_fname);
		t_end= std::time(nullptr);		
		std::cout<<"Umbrella run done at "<<std::asctime(std::localtime(&t_end))<<std::endl;		
	}	

	if (pt){
		std::cout<<"Starting parallel tempering umbrella run at "<<std::asctime(std::localtime(&t_start))<<std::endl;
		run_umbrella_parallel_tempering(input_fname);
		t_end = std::time(nullptr);		
		std::cout<<"Parallel tempering run started at "<<std::asctime(std::localtime(&t_start))<<std::endl;	
		std::cout<<"Parallel tempering run done at "<<std::asctime(std::localtime(&t_end))<<std::endl;
					
	}

	if (wham){
		system(( "python wham.py "+input_fname ).c_str());
	}

	if (convert){
		std::cout<<"Starting conversion at "<<std::asctime(std::localtime(&t_start))<<std::endl;
		std::cout<<om_snapshots_fname<<std::endl;
		std::cout<<conversion_foldername<<std::endl;
		if (b==-1){

			//@@@@ convert_om_to_lammps_trajectory(om_snapshots_fname, conversion_foldername);
			convert_om_to_lammps_snapshots(om_snapshots_fname, conversion_foldername+"lammps_");		  
			 convert_om_to_dat_snapshots(om_snapshots_fname, conversion_foldername+"dat_");
			 convert_om_to_VTK_snapshots(om_snapshots_fname, conversion_foldername+"vtk_", "vtk_colors.conf"); 	

			convert_om_to_lammps_snapshots_double_edge(om_snapshots_fname, conversion_foldername+"double_edge_lammps_");			 	


		}	
		else{
			if (b>=0){
				convert_om_to_dat_snapshots(om_snapshots_fname, conversion_foldername+"dat_", b);
			}
			if (b==-2){
				//convert_om_to_dat_snapshots(om_snapshots_fname, conversion_foldername+"dat_");
				system(("mkdir -p "+conversion_foldername+"/undeformed_frame/").c_str()); 
				convert_om_to_dat_snapshots_undeformed_frame(om_snapshots_fname, conversion_foldername+"/undeformed_frame/dat_");
			}
		}
		t_end = std::time(nullptr);		
		std::cout<<"Conversion run done at "<<std::asctime(std::localtime(&t_end))<<std::endl;	
	}	

	if (debug){
		run_debug();
	}

	if (generate){
		convert_dat_to_om(datfiles+"vertices.dat", datfiles+"faces.dat", datfiles+"edges.dat", structurename);
	}

	if (WHAM){
		system(( "python multi_wham.py "+input_fname ).c_str());
	}

	if (analyze){
		analyze_all(input_fname);
	}

	if (topology_analyze){
		topology_analyze_all(input_fname);
	}

	if (quench){
		std::cout<<"Starting quench run at "<<std::asctime(std::localtime(&t_start))<<std::endl;
		run_quench(input_fname);
		t_end = std::time(nullptr);		
		std::cout<<"Quench run done at "<<std::asctime(std::localtime(&t_end))<<std::endl;		
	}	

	if (wall_squish){
		std::cout<<"Starting wall squish run at "<<std::asctime(std::localtime(&t_start))<<std::endl;
		run_wall(input_fname);
		t_end = std::time(nullptr);		
		std::cout<<"|Wall squish run done at "<<std::asctime(std::localtime(&t_end))<<std::endl;		
	}
	if (thermodynamic_integration){
		std::cout<<"TI run at "<<std::asctime(std::localtime(&t_start))<<std::endl;
		run_thermodynamic_integration(input_fname);
		t_end = std::time(nullptr);		
		std::cout<<"TI run done at "<<std::asctime(std::localtime(&t_end))<<std::endl;		
	}		

	if (append_snapshot_flag){
		std::cout<<"Appending snapshot at pos."<<b<<" from "<<input_fname<<" to "<<om_snapshots_fname<<std::endl;
		long long newpos=append_snapshot(input_fname, b, om_snapshots_fname);
		std::cout<<"Appended to position "<<newpos<<std::endl;
	}

	return 0;
}

