#include <iostream>
#include <stdlib.h>
#include <math.h> 
#include<limits>
#include <fstream>
#include <unistd.h>
#include <getopt.h>
#include <ctime>
#include <cmath>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include<omp.h>
// -------------------- OpenMesh
#include <OpenMesh/Core/System/omstream.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Apps/Assembly/topology.hh>
#include <OpenMesh/Apps/Assembly/geometry.hh>
#include <OpenMesh/Apps/Assembly/energy.hh>
#include <OpenMesh/Apps/Assembly/monte_carlo.hh>
#include <OpenMesh/Apps/Assembly/monte_carlo_wall.hh>
#include <OpenMesh/Apps/Assembly/monte_carlo_external_potential.hh>
#include <OpenMesh/Apps/Assembly/thermodynamic_integration_move.hh>
#include <OpenMesh/Apps/Assembly/excluders.hh>
#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>
#include <OpenMesh/Apps/Assembly/IO.hh>
#include <OpenMesh/Apps/Assembly/run.hh>
#include <OpenMesh/Apps/Assembly/random.hh>
#include <OpenMesh/Apps/Assembly/analyze.hh>
#include <OpenMesh/Apps/Assembly/subset_division.hh>

#include <OpenMesh/Apps/Assembly/json.hh>

// for convenience
using json = nlohmann::json;



std::vector<int> sort_indexes(const std::vector<int> &v) {

  // initialize original index locations
  std::vector<int> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v](int i1, int i2) {return v[i1] > v[i2];});

  return idx;
}


void init_mesh(MyMesh & mesh, std::string init_filename, std::string input_filename, long long byte_position){
  	mesh.request_face_status();
  	mesh.request_edge_status();
  	mesh.request_vertex_status();    
  	mesh.request_halfedge_status();
	mesh.request_face_normals();

	//std::string om_filename = "/home/btyukodi/assembly_openmesh/OpenMesh-9.0/build/outfileom4.om";
	//std::string input_filename = "/home/btyukodi/assembly_openmesh/sandbox/input_carlos2.json";
	//long long byte_position = 145;
	init_mesh_from_om(mesh, init_filename, byte_position, input_filename);


}


//will need to add umbrella window
struct AcceptanceMonitor{

    std::vector<double> n_acc;
    std::vector<double> n_tot;
    std::vector<double> r_acc;	
    int nmoves;
    AcceptanceMonitor(int nmoves=13){
    	for (int i=0; i<nmoves; i++){
    		this->n_acc.push_back(0);
    		this->n_tot.push_back(0);
    		this->r_acc.push_back(0);
    	}

    }
};

static AcceptanceMonitor default_AM = AcceptanceMonitor();

void propagate(MyMesh & mesh, long t_init, long t_final, RunParameters rp, std::ostream& outfile_om,  std::ostream& outfile_data, UmbrellaWindow & uw, std::mt19937 & eng, AcceptanceMonitor & AM = default_AM){

    int nmoves = 13;//10;

	/*for (int iac=0; iac<(nmoves+1); iac++){
		n_acc[iac]=0;
		n_tot[iac]=0;
	}    */

	bool accepted;
	int which_move, which_type, which_rotation, number_of_rotational_configs;
	double k_insertion, k_fusion, p_propose, r;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> , std::vector<int>> wedges;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > fusion_vectors, type2_fission_vectors, holes;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> >	fission_vectors;	
	std::vector <MyMesh::VertexHandle> w1, v1;  
	std::vector<MyMesh::HalfedgeHandle> boundary_halfedges, h1;	  
	std::vector<MyMesh::FaceHandle> removable_faces_geom, removable_faces;	
	std::tuple<std::vector<MyMesh::HalfedgeHandle>, std::vector<MyMesh::HalfedgeHandle> > fusion_halfedges;	
	std::vector<MyMesh::EdgeHandle> fission_edges;	

	std::vector<MyMesh::VertexHandle*> handle_tracking_v;
		std::vector<MyMesh::FaceHandle*> handle_tracking_f;
	std::vector<MyMesh::HalfedgeHandle*> handle_tracking_h;	

	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto proto_face_props = OpenMesh::FProp<FaceProp>((*mesh_props).prototypes_mesh, "face_props");		
	float kT = (*mesh_props).kT;
	double l_fuse = (*mesh_props).l_fuse;	

	//std::ofstream outfile("outfileom5.om", std::ios::out | std::ios::binary);
	// std::ofstream datafile("data_out5.dat");
	long long fpos;
	long long opos;
	int n_boundary_edges;
	int n_type1, n_type2;

//dump_data(mesh, -1, outfile_data, fpos);
	for (long t=t_init; t<t_final; t++){
		if (mesh.n_faces() > rp.Nmax){
			break;
		}



		if (t % rp.dtsave ==0){		
			std::cout<<"t = "<<t<<" n_faces="<<mesh.n_faces()<<std::endl;	
			n_boundary_edges = count_boundary_edges(mesh);
			n_type1 = count_subunits_of_type(mesh, 0);
			n_type2=count_subunits_of_type(mesh, 1);
			
			//skip full dump for the first dtskip timesteps
			
			if ( (t>rp.dtskip) && (t % rp.dtsave_snapshot==0)){
				fpos = dump_om(mesh, t, outfile_om);
			}
			else{
				//save the last snapshot if it will break
				if (rp.stop_at_closure && (n_boundary_edges==0)){
						fpos = dump_om(mesh, t, outfile_om);	
				}
				else{				
					fpos=-1;
				}
			}



			for (int iac=0; iac<(nmoves); iac++){
				AM.r_acc[iac] = AM.n_acc[iac]/AM.n_tot[iac];
				//std::cout<<"r_acc "<< AM.r_acc[iac]<<std::endl;
				//std::cout<<"n_acc "<< AM.n_acc[iac]<<std::endl;
				//std::cout<<"n_tot "<< AM.n_tot[iac]<<std::endl;								
				AM.n_acc[iac]=0;
				AM.n_tot[iac]=0;
			}

			dump_data(mesh, t, outfile_data, fpos, AM.r_acc);
			opos = outfile_data.tellp();
  			outfile_data.seekp (opos-1);			
			outfile_data<<"\t"<< n_type1<<"\t"<<n_type2<<"\t" << n_boundary_edges <<std::endl;
			//dump_om(mesh, "omtest"+std::to_string(t)+".om");


			//don't count boundary edges at every timestep
			if (rp.stop_at_closure && (n_boundary_edges==0)){
					break;				
			}				
		}

		if (t % rp.checkpoint_freq ==0){
			create_checkpoint(mesh, uw, t, rp.data_folder+"checkpoint.ckp");			
		}

		
		/*
		for (int n=0; n<(int)mesh.n_vertices(); n++){
			AM.n_acc[0]+=attempt_move(mesh, eng);
			AM.n_tot[0]++;
		}
		*/
			/*std::cout<<"n_tot0 "<< n_tot[1]<<std::endl;
			std::cout<<"n_acc0 "<< n_acc[1]<<std::endl;			
			std::cout<<"n_acc0/n_tot0 "<< n_acc[1]/n_tot[1]<<std::endl;*/		


		//!! IMPORTANT: in each t iteration, only one move should be selected because the garbage collector
		// is only called once, at the end of each iteration. n_faces, n_vertices, etc won't get updated unless
		// the garbage collector is called and we need them to be up to date for the bias potential and also for attempt_move()
		//which_move = rand() % nmoves;
		which_move = randint(nmoves, eng);		
		//std::cout<<"which move "<<which_move<<std::endl;
		accepted=false;
		switch (which_move){
			case 0:
				k_fusion = (*mesh_props).k_fusion; 
				fusion_vectors = get_type1_fusion_triplets(mesh, l_fuse);	
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size();
				if ( (p_propose>1) && (!rp.adaptive_rates) ) std::cout<<"Warning! Decrease k_fusion rate! - type1_fusion"<<std::endl;
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion rate... - type1_fusion"<<std::endl;
					(*mesh_props).k_fusion*=0.5;
					break;//-------- don't do the move if violates detailed balance
				}				
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);
				if (r<p_propose){	
					accepted = attempt_type1_fusion(mesh, eng);
					AM.n_acc[which_move+1]+=accepted;					
				}
				//if (accepted) std::cout<<"type1 fusion"<<std::endl;
				AM.n_tot[which_move+1]++;
				break;

			case 1:
				k_fusion = (*mesh_props).k_fusion; 
				fission_vectors = get_type1_fission_pairs(mesh);	
    			v1 = std::get<0>(fission_vectors);					
				p_propose = k_fusion * v1.size();
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion rate! - type1_fission"<<std::endl;		
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion rate... - type1_fission"<<std::endl;
					(*mesh_props).k_fusion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}							
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	
					accepted = attempt_type1_fission(mesh, eng);	
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type1 fission"<<std::endl;
				AM.n_tot[which_move+1]++;
				break;

			case 2:		
				k_fusion = (*mesh_props).k_fusion2; 
				fusion_vectors = get_type2_fusion_triplets(mesh, l_fuse);	//each counted TWICE
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size()*0.5;  // x 1/2 to correct for double counting
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fusion"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion2 rate... - type2_fusion"<<std::endl;
					(*mesh_props).k_fusion2*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}								
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_type2_fusion(mesh, eng);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type2 fusion "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;	
				break;
			case 3:			
				k_fusion = (*mesh_props).k_fusion2; 
				type2_fission_vectors = get_type2_fission_triplets(mesh);	//each counted TWICE
    			v1 = std::get<0>(type2_fission_vectors);					
				p_propose = k_fusion * v1.size()*0.5; // x 1/2 to correct for double counting
				if ( (p_propose>1) && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fission"<<std::endl;			
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion2 rate... - type2_fission"<<std::endl;
					(*mesh_props).k_fusion2*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}						
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_type2_fission(mesh, eng);	
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type2 fission "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;
				break;

			case 4:	
				//if (mesh.n_faces()>n_umb_sim ) break;

				//pick subunit type from prototype subunits
				//which_type = rand() % (*mesh_props).prototype_faces.size();
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				number_of_rotational_configs = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].number_of_rotational_configs;
				//which_rotation = rand() % number_of_rotational_configs;
				which_rotation = randint( number_of_rotational_configs, eng );				
				//----
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				boundary_halfedges = find_boundary_halfedges(mesh);			
				p_propose = k_insertion * boundary_halfedges.size() * number_of_rotational_configs;
				if ( (p_propose>1) && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_insertion rate! - insertion"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - insertion"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}								
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	
					accepted = attempt_insertion(mesh, which_type, which_rotation, uw, eng);	
					AM.n_acc[which_move+1]+=accepted;
				}  
				//if (accepted) std::cout<<"insertion type "<< which_type <<" "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;				
				break;

			case 5:	
				//if (mesh.n_faces()>n_umb_sim ) break;			
				//pick subunit type from prototype subunits
				//which_type = rand() % (*mesh_props).prototype_faces.size();
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				number_of_rotational_configs = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].number_of_rotational_configs;
				//which_rotation = rand() % number_of_rotational_configs;
				which_rotation = randint( number_of_rotational_configs, eng );				
				//---
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				wedges = get_open_wedge_triplets(mesh);		
				w1 = std::get<0>(wedges);		
				p_propose = k_insertion * w1.size() * number_of_rotational_configs;
				if ((p_propose>1) && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_insertion rate! - wedge_insertion"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - wedge_insertion"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_wedge_insertion(mesh, which_type, which_rotation, uw, eng);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"wedge insertion type "<< which_type <<" "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;			
				break;
			case 6:
				//pick subunit type from prototype subunits
				//which_type = rand() % (*mesh_props).prototype_faces.size();
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				removable_faces = get_simply_removable_faces(mesh, which_type);

				p_propose = k_insertion * removable_faces.size();
				if ( (p_propose>1) && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_insertion rate! - removal"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - removal"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	

					accepted = attempt_removal(mesh, which_type, uw, eng);
					AM.n_acc[which_move+1]+=accepted;
				}
				AM.n_tot[which_move+1]++;
				break;
			case 7:
				//pick subunit type from prototype subunits
				//which_type = rand() % (*mesh_props).prototype_faces.size();
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				removable_faces = get_wedge_removable_faces(mesh, which_type); 

				p_propose = k_insertion * removable_faces.size();
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_insertion rate! - wedge_removal"<<std::endl;		
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - wedge_removal"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	
					accepted = attempt_wedge_removal(mesh, which_type, uw, eng);
					AM.n_acc[which_move+1]+=accepted;
				}
				AM.n_tot[which_move+1]++;
				break;

			case 8:
				k_fusion = (*mesh_props).k_fusion_edge;
				fusion_halfedges = get_halfedge_fusion_pairs(mesh, l_fuse); //each counted TWICE
    			h1 = std::get<0>(fusion_halfedges);
    			p_propose = k_fusion * h1.size()*0.5; // x 1/2 to correct for double counting
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fusion"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion_edge rate... - edge_fusion"<<std::endl;
					(*mesh_props).k_fusion_edge*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}	
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){					
					accepted = attempt_edge_fusion(mesh, eng);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"EDGE fusion " <<" "<<p_propose<<std::endl;
				AM.n_tot[which_move+1]++;
				break;

			case 9:
				k_fusion = (*mesh_props).k_fusion_edge;
				fission_edges = get_fission_edges(mesh);
    			p_propose = k_fusion * fission_edges.size();	
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fission"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion_edge rate... - edge_fission"<<std::endl;
					(*mesh_props).k_fusion_edge*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}					
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){    									
					accepted = attempt_edge_fission(mesh, eng);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"EDGE fission " <<" "<<p_propose<<std::endl;
				AM.n_tot[which_move+1]++;
				break;

			case 10:
				for (int n=0; n<(int)mesh.n_vertices(); n++){
					AM.n_acc[0]+=attempt_move(mesh, eng);
					AM.n_tot[0]++;
				}	
				break;

			case 11:
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				number_of_rotational_configs = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].number_of_rotational_configs;
				//which_rotation = rand() % number_of_rotational_configs;
				which_rotation = randint( number_of_rotational_configs, eng );				
				//---
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				holes = get_hole_triplets(mesh);		
				w1 = std::get<0>(holes);		
				p_propose = k_insertion * w1.size() * number_of_rotational_configs * 1.0/3.0; //<--- 1/3 because each hole is counted 3x
				if ((p_propose>1) && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_insertion rate! - wedge_insertion"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - hole_insertion"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_hole_insertion(mesh, which_type, which_rotation, uw, eng);
					AM.n_acc[which_move+1]+=accepted;
				}
				if (accepted) std::cout<<"hole insertion type "<< which_type <<" "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;			
				break;	

			case 12:
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				removable_faces = get_hole_removable_faces(mesh, which_type); 

				p_propose = k_insertion * removable_faces.size();
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_insertion rate! - hole_removal"<<std::endl;		
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - hole_removal"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	
					accepted = attempt_hole_removal(mesh, which_type, uw, eng);
					AM.n_acc[which_move+1]+=accepted;
				}
				AM.n_tot[which_move+1]++;
				break;				
		}

		//get all face's neighbor list handles to save them through the garbage collector
    	handle_tracking_f.clear();
    	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
    		for (auto & neighbor: face_props[*fit].neighbor_list){
    			handle_tracking_f.push_back( &neighbor);
    		}
    	}
		mesh.garbage_collection<std::vector<MyMesh::VertexHandle*>, std::vector<MyMesh::HalfedgeHandle*>, std::vector<MyMesh::FaceHandle*> >(handle_tracking_v, handle_tracking_h, handle_tracking_f);

		if (t % 150 ==0){
			update_full_neighbor_list(mesh);
			//Could check overlaps for the full capsid here to verify
		}

		uw.add_value(mesh.n_faces());
	}


	return;
}

void run_dynamical(std::string input_file){

	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);   

    //initialize RNG from ensemble seed
    //std::srand(rp.ensemble*73 + 17);
    std::mt19937 eng(rp.ensemble*73 + 17);

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);

    //should I create data_folder here?
    system(("mkdir -p "+rp.data_folder).c_str());    

    std::ofstream outfile_om(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);
    std::ofstream outfile_data(rp.data_folder+"data_log.dat");	
    outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<std::endl;


    //burnin time to equilibrate initial structure
    for (long t=0; t<rp.dt_burnin; t++){
		for (int n=0; n<(int)mesh.n_vertices(); n++){
			attempt_move(mesh, eng);
		}	
	}    

    long t_init=0;
    long t_final=rp.timesteps;
    UmbrellaWindow uw_fake(0.0, 0.0);
    propagate(mesh, t_init, t_final, rp, outfile_om,  outfile_data, uw_fake, eng);
    outfile_om.close();
    outfile_data.close();


    //do conversions if requested

    //create a conversions directory
    system(("mkdir -p "+rp.data_folder+"conversions").c_str());
    if (rp.convert_to_lammps_trajectory)
		convert_om_to_lammps_trajectory(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/");
    if (rp.convert_to_lammps)
		convert_om_to_lammps_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_");
	if (rp.convert_to_vtk)
		convert_om_to_VTK_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/vtk_", rp.path_to_source+"vtk_colors.conf");    
	if (rp.convert_to_dat)
		convert_om_to_dat_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/dat_");
	return;
}


void run_umbrella(std::string input_file){

	//create_single_subunit_om("single_subunit_trumpet.om", 0, 0, 1, 2);

	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);   

    //initialize RNG from ensemble seed
    //std::srand(rp.ensemble*73 + 17);
    std::mt19937 eng(rp.ensemble*73 + 17);    

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);

    //should I create data_folder here?
    system(("mkdir -p "+rp.data_folder).c_str());    

    std::ofstream outfile_om;
    std::ofstream outfile_data;	
    std::ofstream outfile_histogram;


    bool continue_simulation = false;
    //if required files for restart exist and continue_if_possible flag is set
    if (rp.continue_if_possible & file_exists(rp.data_folder+"data_log.dat") & file_exists(rp.data_folder+"H_N.dat") & file_exists(rp.data_folder+"checkpoint.ckp")){
	    outfile_data.open(rp.data_folder+"data_log.dat", std::ios::app);    
	    outfile_om.open(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary | std::ios::app);    
	    outfile_histogram.open(rp.data_folder+"H_N.dat", std::ios::app); 
	    continue_simulation = true;   
	}
	else{
	    outfile_data.open(rp.data_folder+"data_log.dat");    
	    outfile_om.open(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);    
	    outfile_histogram.open(rp.data_folder+"H_N.dat");  
    	outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<std::endl;
	}
    std::cout<<"Continue? "<<rp.continue_if_possible<<std::endl;


    long t_init=0;
    long t_final=0;//=rp.timesteps;

    //histogram; counts up to 1000 subunits
    std::vector<int> H_N(1000, 0);

    UmbrellaWindow uw(rp.spring_const, mesh.n_faces());

    if (continue_simulation){
    	load_from_checkpoint(mesh, uw, t_final, rp.data_folder+"checkpoint.ckp", input_file);
    }


    while (t_final<rp.timesteps){
    	t_init = t_final;    	
    	t_final = t_init + rp.nsteps_umbrella;
    	//std::cout<<"--- tinit, tfinal "<<t_init<<" "<<t_final<<std::endl;
    	propagate(mesh, t_init, t_final, rp, outfile_om,  outfile_data, uw, eng);    
    	//shift the umbrella window
    	uw.N0+=rp.dN_umbrella;
    	std::cout<<"umbrella N0="<<uw.N0<<std::endl;
    	for (int N: uw.values){
    		H_N[N]++;
    	}
    	for (const auto &hi : H_N) outfile_histogram << hi <<"\t";
    	outfile_histogram << std::endl;

    	uw.values.clear();
    	std::fill(H_N.begin(), H_N.end(), 0);
    }

    outfile_om.close();
    outfile_data.close();
	outfile_histogram.close();    



	/*****checkpoint test*****/
	/*MyMesh mesh_ckp;
	long timestep;
	UmbrellaWindow uwc(0,0);
	load_from_checkpoint(mesh_ckp, uwc, timestep, rp.data_folder+"checkpoint.ckp", input_file);
	std::cout<<"Checkpoint read...."<<std::endl;
	std::cout<<uwc.spring_const<<std::endl;
	std::cout<<uwc.N0<<std::endl;
	std::cout<<timestep<<std::endl;*/
	/**********/
    //do conversions if requested

    //create a conversions directory
    system(("mkdir -p "+rp.data_folder+"conversions").c_str());
    if (rp.convert_to_lammps_trajectory)
		convert_om_to_lammps_trajectory(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/");
    if (rp.convert_to_lammps)
		convert_om_to_lammps_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_");
	if (rp.convert_to_vtk)
		convert_om_to_VTK_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/vtk_", rp.path_to_source+"vtk_colors.conf");    
	if (rp.convert_to_dat)
		convert_om_to_dat_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/dat_");
	//do the WHAM
	//system(("cd "+rp.path_to_source).c_str());
	//system(("cd "+rp.path_to_source + "&& python wham.py "+input_file ).c_str());

	return;	

}


//little helper to switch values of two double variables
/*void swap_uw(UmbrellaWindow & a, UmbrellaWindow & b){
	//UmbrellaWindow tmp=a;
	double spring_const;
	double N0;
	std::vector<int> values;

	spring_const = a.spring_const;
	N0 = a.N0;
	values = a.values;

	a.spring_const = b.spring_const;
	a.N0 = b.N0;
	a.values = b.values;

	b.spring_const = spring_const;
	b.N0 = N0;
	b.values = values;

	return;
}
*/
//work in progress
void run_umbrella_parallel_tempering(std::string PT_input_file){


	
    std::vector<std::string> input_files;


    //read/set input_files here from PT_input_file

    json pt_params;
    std::ifstream input_file(PT_input_file);
    input_file >> pt_params;
    input_file.close();

    pt_params = pt_params["parameters"]["parall_tempering_parameters"];
    long dt_exchange = pt_params["dt_exchange"];
    input_files = pt_params["input_files"].get<std::vector<std::string>>();
    std::cout<<" Exchanging between "<<std::endl;
    for (auto inp: input_files){
    	std::cout<<inp<<std::endl;
    }

    unsigned int nruns = input_files.size();


	std::vector<MyMesh> mesh(nruns);	

    std::vector<RunParameters> rp(nruns);

    std::vector<AcceptanceMonitor> AM(nruns);


    std::ofstream outfile_om[nruns];
    std::ofstream outfile_data[nruns];    
    std::ofstream outfile_histogram[nruns];   
    std::ofstream outfile_k_N0[nruns];

	std::vector<std::vector<int>> H_N( nruns , std::vector<int> (1000, 0));     
    std::vector<UmbrellaWindow> uw;//(nruns);    

    //should "thermalize" all rng engines
	std::vector< std::mt19937 > engines;  
	std::mt19937 exch_engine(std::mt19937(rp[0].ensemble*73 + 17) );               

    for (unsigned int i=0; i<nruns; i++){

    	read_set_run_params(rp[i], input_files[i]);  

		if (rp[i].nsteps_umbrella % dt_exchange != 0 ){
			std::cout<<"simulation #"<<i<<" : nsteps_umbrella % dt_exchange != 0. Exiting..."<<std::endl;
			std::cout<<"This is necessary for proper window shifting. Make sure that nsteps_umbrella % dt_exchange=0 for all simulations."<<std::endl;
			exit(1);
		}    	 

		if (rp[i].dt_burnin % dt_exchange != 0 ){
			std::cout<<"simulation #"<<i<<" : dt_burnin % dt_exchange != 0. Exiting..."<<std::endl;
			std::cout<<"This is necessary for proper window shifting and sampling. Make sure that dt_burnin % dt_exchange=0 for all simulations."<<std::endl;
			exit(1);
		}

		if (rp[i].dt_burnin >= rp[i].nsteps_umbrella ){
			std::cout<<"simulation #"<<i<<" : dt_burnin >= nsteps_umbrella. Exiting..."<<std::endl;
			std::cout<<"Make sure that dt_burnin < nsteps_umbrella for all simulations."<<std::endl;
			exit(1);
		}		
	    //initialize RNG from ensemble seed
	    //std::srand(rp[i].ensemble*73 + 17);
		engines.push_back( std::mt19937(rp[i].ensemble*73 + 17) );	    

		init_mesh(mesh[i], rp[i].init_file, input_files[i], rp[i].init_file_pos);

	    //should I create data_folder here?
	    system(("mkdir -p "+rp[i].data_folder).c_str());    

/*    	outfile_om[i].open(rp[i].data_folder+"snapshots.om", std::ios::out | std::ios::binary);
    	outfile_data[i].open(rp[i].data_folder+"data_log.dat");	
    	outfile_histogram[i].open(rp[i].data_folder+"H_N.dat");    
    	outfile_data[i]<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<std::endl;
*/
    	uw.push_back( UmbrellaWindow(rp[i].spring_const, mesh[i].n_faces()) ); 
    	AM.push_back(AcceptanceMonitor());   	
    }

    long t_init=0;
    long t_final=0;//=rp.timesteps;
    /****** continue? ********/
    bool continue_simulation = true;
    for (unsigned int i=0; i<nruns; i++){
    	if (rp[i].continue_if_possible & file_exists(rp[i].data_folder+"data_log.dat") & file_exists(rp[i].data_folder+"H_N.dat") & file_exists(rp[i].data_folder+"checkpoint.ckp")){
    		//do nothing if files are found; only continue if all files of all systems are there
    	}
    	else{
    		continue_simulation=false;
    	}
    }
    if (continue_simulation){
    	std::cout<<"Continuing.... "<<std::endl;
    	for (unsigned int i=0; i<nruns; i++){
	    	outfile_om[i].open(rp[i].data_folder+"snapshots.om", std::ios::out | std::ios::binary | std::ios::app);
	    	outfile_data[i].open(rp[i].data_folder+"data_log.dat", std::ios::app);	
	    	outfile_histogram[i].open(rp[i].data_folder+"H_N.dat", std::ios::app); 
	    	outfile_k_N0[i].open(rp[i].data_folder+"k_N0.dat", std::ios::app);
	    	//will pick the last t_final here; they should be all equal
	    	load_from_checkpoint(mesh[i], uw[i], t_final, rp[i].data_folder+"checkpoint.ckp", input_files[i]);
    	}
    }
    else{
	    std::cout<<"New simulation.... "<<std::endl;     	
    	for (unsigned int i=0; i<nruns; i++){    	   	
	    	outfile_om[i].open(rp[i].data_folder+"snapshots.om", std::ios::out | std::ios::binary);
	    	outfile_data[i].open(rp[i].data_folder+"data_log.dat");	
	    	outfile_histogram[i].open(rp[i].data_folder+"H_N.dat"); 
	    	outfile_k_N0[i].open(rp[i].data_folder+"k_N0.dat");	    	   
	    	outfile_data[i]<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<std::endl;
	    }
    }
    /*************/


//    auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
    std::vector<OpenMesh::MProp<MeshProp>> mesh_props;// = OpenMesh::MProp<MeshProp>(mesh[0], "mesh_props");
    std::vector<OpenMesh::FProp<FaceProp>> face_props;    
    for (unsigned int i=0; i<nruns; i++){
    	mesh_props.push_back( OpenMesh::MProp<MeshProp>(mesh[i], "mesh_props") );
    	face_props.push_back( OpenMesh::FProp<FaceProp>(mesh[i], "face_props") );    	
    }



    //histogram; counts up to 1000 subunits
    //std::vector<int> H_N(1000, 0);



    unsigned int i;
    //pick the smallest of all timesteps
    long t_min=rp[0].timesteps;
    for (i=0; i<nruns; i++){
    	if (rp[i].timesteps < t_min){
    		t_min = rp[i].timesteps;
    	}
    }

    int ix1, ix2, n1, n2, fpos;
    double r, p, E1, E2, kT1, kT2, wtf, a, b;
    double betaUi, betaUj, betaUi_swapped, betaUj_swapped;
    bool parallel;
    MyMesh * mesh_tmp;
    std::vector<MyMesh*> mesh_ptr;
	std::vector<OpenMesh::MProp<MeshProp>*> mesh_props_ptr;    
    for (i=0; i<nruns; i++){    
    	mesh_ptr.push_back(&mesh[i]);
    	mesh_props_ptr.push_back(&mesh_props[i]);
    }

    auto start = std::chrono::system_clock::now();
    // Some computation here
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::vector<int> n_faces;
    std::vector<int> sorted_i;

    std::ofstream swapfile;
    std::string pt_datafolder = pt_params["data_folder"];
    pt_datafolder=getenv_str("tmp_dir")+pt_datafolder;
    system(("mkdir -p "+pt_datafolder).c_str()); 
    swapfile.open(pt_datafolder+"swap_file.dat");
    omlog();
    //main loop
    while (t_final<t_min){


    	t_init = t_final;    	
    	t_final = t_init + dt_exchange;//rp.nsteps_umbrella; //<--- needs to change to  min( dt_exchange, rp[i].nsteps_umbrella ); timesteps will be the min(rp[i].timesteps)
    	//std::cout<<"--- tinit, tfinal "<<t_init<<" "<<t_final<<std::endl;
    	parallel=true;//(*mesh_ptr[0]).n_faces()>30;
    	//should loop over a sorted array to make sure largest n_faces() meshes get to separate threads
    	//build the vector to be sorted first
    	n_faces.clear();
    	for (i=0; i<nruns; i++){
    		n_faces.push_back( (*mesh_ptr[i]).n_faces() );
    	}

    	//----
    	//#pragma omp parallel for private(i) schedule(static) if(parallel)
/*    	#pragma omp parallel for private(i) if(parallel)    	
    	for (i=0; i<nruns; i++){
    		propagate(*mesh_ptr[i], t_init, t_final, rp[i], outfile_om[i],  outfile_data[i], uw[i], engines[i]); 
    	}
*/
    	sorted_i = sort_indexes(n_faces);
    	/*std::cout<<"----- sorted_i------"<<std::endl;    	
    	for (auto j: sorted_i){

    		std::cout<<j<<std::endl;
    	}*/

    	#pragma omp parallel for private(i) schedule(static, 1) if(parallel)
    	for (i=0; i<nruns; i++){
    		propagate(*mesh_ptr[sorted_i[i]], t_init, t_final, rp[sorted_i[i]], outfile_om[sorted_i[i]],  outfile_data[sorted_i[i]], uw[sorted_i[i]], engines[sorted_i[i]], AM[sorted_i[i]]); 
    	} 

    	//std::cout<<"exchanging..."<<std::endl;
    	//EXCHANGE ATTEMPT HERE ----
    	ix1 = randint(nruns, exch_engine);//int(rand() % nruns);
    	ix2 = randint(nruns, exch_engine);//int(rand() % nruns);

    	n1 = (*mesh_ptr[ix1]).n_faces();
    	n2 = (*mesh_ptr[ix2]).n_faces();
    	betaUi = ( full_energy(*mesh_ptr[ix1]) + uw[ix1].bias_potential(n1) - full_mu_N(*mesh_ptr[ix1]) ) / (**mesh_props_ptr[ix1]).kT;
    	betaUj = ( full_energy(*mesh_ptr[ix2]) + uw[ix2].bias_potential(n2) - full_mu_N(*mesh_ptr[ix2]) ) / (**mesh_props_ptr[ix2]).kT;


		swap_configurations(*mesh_ptr[ix1], *mesh_ptr[ix2]); 
		std::swap(mesh_ptr[ix1], mesh_ptr[ix2]);
		std::swap(mesh_props_ptr[ix1], mesh_props_ptr[ix2]);

		//swap_uw(uw[ix1], uw[ix2]);
    	n1 = (*mesh_ptr[ix1]).n_faces();
    	n2 = (*mesh_ptr[ix2]).n_faces();
    	betaUi_swapped = ( full_energy(*mesh_ptr[ix1]) + uw[ix1].bias_potential(n1) - full_mu_N(*mesh_ptr[ix1]) ) / (**mesh_props_ptr[ix1]).kT;
    	betaUj_swapped = ( full_energy(*mesh_ptr[ix2]) + uw[ix2].bias_potential(n2) - full_mu_N(*mesh_ptr[ix2]) ) / (**mesh_props_ptr[ix2]).kT;

    	p = exp( -betaUi_swapped - betaUj_swapped + betaUi + betaUj );

		//r = rand()/(RAND_MAX + 1.0);	
		r = randdouble(exch_engine);//
		if (r<p){
			//std::cout<<"swapped "<<ix1<<" "<<ix2<<" "<<std::setprecision(10)<<p<<" "<<betaUi<<" "<<betaUi_swapped <<" "<<-betaUi_swapped - betaUj_swapped + betaUi + betaUj<<std::endl;
			swapfile<<t_final<<" "<<ix1<<" "<<ix2<<std::endl;
		}
		else{
			//std::cout<<"swap rejected, swapping back "<<ix1<<" "<<ix2<<std::endl;
			swap_configurations(*mesh_ptr[ix1], *mesh_ptr[ix2]); 
			std::swap(mesh_ptr[ix1], mesh_ptr[ix2]);
			std::swap(mesh_props_ptr[ix1], mesh_props_ptr[ix2]);

			//swap_uw(uw[ix1], uw[ix2]);			
			//switch capsids
			//mesh_tmp
			//std::cout<<std::setprecision(10)<<"pre "<<full_binding_energy(mesh[ix1])<<"  "<<full_binding_energy(mesh[ix2])<<std::endl;
			//std::cout<<std::setprecision(10)<<"pre "<<full_mu_N(mesh[ix1])<<"  "<<full_mu_N(mesh[ix2])<<std::endl;

			//std::cout<<"----> swap "<<ix1<<" "<<ix2<<std::endl;
			//std::cout<<"E_b "<< std::setprecision(40)<<	full_binding_energy(mesh[ix1])<<"  "<<full_binding_energy(mesh[ix2])<<std::endl;
			//std::cout<<"E_elast "<< std::setprecision(40)<<	full_elastic_energy(mesh[ix1])<<"  "<<full_elastic_energy(mesh[ix2])<<std::endl;
			//std::cout<<"mu N "<< std::setprecision(40)<<"post "<<full_mu_N(mesh[ix1])<<"  "<<full_mu_N(mesh[ix2])<<std::endl;						
			//swap_configurations(mesh[ix1], mesh[ix2]);
			//std::cout<<"   sw E_b "<< std::setprecision(20)<<	full_binding_energy(mesh[ix1])<<"  "<<full_binding_energy(mesh[ix2])<<std::endl;			
			//swap_configurations(mesh[ix1], mesh[ix2]);			
//
			//std::cout<<"<-------- "<<std::endl;
			//std::cout<<"E_b "<< std::setprecision(40)<<	full_binding_energy(mesh[ix1])<<"  "<<full_binding_energy(mesh[ix2])<<std::endl;
			//std::cout<<"E_elast "<< std::setprecision(40)<<	full_elastic_energy(mesh[ix1])<<"  "<<full_elastic_energy(mesh[ix2])<<std::endl;
			//std::cout<<"mu N "<< std::setprecision(40)<<"post "<<full_mu_N(mesh[ix1])<<"  "<<full_mu_N(mesh[ix2])<<std::endl;
			
			//std::cout<<"swapped "<<ix1<<" "<<ix2<<std::endl;
			//std::cout<<std::setprecision(10)<<"post "<<full_binding_energy(mesh[ix1])<<"  "<<full_binding_energy(mesh[ix2])<<std::endl;			//switch temperatures too; better define a switch_capsid function
			//std::cout<<std::setprecision(10)<<"post "<<full_mu_N(mesh[ix1])<<"  "<<full_mu_N(mesh[ix2])<<std::endl;	
			//this can get more complicated when tempering for e_b, mu, etc. because then prototype subunits have to be changed too
		}

    	//shift the umbrella windows;
    	//individually?... depending on each nsteps_umbrella?

     	for (i=0; i<nruns; i++){   
     		if (t_final % rp[i].nsteps_umbrella == 0){	
				outfile_k_N0[i]<<uw[i].N0<<"\t"<< uw[i].spring_const <<std::endl;


		    	uw[i].N0+=rp[i].dN_umbrella;
		    	//std::cout<<"umbrella N0="<<uw.N0<<std::endl;
		    	for (int N: uw[i].values){
		    		H_N[i][N]++;
		    	}
		    	for (const auto &hi : H_N[i]) outfile_histogram[i] << hi <<"\t";
		    	outfile_histogram[i] << std::endl;

		    	uw[i].values.clear();
	    		std::fill(H_N[i].begin(), H_N[i].end(), 0);

	    		if (i==0){
	    			end = std::chrono::system_clock::now();
    				std::cout<< "elapsed time: " << elapsed_seconds.count() <<" "<<(*mesh_ptr[0]).n_faces()<< std::endl;
    				elapsed_seconds = end-start;
    				//start=end;
	    		}
	    	}

	    	//clear umbrella window values after burnin time is reached; that's when sampling begins
	    	if ((t_final - rp[i].dt_burnin) % rp[i].nsteps_umbrella == 0 ){
	    		uw[i].values.clear();
	    	}

	    	//need dt dump here, not nstep umbrella
	    	/*if (t_final % rp[i].nsteps_umbrella == 0){	
		    	if (t_final>rp[i].dtskip){
					fpos = dump_om(*mesh_ptr[i], t_final, outfile_om[i]);
				}
				else{
					fpos=-1;
				}

				dump_data(*mesh_ptr[i], t_final, outfile_data[i], fpos);
			}*/
    	}
    	//end it if time limit is passed
    	end = std::chrono::system_clock::now();
    	elapsed_seconds = end-start;
    	if (elapsed_seconds.count()>=rp[0].max_seconds){
    		std::cout<<"Time limit reached, cancelling, saving...."<<std::endl;
    		break;
    	}
    }

    for (i=0; i<nruns; i++){   
    	outfile_om[i].close();
    	outfile_data[i].close();
		outfile_histogram[i].close();
		outfile_k_N0[i].close();    
	}
	swapfile.close();
    //do conversions if requested

    for (i=0; i<nruns; i++){   
	    //create a conversions directory
	    system(("mkdir -p "+rp[i].data_folder+"conversions").c_str());
	    if (rp[i].convert_to_lammps_trajectory)
			convert_om_to_lammps_trajectory(rp[i].data_folder+"snapshots.om", rp[i].data_folder+"conversions/");
	    if (rp[i].convert_to_lammps)
			convert_om_to_lammps_snapshots(rp[i].data_folder+"snapshots.om", rp[i].data_folder+"conversions/lammps_");
		if (rp[i].convert_to_vtk)
			convert_om_to_VTK_snapshots(rp[i].data_folder+"snapshots.om", rp[i].data_folder+"conversions/vtk_", rp[i].path_to_source+"vtk_colors.conf");    
		if (rp[i].convert_to_dat)
			convert_om_to_dat_snapshots(rp[i].data_folder+"snapshots.om", rp[i].data_folder+"conversions/dat_");
	}
	//be careful when parallelizing, there are environmental variables set
	/*------------------ For now, do WHAMs as separate jobs
	#pragma omp parallel for private(i) schedule(static, 1) if(parallel)
    for (i=0; i<nruns; i++){   	
		//do the WHAM
		//system(("cd "+rp[i].path_to_source).c_str());
		system(( "cd "+rp[i].path_to_source + "&& python wham.py "+input_files[i] ).c_str());
	}
	

	//DO THE JOINT WHAM HERE
	system(( "cd "+rp[0].path_to_source + "&& python multi_wham.py "+PT_input_file ).c_str());
	------------------*/
	return;	
}





void propagate_wall(MyMesh & mesh, long t_init, long t_final, RunParameters rp, std::ostream& outfile_om,  std::ostream& outfile_data, std::ostream& outfile_defects ,std::mt19937 & eng, Wall & wall, double dz_wall_per_timestep, bool all_moves, AcceptanceMonitor & AM = default_AM){

    int nmoves = 7;
	bool accepted;
	int which_move, which_type, which_rotation, number_of_rotational_configs;
	double k_insertion, k_fusion, p_propose, r;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> , std::vector<int>> wedges;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > fusion_vectors, type2_fission_vectors;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> >	fission_vectors;	
	std::vector <MyMesh::VertexHandle> w1, v1;  
	std::vector<MyMesh::HalfedgeHandle> boundary_halfedges, h1;	  
	std::vector<MyMesh::FaceHandle> removable_faces_geom, removable_faces;	
	std::tuple<std::vector<MyMesh::HalfedgeHandle>, std::vector<MyMesh::HalfedgeHandle> > fusion_halfedges;	
	std::vector<MyMesh::EdgeHandle> fission_edges;	

	std::vector<MyMesh::VertexHandle*> handle_tracking_v;
		std::vector<MyMesh::FaceHandle*> handle_tracking_f;
	std::vector<MyMesh::HalfedgeHandle*> handle_tracking_h;	

	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto proto_face_props = OpenMesh::FProp<FaceProp>((*mesh_props).prototypes_mesh, "face_props");		
	float kT = (*mesh_props).kT;
	double l_fuse = (*mesh_props).l_fuse;	

	//std::ofstream outfile("outfileom5.om", std::ios::out | std::ios::binary);
	// std::ofstream datafile("data_out5.dat");
	long long fpos;
	long long opos;

	std::tuple< double, double > asphericities;

//dump_data(mesh, -1, outfile_data, fpos);
	for (long t=t_init; t<t_final; t++){
		if ((t-rp.dtskip) % rp.dtsave ==0){		
			std::cout<<"t = "<<t<<" n_faces="<<mesh.n_faces()<<std::endl;	
			
			//skip full dump for the first dtskip timesteps
			
			if ( (t>=rp.dtskip) && ((t-rp.dtskip) % rp.dtsave_snapshot==0)){
					fpos = dump_om(mesh, t, outfile_om);

			}
			else{
				fpos=-1;
			}

			for (int iac=0; iac<(nmoves); iac++){
				AM.r_acc[iac] = AM.n_acc[iac]/AM.n_tot[iac];
				//std::cout<<"r_acc "<< AM.r_acc[iac]<<std::endl;
				//std::cout<<"n_acc "<< AM.n_acc[iac]<<std::endl;
				//std::cout<<"n_tot "<< AM.n_tot[iac]<<std::endl;								
				AM.n_acc[iac]=0;
				AM.n_tot[iac]=0;
			}

			dump_data(mesh, t, outfile_data, fpos, AM.r_acc);
			opos = outfile_data.tellp();
  			outfile_data.seekp (opos-1);
  			//outfile.write (" sam",4);
			//outfile_data<<"\b";
			outfile_data<<"\t"<<wall.wall_energy_1(mesh)<<"\t"<<wall.wall_energy_2(mesh)<<"\t"<<full_stretch_energy(mesh)<<"\t"<<full_bending_energy(mesh)<<"\t"<<full_binding_energy(mesh)<<"\t"<<wall.wall_force_1(mesh)<<"\t"<<wall.wall_force_2(mesh)<<"\t"<<wall.z_position_1<<"\t"<<wall.z_position_2;
			asphericities = compute_asphericity(mesh);
			outfile_data<<"\t"<< count_boundary_edges(mesh) <<"\t"<< std::get<0>(asphericities)<<"\t"<<std::get<1>(asphericities)<<std::endl;
			defect_analyze_snapshot(mesh, t, outfile_defects);
			
			//dump_om(mesh, "omtest"+std::to_string(t)+".om");
		}

		if ((t-rp.dtskip) % rp.checkpoint_freq ==0){
			create_checkpoint_wall(mesh, wall, t, rp.data_folder+"checkpoint.ckp");			
		}
		/*for (int n=0; n<(int)mesh.n_vertices(); n++){
			attempt_move_wall(mesh, eng, wall);
		}*/
		

		//!! IMPORTANT: in each t iteration, only one move should be selected because the garbage collector
		// is only called once, at the end of each iteration. n_faces, n_vertices, etc won't get updated unless
		// the garbage collector is called and we need them to be up to date for the bias potential and also for attempt_move()
		//which_move = rand() % nmoves;
		which_move = randint(nmoves, eng);		
		//std::cout<<"which move "<<which_move<<std::endl;
		accepted=false;
		if (!all_moves) which_move = 6; //only vertex moves
		switch (which_move){
			case 0:
				k_fusion = (*mesh_props).k_fusion; 
				fusion_vectors = get_type1_fusion_triplets(mesh, l_fuse);	
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size();
				if ( (p_propose>1) && (!rp.adaptive_rates) ) std::cout<<"Warning! Decrease k_fusion rate! - type1_fusion"<<std::endl;
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion rate... - type1_fusion"<<std::endl;
					(*mesh_props).k_fusion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}				
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);
				if (r<p_propose){	
					accepted = attempt_type1_fusion_wall(mesh, eng, wall);
					AM.n_acc[which_move+1]+=accepted;
				}
				AM.n_tot[which_move+1]++;
				//if (accepted) std::cout<<"type1 fusion"<<std::endl;
				break;

			case 1:
				k_fusion = (*mesh_props).k_fusion; 
				fission_vectors = get_type1_fission_pairs(mesh);	
    			v1 = std::get<0>(fission_vectors);					
				p_propose = k_fusion * v1.size();
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion rate! - type1_fission"<<std::endl;		
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion rate... - type1_fission"<<std::endl;
					(*mesh_props).k_fusion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}							
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	
					accepted = attempt_type1_fission_wall(mesh, eng, wall);	
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type1 fission"<<std::endl;
				AM.n_tot[which_move+1]++;
				break;

			case 2:		
				k_fusion = (*mesh_props).k_fusion2; 
				fusion_vectors = get_type2_fusion_triplets(mesh, l_fuse);	
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size()*0.5;
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fusion"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion2 rate... - type2_fusion"<<std::endl;
					(*mesh_props).k_fusion2*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}								
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_type2_fusion_wall(mesh, eng, wall);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type2 fusion "<<p_propose<<std::endl;		
				AM.n_tot[which_move+1]++;
				break;
			case 3:			
				k_fusion = (*mesh_props).k_fusion2; 
				type2_fission_vectors = get_type2_fission_triplets(mesh);	
    			v1 = std::get<0>(type2_fission_vectors);					
				p_propose = k_fusion * v1.size()*0.5;
				if ( (p_propose>1) && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fission"<<std::endl;			
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion2 rate... - type2_fission"<<std::endl;
					(*mesh_props).k_fusion2*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}						
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_type2_fission_wall(mesh, eng, wall);	
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type2 fission "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;
				break;



			case 4:
				k_fusion = (*mesh_props).k_fusion_edge;
				fusion_halfedges = get_halfedge_fusion_pairs(mesh, l_fuse);
    			h1 = std::get<0>(fusion_halfedges);
    			p_propose = k_fusion * h1.size()*0.5;
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fusion"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion_edge rate... - edge_fusion"<<std::endl;
					(*mesh_props).k_fusion_edge*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}	
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){					
					accepted = attempt_edge_fusion_wall(mesh, eng, wall);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"EDGE fusion " <<" "<<p_propose<<std::endl;
				AM.n_tot[which_move+1]++;
				break;

			case 5:
				k_fusion = (*mesh_props).k_fusion_edge;
				fission_edges = get_fission_edges(mesh);
    			p_propose = k_fusion * fission_edges.size();	
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fission"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion_edge rate... - edge_fusion"<<std::endl;
					(*mesh_props).k_fusion_edge*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}					
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){    									
					accepted = attempt_edge_fission_wall(mesh, eng, wall);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"EDGE fission " <<" "<<p_propose<<std::endl;
				AM.n_tot[which_move+1]++;
				break;

			case 6:
				for (int n=0; n<(int)mesh.n_vertices(); n++){
					AM.n_acc[0]+=attempt_move_wall(mesh, eng, wall);
					AM.n_tot[0]++;
				}
				break;

		}

		//get all face's neighbor list handles to save them through the garbage collector
    	handle_tracking_f.clear();
    	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
    		for (auto & neighbor: face_props[*fit].neighbor_list){
    			handle_tracking_f.push_back( &neighbor);
    		}
    	}
		mesh.garbage_collection<std::vector<MyMesh::VertexHandle*>, std::vector<MyMesh::HalfedgeHandle*>, std::vector<MyMesh::FaceHandle*> >(handle_tracking_v, handle_tracking_h, handle_tracking_f);

		if (t % 150 ==0){
			update_full_neighbor_list(mesh);
			//Could check overlaps for the full capsid here to verify
		}

	wall.move_1(-dz_wall_per_timestep);
	wall.move_2(dz_wall_per_timestep);	
	}


	return;
}


void run_wall(std::string input_file){

	int n_elastic_half_periods=20;

    auto start = std::chrono::system_clock::now();
    // Some computation here
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;


	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);   

    //initialize RNG from ensemble seed
    //std::srand(rp.ensemble*73 + 17);
    std::mt19937 eng(rp.ensemble*73 + 17);

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);

    //should I create data_folder here?
    system(("mkdir -p "+rp.data_folder).c_str());    

    /*std::ofstream outfile_om(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);
    std::ofstream outfile_data(rp.data_folder+"data_log.dat");	
    outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<"\t"<<"E_wall_1"<<"\t"<<"E_wall_2"<<"\tE_stretch\tE_bend\tE_bind\tF_wall_1\tF_wall_2\tz_wall_1\tz_wall_2"<<std::endl;
	*/
    std::ofstream outfile_om;
    std::ofstream outfile_data;	
    std::ofstream outfile_defects;	    

    bool continue_simulation = false;
    //if required files for restart exist and continue_if_possible flag is set
    if (rp.continue_if_possible & file_exists(rp.data_folder+"data_log.dat") & file_exists(rp.data_folder+"checkpoint.ckp")){
	    outfile_data.open(rp.data_folder+"data_log.dat", std::ios::ate | std::ios::in);  
	    outfile_defects.open(rp.data_folder+"topology_defect.dat", std::ios::ate | std::ios::in);   	      
	    outfile_om.open(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary | std::ios::app);    
	    continue_simulation = true;   
	    std::cout<<"Continuing existing simulation... "<<std::endl;
	}
	else{
	    outfile_data.open(rp.data_folder+"data_log.dat");  
	    outfile_defects.open(rp.data_folder+"topology_defect.dat"); 	      
	    outfile_om.open(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);    
    	outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<"\t"<<"E_wall_1"<<"\t"<<"E_wall_2"<<"\tE_stretch\tE_bend\tE_bind\tF_wall_1\tF_wall_2\tz_wall_1\tz_wall_2"<<std::endl;
	    std::cout<<"Starting new simulation... "<<std::endl;
	}
    std::cout<<"Continue? "<<continue_simulation<<std::endl;


    long t_init=0;
    long t_final=rp.timesteps;
    //UmbrellaWindow uw_fake(0.0, 0.0);
    //shift the structure to COM=(0,0,0) and give it a random orientation
    std::cout<<"Energy before "<<full_elastic_energy(mesh)<<std::endl;

    //new simulation; shift to center, give random rotation
    double R=0.0;
    if (!continue_simulation){
	    zero_com_shift_mesh(mesh);
	    //get approximate radius to set the walls:
	    
	    MyMesh::Point p;
		for(MyMesh::VertexIter vit = mesh.vertices_sbegin(); vit != mesh.vertices_end(); ++vit) {
			p = mesh.point(*vit);
			R+=p.norm();
		}
		R/=mesh.n_vertices();
		std::cout<<"R="<<R<<std::endl;

		random_rotate_mesh(mesh, eng);
	   /* for (int t=0; t<1000; t++){
	    	random_rotate_mesh(mesh, eng);
	    	dump_om(mesh, t, outfile_om);
		}*/
	    std::cout<<"Energy after "<<full_elastic_energy(mesh)<<std::endl;    
	}
    ///Wall(double hardness, double z_position_1, double z_position_2, double amplitude)
    //(100.0, 5.8, 0.0, 10.0);
    Wall wall(rp.wall_hardness, 1.1*R, -1.1*R, rp.wall_amplitude);


    if (continue_simulation){
    	load_from_checkpoint_wall(mesh, wall, t_init, rp.data_folder+"checkpoint.ckp", input_file);
    	//t_init+=1;
    }

    //!!!!!!!!!!!!!!!!!!!
    //should do a random rotation here of the structure

    if (!continue_simulation){
	    //thermalize before loading; if necessary, add code to optimize d_max and l_fuse by the acceptance rate; hopefully not necessary for small FvK number
	    propagate_wall(mesh, 0, rp.dt_burnin, rp, outfile_om,  outfile_data,outfile_defects , eng, wall, 0.0, false);
	    std::cout<<"Energy thermalized "<<full_elastic_energy(mesh)<<" t_init="<< t_init<<std::endl;      
    /*for (int t=0; t<rp.dt_burnin; t++){
		for (int n=0; n<(int)mesh.n_vertices(); n++){
			attempt_move_wall(mesh, eng, wall);
		}
	}*/

	    //wall = Wall(hardness=10, z_position = 0.2+diameter, amplitude=10.0)
		std::cout<<"pre "<<  rp.wall_period <<std::endl;
	}


    int n_half_periods = (int)((t_final-t_init) / (rp.wall_period/2 ));
    std::cout<<"n half periods "<<n_half_periods<<std::endl;
    bool crack_allowed=false;

    //start loading elastically
    //increase wall speed to match the one with all 7 moves enabled
    long full_wall_period=rp.wall_period; //save it in a variable to avoid division+multiplication error
    double full_wall_speed=rp.wall_speed;
    rp.wall_period/=7; //10
    rp.wall_speed*=7;  //10

//PUT IT BACK!!!!!!!
    if (!continue_simulation){
	    if (n_half_periods>0){
		    for (int halfperiod=0; halfperiod<n_elastic_half_periods; halfperiod++){
		    	propagate_wall(mesh, t_init+rp.dt_burnin+halfperiod*rp.wall_period/2, t_init+rp.dt_burnin+(halfperiod+1)*rp.wall_period/2, rp, outfile_om,  outfile_data,outfile_defects , eng, wall, rp.wall_speed, crack_allowed);
		    	rp.wall_speed*=-1.0;
			}
		}
	}
	
	//steady load if time too small to do one half period
	//else{
	//	propagate_wall(mesh, t_init+rp.dt_burnin, rp.dt_burnin+t_final, rp, outfile_om,  outfile_data, eng, wall, rp.wall_speed, crack_allowed);
	//}

	crack_allowed=true;
    //rp.wall_period*=7; //10;
    //rp.wall_speed/=7; //10;	
    rp.wall_period=full_wall_period; //10;
    rp.wall_speed=full_wall_speed; //10;
    std::cout<<"WALL "<<rp.wall_period<<" "<<rp.wall_speed<<std::endl;	    
    //start loading with crack
    if (!continue_simulation){
	    if (n_half_periods>0){
		    for (int halfperiod=n_elastic_half_periods; halfperiod<2*n_half_periods; halfperiod++){
	    	//one extra half period for first crack only
	    	//for (int halfperiod=n_half_periods; halfperiod<n_half_periods+1; halfperiod++){
		    	propagate_wall(mesh, t_init+rp.dt_burnin+halfperiod*rp.wall_period/2, t_init+rp.dt_burnin+(halfperiod+1)*rp.wall_period/2, rp, outfile_om,  outfile_data, outfile_defects ,eng, wall, rp.wall_speed, crack_allowed);
		    	rp.wall_speed*=-1.0;


    			end = std::chrono::system_clock::now();
				std::cout<< "elapsed time: " << elapsed_seconds.count() << std::endl;
				elapsed_seconds = end-start;		 

		    	if (elapsed_seconds.count()>=rp.max_seconds){
		    		std::cout<<"Time limit reached, cancelling, saving...."<<std::endl;
		    		break;
		    	}				   	
			}
		}
	}
	else{
		//remaining time till the next halfperiod end
		long t_remainder=((t_init-rp.dt_burnin)/(rp.wall_period/2)+1)*(rp.wall_period/2)+rp.dt_burnin-t_init;
		int halfperiod_count=((t_init-rp.dt_burnin)/(rp.wall_period/2)+1);
		if (halfperiod_count%2==0){
			rp.wall_speed*=-1.0;
		}
		std::cout<<"t_remainder: "<< t_remainder<<std::endl;
		std::cout<<"t_init: "<< t_init<<std::endl;		
	    if (n_half_periods>0){
		    for (int halfperiod=n_elastic_half_periods; halfperiod<2*n_half_periods; halfperiod++){
	    	//one extra half period for first crack only
	    	//for (int halfperiod=n_half_periods; halfperiod<n_half_periods+1; halfperiod++){
		    	propagate_wall(mesh, t_init, t_init+t_remainder, rp, outfile_om,  outfile_data,outfile_defects , eng, wall, rp.wall_speed, crack_allowed);

		    	//propagate_wall(mesh, t_init+rp.dt_burnin+halfperiod*rp.wall_period/2, t_init+rp.dt_burnin+(halfperiod+1)*rp.wall_period/2, rp, outfile_om,  outfile_data, eng, wall, rp.wall_speed, crack_allowed);
		    	rp.wall_speed*=-1.0;
		    	//t_init+=rp.wall_period/2;
		    	t_init+=t_remainder;
		    	t_remainder=rp.wall_period/2;



	    		end = std::chrono::system_clock::now();
				std::cout<< "elapsed time: " << elapsed_seconds.count() <<std::endl;
				elapsed_seconds = end-start;

		    	if (elapsed_seconds.count()>=rp.max_seconds){
		    		std::cout<<"Time limit reached, cancelling, saving...."<<std::endl;
		    		break;
		    	}				
			}
		}		
	}
	//steady load if time too small to do one half period
	//else{
	//	propagate_wall(mesh, t_init+rp.dt_burnin, rp.dt_burnin+t_final, rp, outfile_om,  outfile_data, eng, wall, rp.wall_speed, crack_allowed);
	//}


    //propagate_wall(mesh, t_init+rp.dt_burnin, t_final+rp.dt_burnin, rp, outfile_om,  outfile_data, eng, wall, rp.wall_speed, true);
    outfile_om.close();
    outfile_data.close();
    outfile_defects.close();


    //do conversions if requested

    //create a conversions directory
    system(("mkdir -p "+rp.data_folder+"conversions").c_str());
    if (rp.convert_to_lammps_trajectory)
		convert_om_to_lammps_trajectory(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/");
    if (rp.convert_to_lammps)
		convert_om_to_lammps_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_");
	if (rp.convert_to_vtk)
		convert_om_to_VTK_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/vtk_", rp.path_to_source+"vtk_colors.conf");    
	if (rp.convert_to_dat)
		convert_om_to_dat_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/dat_");
	return;
}


void rand_orientation_test(){

	MyMesh::Point p_from = MyMesh::Point(0,0,1);
	MyMesh::Point axis;
	MyMesh::Point p_to, p;
	std::mt19937 eng(73 + 17);
	double angle;
	for (int i=0; i<1000; i++){
		p_to = randvec(eng);
		p_to = p_to.normalize();
		p_to*=p_from.norm();

		axis = p_from.cross(p_to);
		axis = axis.normalize();
		p = p_from;
		angle = OpenMesh::angle( p_from.dot(p_to), (p_from.cross(p_to)).norm() );
		rotatevec(p, axis, angle);

		//std::cout<<p[0]<<" "<<p0[0]<<std::endl;
		std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
		std::cout<<p_to[0]<<" "<<p_to[1]<<" "<<p_to[2]<<std::endl;		
	}


}


//Warning: does not update any neighbor list to make it faster
//if overlaps become a problem, this needs to be rethought and properly separated in a propagate - run pair of functions as all the others
void run_quench(std::string input_file){
	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);   

    //initialize RNG from ensemble seed
    //std::srand(rp.ensemble*73 + 17);
    std::mt19937 eng(rp.ensemble*73 + 17);

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);

    //should I create data_folder here?
    system(("mkdir -p "+rp.data_folder).c_str());    

    std::ofstream outfile_om(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);
    std::ofstream outfile_data(rp.data_folder+"data_log.dat");	
    outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full\tkT\tMQS"<<std::endl;


    long t_init=0;
    long t_final=rp.timesteps;
    long long fpos;
    long long opos;


	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");


	//before everything, the structure needs to be heated and equilibrated because we start from a cylinder
	//and that's far from the equilibrium shape of the trumpet
	//so give a T>>1 in the input json file
	double T_init = (*mesh_props).kT;
	double T_final = 1e-5;
	double T_old, T_new;
	double acc = 0.0;
	double tot=0.0;

	double quench_factor=pow(T_init/T_final, 1.0/rp.timesteps);
	double dT = (T_init - T_final)/(rp.timesteps);
	std::cout<<"dT "<<dT<<std::endl;
	std::cout<<"1. d_max "<<(*mesh_props).d_max<<std::endl;
	//bring it up to large temperature
	for (long t=0; t<rp.timesteps; t++){
		if (t % rp.dtsave ==0 ){		
			std::cout<<"t = "<<t<<" n_faces="<<mesh.n_faces()<<std::endl;	
			std::cout<<"2. d_max "<<(*mesh_props).d_max<<std::endl;
			//skip full dump for the first dtskip timesteps
			
			if (t>rp.dtskip && (t % (100*rp.dtsave) ==0)){
				fpos = dump_om(mesh, t, outfile_om);
			}
			else{
				fpos=-1;
			}

			dump_data(mesh, t, outfile_data, fpos);
			opos = outfile_data.tellp();
  			outfile_data.seekp (opos-1);
			outfile_data<<"\t"<< (*mesh_props).kT;// <<std::endl;
			outfile_data<<"\t"<<acc/tot;
			outfile_data<<"\t"<<full_stretch_energy(mesh);
			outfile_data<<"\t"<<compute_mean_quadratic_strain(mesh)<<std::endl;			
			acc=0;
			tot=0;

		}

		for (int n=0; n<(int)mesh.n_vertices(); n++){
			acc+=attempt_move(mesh, eng);
			tot++;
		}

	}	


	std::cout<<"T init "<<T_init<<std::endl;
	std::cout<<" (*mesh_props).kT "<<(*mesh_props).kT<<std::endl;
	std::cout<<" (*mesh_props).d_max "<<(*mesh_props).d_max<<std::endl;		
	//quench with thermal MC
	int nsub=0;
	double kT = (*mesh_props).kT;
	double d_max = (*mesh_props).d_max;
	//(*mesh_props).kT = T_test;
	acc=0.0;
	tot=0.0;
	for (long t=rp.timesteps; t<2*rp.timesteps; t++){
	//for (long t=0; t<rp.timesteps; t++){	
		if (t % rp.dtsave ==0){		
			std::cout<<"t = "<<t<<" n_faces="<<mesh.n_faces()<<std::endl;	
			std::cout<<"3. d_max "<<(*mesh_props).d_max<<std::endl;

			//skip full dump for the first dtskip timesteps
			
			if (t>rp.dtskip && (t % (100*rp.dtsave) ==0)){
				fpos = dump_om(mesh, t, outfile_om);
			}
			else{
				fpos=-1;
			}

			dump_data(mesh, t, outfile_data, fpos);
			opos = outfile_data.tellp();
  			outfile_data.seekp (opos-1);
			outfile_data<<"\t"<< (*mesh_props).kT;// <<std::endl;
			outfile_data<<"\t"<<acc/tot;
			outfile_data<<"\t"<<full_stretch_energy(mesh);//<<std::endl;				
			outfile_data<<"\t"<<compute_mean_quadratic_strain(mesh)<<std::endl;				
			acc=0;
			tot=0;
		}



		for (int n=0; n<(int)mesh.n_vertices(); n++){
			acc+=attempt_move(mesh, eng);
			tot++;
		}
		//(*mesh_props).kT/=quench_factor;
		//(*mesh_props).d_max/=sqrt(quench_factor);
		T_old = kT;
		//T_test = (*mesh_props).kT;
		//T_test-=dT;
		//kT-=dT;
		kT/=quench_factor;
		(*mesh_props).kT =kT; //T_test;//(*mesh_props).kT - dT; nsub++; //T_test-=dT;
		T_new = kT;

		d_max/=sqrt( quench_factor);
		(*mesh_props).d_max=d_max;
		//T_test-=dT;
		//std::cout<<T_test<<" "<<(*mesh_props).kT<<std::endl;

	}


	std::cout<<"nsubtr "<<nsub<<std::endl;
	std::cout<<"Tnew "<<T_new<<std::endl;
	std::cout<<"Told "<<T_old<<std::endl;
	std::cout<<"kT "<<kT<<std::endl;	
	std::cout<<" (*mesh_props).kT "<<(*mesh_props).kT<<std::endl;
	std::cout<<" (*mesh_props).d_max "<<(*mesh_props).d_max<<std::endl;	
	double d_max_start = (*mesh_props).d_max;
	std::cout<<"d_max_start "<<d_max_start<<std::endl;
	//from here on do athermal move only, but still decrease d_max to allow for more accurate convergence
	double ddmax = d_max_start/(2*rp.timesteps);
	acc=0.0;
	tot=0.0;	
	for (long t=2*rp.timesteps; t<4*rp.timesteps; t++){
		if (t % rp.dtsave ==0){		
			std::cout<<"t = "<<t<<" n_faces="<<mesh.n_faces()<<std::endl;	
			std::cout<<"4. d_max "<<(*mesh_props).d_max<<std::endl;
			//skip full dump for the first dtskip timesteps
			
			if (t>rp.dtskip && (t % (100*rp.dtsave) ==0) ){
				fpos = dump_om(mesh, t, outfile_om);
			}
			else{
				fpos=-1;
			}

			dump_data(mesh, t, outfile_data, fpos);
			opos = outfile_data.tellp();
  			outfile_data.seekp (opos-1);
			outfile_data<<"\t"<< (*mesh_props).kT;// <<std::endl;
			outfile_data<<"\t"<<acc/tot;
			outfile_data<<"\t"<<full_stretch_energy(mesh);//<<std::endl;	
			outfile_data<<"\t"<<compute_mean_quadratic_strain(mesh)<<std::endl;							
			acc=0;
			tot=0;
		}


		for (int n=0; n<(int)mesh.n_vertices(); n++){
			acc+=attempt_move_athermal(mesh, eng);
			tot++;
		}

		//(*mesh_props).d_max-=ddmax;
		if ( (acc/tot)<0.2){
				d_max/=1.1;
				(*mesh_props).d_max=d_max;
		}
		else{
				d_max*=1.1;
				(*mesh_props).d_max=d_max;
		}
		//d_max/=sqrt( quench_factor);
		//(*mesh_props).d_max=d_max;
		
	}	
    //propagate(mesh, t_init, t_final, rp, outfile_om,  outfile_data, uw_fake, eng);

	//finally, decrease d_max exponentially
/*	d_max_start = (*mesh_props).d_max;
	double q = 1.001;
	for (long t=3*rp.timesteps; t<4*rp.timesteps; t++){
		if (t % rp.dtsave ==0){		
			std::cout<<"t = "<<t<<" n_faces="<<mesh.n_faces()<<std::endl;	
			std::cout<<"4. d_max "<<(*mesh_props).d_max<<std::endl;
			//skip full dump for the first dtskip timesteps
			
			if (t>rp.dtskip && (t % (100*rp.dtsave) ==0) ){
				fpos = dump_om(mesh, t, outfile_om);
			}
			else{
				fpos=-1;
			}

			dump_data(mesh, t, outfile_data, fpos);
			opos = outfile_data.tellp();
  			outfile_data.seekp (opos-1);
			outfile_data<<"\t"<< (*mesh_props).kT <<std::endl;

		}


		for (int n=0; n<(int)mesh.n_vertices(); n++){
			attempt_move_athermal(mesh, eng);
		}

		(*mesh_props).d_max/=q;
	}	
*/


    outfile_om.close();
    outfile_data.close();	


    std::ofstream outfile_om_last(rp.data_folder+"last_snapshots.om", std::ios::out | std::ios::binary); 
    dump_om(mesh, 0, outfile_om_last);
    outfile_om_last.close();   


    //create a conversions directory
    system(("mkdir -p "+rp.data_folder+"conversions").c_str());
    if (rp.convert_to_lammps_trajectory)
		convert_om_to_lammps_trajectory(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/");
    if (rp.convert_to_lammps)
		convert_om_to_lammps_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_");
	if (rp.convert_to_vtk)
		convert_om_to_VTK_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/vtk_", rp.path_to_source+"vtk_colors.conf");    
	if (rp.convert_to_dat)
		convert_om_to_dat_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/dat_");
	return;


}



//dump mesh in a separate .om file, no timestep
int dump_off(MyMesh mesh, std::string off_filename){
	MyMesh simp_mesh = get_stripped_mesh(mesh);
		
	try
	  {
	    if ( !OpenMesh::IO::write_mesh(simp_mesh, off_filename) )
	    {
	      std::cerr << "Cannot write mesh to stream" << std::endl;
	      return 0;
	    }
	  }
	  catch( std::exception& x )
	  {
	    std::cerr << x.what() << std::endl;
	    return 0;
	  }	

	//add separator in case keys are lost
	//probably no nead; the reader should know where the object ends  
	//os<<"-----------------";
	return 1;
}



void run_umbrella_cycling_referee(std::string input_file){

	//create_single_subunit_om("single_subunit_trumpet.om", 0, 0, 1, 2);

	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);   

    //initialize RNG from ensemble seed
    //std::srand(rp.ensemble*73 + 17);
    std::mt19937 eng(rp.ensemble*73 + 17);    

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);

    //should I create data_folder here?
    system(("mkdir -p "+rp.data_folder).c_str());    

    std::ofstream outfile_om;
    std::ofstream outfile_data;	
    std::ofstream outfile_histogram;


    bool continue_simulation = false;
    //if required files for restart exist and continue_if_possible flag is set
    if (rp.continue_if_possible & file_exists(rp.data_folder+"data_log.dat") & file_exists(rp.data_folder+"H_N.dat") & file_exists(rp.data_folder+"checkpoint.ckp")){
	    outfile_data.open(rp.data_folder+"data_log.dat", std::ios::app);    
	    outfile_om.open(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary | std::ios::app);    
	    outfile_histogram.open(rp.data_folder+"H_N.dat", std::ios::app); 
	    continue_simulation = true;   
	}
	else{
	    outfile_data.open(rp.data_folder+"data_log.dat");    
	    outfile_om.open(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);    
	    outfile_histogram.open(rp.data_folder+"H_N.dat");  
    	outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<std::endl;
	}
    std::cout<<"Continue? "<<rp.continue_if_possible<<std::endl;


    long t_init=0;
    long t_final=0;//=rp.timesteps;

    //histogram; counts up to 1000 subunits
    std::vector<int> H_N(1000, 0);

    UmbrellaWindow uw(rp.spring_const, mesh.n_faces());

    if (continue_simulation){
    	load_from_checkpoint(mesh, uw, t_final, rp.data_folder+"checkpoint.ckp", input_file);
    }


    while (t_final<rp.timesteps){
    	t_init = t_final;    	
    	t_final = t_init + rp.nsteps_umbrella;
    	//std::cout<<"--- tinit, tfinal "<<t_init<<" "<<t_final<<std::endl;
    	propagate(mesh, t_init, t_final, rp, outfile_om,  outfile_data, uw, eng);    
    	//shift the umbrella window
    	uw.N0+=rp.dN_umbrella;
    	std::cout<<"umbrella N0="<<uw.N0<<std::endl;
    	for (int N: uw.values){
    		H_N[N]++;
    	}
    	for (const auto &hi : H_N) outfile_histogram << hi <<"\t";
    	outfile_histogram << std::endl;

    	uw.values.clear();
    	std::fill(H_N.begin(), H_N.end(), 0);

    	//will start from N0=137. Go up and down 2 rings, dN=52
    	if (uw.N0==(137+78)){
    		rp.dN_umbrella*=-1;
    	}
    	if (uw.N0==(137-78)){
    		rp.dN_umbrella*=-1;
    	}

    }

    outfile_om.close();
    outfile_data.close();
	outfile_histogram.close();    



	/*****checkpoint test*****/
	/*MyMesh mesh_ckp;
	long timestep;
	UmbrellaWindow uwc(0,0);
	load_from_checkpoint(mesh_ckp, uwc, timestep, rp.data_folder+"checkpoint.ckp", input_file);
	std::cout<<"Checkpoint read...."<<std::endl;
	std::cout<<uwc.spring_const<<std::endl;
	std::cout<<uwc.N0<<std::endl;
	std::cout<<timestep<<std::endl;*/
	/**********/
    //do conversions if requested

    //create a conversions directory
    system(("mkdir -p "+rp.data_folder+"conversions").c_str());
    if (rp.convert_to_lammps_trajectory)
		convert_om_to_lammps_trajectory(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/");
    if (rp.convert_to_lammps)
		convert_om_to_lammps_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_");
	if (rp.convert_to_vtk)
		convert_om_to_VTK_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/vtk_", rp.path_to_source+"vtk_colors.conf");    
	if (rp.convert_to_dat)
		convert_om_to_dat_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/dat_");
	//do the WHAM
	//system(("cd "+rp.path_to_source).c_str());
	//system(("cd "+rp.path_to_source + "&& python wham.py "+input_file ).c_str());

	return;	

}


void propagate_external_potential(MyMesh & mesh, long t_init, long t_final, RunParameters rp, std::ostream& outfile_om,  std::ostream& outfile_data, UmbrellaWindow & uw, std::mt19937 & eng, ExternalPotential & externalpotential, AcceptanceMonitor & AM = default_AM){

    int nmoves = 13;  //13;//10;

	/*for (int iac=0; iac<(nmoves+1); iac++){
		n_acc[iac]=0;
		n_tot[iac]=0;
	}    */

    //time limit on propagator only, not the whole run function
	auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;

	bool accepted;
	int which_move, which_type, which_rotation, number_of_rotational_configs;
	double k_insertion, k_fusion, p_propose, r;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> , std::vector<int>> wedges;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > fusion_vectors, type2_fission_vectors, holes;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> >	fission_vectors;	
	std::vector <MyMesh::VertexHandle> w1, v1;  
	std::vector<MyMesh::HalfedgeHandle> boundary_halfedges, h1;	  
	std::vector<MyMesh::FaceHandle> removable_faces_geom, removable_faces;	
	std::tuple<std::vector<MyMesh::HalfedgeHandle>, std::vector<MyMesh::HalfedgeHandle> > fusion_halfedges;	
	std::vector<MyMesh::EdgeHandle> fission_edges;	

	std::vector<MyMesh::VertexHandle*> handle_tracking_v;
		std::vector<MyMesh::FaceHandle*> handle_tracking_f;
	std::vector<MyMesh::HalfedgeHandle*> handle_tracking_h;	

	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto proto_face_props = OpenMesh::FProp<FaceProp>((*mesh_props).prototypes_mesh, "face_props");		
	float kT = (*mesh_props).kT;
	double l_fuse = (*mesh_props).l_fuse;	
	float k_vertex_move = (*mesh_props).k_vertex_move;

	//std::ofstream outfile("outfileom5.om", std::ios::out | std::ios::binary);
	// std::ofstream datafile("data_out5.dat");
	long long fpos;
	long long opos;
	int n_boundary_edges;
	int n_type1, n_type2;
	int ntypes;

	int n_proposals_rejected=0;
	int n_proposals_accepted=0;



//dump_data(mesh, -1, outfile_data, fpos);
	for (long t=t_init; t<t_final; t++){
		if (mesh.n_faces() > rp.Nmax){
            fpos = dump_om(mesh, t, outfile_om);
            dump_data(mesh, t, outfile_data, fpos, AM.r_acc);	
			n_boundary_edges = count_boundary_edges(mesh);
			n_type1 = count_subunits_of_type(mesh, 0);
			n_type2=count_subunits_of_type(mesh, 1);            
			opos = outfile_data.tellp();
  			outfile_data.seekp (opos-1);			
			outfile_data<<"\t"<< n_type1<<"\t"<<n_type2<<"\t" << n_boundary_edges <<std::endl;            
			create_checkpoint(mesh, uw, t, rp.data_folder+"checkpoint.ckp");            		
			break;
		}



		if (t % rp.dtsave ==0){		
			std::cout<<"t = "<<t<<" n_faces="<<mesh.n_faces()<<std::endl;	
			//std::cout<<"t = "<<t<<" n_rej="<<n_proposals_rejected<<" n_acc="<<n_proposals_accepted<<" ratio="<<(float)n_proposals_accepted/(n_proposals_rejected+n_proposals_accepted)<<std::endl;	
			std::cout<<"vertex move acc: t = "<<t<<" n_faces="<<mesh.n_faces()<<" n_tot="<<AM.n_tot[0]<<" n_acc="<<AM.n_acc[0]<<" ratio="<<(float)AM.n_acc[0]/AM.n_tot[0]<<std::endl;	

			n_boundary_edges = count_boundary_edges(mesh);
			n_type1 = count_subunits_of_type(mesh, 0);
			n_type2=count_subunits_of_type(mesh, 1);
			
			//skip full dump for the first dtskip timesteps
			
			if ( (t>rp.dtskip) && (t % rp.dtsave_snapshot==0)){
				fpos = dump_om(mesh, t, outfile_om);
				
			}
			else{
				//save the last snapshot if it will break
				if (rp.stop_at_closure && (n_boundary_edges==0)){
						fpos = dump_om(mesh, t, outfile_om);	
				}
				else{				
					fpos=-1;
				}
			}



			for (int iac=0; iac<(nmoves); iac++){
				AM.r_acc[iac] = AM.n_acc[iac]/AM.n_tot[iac];
				//std::cout<<"r_acc "<< AM.r_acc[iac]<<std::endl;
				//std::cout<<"n_acc "<< AM.n_acc[iac]<<std::endl;
				//std::cout<<"n_tot "<< AM.n_tot[iac]<<std::endl;								
				AM.n_acc[iac]=0;
				AM.n_tot[iac]=0;
			}

			dump_data(mesh, t, outfile_data, fpos, AM.r_acc);
			opos = outfile_data.tellp();
  			outfile_data.seekp (opos-1);			
			outfile_data<<"\t"<< n_type1<<"\t"<<n_type2<<"\t" << n_boundary_edges <<std::endl;
			//dump_om(mesh, "omtest"+std::to_string(t)+".om");


			//don't count boundary edges at every timestep
			if (rp.stop_at_closure && (n_boundary_edges==0)){
					break;				
			}


			//if running for too long, stop
			end = std::chrono::system_clock::now();
			elapsed_seconds = end-start;
	    	if (elapsed_seconds.count()>=rp.max_seconds){
	    		std::cout<<"Time limit reached, cancelling, saving...."<<std::endl;
	    		break;
	    	}			

		}

		if (t % rp.checkpoint_freq ==0){
			create_checkpoint(mesh, uw, t, rp.data_folder+"checkpoint.ckp");			
		}

		
		/*
		for (int n=0; n<(int)mesh.n_vertices(); n++){
			AM.n_acc[0]+=attempt_move(mesh, eng);
			AM.n_tot[0]++;
		}
		*/
			/*std::cout<<"n_tot0 "<< n_tot[1]<<std::endl;
			std::cout<<"n_acc0 "<< n_acc[1]<<std::endl;			
			std::cout<<"n_acc0/n_tot0 "<< n_acc[1]/n_tot[1]<<std::endl;*/		


		//!! IMPORTANT: in each t iteration, only one move should be selected because the garbage collector
		// is only called once, at the end of each iteration. n_faces, n_vertices, etc won't get updated unless
		// the garbage collector is called and we need them to be up to date for the bias potential and also for attempt_move()
		//which_move = rand() % nmoves;
		which_move = randint(nmoves, eng);		
		//std::cout<<"which move "<<which_move<<std::endl;
		accepted=false;
		switch (which_move){
			case 0:
				k_fusion = (*mesh_props).k_fusion; 
				fusion_vectors = get_type1_fusion_triplets(mesh, l_fuse);	
	    		v1 = std::get<0>(fusion_vectors);				
				p_propose = k_fusion * v1.size();
				if ( (p_propose>1) && (!rp.adaptive_rates) ) {
					std::cout<<"Warning! Decrease k_fusion rate! - type1_fusion"<<std::endl;
					dump_error("case 0 ", rp.data_folder+"/error.txt");
				}
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion rate... - type1_fusion"<<std::endl;
					(*mesh_props).k_fusion*=0.5;
					break;//-------- don't do the move if violates detailed balance
				}				
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);
				if (r<p_propose){	
					accepted = attempt_type1_fusion_external_potential(mesh, eng, externalpotential);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type1 fusion"<<std::endl;
				AM.n_tot[which_move+1]++;

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}
				break;

			case 1:
				k_fusion = (*mesh_props).k_fusion; 
				fission_vectors = get_type1_fission_pairs(mesh);	
    			v1 = std::get<0>(fission_vectors);					
				p_propose = k_fusion * v1.size();
				if ( (p_propose>1)  && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_fusion rate! - type1_fission"<<std::endl;
					dump_error("case 1 ", rp.data_folder+"/error.txt");
				}		
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion rate... - type1_fission"<<std::endl;
					(*mesh_props).k_fusion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}							
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	
					accepted = attempt_type1_fission_external_potential(mesh, eng, externalpotential);	
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type1 fission"<<std::endl;
				AM.n_tot[which_move+1]++;

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}				
				break;

			case 2:		
				k_fusion = (*mesh_props).k_fusion2; 
				fusion_vectors = get_type2_fusion_triplets(mesh, l_fuse);	//each counted TWICE
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size()*0.5;  // x 1/2 to correct for double counting
				if ( (p_propose>1)  && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fusion"<<std::endl;
					dump_error("case 2 ", rp.data_folder+"/error.txt");
				}	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion2 rate... - type2_fusion"<<std::endl;
					(*mesh_props).k_fusion2*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}								
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_type2_fusion_external_potential(mesh, eng, externalpotential);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type2 fusion "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;	

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}				
				break;
			case 3:			
				k_fusion = (*mesh_props).k_fusion2; 
				type2_fission_vectors = get_type2_fission_triplets(mesh);	//each counted TWICE
    			v1 = std::get<0>(type2_fission_vectors);					
				p_propose = k_fusion * v1.size()*0.5; // x 1/2 to correct for double counting
				if ( (p_propose>1) && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fission"<<std::endl;
					dump_error("case 3 ", rp.data_folder+"/error.txt");
				}			
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion2 rate... - type2_fission"<<std::endl;
					(*mesh_props).k_fusion2*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}						
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_type2_fission_external_potential(mesh, eng, externalpotential);	
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type2 fission "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}				
				break;

			case 4:	
				//if (mesh.n_faces()>n_umb_sim ) break;

				//pick subunit type from prototype subunits
				//which_type = rand() % (*mesh_props).prototype_faces.size();
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				number_of_rotational_configs = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].number_of_rotational_configs;
				//which_rotation = rand() % number_of_rotational_configs;
				which_rotation = randint( number_of_rotational_configs, eng );		
				ntypes=(*mesh_props).prototype_faces.size();		
				//----
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				boundary_halfedges = find_boundary_halfedges(mesh);			
				p_propose = k_insertion * boundary_halfedges.size() * number_of_rotational_configs * ntypes;
				if ( (p_propose>1) && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_insertion rate! - insertion"<<std::endl;
					dump_error("case 4 ", rp.data_folder+"/error.txt");
				}	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - insertion"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}								
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	
					accepted = attempt_insertion_external_potential(mesh, which_type, which_rotation, uw, eng, externalpotential);	
					AM.n_acc[which_move+1]+=accepted;
				}  
				//if (accepted) std::cout<<"insertion type "<< which_type <<" "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;	

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}							
				break;

			case 5:	
				//if (mesh.n_faces()>n_umb_sim ) break;			
				//pick subunit type from prototype subunits
				//which_type = rand() % (*mesh_props).prototype_faces.size();
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				number_of_rotational_configs = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].number_of_rotational_configs;
				//which_rotation = rand() % number_of_rotational_configs;
				which_rotation = randint( number_of_rotational_configs, eng );
				ntypes=(*mesh_props).prototype_faces.size();	
				//std::cout<<"ntypes "<<ntypes<<std::endl;			
				//---
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				wedges = get_open_wedge_triplets(mesh);		
				w1 = std::get<0>(wedges);		
				p_propose = k_insertion * w1.size() * number_of_rotational_configs * ntypes;
				if ((p_propose>1) && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_insertion rate! - wedge_insertion"<<std::endl;
					dump_error("case 5 ", rp.data_folder+"/error.txt");
				}	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - wedge_insertion"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_wedge_insertion_external_potential(mesh, which_type, which_rotation, uw, eng, externalpotential);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"wedge insertion type "<< which_type <<" "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;		

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}					
				break;
			case 6:
				//pick subunit type from prototype subunits
				//which_type = rand() % (*mesh_props).prototype_faces.size();
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				removable_faces = get_simply_removable_faces(mesh, which_type);
				ntypes=(*mesh_props).prototype_faces.size();

				p_propose = k_insertion * removable_faces.size() * ntypes;
				if ( (p_propose>1) && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_insertion rate! - removal"<<std::endl;	
					dump_error("case 6 ", rp.data_folder+"/error.txt");
				}
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - removal"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	

					accepted = attempt_removal_external_potential(mesh, which_type, uw, eng, externalpotential);
					AM.n_acc[which_move+1]+=accepted;
				}
				AM.n_tot[which_move+1]++;

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}				
				break;
			case 7:
				//pick subunit type from prototype subunits
				//which_type = rand() % (*mesh_props).prototype_faces.size();
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				removable_faces = get_wedge_removable_faces(mesh, which_type); 
				ntypes=(*mesh_props).prototype_faces.size();

				p_propose = k_insertion * removable_faces.size() * ntypes;
				if ( (p_propose>1)  && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_insertion rate! - wedge_removal"<<std::endl;
					dump_error("case 7 ", rp.data_folder+"/error.txt");
				}		
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - wedge_removal"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	
					accepted = attempt_wedge_removal_external_potential(mesh, which_type, uw, eng, externalpotential);
					AM.n_acc[which_move+1]+=accepted;
				}
				AM.n_tot[which_move+1]++;

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}				
				break;

			case 8:
				k_fusion = (*mesh_props).k_fusion_edge;
				fusion_halfedges = get_halfedge_fusion_pairs(mesh, l_fuse); //each counted TWICE
    			h1 = std::get<0>(fusion_halfedges);
    			p_propose = k_fusion * h1.size()*0.5; // x 1/2 to correct for double counting
				if ( (p_propose>1)  && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fusion"<<std::endl;
					dump_error("case 8 ", rp.data_folder+"/error.txt");
				}	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion_edge rate... - edge_fusion"<<std::endl;
					(*mesh_props).k_fusion_edge*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}	
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){					
					accepted = attempt_edge_fusion_external_potential(mesh, eng, externalpotential);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"EDGE fusion " <<" "<<p_propose<<std::endl;
				AM.n_tot[which_move+1]++;

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}				
				break;

			case 9:
				k_fusion = (*mesh_props).k_fusion_edge;
				fission_edges = get_fission_edges(mesh);
    			p_propose = k_fusion * fission_edges.size();	
				if ( (p_propose>1)  && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fission"<<std::endl;	
					dump_error("case 9 ", rp.data_folder+"/error.txt");
				}
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion_edge rate... - edge_fission"<<std::endl;
					(*mesh_props).k_fusion_edge*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}					
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){    									
					accepted = attempt_edge_fission_external_potential(mesh, eng, externalpotential);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"EDGE fission " <<" "<<p_propose<<std::endl;
				AM.n_tot[which_move+1]++;

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}				
				break;

			case 10:
				p_propose = k_vertex_move;
				r = randdouble(eng);				
				if (r<p_propose){ 
					for (int n=0; n<(int)mesh.n_vertices(); n++){ //<--- need this loop. Otherwise vertex moves are too rare, things don't assemble properly
						AM.n_acc[0]+=attempt_move_external_potential(mesh, eng, externalpotential);
						AM.n_tot[0]++;
					}
				}
				/*p_propose = k_vertex_move*mesh.n_vertices();
				if ((p_propose>1) && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_vertex_move rate! - vertex_move"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing vertex move rate"<<std::endl;
					(*mesh_props).k_vertex_move*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				r = randdouble(eng);				
				if (r<p_propose){ 
					AM.n_acc[0]+=attempt_move_external_potential(mesh, eng, externalpotential);
				}
				AM.n_tot[0]++;	
				*/
				/*if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}	*/						

				break;

			case 11:
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				number_of_rotational_configs = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].number_of_rotational_configs;
				//which_rotation = rand() % number_of_rotational_configs;
				which_rotation = randint( number_of_rotational_configs, eng );	
				ntypes=(*mesh_props).prototype_faces.size();			
				//---
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				holes = get_hole_triplets(mesh);		
				w1 = std::get<0>(holes);		
				p_propose = k_insertion * w1.size() * number_of_rotational_configs * 1.0/3.0 * ntypes; //<--- 1/3 because each hole is counted 3x
				if ((p_propose>1) && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_insertion rate! - wedge_insertion"<<std::endl;
					dump_error("case 11 ", rp.data_folder+"/error.txt");
				}	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - hole_insertion"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_hole_insertion_external_potential(mesh, which_type, which_rotation, uw, eng, externalpotential);
					AM.n_acc[which_move+1]+=accepted;
				}
				if (accepted) std::cout<<"hole insertion type "<< which_type <<" "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;	

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}							
				break;	

			case 12:
				which_type = randint( (*mesh_props).prototype_faces.size(), eng );				
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				removable_faces = get_hole_removable_faces(mesh, which_type); 
				ntypes=(*mesh_props).prototype_faces.size();
				p_propose = k_insertion * removable_faces.size() * ntypes;
				if ( (p_propose>1)  && (!rp.adaptive_rates)) {
					std::cout<<"Warning! Decrease k_insertion rate! - hole_removal"<<std::endl;
					dump_error("case 12 ", rp.data_folder+"/error.txt");

				}		
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_insertion rate... - hole_removal"<<std::endl;
					proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}									
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	
					accepted = attempt_hole_removal_external_potential(mesh, which_type, uw, eng, externalpotential);
					AM.n_acc[which_move+1]+=accepted;
				}
				AM.n_tot[which_move+1]++;

				if (r<p_propose){
					n_proposals_accepted++;
				}
				else{
					n_proposals_rejected++;
				}					
				break;				
		}

		//get all face's neighbor list handles to save them through the garbage collector
    	handle_tracking_f.clear();
    	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
    		for (auto & neighbor: face_props[*fit].neighbor_list){
    			handle_tracking_f.push_back( &neighbor);
    		}
    	}
		mesh.garbage_collection<std::vector<MyMesh::VertexHandle*>, std::vector<MyMesh::HalfedgeHandle*>, std::vector<MyMesh::FaceHandle*> >(handle_tracking_v, handle_tracking_h, handle_tracking_f);

		if (t % 150 ==0){
			update_full_neighbor_list(mesh);
			//Could check overlaps for the full capsid here to verify
		}

	//!!	uw.add_value(mesh.n_faces());
	}


	return;
}

void run_dynamical_external_potential(std::string input_file){

	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);  

    //make a copy of the input json to the data folder


    								//radius, hardness 
    //https://stackoverflow.com/questions/36218655/c-declare-derived-class-object-inside-of-if-else-and-use-it-outside
    //need dynamic declaration
    //if (rp.corepotential){
    //	CorePotential externalpotential(5.0, 3.0);
    //}
    //else{
    //	MorsePotential externalpotential(0.35, 1.0, 10.0);
	//}

	ExternalPotential * externalpotential;
	if (rp.externalpotential=="DummyPotential"){
		externalpotential = new DummyPotential();
	}

	if (rp.externalpotential=="MorsePotential"){
		externalpotential = new MorsePotential(rp.extparam1, rp.extparam2, rp.extparam3, rp.extparam4, rp.extparam5, rp.extparam6); //(0.35, 1.0, 10.0);
	}	

    //initialize RNG from ensemble seed
    //std::srand(rp.ensemble*73 + 17);
    std::mt19937 eng(rp.ensemble*73 + 17);

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);

	//shift to the surface of sphere
	/*MyMesh::Point shift = MyMesh::Point(0,0, rp.extparam3);
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		mesh.set_point(*v, mesh.point(*v) + shift );
	}*/

	/*std::vector<double> xinit={-0.5, 0.0, 0.5};
	std::vector<double> yinit={-sqrt(3.0)/6.0, sqrt(3.0)/3.0, -sqrt(3.0)/6.0};
	std::vector<double> zinit={-0.5, 0.0, 0.5};		
	*/


												//height of a tetrahedron + 1/2 cylinder length

	MyMesh::Point shift = MyMesh::Point(0,0, -sqrt(rp.extparam3*rp.extparam3 - 1.0/3.0)-0.5*rp.extparam6);
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		mesh.set_point(*v, mesh.point(*v) + shift );	
	}


    //should I create data_folder here?
    system(("mkdir -p "+rp.data_folder).c_str());    

    std::cout<<"copying... "<<("cp "+input_file+" "+rp.data_folder).c_str()<<std::endl;
    system(("cp "+input_file+" "+rp.data_folder).c_str());     

    std::ofstream outfile_om(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);
    std::ofstream outfile_data(rp.data_folder+"data_log.dat");	
    outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<std::endl;


    //burnin time to equilibrate initial structure
    for (long t=0; t<rp.dt_burnin; t++){
		for (int n=0; n<(int)mesh.n_vertices(); n++){
			attempt_move_external_potential(mesh, eng, *externalpotential);
		}	
	}    

    long t_init=0;
    long t_final=rp.timesteps;
    UmbrellaWindow uw_fake(0.0, 0.0);
    propagate_external_potential(mesh, t_init, t_final, rp, outfile_om,  outfile_data, uw_fake, eng, *externalpotential);
    outfile_om.close();
    outfile_data.close();

    //dump last snapshot and do some analysis
    system(("mkdir -p "+rp.data_folder+"/final_structure").c_str());
    process_final_structure(mesh, rp.data_folder+"/final_structure");

    //do conversions if requested

    //create a conversions directory
    system(("mkdir -p "+rp.data_folder+"conversions").c_str());
    if (rp.convert_to_lammps_trajectory)
		convert_om_to_lammps_trajectory(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/");
    if (rp.convert_to_lammps){
		convert_om_to_lammps_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_");
		convert_om_to_lammps_snapshots_double_edge(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_double_edge_");	
	}
	if (rp.convert_to_vtk)
		convert_om_to_VTK_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/vtk_", rp.path_to_source+"vtk_colors.conf");    
	if (rp.convert_to_dat)
		convert_om_to_dat_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/dat_");
	return;

	delete externalpotential;
}



void run_dynamical_adhesion_potential(std::string input_file){

	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);  

    //make a copy of the input json to the data folder


    								//radius, hardness 
    //https://stackoverflow.com/questions/36218655/c-declare-derived-class-object-inside-of-if-else-and-use-it-outside
    //need dynamic declaration
    //if (rp.corepotential){
    //	CorePotential externalpotential(5.0, 3.0);
    //}
    //else{
    //	MorsePotential externalpotential(0.35, 1.0, 10.0);
	//}



    //initialize RNG from ensemble seed
    //std::srand(rp.ensemble*73 + 17);
    std::mt19937 eng(rp.ensemble*73 + 17);

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);

	zero_com_shift_mesh(mesh);
    MyMesh::Point p;
    double R=0.0;
	for(MyMesh::VertexIter vit = mesh.vertices_sbegin(); vit != mesh.vertices_end(); ++vit) {
		p = mesh.point(*vit);
		R+=p.norm();
	}
	R/=mesh.n_vertices();
	std::cout<<"R="<<R<<std::endl;

	random_rotate_mesh(mesh, eng);



	ExternalPotential * externalpotential;
	if (rp.externalpotential!="WallPotentialWithAdhesion"){
		std::cout<<"WARNING!!! Only use this runner with WallPotentialWithAdhesion!"<<std::endl;
	}

	else{
		externalpotential = new WallPotentialWithAdhesion(rp.extparam1, 10000.1*R, -1.1*R, rp.extparam2, rp.extparam3, rp.extparam4, 0.0); //(0.35, 1.0, 10.0);
	}	

    //should I create data_folder here?
    system(("mkdir -p "+rp.data_folder).c_str());    

    std::cout<<"copying... "<<("cp "+input_file+" "+rp.data_folder).c_str()<<std::endl;
    system(("cp "+input_file+" "+rp.data_folder).c_str());     

    std::ofstream outfile_om(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);
    std::ofstream outfile_data(rp.data_folder+"data_log.dat");	
    outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<std::endl;


    //burnin time to equilibrate initial structure
    for (long t=0; t<rp.dt_burnin; t++){
		for (int n=0; n<(int)mesh.n_vertices(); n++){
			attempt_move_external_potential(mesh, eng, *externalpotential);
		}	
		std::cout<<" Substrate energy "<< (*externalpotential).energy(mesh) <<std::endl;
	}    

    long t_init=0;
    long t_final=rp.timesteps;
    std::cout<<"tinit, tfinal "<<t_init<<" "<<t_final<<std::endl;
    UmbrellaWindow uw_fake(0.0, 0.0);
    //should use a different propagator to forbid insertions/removals
    propagate_external_potential(mesh, t_init, t_final, rp, outfile_om,  outfile_data, uw_fake, eng, *externalpotential);
    outfile_om.close();
    outfile_data.close();

    //dump last snapshot and do some analysis
    system(("mkdir -p "+rp.data_folder+"/final_structure").c_str());
    process_final_structure(mesh, rp.data_folder+"/final_structure");

    //do conversions if requested

    //create a conversions directory
    system(("mkdir -p "+rp.data_folder+"conversions").c_str());
    if (rp.convert_to_lammps_trajectory)
		convert_om_to_lammps_trajectory(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/");
    if (rp.convert_to_lammps){
		convert_om_to_lammps_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_");
		convert_om_to_lammps_snapshots_double_edge(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_double_edge_");	
	}
	if (rp.convert_to_vtk)
		convert_om_to_VTK_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/vtk_", rp.path_to_source+"vtk_colors.conf");    
	if (rp.convert_to_dat)
		convert_om_to_dat_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/dat_");
	return;

	delete externalpotential;
}


void propagate_wall_potential_constant_load(MyMesh & mesh, long t_init, long t_final, RunParameters rp, std::ostream& outfile_om,  std::ostream& outfile_data, std::mt19937 & eng, WallPotential & wall, bool all_moves, AcceptanceMonitor & AM = default_AM){

    int nmoves = 8;
	bool accepted;
	int which_move, which_type, which_rotation, number_of_rotational_configs;
	double k_insertion, k_fusion, p_propose, r;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> , std::vector<int>> wedges;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > fusion_vectors, type2_fission_vectors;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> >	fission_vectors;	
	std::vector <MyMesh::VertexHandle> w1, v1;  
	std::vector<MyMesh::HalfedgeHandle> boundary_halfedges, h1;	  
	std::vector<MyMesh::FaceHandle> removable_faces_geom, removable_faces;	
	std::tuple<std::vector<MyMesh::HalfedgeHandle>, std::vector<MyMesh::HalfedgeHandle> > fusion_halfedges;	
	std::vector<MyMesh::EdgeHandle> fission_edges;	

	std::vector<MyMesh::VertexHandle*> handle_tracking_v;
		std::vector<MyMesh::FaceHandle*> handle_tracking_f;
	std::vector<MyMesh::HalfedgeHandle*> handle_tracking_h;	

	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto proto_face_props = OpenMesh::FProp<FaceProp>((*mesh_props).prototypes_mesh, "face_props");		
	float kT = (*mesh_props).kT;
	double l_fuse = (*mesh_props).l_fuse;	

	//std::ofstream outfile("outfileom5.om", std::ios::out | std::ios::binary);
	// std::ofstream datafile("data_out5.dat");
	long long fpos;
	long long opos;

	std::tuple< double, double > asphericities;
	double E_before, E_after, dE, dz,p;

//dump_data(mesh, -1, outfile_data, fpos);
	for (long t=t_init; t<t_final; t++){
		if (t % rp.dtsave ==0){		
			std::cout<<"t = "<<t<<" n_faces="<<mesh.n_faces()<<std::endl;	
			
			//skip full dump for the first dtskip timesteps
			
			if ( (t>rp.dtskip) && (t % rp.dtsave_snapshot==0)){
					fpos = dump_om(mesh, t, outfile_om);

			}
			else{
				fpos=-1;
			}

			for (int iac=0; iac<=(nmoves); iac++){
				AM.r_acc[iac] = AM.n_acc[iac]/AM.n_tot[iac];
				//std::cout<<"r_acc "<< AM.r_acc[iac]<<std::endl;
				//std::cout<<"n_acc "<< AM.n_acc[iac]<<std::endl;
				//std::cout<<"n_tot "<< AM.n_tot[iac]<<std::endl;								
				AM.n_acc[iac]=0;
				AM.n_tot[iac]=0;
			}

			dump_data(mesh, t, outfile_data, fpos, AM.r_acc);
			opos = outfile_data.tellp();
  			outfile_data.seekp (opos-1);
  			//outfile.write (" sam",4);
			//outfile_data<<"\b";
			outfile_data<<"\t"<<wall.wall_energy_1(mesh)<<"\t"<<wall.wall_energy_2(mesh)<<"\t"<<full_stretch_energy(mesh)<<"\t"<<full_bending_energy(mesh)<<"\t"<<full_binding_energy(mesh)<<"\t"<<wall.wall_force_1(mesh)<<"\t"<<wall.wall_force_2(mesh)<<"\t"<<wall.z_position_1<<"\t"<<wall.z_position_2;
			asphericities = compute_asphericity(mesh);
			outfile_data<<"\t"<< count_boundary_edges(mesh) <<"\t"<< std::get<0>(asphericities)<<"\t"<<std::get<1>(asphericities)<<"\t"<<count_boundary_edges(mesh)<<"\t"<< get_max_height(mesh) <<std::endl;
			
			//dump_om(mesh, "omtest"+std::to_string(t)+".om");
		}

		
		/*for (int n=0; n<(int)mesh.n_vertices(); n++){
			attempt_move_wall(mesh, eng, wall);
		}*/
		

		//!! IMPORTANT: in each t iteration, only one move should be selected because the garbage collector
		// is only called once, at the end of each iteration. n_faces, n_vertices, etc won't get updated unless
		// the garbage collector is called and we need them to be up to date for the bias potential and also for attempt_move()
		//which_move = rand() % nmoves;
		which_move = randint(nmoves, eng);		
		//std::cout<<"which move "<<which_move<<std::endl;
		accepted=false;
		if (!all_moves) which_move = (6+randint(2,eng)); //only vertex moves or wall equilibration
		switch (which_move){
			case 0:
				k_fusion = (*mesh_props).k_fusion; 
				fusion_vectors = get_type1_fusion_triplets(mesh, l_fuse);	
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size();
				if ( (p_propose>1) && (!rp.adaptive_rates) ) std::cout<<"Warning! Decrease k_fusion rate! - type1_fusion"<<std::endl;
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion rate... - type1_fusion"<<std::endl;
					(*mesh_props).k_fusion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}				
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);
				if (r<p_propose){	
					accepted = attempt_type1_fusion_external_potential(mesh, eng, wall);
					AM.n_acc[which_move+1]+=accepted;
				}
				AM.n_tot[which_move+1]++;
				//if (accepted) std::cout<<"type1 fusion"<<std::endl;
				break;

			case 1:
				k_fusion = (*mesh_props).k_fusion; 
				fission_vectors = get_type1_fission_pairs(mesh);	
    			v1 = std::get<0>(fission_vectors);					
				p_propose = k_fusion * v1.size();
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion rate! - type1_fission"<<std::endl;		
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion rate... - type1_fission"<<std::endl;
					(*mesh_props).k_fusion*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}							
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){	
					accepted = attempt_type1_fission_external_potential(mesh, eng, wall);	
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type1 fission"<<std::endl;
				AM.n_tot[which_move+1]++;
				break;

			case 2:		
				k_fusion = (*mesh_props).k_fusion2; 
				fusion_vectors = get_type2_fusion_triplets(mesh, l_fuse);	
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size()*0.5;
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fusion"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion2 rate... - type2_fusion"<<std::endl;
					(*mesh_props).k_fusion2*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}								
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_type2_fusion_external_potential(mesh, eng, wall);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type2 fusion "<<p_propose<<std::endl;		
				AM.n_tot[which_move+1]++;
				break;
			case 3:			
				k_fusion = (*mesh_props).k_fusion2; 
				type2_fission_vectors = get_type2_fission_triplets(mesh);	
    			v1 = std::get<0>(type2_fission_vectors);					
				p_propose = k_fusion * v1.size()*0.5;
				if ( (p_propose>1) && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fission"<<std::endl;			
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion2 rate... - type2_fission"<<std::endl;
					(*mesh_props).k_fusion2*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}						
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){
					accepted = attempt_type2_fission_external_potential(mesh, eng, wall);	
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"type2 fission "<<p_propose<<std::endl;	
				AM.n_tot[which_move+1]++;
				break;



			case 4:
				k_fusion = (*mesh_props).k_fusion_edge;
				fusion_halfedges = get_halfedge_fusion_pairs(mesh, l_fuse);
    			h1 = std::get<0>(fusion_halfedges);
    			p_propose = k_fusion * h1.size()*0.5;
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fusion"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion_edge rate... - edge_fusion"<<std::endl;
					(*mesh_props).k_fusion_edge*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}	
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){					
					accepted = attempt_edge_fusion_external_potential(mesh, eng, wall);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"EDGE fusion " <<" "<<p_propose<<std::endl;
				AM.n_tot[which_move+1]++;
				break;

			case 5:
				k_fusion = (*mesh_props).k_fusion_edge;
				fission_edges = get_fission_edges(mesh);
    			p_propose = k_fusion * fission_edges.size();	
				if ( (p_propose>1)  && (!rp.adaptive_rates)) std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fission"<<std::endl;	
				if ( (p_propose>1) && (rp.adaptive_rates) ) {
					std::cout<<"Decreasing k_fusion_edge rate... - edge_fusion"<<std::endl;
					(*mesh_props).k_fusion_edge*=0.5;
					break;//-------- don't do the move if violates detailed balance					
				}					
				//r = rand()/(RAND_MAX + 1.0);
				r = randdouble(eng);				
				if (r<p_propose){    									
					accepted = attempt_edge_fission_external_potential(mesh, eng, wall);
					AM.n_acc[which_move+1]+=accepted;
				}
				//if (accepted) std::cout<<"EDGE fission " <<" "<<p_propose<<std::endl;
				AM.n_tot[which_move+1]++;
				break;

			case 6:
				for (int n=0; n<(int)mesh.n_vertices(); n++){
					AM.n_acc[0]+=attempt_move_external_potential(mesh, eng, wall);
					AM.n_tot[0]++;
				}
				break;

			case 7:
				//Upper wall move MC here (possibly lower wall too)
				p_propose = wall.move_rate;
				r = randdouble(eng);				
				if (r<p_propose){  
					//attempt moving the wall; expensive operation, has to compute full capsid-wall interaction energy
					//could be optimized by
					//1) cutoff on the potential
					//2) only computing upper/lower wall potentials because only those change
					E_before=wall.energy(mesh);
					dz = randdouble(-wall.dz_max, wall.dz_max, eng);
					wall.move_1(dz);
					wall.move_2(-dz);					
					E_after=wall.energy(mesh)+dz*wall.weight;
					dE = E_after-E_before;
					p = exp(-dE/kT);
					r = randdouble(eng);
					if ( r<p){
						//accepted
						AM.n_acc[which_move+1]++;
						
					}
					else{	
						wall.move_1(-dz);
						wall.move_2(dz);						
					}				
				}
				AM.n_tot[which_move+1]++;
				break;

		}

		//get all face's neighbor list handles to save them through the garbage collector
    	handle_tracking_f.clear();
    	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
    		for (auto & neighbor: face_props[*fit].neighbor_list){
    			handle_tracking_f.push_back( &neighbor);
    		}
    	}
		mesh.garbage_collection<std::vector<MyMesh::VertexHandle*>, std::vector<MyMesh::HalfedgeHandle*>, std::vector<MyMesh::FaceHandle*> >(handle_tracking_v, handle_tracking_h, handle_tracking_f);

		if (t % 150 ==0){
			update_full_neighbor_list(mesh);
			//Could check overlaps for the full capsid here to verify
		}

	//wall.move_1(-dz_wall_per_timestep);
	//wall.move_2(dz_wall_per_timestep);	
	}


	return;
}

void run_wall_potential_constant_load(std::string input_file){

	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);  
    								//radius, hardness 
    //https://stackoverflow.com/questions/36218655/c-declare-derived-class-object-inside-of-if-else-and-use-it-outside
    //need dynamic declaration
    //if (rp.corepotential){
    //	CorePotential externalpotential(5.0, 3.0);
    //}
    //else{
    //	MorsePotential externalpotential(0.35, 1.0, 10.0);
	//}





    //initialize RNG from ensemble seed
    //std::srand(rp.ensemble*73 + 17);
    std::mt19937 eng(rp.ensemble*73 + 17);

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);


    //should I create data_folder here?
    system(("mkdir -p "+rp.data_folder).c_str());    

    std::ofstream outfile_om(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);
    std::ofstream outfile_data(rp.data_folder+"data_log.dat");	
    outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<std::endl;



    std::cout<<"Energy before "<<full_elastic_energy(mesh)<<std::endl;
    zero_com_shift_mesh(mesh);
    //get approximate radius to set the walls:
    double R=0.0;
    MyMesh::Point p;
	for(MyMesh::VertexIter vit = mesh.vertices_sbegin(); vit != mesh.vertices_end(); ++vit) {
		p = mesh.point(*vit);
		R+=p.norm();
	}
	R/=mesh.n_vertices();
	std::cout<<"R="<<R<<std::endl;

	random_rotate_mesh(mesh, eng);
   /* for (int t=0; t<1000; t++){
    	random_rotate_mesh(mesh, eng);
    	dump_om(mesh, t, outfile_om);
	}*/
    std::cout<<"Energy after "<<full_elastic_energy(mesh)<<std::endl;  


	//WallPotential(double hardness, double z_position_1, double z_position_2, double amplitude, double move_rate, double dz_max, double weight)
	WallPotentialWithAdhesion externalpotential(rp.extparam1, 11*R, -11*R, rp.extparam2, rp.extparam3, rp.extparam4, rp.extparam5);


	double wall_move_rate=rp.extparam3;
	rp.extparam3=0.0;
    //burnin time to equilibrate initial structure
    for (long t=0; t<rp.dt_burnin; t++){
		for (int n=0; n<(int)mesh.n_vertices(); n++){
			attempt_move_external_potential(mesh, eng, externalpotential);
		}	
	}    

	//1. push in for dt_burnin timesteps (I don't want to use a separate time parameter for this)
	rp.extparam3=wall_move_rate;
	long t_init=0;
    long t_final=rp.dt_burnin;
    rp.dt_burnin=0; //RESET to save during compression too
	bool crack_allowed=false;

	propagate_wall_potential_constant_load(mesh, t_init, t_final, rp, outfile_om,  outfile_data, eng, externalpotential, crack_allowed);
	 

	//2. hold
	t_init=t_final;
	long dwell_time=rp.extparam6;
    t_final=t_init+dwell_time;
	crack_allowed=true;	
	propagate_wall_potential_constant_load(mesh, t_init, t_final, rp, outfile_om,  outfile_data, eng, externalpotential, crack_allowed);

	//3. remove wall
	t_init=t_final;
	t_final=t_init+rp.timesteps;
	externalpotential.weight=0;
	externalpotential.move_rate=0;
	externalpotential.z_position_1 = 100.0*R;
	//externalpotential.z_position_2 = -100.0*R;
	propagate_wall_potential_constant_load(mesh, t_init, t_final, rp, outfile_om,  outfile_data, eng, externalpotential, crack_allowed);


    //bool crack_allowed=false;
    //propagate_wall_potential_constant_load(mesh, t_init, t_final, rp, outfile_om,  outfile_data, eng, externalpotential, crack_allowed);
    outfile_om.close();
    outfile_data.close();


    //do conversions if requested

    //create a conversions directory
    system(("mkdir -p "+rp.data_folder+"conversions").c_str());
    if (rp.convert_to_lammps_trajectory)
		convert_om_to_lammps_trajectory(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/");
    if (rp.convert_to_lammps)
		convert_om_to_lammps_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_");
	if (rp.convert_to_vtk)
		convert_om_to_VTK_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/vtk_", rp.path_to_source+"vtk_colors.conf");    
	if (rp.convert_to_dat)
		convert_om_to_dat_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/dat_");
	return;

}






//Warning: does not update any neighbor list to make it faster
//if overlaps become a problem, this needs to be rethought and properly separated in a propagate - run pair of functions as all the others
void run_thermodynamic_integration(std::string input_file){

	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);  

    //these will be read from input file
    double k_einstein=rp.k_einstein;



    //initialize RNG from ensemble seed
    //std::srand(rp.ensemble*73 + 17);
    std::mt19937 eng(rp.ensemble*73 + 17);

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);


	//
	//fix the Einstein spring centers
	auto vertex_props = OpenMesh::VProp<VertexProp>(mesh, "vertex_props");
    for(MyMesh::VertexIter it = mesh.vertices_sbegin(); it != mesh.vertices_end(); ++it) {
    	vertex_props[*it].r0=mesh.point(*it);
    	vertex_props[*it].vertex_type=0;
    }
    // ============ define constrained points; ideally would find the largest area triangle. But that's O(n_vertices^3)
    //v_fix1 --> random
    //v_fix2 --> farthest from v_fix1
    //v_fix3 --> farthest from both
    int v1r = randint( mesh.n_vertices(), eng);	
	MyMesh::VertexHandle v_fix1 = mesh.vertex_handle(v1r);
	MyMesh::VertexHandle v_fix2, v_fix3;
	double dist1=0.0;
	double dist2=0.0;
	double dist;
	v_fix2=mesh.vertex_handle(0);
    for(MyMesh::VertexIter it = mesh.vertices_sbegin(); it != mesh.vertices_end(); ++it) {
    	dist = (mesh.point(*it) - mesh.point(v_fix1)).norm();
    	if (dist>dist1){
    		dist1=dist;
    		v_fix2=(*it);
    	}
    	
    }	
    v_fix3=mesh.vertex_handle(0);
    for(MyMesh::VertexIter it = mesh.vertices_sbegin(); it != mesh.vertices_end(); ++it) {
    	dist = (mesh.point(*it) - mesh.point(v_fix1)).norm() + (mesh.point(*it) - mesh.point(v_fix2)).norm();
    	if (dist>dist2){
    		dist2=dist;
    		v_fix3=(*it);
    	}
    	
    }	    

    MyMesh::Point e1=(mesh.point(v_fix1) - mesh.point(v_fix2)).normalize();
    MyMesh::Point r13=(mesh.point(v_fix3)-mesh.point(v_fix1));
    MyMesh::Point e2 = (  r13 - r13.dot(e1)*e1   ).normalize();

    //for lammps visualization 
    vertex_props[v_fix1].vertex_type=1;
    vertex_props[v_fix2].vertex_type=2;
    vertex_props[v_fix3].vertex_type=3;

    std::cout<<"restricted vertices: "<<v_fix1<<" "<<v_fix2<<" "<<v_fix3<<std::endl;
    //===================================================================

    //should I create data_folder here?
    system(("mkdir -p "+rp.data_folder).c_str());    

    std::ofstream outfile_om(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);
    std::ofstream outfile_data(rp.data_folder+"data_log.dat");	
    outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<std::endl;


    std::ofstream outfile_thermodynamic_integration(rp.data_folder+"thermodynamic_integration.dat");	



    double lambda_ti;
    double dlambda_ti=1e-2;
    long dtsample=100;
    long nsamples_per_lambda=5000;
    double U0, U1;
    double dU_avg;

    //burnin time to equilibrate initial structure
    lambda_ti=0.0;


//================== k sweep for test ================== 
/*
	for (int ik=1; ik<=10; ik++){
		k_einstein=exp(ik);
	    for (long ilambda=0; ilambda<=(long)(1.01/dlambda_ti); ilambda++){
	    	lambda_ti=ilambda*dlambda_ti;

	    	//equilibrate after adjusting lambda_ti
			for (long t=0; t<rp.dt_burnin; t++){
				for (int n=0; n<(int)mesh.n_vertices(); n++){
					attempt_move_thermodynamic_integration(mesh, eng, v_fix1, v_fix2, v_fix3, e1, e2, k_einstein, lambda_ti);
				}
			}


			//start sampling
			dU_avg=0.0;
			for (long t=0; t<dtsample*nsamples_per_lambda; t++){
				for (int n=0; n<(int)mesh.n_vertices(); n++){
					attempt_move_thermodynamic_integration(mesh, eng, v_fix1, v_fix2, v_fix3, e1, e2, k_einstein, lambda_ti);
				}

				if (t % dtsample==0){
					//sample
					U0=einstein_solid_energy(mesh, k_einstein);
					U1=full_energy(mesh);		
					dU_avg+=(U0-U1)/nsamples_per_lambda;		
				}	
			}
			outfile_thermodynamic_integration<<lambda_ti<<"\t"<<dU_avg<<"\t"<<k_einstein<<std::endl;

		}
	}
*/

//============================================    

    double n_acc, n_tot;

    for (long ilambda=0; ilambda<=(long)(1.01/dlambda_ti); ilambda++){
    	lambda_ti=ilambda*dlambda_ti;

    	//equilibrate after adjusting lambda_ti
		for (long t=0; t<rp.dt_burnin; t++){
			for (int n=0; n<(int)mesh.n_vertices(); n++){
				attempt_move_thermodynamic_integration(mesh, eng, v_fix1, v_fix2, v_fix3, e1, e2, k_einstein, lambda_ti);
			}
		}


		//start sampling
		dU_avg=0.0;
		n_acc=0.0;
		n_tot=0.0;
		for (long t=0; t<dtsample*nsamples_per_lambda; t++){
			for (int n=0; n<(int)mesh.n_vertices(); n++){
				n_acc+=attempt_move_thermodynamic_integration(mesh, eng, v_fix1, v_fix2, v_fix3, e1, e2, k_einstein, lambda_ti);
				n_tot++;
			}

			if (t % dtsample==0){
				//sample
				U0=einstein_solid_energy(mesh, k_einstein);
				U1=full_energy(mesh);		
				dU_avg+=(U0-U1)/nsamples_per_lambda;		
			}	

			if ( (t>rp.dtskip) && (t % rp.dtsave_snapshot==0)){
				dump_om(mesh, ilambda*dtsample*nsamples_per_lambda + t, outfile_om);
				std::cout<<t/(double)rp.dt_burnin<<std::endl;
			}
		}
		outfile_thermodynamic_integration<<lambda_ti<<"\t"<<dU_avg<<"\t"<<k_einstein<<"\t"<<n_acc/n_tot<<std::endl;

	}


/*    for (long t=0; t<rp.dt_burnin; t++){
		for (int n=0; n<(int)mesh.n_vertices(); n++){
			attempt_move_thermodynamic_integration(mesh, eng, v_fix1, v_fix2, v_fix3, e1, e2, k_einstein, lambda_ti);
		}	

		if ( (t>rp.dtskip) && (t % rp.dtsave_snapshot==0)){
			dump_om(mesh, t, outfile_om);
			std::cout<<t/(double)rp.dt_burnin<<std::endl;
		}


		//DO NOT HARDCODE. USE a rp.dt_sample parameter instead 
		if (t%100==0){
			U0=einstein_solid_energy(mesh, k_einstein);
			U1=full_energy(mesh);

			std::cout<<"U0, U1 "<<U0<<" "<<U1<<std::endl;
		}

	} 

*/   

    long t_init=0;
    long t_final=rp.timesteps;





    outfile_om.close();
    outfile_data.close();
	outfile_thermodynamic_integration.close();

    //do conversions if requested

    //create a conversions directory
    system(("mkdir -p "+rp.data_folder+"conversions").c_str());
    if (rp.convert_to_lammps_trajectory)
		convert_om_to_lammps_trajectory(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/");
    if (rp.convert_to_lammps)
		convert_om_to_lammps_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/lammps_");
	if (rp.convert_to_vtk)
		convert_om_to_VTK_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/vtk_", rp.path_to_source+"vtk_colors.conf");    
	if (rp.convert_to_dat)
		convert_om_to_dat_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/dat_");
	return;

}






void run_debug(){

///---	make_capsid_SAT_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/SAT/bmn_colorings/bmn_coloring_0.dat" ,"/home/btyukodi/SAT/bmn_colorings/0/");
///---	make_capsid_SAT_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/SAT/bmn_colorings/bmn_coloring_1.dat" ,"/home/btyukodi/SAT/bmn_colorings/1/");
///---	make_capsid_SAT_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/SAT/bmn_colorings/bmn_coloring_2.dat" ,"/home/btyukodi/SAT/bmn_colorings/2/");
///---	make_capsid_SAT_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/SAT/bmn_colorings/bmn_coloring_3.dat" ,"/home/btyukodi/SAT/bmn_colorings/3/");
///---	make_capsid_SAT_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/SAT/bmn_colorings/bmn_coloring_4.dat" ,"/home/btyukodi/SAT/bmn_colorings/4/");
///---	make_capsid_SAT_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/SAT/bmn_colorings/bmn_coloring_5.dat" ,"/home/btyukodi/SAT/bmn_colorings/5/");
///---	make_capsid_SAT_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/SAT/bmn_colorings/bmn_coloring_6.dat" ,"/home/btyukodi/SAT/bmn_colorings/6/");
///---	make_capsid_SAT_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/SAT/bmn_colorings/bmn_coloring_7.dat" ,"/home/btyukodi/SAT/bmn_colorings/7/");
///---	make_capsid_SAT_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/SAT/bmn_colorings/bmn_coloring_8.dat" ,"/home/btyukodi/SAT/bmn_colorings/8/");
///---	make_capsid_SAT_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/SAT/bmn_colorings/bmn_coloring_9.dat" ,"/home/btyukodi/SAT/bmn_colorings/9/");
///---
///---return;

///////---	print_mesh_SAT_input("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT1/icoT1.om", "/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/SAT_inputs/T1/");
///////---	print_mesh_SAT_input("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT3/icoT3.om", "/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/SAT_inputs/T3/");
///////---	print_mesh_SAT_input("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT4/icoT4.om", "/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/SAT_inputs/T4/");
///////---	print_mesh_SAT_input("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT7/icoT7.om", "/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/SAT_inputs/T7/");
///////---	print_mesh_SAT_input("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT9/icoT9.om", "/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/SAT_inputs/T9/");
///////---
///////---	return;

//	run_wall_potential_constant_load("/home/btyukodi/testdatafolder/wall_potential_with_adhesion/capsid_input_wall_potential.json");
//	return;

/*	count_disclination_spots("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/532/T7/sample.json");
	count_disclination_spots("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/332/T7/sample.json");
	count_disclination_spots("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/222/T7/sample.json");
	count_disclination_spots("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/33/T7/sample.json");
	count_disclination_spots("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/0/T7/sample.json");


	count_disclination_spots("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/532/T9/sample.json");
	count_disclination_spots("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/332/T9/sample.json");
	count_disclination_spots("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/222/T9/sample.json");
	count_disclination_spots("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/33/T9/sample.json");
	count_disclination_spots("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/0/T9/sample.json");
*/

/***************
	std::string sfile="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/Daichi/";
	std::vector<int> T0 = {1, 4, 7, 9, 12, 13, 16, 19, 21, 25};//, 27, 28, 31, 36, 37, 39, 43, 48, 49, 52};
	std::vector<int> daichi_sym={532, 522, 332, 322, 55,222,33,22,0};
	for (auto t: T0){
		for (auto daichi: daichi_sym){
			count_disclination_spots(sfile+std::to_string(daichi)+"/T"+std::to_string(t)+"/sample.json");
		}
	}
	return;
********************/


//////	std::string root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/";
//////	std::string source_root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/";
//////	std::string design_file_root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/";
//////
//////	std::vector<int> T2 = {1, 4, 7, 9, 12, 13, 16, 19, 21, 25, 27, 28, 31, 36, 37, 39, 43, 48, 49, 52};//{1, 4, 9, 16};
//////	std::vector<int> daichi_symm={532, 522, 332, 322, 55,222,33,22,0};
//////	for (auto t: T2){
//////		for (auto daichi: daichi_symm){
//////			make_capsid_Sigl_specific(source_root+"/icoT"+std::to_string(t)+"/icoT"+std::to_string(t)+".om", source_root+"/Sigl/T"+std::to_string(t)+"/");
//////			make_capsid_daichi_specific(source_root+"/Sigl/T"+std::to_string(t)+"/Sigl_specific.om", design_file_root+std::to_string(daichi)+".dat", root+"Daichi/"+std::to_string(daichi)+"/T"+std::to_string(t)+"/");
//////			//make_capsid_daichi_specific(source_root+"/T"+std::to_string(t)+"/ico_greg.om", design_file_root+std::to_string(daichi)+".dat", root+std::to_string(daichi)+"/T"+std::to_string(t)+"/");
//////
//////		}
//////	}
//////	return;


///	std::string root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/";
///	std::string source_root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/";
///
///
///	std::vector<int> T2 = {7,12, 13, 16, 19, 21, 36};//{1, 4, 7, 9, 12, 13, 16, 19, 21, 25, 27, 28, 31, 36, 37, 39, 43, 48, 49, 52};//{1, 4, 9, 16};
///	for (auto t: T2){
///			make_capsid_Sigl_specific(source_root+"/icoT"+std::to_string(t)+"/icoT"+std::to_string(t)+".om", source_root+"/Sigl/T"+std::to_string(t)+"/");
///			make_capsid_botond_specific(source_root+"/Sigl/T"+std::to_string(t)+"/Sigl_specific.om",  root+"Botond/T"+std::to_string(t)+"/");
///			//make_capsid_daichi_specific(source_root+"/T"+std::to_string(t)+"/ico_greg.om", design_file_root+std::to_string(daichi)+".dat", root+std::to_string(daichi)+"/T"+std::to_string(t)+"/");
///
///	}
///	return;


/************************
	std::string root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Botond_design/equal_dihedrals/";
	std::string source_root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/greg_icosahedra/equal_dihedrals/";

	std::vector<int> T2 = {1, 4, 7, 9, 12, 13, 16};//{1, 4, 9, 16};
	for (auto t: T2){
			make_capsid_Sigl_specific(source_root+"/T"+std::to_string(t)+"/ico_greg.om", source_root+"/Sigl/T"+std::to_string(t)+"/");
			make_capsid_botond_specific(source_root+"/Sigl/T"+std::to_string(t)+"/Sigl_specific.om", root+"/T"+std::to_string(t)+"/");
			//make_capsid_daichi_specific(source_root+"/T"+std::to_string(t)+"/ico_greg.om", design_file_root+std::to_string(daichi)+".dat", root+std::to_string(daichi)+"/T"+std::to_string(t)+"/");
	}

	return;
************************/

///
//run_dynamical_adhesion_potential("/home/btyukodi/testdatafolder/adhesion/adhesiontest.json");
//compute_curvature("/home/btyukodi/testdatafolder/adhesion/adhesion/data/data/snapshots.om", "/home/btyukodi/testdatafolder/adhesion/adhesion/data/data/curvature/");

//color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT9/icoT9.om", "/home/btyukodi/testdatafolder/facet_numbering/T9/");
//	color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT3/icoT3.om", "/home/btyukodi/testdatafolder/facet_numbering/T3/");
//
//	color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT4/icoT4.om", "/home/btyukodi/testdatafolder/facet_numbering/T4/");
//
//	color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT7/icoT7.om", "/home/btyukodi/testdatafolder/facet_numbering/T7/");


/*	color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT4/icoT4.om", "/home/btyukodi/testdatafolder/facet_numbering/T4/");


	color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT13/icoT13.om", "/home/btyukodi/testdatafolder/facet_numbering/T13/");

	color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT7/icoT7.om", "/home/btyukodi/testdatafolder/facet_numbering/T7/");

	color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT9/icoT9.om", "/home/btyukodi/testdatafolder/facet_numbering/T9/");


	color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT48/icoT48.om", "/home/btyukodi/testdatafolder/facet_numbering/T48/");

	color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT16/icoT16.om", "/home/btyukodi/testdatafolder/facet_numbering/T16/");
*/
//		color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT16/icoT16.om", "/home/btyukodi/testdatafolder/facet_numbering/T16/");
//
//		color_facets("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT12/icoT12.om", "/home/btyukodi/testdatafolder/facet_numbering/T12/");
//
//make_capsid_Sigl_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT48/icoT48.om", "/home/btyukodi/testdatafolder/specific_coloring/T48/");



//int ncalls=0;
//MyMesh mesh;
//Rotate(5, /*mesh,*/ ncalls);
//
//return;

//	make_capsid_Sigl_specific("/home/btyukodi/testdatafolder/specific_coloring/icoT7.om", "/home/btyukodi/testdatafolder/specific_coloring/T7");
//	make_capsid_Sigl_specific("/home/btyukodi/testdatafolder/specific_coloring/icoT4.om", "/home/btyukodi/testdatafolder/specific_coloring/T4");
//	make_capsid_Sigl_specific("/home/btyukodi/testdatafolder/specific_coloring/icoT9.om", "/home/btyukodi/testdatafolder/specific_coloring/T9");
//	make_capsid_Sigl_specific("/home/btyukodi/testdatafolder/specific_coloring/icoT13.om", "/home/btyukodi/testdatafolder/specific_coloring/T13");
//	make_capsid_Sigl_specific("/home/btyukodi/testdatafolder/specific_coloring/icoT19.om", "/home/btyukodi/testdatafolder/specific_coloring/T19");	
//	make_capsid_Sigl_specific("/home/btyukodi/testdatafolder/specific_coloring/icoT19.om", "/home/btyukodi/testdatafolder/specific_coloring/T19");		
//	make_capsid_Sigl_specific("/home/btyukodi/testdatafolder/specific_coloring/icoT3.om", "/home/btyukodi/testdatafolder/specific_coloring/T3");	


//    push_vertices_to_sphere("/home/btyukodi/testdatafolder/specific_coloring/faceting/data/last_snapshots.om", "whatever");


/******************************************************************************************
	std::vector<int> T = {1,3,4,7,9,12, 13,16};

	std::string root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/greg_icosahedra/equal_dihedrals/";
	std::string source_root;
	for (auto t: T){
		//std::cout<<"/home/btyukodi/greg/equal_dihedrals/T"+std::to_string(t)+"/ico_greg.om"<<std::endl;

		//print_capsid_data(root+"/T"+std::to_string(t)+"/ico_greg.om", "/home/btyukodi/greg/equal_dihedrals/T"+std::to_string(t)+"/");
		//make_capsid_Sigl_specific(root+"/T"+std::to_string(t)+"/ico_greg.om", root+"/Sigl/T"+std::to_string(t)+"/");
		//make_capsid_full_specific(root+"/T"+std::to_string(t)+"/ico_greg.om", root+"/full/T"+std::to_string(t)+"/");
		create_simple_seeds(root+"/T"+std::to_string(t)+"/ico_greg.om", "/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/sub_CK_merged_design/equal_dihedrals/Ntypes0/T"+std::to_string(t)+"/");


	}
return;
/*	root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/greg_icosahedra/equal_lengths/";
	for (auto t: T){
		make_capsid_Sigl_specific(root+"/T"+std::to_string(t)+"/ico_greg.om", root+"/Sigl/T"+std::to_string(t)+"/");
		make_capsid_full_specific(root+"/T"+std::to_string(t)+"/ico_greg.om", root+"/full/T"+std::to_string(t)+"/");

	}
*/	
/****************************************************************************************** */

/*
	convert_om_to_VTK_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/ugly/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/ugly/vtk_","/home/btyukodi/testdatafolder/carlos/colorListT7_rgb.txt");
    convert_om_to_VTK_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/good/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/good/vtk_","/home/btyukodi/testdatafolder/carlos/colorListT7_rgb.txt");
    convert_om_to_VTK_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/bad/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/bad/vtk_","/home/btyukodi/testdatafolder/carlos/colorListT7_rgb.txt");
    convert_om_to_VTK_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G3/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G3/vtk_","/home/btyukodi/testdatafolder/carlos/colorList_G_T3_rgb.txt");
    convert_om_to_VTK_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/P3/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/P3/vtk_","/home/btyukodi/testdatafolder/carlos/colorList_P_T3_rgb.txt");
    convert_om_to_VTK_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/D3/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/D3/vtk_","/home/btyukodi/testdatafolder/carlos/colorList_D_T3_rgb.txt");



	convert_om_to_lammps_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/ugly/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/ugly/lammps_");
    convert_om_to_lammps_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/good/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/good/lammps_");
    convert_om_to_lammps_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/bad/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G7/bad/lammps_");
    convert_om_to_lammps_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G3/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/G3/lammps_");
    convert_om_to_lammps_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/P3/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/P3/lammps_");
    convert_om_to_lammps_snapshots("/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/D3/snapshots.om","/home/btyukodi/testdatafolder/carlos/omfiles/omfiles/D3/lammps_");
*/

/***************************************/
	std::string root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/";
	std::string source_root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/greg_icosahedra/equal_dihedrals/";
	std::string design_file_root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/";

	std::vector<int> T2 = {1, 4, 7, 9, 12, 13, 16};//{1, 4, 9, 16};
	std::vector<int> daichi_symm={532, 522, 332, 322, 55,222,33,22,0};
	for (auto t: T2){
		for (auto daichi: daichi_symm){
			make_capsid_Sigl_specific(source_root+"/T"+std::to_string(t)+"/ico_greg.om", source_root+"/Sigl/T"+std::to_string(t)+"/");
			make_capsid_daichi_specific(source_root+"/Sigl/T"+std::to_string(t)+"/Sigl_specific.om", design_file_root+std::to_string(daichi)+".dat", root+std::to_string(daichi)+"/T"+std::to_string(t)+"/");
			//make_capsid_daichi_specific(source_root+"/T"+std::to_string(t)+"/ico_greg.om", design_file_root+std::to_string(daichi)+".dat", root+std::to_string(daichi)+"/T"+std::to_string(t)+"/");

		}
	}
return;
/*******************************////////////////////

/*	root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_lengths/";
	source_root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/greg_icosahedra/equal_lengths/";
	for (auto t: T2){
		for (auto daichi: daichi_symm){
			make_capsid_daichi_specific(source_root+"/Sigl/T"+std::to_string(t)+"/Sigl_specific.om", design_file_root+std::to_string(daichi)+".dat", root+std::to_string(daichi)+"/T"+std::to_string(t)+"/");
		}
	}	
/**********************************************/


//@@@@@	root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/sub_CK_merged_design/equal_dihedrals/";
//@@@@@	source_root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/equal_dihedrals/";
//@@@@@
//@@@@@
//@@@@@	std::vector<int> T = {1, 4, 7, 9, 12, 13, 16};//{4,9,16};
//@@@@@	int n_ck_types;
//@@@@@	for (auto t: T){
//@@@@@
//@@@@@		if (t%3==0){
//@@@@@			n_ck_types=t/3;
//@@@@@		}
//@@@@@		else{
//@@@@@			n_ck_types=ceil(t/3.0);
//@@@@@		}
//@@@@@		std::cout<<"T "<<t<<" nCK"<< n_ck_types<<std::endl;
//@@@@@		for (int k=1; k<n_ck_types; k++){
//@@@@@			//make_capsid_merged_specific(source_root+"/Sigl/T"+std::to_string(t)+"/Sigl_specific.om", root+"/Ntypes"+std::to_string(k)+"/T"+std::to_string(t)+"/", n_ck_types, k);
//@@@@@			make_capsid_merged_specific(source_root+"/532/T"+std::to_string(t)+"/full_capsid.om", root+"/Ntypes"+std::to_string(k)+"/T"+std::to_string(t)+"/", n_ck_types, k);
//@@@@@
//@@@@@		}
//@@@@@	}
//@@@@@
//@@@@@
//@@@@@	return;

//std::string design_file_root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/Daichi_design/";

//make_capsid_daichi_specific("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/greg_icosahedra/equal_dihedrals/Sigl/T16/Sigl_specific.om", design_file_root+"532.dat", "/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/");


	//std::vector<int> T = {1,3,4,7,9,13,16};
	//int t=1;
	//root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/greg_icosahedra/equal_dihedrals/";
	
	//std::vector<int> merge_pattern={4, 1, 1, 1, 1, 0, 0};//{4,0,0,4,0,0}; //{1,2,3,4,0,0}; //{0, 0, 1, 2, 3, 4}; //
	//std::vector<std::vector<int> > merge_pattern{{0,0}} ;//{ {1,2}};//,{1,3 } , {0,4} };

		//std::cout<<"/home/btyukodi/greg/equal_dihedrals/T"+std::to_string(t)+"/ico_greg.om"<<std::endl;

		//print_capsid_data(root+"/T"+std::to_string(t)+"/ico_greg.om", "/home/btyukodi/greg/equal_dihedrals/T"+std::to_string(t)+"/");
		//make_capsid_Sigl_specific_merge("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/icoT25/icoT25.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);
		
	//std::vector<std::vector<int> > merge_pattern;	
/*	std::vector<std::vector<std::vector<int> > > all_patterns;

	int n=3;
	int k=2;
	all_patterns = get_partitions( n, k); 


    for (auto subset: all_patterns){
        std::cout<<"-------------"<<std::endl;
        for (auto group: subset){
            for (auto element: group){
                std::cout<<" "<<element;
            }
            std::cout<<std::endl;
        }
        
    }

	//for (auto merge_pattern: all_patterns){
		std::vector<std::vector<int> > merge_pattern=all_patterns[0];
		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/Ttest/", merge_pattern);

	//}
*/
/*
		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);


		merge_pattern = { {1,2,4}, {0,3}};
		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);

		merge_pattern = { {1,3}, {2,4}};

		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);

		merge_pattern = { {0,3}, {0,4}};
		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);

		merge_pattern = { {2,3}, {0,4}};
		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);

	
		merge_pattern = { {2,3}, {1,4}};
		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);


		merge_pattern = { {2,4}, {1,2}};
		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);



		merge_pattern = { {2,4}};
		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);



		merge_pattern = { {1,2}};
		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);



		merge_pattern = { {1,3}};
		make_capsid_merged_specific("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om", "/home/btyukodi/testdatafolder/specific_coloring_merge/", merge_pattern);
*/


//	std::cout<<"COPYMESH TEST"<<std::endl;
//	copy_mesh_test("/home/btyukodi/testdatafolder/specific_coloring_merge/T25_daichi/full_capsid.om");
	/*root="/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/greg_icosahedra/equal_lengths/";
	for (auto t: T){
		make_capsid_Sigl_specific(root+"/T"+std::to_string(t)+"/ico_greg.om", root+"/Sigl/T"+std::to_string(t)+"/");
		make_capsid_full_specific(root+"/T"+std::to_string(t)+"/ico_greg.om", root+"/full/T"+std::to_string(t)+"/");

	}*/
/*	std::vector<int> T = {1,3,4,7,9,13,16};
	for (auto t: T){
		std::cout<<"/home/btyukodi/greg/equal_dihedrals/T"+std::to_string(t)+"/ico_greg.om"<<std::endl;
		print_capsid_data("/home/btyukodi/greg/equal_dihedrals/T"+std::to_string(t)+"/ico_greg.om", "/home/btyukodi/greg/equal_dihedrals/T"+std::to_string(t)+"/");

	}

	for (auto t: T){
		print_capsid_data("/home/btyukodi/greg/equal_lengths/T"+std::to_string(t)+"/ico_greg.om", "/home/btyukodi/greg/equal_lengths/T"+std::to_string(t)+"/");

	}	
*/
    //print_capsid_data("/home/btyukodi/testdatafolder/specific_coloring/Greg_icosahedra/greg.om", "whatever");


//	std::cout<<"------------------------ T7 -------------------------------"<<std::endl;
//	print_capsid_data("/home/btyukodi/testdatafolder/specific_coloring/faceted_icoT7.om", "/home/btyukodi/testdatafolder/specific_coloring/T7");

/*	std::cout<<"------------------------ T4 -------------------------------"<<std::endl;	
	print_capsid_data("/home/btyukodi/testdatafolder/specific_coloring/faceted_icoT4.om", "/home/btyukodi/testdatafolder/specific_coloring/T4");

	std::cout<<"------------------------ T9 -------------------------------"<<std::endl;	
	print_capsid_data("/home/btyukodi/testdatafolder/specific_coloring/faceted_icoT9.om", "/home/btyukodi/testdatafolder/specific_coloring/T9");

	std::cout<<"------------------------ T13 -------------------------------"<<std::endl;	
	print_capsid_data("/home/btyukodi/testdatafolder/specific_coloring/faceted_icoT13.om", "/home/btyukodi/testdatafolder/specific_coloring/T13");

	std::cout<<"------------------------ T19 -------------------------------"<<std::endl;	
	print_capsid_data("/home/btyukodi/testdatafolder/specific_coloring/faceted_icoT19.om", "/home/btyukodi/testdatafolder/specific_coloring/T19");		
*/
//	color_faces_Sigl_specific_bak("/home/btyukodi/testdatafolder/specific_coloring/icoT19.om", "/home/btyukodi/testdatafolder/specific_coloring/T19");



	//run_thermodynamic_integration("/home/btyukodi/testdatafolder/thermodynamic_integration/sweep/10.json");
	//make_capsid_full_specific("/home/btyukodi/testdatafolder/full_specific/capsid/icoT9.om", "/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/init_structures/icosahedra/T9_full_specific/");
	//run_umbrella_cycling_referee("/home/btyukodi/testdatafolder/referee_cycling/cycling_60d22ba2946df6e5801a308a.json");
	//run_dynamical_external_potential("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/capsid_input_external_potential.json");
	//run_wall_potential_constant_load("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/capsid_input_wall_potential.json");
	//run_dynamical("/home/btyukodi/testdatafolder/core/capsid_input.json");
	//create_single_subunit_om("/home/btyukodi/testdatafolder/init/single_subunit_plane_triangle_non_removable_fixed.om", -2, 0, 0, 0);

	//create_hexamer_subunit_om("/home/btyukodi/testdatafolder/init/hexamer_subunit_plane_triangle_non_removable_fixed.om", -2, 0, 0, 0);
	//convert_om_to_dat_snapshots("/home/btyukodi/testdatafolder/carlos4/snapshots.om", "/home/btyukodi/testdatafolder/carlos4/conversions/", 33165108);

/*	UmbrellaWindow uw_fake(0.0, 0.0);
	std::string input_file = "/home/btyukodi/testdatafolder/spherocylinder/hole_insertion_debug/23.json";
	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);   

    //initialize RNG from ensemble seed
    //std::srand(rp.ensemble*73 + 17);
    std::mt19937 eng(rp.ensemble*73 + 17);

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);

 	dump_off(mesh, "/home/btyukodi/testdatafolder/spherocylinder/hole_insertion_debug/dump.off");

	std::vector<MyMesh::VertexHandle> v1h, v2h, v3h;	
    auto hole_vectors = get_hole_triplets(mesh);
    v1h = std::get<0>(hole_vectors);
    v2h = std::get<1>(hole_vectors);
    v3h = std::get<2>(hole_vectors);   

    std::cout<<"============== Hole triplets ===================="<<std::endl;
    for (int i=0; i<v1h.size(); i++){
    	std::cout<<"hole #"<<i<<" "<<v1h[i]<<" "<<v2h[i]<<" "<<v3h[i]<<std::endl;

    }
    std::cout<<"================================================="<<std::endl;


	std::vector<MyMesh::FaceHandle> type1rem = get_hole_removable_faces( mesh, 0);
    std::cout<<"============== Removable type1 holes ===================="<<std::endl;
    for (int i=0; i<type1rem.size(); i++){
    	std::cout<<"face #"<<i<<" "<<type1rem[i]<<std::endl;

    }
    std::cout<<"================================================="<<std::endl;  


	std::vector<MyMesh::FaceHandle> type2rem = get_hole_removable_faces( mesh, 1);
    std::cout<<"============== Removable type2 holes ===================="<<std::endl;
    for (int i=0; i<type2rem.size(); i++){
    	std::cout<<"face #"<<i<<" "<<type2rem[i]<<std::endl;

    }
    std::cout<<"================================================="<<std::endl;      
    std::cout<<"Energy: "<<full_energy(mesh)<<std::endl;
    while (!attempt_hole_insertion(mesh, 1, 1, uw_fake, eng)){
    	std::cout<<"attempting... "<<std::endl;
    }   
 	dump_off(mesh, "/home/btyukodi/testdatafolder/spherocylinder/hole_insertion_debug/dump_after.off");
 	std::cout<<"Energy: "<<full_energy(mesh)<<std::endl;

    while (!attempt_hole_removal(mesh, 1, uw_fake, eng)){
    	std::cout<<"attempting... "<<std::endl;
    }   


		mesh.garbage_collection();

 	dump_off(mesh, "/home/btyukodi/testdatafolder/spherocylinder/hole_insertion_debug/dump_after2.off");
 	std::cout<<"Energy: "<<full_energy(mesh)<<std::endl; 	
 */
	//rand_orientation_test();
	//run_wall("/home/btyukodi/testdatafolder/wall/capsid_input.json");
	//run_umbrella("/home/btyukodi/testdatafolder/checkpoint/capsid_input.json");
	//run_quench("/home/btyukodi/testdatafolder/quench_test/trumpet_quench.json");
	//run_dynamical("/home/btyukodi/testdatafolder/carlos3/minimal_failure/minimal_failure/81.json");

/*
	convert_dat_to_om("/home/btyukodi/assembly_openmesh/sandbox/carlos/vertices.dat",
		"/home/btyukodi/assembly_openmesh/sandbox/carlos/faces.dat",
		"/home/btyukodi/assembly_openmesh/sandbox/carlos/edges.dat",
		"/home/btyukodi/assembly_openmesh/sandbox/carlos/minimal.om");
*/
/*
	MyMesh mesh;
  	mesh.request_face_status();
  	mesh.request_edge_status();
  	mesh.request_vertex_status();    
  	mesh.request_halfedge_status();
	mesh.request_face_normals();

	auto face_props_bare = OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_bare = OpenMesh::HProp<int>(mesh, "edge_type");	

	std::cout<<"run_debug()"<<std::endl;
	std::ifstream vertices("/home/btyukodi/assembly_openmesh/sandbox/carlos/vertices.dat");
	std::ifstream faces("/home/btyukodi/assembly_openmesh/sandbox/carlos/faces.dat");
	std::ifstream edges("/home/btyukodi/assembly_openmesh/sandbox/carlos/edges.dat");	

	double x, y, z;
	long ix, ix1, ix2, ix3;
	int face_type, edge_type;
	long t=0;
	std::vector<MyMesh::VertexHandle> vhandles;
	while(vertices >> ix){
		vertices >> x;
		vertices >> y;
		vertices >> z;

		std::cout<<ix<<" "<<x<<" "<<y<<" "<<z<<std::endl;
		vhandles.push_back(mesh.add_vertex(MyMesh::Point(x,y,z)));
	}
	vertices.close();

	std::vector<MyMesh::VertexHandle> vfhandles;	
	MyMesh::FaceHandle fh;
	while(faces >> ix){
		vfhandles.clear();
		faces >> ix1;
		faces >> ix2;
		faces >> ix3;
		faces >> face_type;
		vfhandles.push_back( vhandles[ix1-1] );
		vfhandles.push_back( vhandles[ix2-1] );
		vfhandles.push_back( vhandles[ix3-1] );				

		fh = mesh.add_face(vfhandles);
		face_props_bare[fh] = face_type;
		std::cout<<ix1<<" "<<ix2<<" "<<ix3<<" "<<face_type<<std::endl;
	}
	faces.close();

	MyMesh::HalfedgeHandle he;
	while (edges >> ix){
		edges >> ix1;
		edges >> ix2;
		edges >> edge_type;
		he = mesh.find_halfedge(vhandles[ix1-1], vhandles[ix2-1]);
		edge_props_bare[he] = edge_type;
	}
	edges.close();

	OpenMesh::FPropHandleT< int > ifprop;	
	OpenMesh::HPropHandleT< int > ihprop;
	mesh.get_property_handle(ifprop, "face_type");
	mesh.get_property_handle(ihprop, "edge_type");
	mesh.property(ifprop).set_persistent(true);
	mesh.property(ihprop).set_persistent(true);

  	mesh.release_face_status();
  	mesh.release_edge_status();
  	mesh.release_vertex_status();    
  	mesh.release_halfedge_status();
	mesh.release_face_normals();	

	std::ofstream os("/home/btyukodi/assembly_openmesh/sandbox/carlos/minimal.om", std::ios::out | std::ios::binary);
	os<<t;			
	try
	  {
	    if ( !OpenMesh::IO::write_mesh(mesh, os, ".OM") )
	    {
	      std::cerr << "Cannot write mesh to stream" << std::endl;
	      return;
	    }
	  }
	  catch( std::exception& x )
	  {
	    std::cerr << x.what() << std::endl;
	    return;
	  }	
	os.close();
	*/
	return;
}
