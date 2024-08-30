#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <regex>
#include <cstdio>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Apps/Assembly/json.hh>
#include <OpenMesh/Apps/Assembly/IO.hh>
#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>
#include <OpenMesh/Apps/Assembly/energy.hh>
#include <OpenMesh/Apps/Assembly/excluders.hh>
#include <sys/stat.h>

// for convenience
using json = nlohmann::json;

/*
get environmental variable name in string format
return empty string if not set
*/
std::string getenv_str(std::string env_name){
	char * env = std::getenv(env_name.c_str());
	std::string env_str="";
	if (env==NULL) return env_str;
	env_str = env;
	return env_str;
}

/*
Read and set mesh parameters from input file
*/
void read_set_mesh_params(MyMesh & mesh, std::string input_filename){

	// read a JSON file
	std::ifstream input_file(input_filename);
	json js;
	input_file >> js;

	json parameters = js["parameters"];
	json mesh_parameters = parameters["mesh_parameters"];

	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	(*mesh_props).kT = mesh_parameters["kT"];
	(*mesh_props).d_max = mesh_parameters["d_max"];
	(*mesh_props).R_add = mesh_parameters["R_add"];	
	(*mesh_props).l_fuse = mesh_parameters["l_fuse"];	
	(*mesh_props).L_neighbor = mesh_parameters["L_neighbor"];
	(*mesh_props).k_fusion = mesh_parameters["k_fusion"]; //0.001;//0.01;//0.01;		
	(*mesh_props).k_fusion2 = mesh_parameters["k_fusion2"];//0.0001;//0.001;	
	(*mesh_props).k_fusion_edge = mesh_parameters["k_fusion_edge"];//0.01;


	if (mesh_parameters.contains("DOS")){
		(*mesh_props).DOS = mesh_parameters["DOS"];	
	}
	else{
		(*mesh_props).DOS = 1.0;
	}
	
	if (mesh_parameters.contains("k_vertex_move")){
		(*mesh_props).k_vertex_move = mesh_parameters["k_vertex_move"];	
	}
	else{
		(*mesh_props).k_vertex_move = 1.0;
	}	

	//loop over all subunit types, make sure they are continuously labeled for the prototype subunits (identified by position in vector)	
	json subunit_parameters = parameters["subunits"];
	std::string subunit_type_str;
	std::vector<int> subunit_types;

	for (json::iterator it = subunit_parameters.begin(); it != subunit_parameters.end(); ++it) {
	  	subunit_type_str = it.key();		
	  	subunit_types.push_back( std::stoi(std::regex_replace(subunit_type_str, std::regex("sub"), ""))-1 );
	}	
	unsigned int max_sub_types = *std::max_element(subunit_types.begin(), subunit_types.end());
	if (max_sub_types != subunit_types.size()-1){
		std::cout<<"Warning! Subunit types have to be counted continuously by the pattern sub1, sub2, ..., subn "<<std::endl;
	}

	//loop over all edge types, make sure they are continuously labeled 1, 2, 3, ..., n for the interaction "matrix" (identified by position in vector)
	json edge_parameters = parameters["edges"];	
	std::string edge_type_str;
	std::vector<int> edge_types;

	for (json::iterator it = edge_parameters.begin(); it != edge_parameters.end(); ++it) {
	  	edge_type_str = it.key();		
	  	edge_types.push_back( std::stoi(std::regex_replace(edge_type_str, std::regex("type"), ""))-1 );
	}	
	unsigned int max_edge_types = *std::max_element(edge_types.begin(), edge_types.end());
	if (max_edge_types != edge_types.size()-1){
		std::cout<<"Warning! Edge types have to be counted continuously by the pattern type1, type2, ..., typen "<<std::endl;
	}


	MyMesh prototypes_mesh;
	std::vector<MyMesh::VertexHandle> vhandles;
	std::vector<MyMesh::FaceHandle> prototype_faces;
	std::vector<MyMesh::HalfedgeHandle> anchor_halfedges;
	MyMesh::FaceHandle proto_face;
	MyMesh::HalfedgeHandle halfedge;	
	std::vector<float> e_b, bend_mod, theta0;	
	auto proto_edge_props = OpenMesh::HProp<HalfedgeProp>(prototypes_mesh, "edge_props");
	auto proto_face_props = OpenMesh::FProp<FaceProp>(prototypes_mesh, "face_props");
	//go over subunit types and add an instance of each to the prototypes mesh
	//for (int sub : subunit_types){
	for (int sub=0; sub<subunit_types.size(); sub++){
		std::cout<<"**sub "<<sub<<std::endl;
		//add and set prototype faces here
		vhandles.clear();
		vhandles.push_back(prototypes_mesh.add_vertex(MyMesh::Point(0,0,0)));
		vhandles.push_back(prototypes_mesh.add_vertex(MyMesh::Point(0,1,0)));
		vhandles.push_back(prototypes_mesh.add_vertex(MyMesh::Point(1,0,0)));	
		proto_face = prototypes_mesh.add_face(vhandles);

		subunit_type_str = "sub"+std::to_string(sub+1);
		proto_face_props[proto_face].face_type = sub;
		proto_face_props[proto_face].mu = subunit_parameters[subunit_type_str]["mu"];
		proto_face_props[proto_face].R_exc = subunit_parameters[subunit_type_str]["R_exc"];	
		proto_face_props[proto_face].k_insertion = subunit_parameters[subunit_type_str]["k_insertion"];	
		proto_face_props[proto_face].number_of_rotational_configs = subunit_parameters[subunit_type_str]["number_of_rotational_configs"];			

		halfedge = prototypes_mesh.halfedge_handle(proto_face);
		//for (int i=0; i<3; i++){
		//go over edges of the subunit type
		//for (std::string sub_edge: subunit_parameters[subunit_type_str]["edges"]){
		for (int se=1; se<=3; se++){
			std::string sub_edge = subunit_parameters[subunit_type_str]["edges"]["edge"+std::to_string(se)];
			//std::cout<<"** "<<sub_edge<<std::endl;

			e_b.clear();
			bend_mod.clear();
			theta0.clear();
			//for -- over all types to ensure order in e_b, bend_mod, theta0 regardless of types order in edges.type1.bending_modulus
			//go over all edge types
			//for (int edge_t: edge_types){ <-- edge_types array is not sorted!!!
			for (int edge_t=0; edge_t<edge_types.size(); edge_t++){
				//std::cout<<"edge_t "<<edge_t<<std::endl;
				//if binding energy provided in the input json
				if (edge_parameters[sub_edge]["e_b"].contains("type"+std::to_string(edge_t+1))){
					e_b.push_back( edge_parameters[sub_edge]["e_b"]["type"+std::to_string(edge_t+1)] );
					bend_mod.push_back( edge_parameters[sub_edge]["bending_modulus"]["type"+std::to_string(edge_t+1)] );
					theta0.push_back( edge_parameters[sub_edge]["theta0"]["type"+std::to_string(edge_t+1)] );		
				}	
				//if binding energy not provided, it is assumed that they don't bind
				else{
					//std::cout<<"binding not provided, assuming infinite repulsion"<<std::endl;
					e_b.push_back( 1000.0 );
					bend_mod.push_back( 1000.0 );
					theta0.push_back( 1000.0 );
				}					
			}


			proto_edge_props[halfedge].l0 = edge_parameters[sub_edge]["l0"];
			proto_edge_props[halfedge].edge_type=std::stoi(std::regex_replace(sub_edge, std::regex("type"), ""))-1;
			proto_edge_props[halfedge].stretch_mod=edge_parameters[sub_edge]["stretch_mod"];
			proto_edge_props[halfedge].e_b = e_b;
			proto_edge_props[halfedge].bend_mod = bend_mod;
			proto_edge_props[halfedge].theta0 = theta0;
			

			halfedge = prototypes_mesh.next_halfedge_handle(halfedge);
		}	

		prototype_faces.push_back(proto_face);
		anchor_halfedges.push_back(halfedge);

		//test if interaction matrices are symmetric
		double val1, val2;
		for (int edget_1: edge_types){		
			for (int edget_2: edge_types){	
				//will throw error if the json is not symmetric at least in the binding edges	
				if (edge_parameters["type"+std::to_string(edget_1+1)]["e_b"].contains("type"+std::to_string(edget_2+1))){

					val1 = edge_parameters["type"+std::to_string(edget_1+1)]["theta0"]["type"+std::to_string(edget_2+1)];
					val2 = edge_parameters["type"+std::to_string(edget_2+1)]["theta0"]["type"+std::to_string(edget_1+1)];
					if (val1 != val2) std::cout<<"Warning, theta0 not symmetric! Will take average "<<val1<<" "<<val2<<std::endl;

					val1 = edge_parameters["type"+std::to_string(edget_1+1)]["e_b"]["type"+std::to_string(edget_2+1)];
					val2 = edge_parameters["type"+std::to_string(edget_2+1)]["e_b"]["type"+std::to_string(edget_1+1)];
					if (val1 != val2) std::cout<<"Warning, e_b not symmetric!  Will take average "<<val1<<" "<<val2<<std::endl;	

					val1 = edge_parameters["type"+std::to_string(edget_1+1)]["bending_modulus"]["type"+std::to_string(edget_2+1)];
					val2 = edge_parameters["type"+std::to_string(edget_2+1)]["bending_modulus"]["type"+std::to_string(edget_1+1)];
					if (val1 != val2) std::cout<<"Warning, bending_modulus not symmetric!  Will take average "<<val1<<" "<<val2<<std::endl;			
				}				
			}
		}

		//std::cout<<"sub"+std::to_string(sub)<<std::endl;		
		//std::cout<<subunit_parameters["sub"+std::to_string(sub+1)]<<std::endl;
	}
	(*mesh_props).prototypes_mesh = prototypes_mesh;
	(*mesh_props).prototype_faces = prototype_faces;	
	(*mesh_props).anchor_halfedges = anchor_halfedges;

	input_file.close();
}


/*
should be a simple simulation param data structure for
data_folder
tmax
dtsave
spring_const
ensemble -> use it for random seed!!
Nperim?..  - no.
alpha?...  - no.
*/
void read_set_run_params(RunParameters & rp, std::string input_filename){
	// read a JSON file
	std::ifstream input_file(input_filename);
	json js;
	input_file >> js;

	json parameters = js["parameters"];
	json run_parameters = parameters["run_parameters"];

	//try to add a root if specified by the environmental variable
	rp.data_folder = run_parameters["data_folder"];
	rp.data_folder=getenv_str("tmp_dir")+rp.data_folder;


	rp.timesteps = run_parameters["timesteps"];
	rp.dtsave = run_parameters["dtsave"];
	rp.dtskip = run_parameters["dtskip"];
	rp.ensemble = run_parameters["ensemble"];
	rp.spring_const = run_parameters["spring_const"];
	rp.nsteps_umbrella = run_parameters["nsteps_umbrella"];

	rp.dN_umbrella= run_parameters["dN_umbrella"];
	rp.init_file = run_parameters["init_file"];
	rp.init_file=getenv_str("init_file_path")+rp.init_file;

	rp.init_file_pos = run_parameters["init_file_pos"];
	rp.adaptive_rates = run_parameters["adaptive_rates"];
	rp.convert_to_lammps = run_parameters["convert_to_lammps"];
	rp.convert_to_vtk = run_parameters["convert_to_vtk"];
	rp.convert_to_lammps_trajectory = run_parameters["convert_to_lammps_trajectory"];
	rp.convert_to_dat = run_parameters["convert_to_dat"];
	rp.dt_burnin = run_parameters["dt_burnin"];

	if (getenv_str("openmesh_src")==""){
		rp.path_to_source = run_parameters["path_to_source"];
	}
	else{
		rp.path_to_source = getenv_str("openmesh_src");
	}
	rp.Nmax = run_parameters["Nmax"];


	if (run_parameters.contains("checkpoint_freq")){
		rp.checkpoint_freq=run_parameters["checkpoint_freq"];
	}
	else{
		rp.checkpoint_freq=99999999;
	}

	if (run_parameters.contains("continue_if_possible")){
		rp.continue_if_possible=run_parameters["continue_if_possible"];
	}
	else{
		rp.continue_if_possible=false;
	}

	if (run_parameters.contains("max_seconds")){
		rp.max_seconds=run_parameters["max_seconds"];
	}
	else{
		rp.max_seconds=99999999999;
	}	


	if (run_parameters.contains("wall_hardness")){
		rp.wall_hardness=run_parameters["wall_hardness"];
	}
	else{
		rp.wall_hardness = 0.0;
	}
	if (run_parameters.contains("wall_amplitude")){
		rp.wall_amplitude=run_parameters["wall_amplitude"];
	}
	else{
		rp.wall_amplitude = 0.0;
	}	
	if (run_parameters.contains("wall_speed")){
		rp.wall_speed=run_parameters["wall_speed"];
	}
	else{
		rp.wall_speed = 0.0;
	}

	if (run_parameters.contains("wall_period")){
		rp.wall_period=run_parameters["wall_period"];
	}
	else{
		rp.wall_period = 99999999999999;
	}	
		
	if (run_parameters.contains("dtsave_snapshot")){
		rp.dtsave_snapshot=run_parameters["dtsave_snapshot"];
	}
	else{
		rp.dtsave_snapshot = run_parameters["dtsave"];
	}	


	if (run_parameters.contains("stop_at_closure")){
		rp.stop_at_closure=run_parameters["stop_at_closure"];
	}
	else{
		rp.stop_at_closure = false;
	}	


	if (run_parameters.contains("externalpotential")){
		rp.externalpotential=run_parameters["externalpotential"];
		rp.extparam1=run_parameters["extparam1"];
		rp.extparam2=run_parameters["extparam2"];
		rp.extparam3=run_parameters["extparam3"];
		rp.extparam4=run_parameters["extparam4"];
		rp.extparam5=run_parameters["extparam5"];
		if (run_parameters.contains("extparam6")){
			rp.extparam6=run_parameters["extparam6"];
		}


	}
	else{
		rp.externalpotential = "DummyPotential";
	}	


	if (run_parameters.contains("k_einstein")){
		rp.k_einstein=run_parameters["k_einstein"];
	}
	else{
		rp.k_einstein = 0.0;
	}		

/*	if (run_parameters.contains("find_optimal_k_einstein")){
		rp.find_optimal_k_einstein=run_parameters["find_optimal_k_einstein"];
	}
	else{
		rp.find_optimal_k_einstein = false;
	}		
*/	
	input_file.close();
	//put seed file here once you figure out how to save persistently
}

//get a mesh with basic properties only that can be saved
//face_type and edge_type only
//mesh passed by copy so that the original mesh is not messed up
//IF THIS BREAKS, remove lines with //%%
MyMesh get_stripped_mesh(MyMesh mesh){
	//drop properties, update with simple properties
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");

	//%%
	auto vertex_props = OpenMesh::VProp<VertexProp>(mesh, "vertex_props");	

	auto face_props_bare = OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_bare = OpenMesh::HProp<int>(mesh, "edge_type");


	//%%
	auto vertex_props_bare = OpenMesh::VProp<int>(mesh, "vertex_type");	

	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
		face_props_bare[*fit] = face_props[*fit].face_type;
	}

	for (MyMesh::HalfedgeIter hit = mesh.halfedges_sbegin(); hit!=mesh.halfedges_end(); ++hit){
		edge_props_bare[*hit] = edge_props[*hit].edge_type;
	}	


	//%%
    for(MyMesh::VertexIter vit = mesh.vertices_sbegin(); vit != mesh.vertices_end(); ++vit) {
    	vertex_props_bare[*vit] = vertex_props[*vit].vertex_type;
	}


	//remove the complex properties
	OpenMesh::FPropHandleT< FaceProp > fprop;
	OpenMesh::MPropHandleT< MeshProp > mprop;
	OpenMesh::HPropHandleT< HalfedgeProp > hprop;

	//%%
	OpenMesh::VPropHandleT< VertexProp > vprop;			

	mesh.get_property_handle(fprop, "face_props");
	mesh.remove_property( fprop );		

	mesh.get_property_handle(hprop, "edge_props");
	mesh.remove_property( hprop );	

	mesh.get_property_handle(mprop, "mesh_props");
	mesh.remove_property( mprop );	

	//%%
	mesh.get_property_handle(vprop, "vertex_props");
	mesh.remove_property( vprop );		

	OpenMesh::FPropHandleT< int > ifprop;	
	OpenMesh::HPropHandleT< int > ihprop;

	//%%
	OpenMesh::VPropHandleT< int > ivprop;

	mesh.get_property_handle(ifprop, "face_type");
	mesh.get_property_handle(ihprop, "edge_type");

	//%%
	mesh.get_property_handle(ivprop, "vertex_type");

	mesh.property(ifprop).set_persistent(true);
	mesh.property(ihprop).set_persistent(true);

	//%%
	mesh.property(ivprop).set_persistent(true);	

  	mesh.release_face_status();
  	mesh.release_edge_status();
  	mesh.release_vertex_status();    
  	mesh.release_halfedge_status();
	mesh.release_face_normals(); 

	return mesh;

}

/*
Initialize a mesh from .om file. Needs the position in the file and parameter filename because
.om file only takes the edge and face types, not the actual e_b, bend_mod, theta0, mu, etc. values
*/
void init_mesh_from_om(MyMesh & mesh, std::string om_file, long long om_byte_position, std::string input_filename){
	//request strip mesh properties

	auto edge_type_prop = OpenMesh::HProp<int>(mesh, "edge_type");
	auto face_type_prop = OpenMesh::FProp<int>(mesh, "face_type");	


	//read the mesh here from .om
	const std::string _ext=".OM";
	OpenMesh::IO::Options ropt, wopt;
	//ropt+=OpenMesh::IO::Options::FaceColor; 	
	std::ifstream omfile(om_file, std::ios::out | std::ios::binary);
	if (!omfile.is_open()) std::cout<<"init_mesh_from_om() could not open om init file"<<std::endl;
	omfile.seekg(om_byte_position, omfile.beg); 
	long timestep;
	omfile >> timestep;
	if (!OpenMesh::IO::read_mesh(mesh, omfile, _ext, ropt)) 			
	{
		std::cerr << "init config read error\n";
    	exit(1);
	}	

	omfile.close();
	read_set_mesh_params(mesh, input_filename);

	//at this point there is a stripped mesh with prototypes

	//loop through the mesh, set face and edge properties

	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");	
	auto proto_edge_props = OpenMesh::HProp<HalfedgeProp>((*mesh_props).prototypes_mesh, "edge_props");
	


	MyMesh::FaceHandle fh;
	MyMesh::HalfedgeHandle anchor_halfedge, anchor_halfedge_proto;
	int face_type, edge_type;
	for (MyMesh::FaceIter face=mesh.faces_begin(); face!=mesh.faces_end(); ++face){
		fh = *face;
		face_type = face_type_prop[fh];
		if (face_type>=0){
			anchor_halfedge_proto = (*mesh_props).anchor_halfedges[face_type];
		}
		else{
			anchor_halfedge_proto = (*mesh_props).anchor_halfedges[-(face_type+2)];			
		}
		//for now, find any edge as anchor
		//NEEDS to match anchor he
		anchor_halfedge = mesh.halfedge_handle(fh);
		for (int ir=0; ir<3; ir++){
			edge_type = edge_type_prop[anchor_halfedge];
			if ( proto_edge_props[anchor_halfedge_proto].edge_type == edge_type){
				//std::cout<<"ROT POS FOUND facetype "<<face_type <<std::endl;
				break;
			}

			anchor_halfedge = mesh.next_halfedge_handle(anchor_halfedge);
		}		

		//try to do abs(face_type) below; then if face_type<0, then it should be non-removable
		//what to do with face_type=0?? 0->-1; 1->-2; ... ; n -> -n-1
		//then not quite absolute value, do a test based on sign
		if (face_type==-1){
			std::cout<<"something's screwed up, found face_type=-1"<<std::endl;
		}
		if (face_type>=0){
			clone_properties_from_prototype(mesh, face_type, fh, anchor_halfedge);
		}
		//negative face_type is reserved for non-removable faces as per Carlos request.
		//in the input files, face types are 1, 2, 3 ... and corresponding non-removable ones -1, -2, -3
		//in the internal representation these transform to 0, 1, 2, ... and -2, -3, -4, ...
		else{
			clone_properties_from_prototype(mesh, -(face_type+2), fh, anchor_halfedge);
			face_props[fh].face_type = face_type;
		}

		mesh.calc_face_centroid(fh, face_props[fh].COM );		
	}
	update_full_neighbor_list(mesh);	



	//remove the stripped mesh properties
	OpenMesh::FPropHandleT< int > fprop;
	OpenMesh::HPropHandleT< int > hprop;

	mesh.get_property_handle(fprop, "face_type");	
	mesh.get_property_handle(hprop, "edge_type");	
	mesh.remove_property( fprop );		
	mesh.remove_property( hprop );	

	std::cout<<"----- Initial mesh -------------"<<std::endl;
	print_mesh(mesh);
	std::cout<<"--------------------------------"<<std::endl;
	return;
}

//dump the into an .om stream
//maybe should add the key stream here too? and the timestep?
long long dump_om(MyMesh mesh, long t, std::ostream& os){

	MyMesh simp_mesh = get_stripped_mesh(mesh);
	auto face_props_bare = OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_bare = OpenMesh::HProp<int>(mesh, "edge_type");
	
	long long current_position = os.tellp();
	//std::cout<<"tellp1 "<<os.tellp()<<std::endl; 
	OpenMesh::FPropHandleT< int > ifprop;	
	OpenMesh::HPropHandleT< int > ihprop;
	simp_mesh.get_property_handle(ifprop, "face_type");
	simp_mesh.get_property_handle(ihprop, "edge_type");
	simp_mesh.property(ifprop).set_persistent(true);
	simp_mesh.property(ihprop).set_persistent(true);

	//write the timestep first, then the mesh 			
	os<<t;			
	try
	  {
	    if ( !OpenMesh::IO::write_mesh(simp_mesh, os, ".OM") )
	    {
	      std::cerr << "Cannot write mesh to stream" << std::endl;
	      return 0;
	    }
	  }
	  catch( std::exception& x )
	  {
	    std::cerr << "EXCEPTION" << x.what() << std::endl;
	    return 0;
	  }	

	//add separator in case keys are lost
	//probably no nead; the reader should know where the object ends  
	//os<<"-----------------";
	return current_position;
}

//dump mesh in a separate .om file, no timestep
int dump_om(MyMesh mesh, std::string om_filename){
	MyMesh simp_mesh = get_stripped_mesh(mesh);
	OpenMesh::FPropHandleT< int > ifprop;	
	OpenMesh::HPropHandleT< int > ihprop;
	mesh.get_property_handle(ifprop, "face_type");
	mesh.get_property_handle(ihprop, "edge_type");
	mesh.property(ifprop).set_persistent(true);
	mesh.property(ihprop).set_persistent(true);
		
	try
	  {
	    if ( !OpenMesh::IO::write_mesh(simp_mesh, om_filename) )
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


//dump mesh in a separate .om file, no timestep
//assumes stripped mesh input
long long dump_om_stripped(MyMesh mesh, long t, std::ostream& os){

	//MyMesh simp_mesh = get_stripped_mesh(mesh);
	auto face_props_bare = OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_bare = OpenMesh::HProp<int>(mesh, "edge_type");
	
	long long current_position = os.tellp();
	//std::cout<<"tellp1 "<<os.tellp()<<std::endl; 
	OpenMesh::FPropHandleT< int > ifprop;	
	OpenMesh::HPropHandleT< int > ihprop;
	mesh.get_property_handle(ifprop, "face_type");
	mesh.get_property_handle(ihprop, "edge_type");
	mesh.property(ifprop).set_persistent(true);
	mesh.property(ihprop).set_persistent(true);

	//write the timestep first, then the mesh 			
	os<<t;			
	try
	  {
	    if ( !OpenMesh::IO::write_mesh(mesh, os, ".OM") )
	    {
	      std::cerr << "Cannot write mesh to stream" << std::endl;
	      return 0;
	    }
	  }
	  catch( std::exception& x )
	  {
	    std::cerr << "EXCEPTION" << x.what() << std::endl;
	    return 0;
	  }	

	//add separator in case keys are lost
	//probably no nead; the reader should know where the object ends  
	//os<<"-----------------";
	return current_position;
}


//dump columns to the data file
int dump_data(MyMesh mesh, long timestep, std::ostream& os, long long key_to_snapshot){

	os.precision(6);	
	os.setf(std::ios::fixed);
	os.setf(std::ios::showpoint);

	double E_el = full_elastic_energy(mesh);
	double E_full = full_energy(mesh);
	os<<timestep<<"\t"<<key_to_snapshot <<"\t" <<mesh.n_faces() <<"\t"<<mesh.n_vertices()<<"\t"<<mesh.n_edges()<<"\t"<< E_el<<"\t"<<E_full<<std::endl;

	return 0;
}

//dump columns to the data file
int dump_data(MyMesh mesh, long timestep, std::ostream& os, long long key_to_snapshot, std::vector<double> r_acc){

	os.precision(6);	
	os.setf(std::ios::fixed);
	os.setf(std::ios::showpoint);

	double E_el = full_elastic_energy(mesh);
	double E_full = full_energy(mesh);
	os<<timestep<<"\t"<<key_to_snapshot <<"\t" <<mesh.n_faces() <<"\t"<<mesh.n_vertices()<<"\t"<<mesh.n_edges()<<"\t"<< E_el<<"\t"<<E_full<<"\t";
	for (double &r: r_acc){
		os<<r<<"\t";
	}

	os<<std::endl;

	return 0;
}

//convert a mesh to lammps string
//and writes it to the two output streams
int dump_lammps(MyMesh mesh,  std::ostream& os,  std::ostream& os_bonds){
	os.precision(6);	
	os.setf(std::ios::fixed);
	os.setf(std::ios::showpoint);
	os_bonds.precision(6);	
	os_bonds.setf(std::ios::fixed);
	os_bonds.setf(std::ios::showpoint);	

	//1. write the lammps trajectory files; note that atom indices change between snapshots due to garbage collection
	//center it; subtract COM
	MyMesh::Point com = MyMesh::Point(0,0,0);
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		com+=mesh.point(*v) / mesh.n_vertices();
	}
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		mesh.set_point(*v, mesh.point(*v) - com );
	}

	os<<"ITEM: NUMBER OF ATOMS"<<std::endl;
	os<<mesh.n_vertices()<<std::endl;
	os<<"ITEM: BOX BOUNDS pp pp pp"<<std::endl;
	os<<"-25 25"<<std::endl;
	os<<"-25 25"<<std::endl;	
	os<<"-25 25"<<std::endl;
	os<<"ITEM: ATOMS id type x y z"<<std::endl;
	int atom_type=1;
	MyMesh::Point p;	
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		p = mesh.point(*v);
		os<<(*v).idx()+1<<" "<<atom_type<<" "<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
	}

	//2. write the bonds file
	os_bonds<<"ITEM: NUMBER OF ENTRIES"<<std::endl;
	os_bonds<<mesh.n_edges()<<std::endl;
	os_bonds<<"ITEM: BOX BOUNDS pp pp pp"<<std::endl;
	os_bonds<<"-25 25"<<std::endl;
	os_bonds<<"-25 25"<<std::endl;
	os_bonds<<"-25 25"<<std::endl;		
	os_bonds<<"ITEM: ENTRIES index c_1[1] c_1[2] c_1[3] c_2[1] c_2[2]"<<std::endl;


	int bond_type=1;
	//!!! add calculation of engpot vor visualization once proper init from .om is solved;
	//until then mesh will only have simple face and edge types, no elastic parameters for energy
	double dist, engpot;
	engpot=0.0;
	MyMesh::HalfedgeHandle hedge;
	MyMesh::VertexHandle v1, v2;
	for (MyMesh::EdgeIter e_it=mesh.edges_sbegin(); e_it!=mesh.edges_end(); ++e_it){
		hedge = mesh.halfedge_handle(*e_it, 0);
		v1 = mesh.from_vertex_handle(hedge);
		v2 = mesh.to_vertex_handle(hedge);
		dist = mesh.calc_edge_length(*e_it);	
		engpot=0.0;
		os_bonds<<(*e_it).idx()+1<<" "<< bond_type <<" "<<v1.idx()+1<<" "<<v2.idx()+1<<" "<<dist<<" "<<engpot<<std::endl;
	}



	return 1;

}	

int dump_lammps_snapshot(MyMesh mesh, std::string outfile_name){
   	std::ofstream os( outfile_name );

	os.precision(6);	
	os.setf(std::ios::fixed);
	os.setf(std::ios::showpoint);

	MyMesh::Point com = MyMesh::Point(0,0,0);
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		com+=mesh.point(*v) / mesh.n_vertices();
	}
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
	//@@@	mesh.set_point(*v, mesh.point(*v) - com );
	}	


	os<<"LAMMPSDescription-Generated by btyukodi"<<std::endl;
	os<<std::endl;
	os<<mesh.n_vertices()<<" atoms"<<std::endl;
	os<<mesh.n_edges()<<" bonds"<<std::endl;
	os<<std::endl;
	os<<"6 atom types"<<std::endl;
	os<<"2 bond types"<<std::endl;
	os<<std::endl;
	os<<"  -12.000    12.000 xlo xhi"<<std::endl;
	os<<"  -12.000    12.000 ylo yhi"<<std::endl;
	os<<"  -12.000    12.000 zlo zhi"<<std::endl;		
	os<<std::endl;
	os<<"Atoms"<<std::endl;
	os<<std::endl;
	int atom_type=1;
	MyMesh::Point p;	
	auto vertex_props_stripped = OpenMesh::VProp<int>(mesh, "vertex_type");
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		p = mesh.point(*v);
		atom_type=1;
		/*if (mesh.valence(*v)==4){
			atom_type=2;
		}
		if (mesh.valence(*v)==5){
			atom_type=3;
		}		
		if (mesh.valence(*v)==7){
			atom_type=4;
		}
		*/
		
		if (!mesh.is_boundary(*v)){
			for (int disclination_type=4; disclination_type<=8; disclination_type++)
			if (mesh.valence(*v)==disclination_type ){
				atom_type=disclination_type-2;
			}
		}		
		
		//atom_type=vertex_props_stripped[*v]+1;

		os<<(*v).idx()+1<<" "<<atom_type<<" "<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
		//std::cout<<"valence "<<mesh.valence(*v)<<std::endl;
	}
	os<<std::endl;
	os<<"Bonds"<<std::endl;
	os<<std::endl;
	//int bond_type=1;
	//for now, use bond_type to color boundary edges to follow cracks
	int bond_type;
	MyMesh::HalfedgeHandle hedge;
	MyMesh::VertexHandle v1, v2;

	auto edge_props_stripped = OpenMesh::HProp<int>(mesh, "edge_type");
	for (MyMesh::EdgeIter e_it=mesh.edges_sbegin(); e_it!=mesh.edges_end(); ++e_it){
		hedge = mesh.halfedge_handle(*e_it, 0);
		v1 = mesh.from_vertex_handle(hedge);
		v2 = mesh.to_vertex_handle(hedge);
		
		if (mesh.is_boundary(*e_it)){
			bond_type=2;
		}
		else{
			bond_type=1;
		}
		//bond_type=edge_props_stripped[hedge]+1;
		os<<(*e_it).idx()+1<<" "<< bond_type <<" "<<v1.idx()+1<<" "<<v2.idx()+1<<std::endl;

	}

	os.close();
	return 1;
}



int dump_lammps_snapshot_double_edge(MyMesh mesh, std::string outfile_name){
   	std::ofstream os( outfile_name );

	os.precision(6);	
	os.setf(std::ios::fixed);
	os.setf(std::ios::showpoint);

   	std::ofstream os_vertex_tmp( outfile_name + "vtmp");

	os_vertex_tmp.precision(6);	
	os_vertex_tmp.setf(std::ios::fixed);
	os_vertex_tmp.setf(std::ios::showpoint);	

	MyMesh::Point com = MyMesh::Point(0,0,0);
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		com+=mesh.point(*v) / mesh.n_vertices();
	}
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		mesh.set_point(*v, mesh.point(*v) - com );
	}	


	int atom_type=1;
	MyMesh::Point p, dp;	


   	std::ofstream os_edge_tmp( outfile_name+"tmp" );

	os_edge_tmp.precision(6);	
	os_edge_tmp.setf(std::ios::fixed);
	os_edge_tmp.setf(std::ios::showpoint);

	os_edge_tmp<<"Bonds"<<std::endl;
	os_edge_tmp<<std::endl;

	auto edge_props_stripped = OpenMesh::HProp<int>(mesh, "edge_type");

	int vertex_id=1;
	double bond_width=0.1;
	int bond_type;
	int atom_types=1;
	int bond_types=999;
	int n_bonds=0;
	int t=1;
	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
  		for (MyMesh::FaceHalfedgeIter fh = mesh.fh_iter(*face); fh.is_valid(); ++fh){
			atom_type=1;

			mesh.calc_face_centroid(face, com );

			p = mesh.point( mesh.from_vertex_handle(*fh)  );

			dp = (com-p).normalize();
			p=p+dp*bond_width;
			os_vertex_tmp<<vertex_id<<" "<<atom_type<<" "<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;

			bond_type=edge_props_stripped[*fh]+1;
			if (t%3 != 0){
				os_edge_tmp<<(*fh).idx()+1<<" "<< bond_type <<" "<<vertex_id<<" "<<vertex_id+1<<std::endl;
				t++;
			}
			else{
				os_edge_tmp<<(*fh).idx()+1<<" "<< bond_type <<" "<<vertex_id<<" "<<vertex_id-2<<std::endl;
				t=1;
			}
			n_bonds++;
			vertex_id++;

		}
	}

	os_edge_tmp.close();
	os_vertex_tmp.close();
	std::ifstream if_edge(outfile_name+"tmp");
	std::ifstream if_vertex(outfile_name+"vtmp");


	os<<"LAMMPSDescription-Generated by btyukodi"<<std::endl;
	os<<std::endl;
	os<<vertex_id-1<<" atoms"<<std::endl;
	os<<n_bonds<<" bonds"<<std::endl;
	os<<std::endl;
	os<< atom_types <<" atom types"<<std::endl;
	os<<bond_types<<" bond types"<<std::endl;
	os<<std::endl;
	os<<"  -12.000    12.000 xlo xhi"<<std::endl;
	os<<"  -12.000    12.000 ylo yhi"<<std::endl;
	os<<"  -12.000    12.000 zlo zhi"<<std::endl;		
	os<<std::endl;
	os<<"Atoms"<<std::endl;
	os<<std::endl;

	os<< if_vertex.rdbuf() << std::endl<<if_edge.rdbuf();
	os.close();
	os_edge_tmp.close();
	os_vertex_tmp.close();
	std::remove((outfile_name+"tmp").c_str());
	std::remove((outfile_name+"vtmp").c_str());
	
//-----------

//	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
//		p = mesh.point(*v);
//		atom_type=1;
//		/*if (mesh.valence(*v)==4){
//			atom_type=2;
//		}
//		if (mesh.valence(*v)==5){
//			atom_type=3;
//		}		
//		if (mesh.valence(*v)==7){
//			atom_type=4;
//		}*/
//
//		if (!mesh.is_boundary(*v)){
//			for (int disclination_type=4; disclination_type<=8; disclination_type++)
//			if (mesh.valence(*v)==disclination_type ){
//				atom_type=disclination_type-2;
//			}
//		}						
//
//		os<<(*v).idx()+1<<" "<<atom_type<<" "<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
//		//std::cout<<"valence "<<mesh.valence(*v)<<std::endl;
//	}
//	os<<std::endl;
//	os<<"Bonds"<<std::endl;
//	os<<std::endl;
//	//int bond_type=1;
//	//for now, use bond_type to color boundary edges to follow cracks
//	int bond_type;
//	MyMesh::HalfedgeHandle hedge;
//	MyMesh::VertexHandle v1, v2;
//
//	auto edge_props_stripped = OpenMesh::HProp<int>(mesh, "edge_type");
//	for (MyMesh::EdgeIter e_it=mesh.edges_sbegin(); e_it!=mesh.edges_end(); ++e_it){
//		hedge = mesh.halfedge_handle(*e_it, 0);
//		v1 = mesh.from_vertex_handle(hedge);
//		v2 = mesh.to_vertex_handle(hedge);
//		
//		if (mesh.is_boundary(*e_it)){
//			bond_type=2;
//		}
//		else{
//			bond_type=1;
//		}
//		//bond_type=edge_props_stripped[hedge]+1;
//		os<<(*e_it).idx()+1<<" "<< bond_type <<" "<<v1.idx()+1<<" "<<v2.idx()+1<<std::endl;
//
//	}
//
//	os.close();
	return 1;
}
	


int dump_lammps_typed_edges_snapshot(MyMesh mesh, std::string outfile_name){
   	std::ofstream os( outfile_name );

	os.precision(6);	
	os.setf(std::ios::fixed);
	os.setf(std::ios::showpoint);

   	std::ofstream os_vertex_tmp( outfile_name + "vtmp");

	os_vertex_tmp.precision(6);	
	os_vertex_tmp.setf(std::ios::fixed);
	os_vertex_tmp.setf(std::ios::showpoint);	

	MyMesh::Point com = MyMesh::Point(0,0,0);
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		com+=mesh.point(*v) / mesh.n_vertices();
	}
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		mesh.set_point(*v, mesh.point(*v) - com );
	}	


	int atom_type=1;
	MyMesh::Point p, dp;	


   	std::ofstream os_edge_tmp( outfile_name+"tmp" );

	os_edge_tmp.precision(6);	
	os_edge_tmp.setf(std::ios::fixed);
	os_edge_tmp.setf(std::ios::showpoint);

	os_edge_tmp<<"Bonds"<<std::endl;
	os_edge_tmp<<std::endl;

	auto edge_props_stripped = OpenMesh::HProp<int>(mesh, "edge_type");

	int vertex_id=1;
	double bond_width=0.0;
	int bond_type;
	int atom_types=1;
	int bond_types=999;
	int n_bonds=0;
	int t=1;
	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
  		for (MyMesh::FaceHalfedgeIter fh = mesh.fh_iter(*face); fh.is_valid(); ++fh){
			atom_type=1;

			mesh.calc_face_centroid(face, com );

			p = mesh.point( mesh.from_vertex_handle(*fh)  );

			dp = (com-p).normalize();
			p=p+dp*bond_width;
			os_vertex_tmp<<vertex_id<<" "<<atom_type<<" "<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;

			bond_type=edge_props_stripped[*fh]+1;
			if (t%3 != 0){
				os_edge_tmp<<(*fh).idx()+1<<" "<< bond_type <<" "<<vertex_id<<" "<<vertex_id+1<<std::endl;
				t++;
			}
			else{
				os_edge_tmp<<(*fh).idx()+1<<" "<< bond_type <<" "<<vertex_id<<" "<<vertex_id-2<<std::endl;
				t=1;
			}
			n_bonds++;
			vertex_id++;

		}
	}

	os_edge_tmp.close();
	os_vertex_tmp.close();
	std::ifstream if_edge(outfile_name+"tmp");
	std::ifstream if_vertex(outfile_name+"vtmp");


	os<<"LAMMPSDescription-Generated by btyukodi"<<std::endl;
	os<<std::endl;
	os<<vertex_id-1<<" atoms"<<std::endl;
	os<<n_bonds<<" bonds"<<std::endl;
	os<<std::endl;
	os<< atom_types <<" atom types"<<std::endl;
	os<<bond_types<<" bond types"<<std::endl;
	os<<std::endl;
	os<<"  -12.000    12.000 xlo xhi"<<std::endl;
	os<<"  -12.000    12.000 ylo yhi"<<std::endl;
	os<<"  -12.000    12.000 zlo zhi"<<std::endl;		
	os<<std::endl;
	os<<"Atoms"<<std::endl;
	os<<std::endl;

	os<< if_vertex.rdbuf() << std::endl<<if_edge.rdbuf();
	os.close();
	os_edge_tmp.close();
	os_vertex_tmp.close();
	std::remove((outfile_name+"tmp").c_str());
	std::remove((outfile_name+"vtmp").c_str());
	
//-----------

//	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
//		p = mesh.point(*v);
//		atom_type=1;
//		/*if (mesh.valence(*v)==4){
//			atom_type=2;
//		}
//		if (mesh.valence(*v)==5){
//			atom_type=3;
//		}		
//		if (mesh.valence(*v)==7){
//			atom_type=4;
//		}*/
//
//		if (!mesh.is_boundary(*v)){
//			for (int disclination_type=4; disclination_type<=8; disclination_type++)
//			if (mesh.valence(*v)==disclination_type ){
//				atom_type=disclination_type-2;
//			}
//		}						
//
//		os<<(*v).idx()+1<<" "<<atom_type<<" "<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
//		//std::cout<<"valence "<<mesh.valence(*v)<<std::endl;
//	}
//	os<<std::endl;
//	os<<"Bonds"<<std::endl;
//	os<<std::endl;
//	//int bond_type=1;
//	//for now, use bond_type to color boundary edges to follow cracks
//	int bond_type;
//	MyMesh::HalfedgeHandle hedge;
//	MyMesh::VertexHandle v1, v2;
//
//	auto edge_props_stripped = OpenMesh::HProp<int>(mesh, "edge_type");
//	for (MyMesh::EdgeIter e_it=mesh.edges_sbegin(); e_it!=mesh.edges_end(); ++e_it){
//		hedge = mesh.halfedge_handle(*e_it, 0);
//		v1 = mesh.from_vertex_handle(hedge);
//		v2 = mesh.to_vertex_handle(hedge);
//		
//		if (mesh.is_boundary(*e_it)){
//			bond_type=2;
//		}
//		else{
//			bond_type=1;
//		}
//		//bond_type=edge_props_stripped[hedge]+1;
//		os<<(*e_it).idx()+1<<" "<< bond_type <<" "<<v1.idx()+1<<" "<<v2.idx()+1<<std::endl;
//
//	}
//
//	os.close();
	return 1;
}	


//convert om file to lammps
//for now, only the file name; once energy is added to lammps file, param file will be needed too
int convert_om_to_lammps_trajectory(std::string om_filename, std::string lammps_base){
	//std::cout<<"lammps traj"<<std::endl;
   	std::ifstream omfile(om_filename, std::ios::out | std::ios::binary);
	if (!omfile.is_open()) std::cout<<"convert_om_to_lammps_trajectory could not open om file"<<std::endl;   	


   	//std::ofstream trajfile(lammps_base + std::regex_replace(om_filename, std::regex("om$"), "traj") );
   	//std::ofstream bondfile(lammps_base + std::regex_replace(om_filename, std::regex("om$"), "bond") ); 
   	std::ofstream trajfile(lammps_base + "lammps_traj" );
   	std::ofstream bondfile(lammps_base + "lammps_bond" );      	  	

	//alernatingly read a timestep (long) and a mesh (.om) objects

	const std::string _ext=".OM";
	OpenMesh::IO::Options ropt, wopt;	
	ropt+=OpenMesh::IO::Options::FaceColor; 	
	long timestep;
	MyMesh mesh;
	while (omfile >> timestep){
   		//std::cout<<"lammps traj"<<  lammps_base + std::regex_replace(om_filename, std::regex("om$"), "traj") <<std::endl;
   		OpenMesh::IO::read_mesh(mesh, omfile, _ext, ropt);
   		//trajfile<<"ITEM: TIME"<<std::endl;
   		//trajfile<<timestep<<std::endl;   		
   		trajfile<<"ITEM: TIMESTEP"<<std::endl;
   		trajfile<<timestep<<std::endl;

   		bondfile<<"ITEM: TIMESTEP"<<std::endl;
   		bondfile<<timestep<<std::endl;   		
   		dump_lammps(mesh,  trajfile,  bondfile);

   		//std::cout<<"timestep " <<timestep<<" "<<mesh.n_faces()<<" "<< omfile.tellg()<<std::endl;

	}

	omfile.close();
	trajfile.close();
	bondfile.close();
	return 1;
}

int convert_om_to_lammps_snapshots(std::string om_filename, std::string lammps_base){

	std::ifstream omfile(om_filename, std::ios::out | std::ios::binary);

 	//alernatingly read a timestep (long) and a mesh (.om) objects

	const std::string _ext=".OM";
	OpenMesh::IO::Options ropt, wopt;	
	ropt+=OpenMesh::IO::Options::FaceColor; 	
	//%
	//ropt+=OpenMesh::IO::Options::VertexColor; 		

	long timestep;
	MyMesh mesh;
	auto face_props_bare = OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_bare = OpenMesh::HProp<int>(mesh, "edge_type");	
	//%%
	auto vertex_props_bare = OpenMesh::VProp<int>(mesh, "vertex_type");				
	while (omfile >> timestep){
   		
   		OpenMesh::IO::read_mesh(mesh, omfile, _ext, ropt);
   		dump_lammps_snapshot(mesh,  lammps_base + std::to_string(timestep) +".dat");

	}
	omfile.close();
	return 1;

}


int convert_om_to_lammps_snapshots_double_edge(std::string om_filename, std::string lammps_base){

	std::ifstream omfile(om_filename, std::ios::out | std::ios::binary);

 	//alernatingly read a timestep (long) and a mesh (.om) objects

	const std::string _ext=".OM";
	OpenMesh::IO::Options ropt, wopt;	
	ropt+=OpenMesh::IO::Options::FaceColor; 	
	long timestep;
	MyMesh mesh;
	auto face_props_bare = OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_bare = OpenMesh::HProp<int>(mesh, "edge_type");		
	//%% WHY ISN'T THIS VProp ??? CHECK!!!!
	auto vertex_props_bare = OpenMesh::HProp<int>(mesh, "vertex_type");		
	while (omfile >> timestep){
   		
   		OpenMesh::IO::read_mesh(mesh, omfile, _ext, ropt);
   		dump_lammps_snapshot_double_edge(mesh,  lammps_base + std::to_string(timestep) +".dat");

	}
	omfile.close();
	return 1;

}


//this should be vtk_snapshot and filename be given
int dump_vtk_snapshot(MyMesh mesh, std::string outfile_name, std::string vtk_conf_file){
   	std::ofstream os( outfile_name );

	os.precision(6);	
	os.setf(std::ios::fixed);
	os.setf(std::ios::showpoint);

	MyMesh::Point com = MyMesh::Point(0,0,0);
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		com+=mesh.point(*v) / mesh.n_vertices();
	}

	//SUBTRACT COM
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		mesh.set_point(*v, mesh.point(*v) - com );
	}

	//std::cout<<"breakp1";
	//std::ifstream colors("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/vtk_colors.conf");
	std::ifstream colors(vtk_conf_file);
	std::vector<double> red, green, blue, alpha;
	double tmp;
	//std::cout<<"colors reading";	

	while (!colors.eof()){		
		colors >> tmp;
		red.push_back(tmp);
		colors >> tmp;
		blue.push_back(tmp);
		colors >>tmp;
		green.push_back(tmp);
		colors >>tmp;
		alpha.push_back(tmp);
	}
	//std::cout<<"colors read";
	colors.close();
	os<<"# vtk DataFile Version 3.1"<<std::endl;
	os<<"# created by btyukodi"<<std::endl;	
	os<<"ASCII"<<std::endl;
	os<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	os<<"POINTS "<< mesh.n_vertices()<<" FLOAT"<<std::endl;	
	MyMesh::Point p;
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		p = mesh.point(*v);
		os<<"\t "<<p[0]<<" \t "<<p[1]<<" \t "<<p[2]<<std::endl;
	}
	os<<"CELLS "<<mesh.n_faces()<<" "<<mesh.n_faces()<<std::endl;

	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
		os<<"3";
		for (MyMesh::FaceVertexIter fv = mesh.fv_iter(*face); fv.is_valid(); ++fv ){
			os<<" "<<(*fv).idx();
		}
		os<<std::endl;
	}
	os<<"CELL_TYPES "<<mesh.n_faces()<<std::endl;
	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
		os<<"5"<<std::endl;
	}		
	os<<"CELL_DATA "<<mesh.n_faces()<<std::endl;
	os<<"COLOR_SCALARS faceColor 4"<<std::endl;	

	bool stripped_mesh;
	if (OpenMesh::hasProperty<OpenMesh::FaceHandle, FaceProp>(mesh, "face_props")) {
		stripped_mesh=false;		
	}
	else{
		stripped_mesh=true;
	}

	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto face_props_stripped =  OpenMesh::FProp<int>(mesh, "face_type");

	int face_type;
	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
		//get face type
		if (stripped_mesh){
			face_type = face_props_stripped[*face];
		}
		else{
			face_type = face_props[*face].face_type;			
		}
	//std::cout<<"strp? "<<stripped_mesh<<std::endl;	
	//std::cout<<"face type  "<<face_type<<std::endl;	
	if (face_type>=0){	
		os<<red[face_type]<<" "<<green[face_type]<<" "<<blue[face_type]<<" "<<alpha[face_type]<<std::endl;
	}
	else{
		os<<"1.0 1.0 0.0 1.0"<<std::endl;
	}

	}		
	os.close();
	return 1;
}


//this should be vtk_snapshot and filename be given
int dump_vtk_snapshot_corrected_but_vtk_conf_needs_column_swap(MyMesh mesh, std::string outfile_name, std::string vtk_conf_file){
  	mesh.request_face_status();
  	mesh.request_edge_status();
  	mesh.request_vertex_status();    
  	mesh.request_halfedge_status();
	mesh.request_face_normals();	
   	std::ofstream os( outfile_name );

	os.precision(6);	
	os.setf(std::ios::fixed);
	os.setf(std::ios::showpoint);

	MyMesh::Point com = MyMesh::Point(0,0,0);
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		com+=mesh.point(*v) / mesh.n_vertices();
	}
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		mesh.set_point(*v, mesh.point(*v) - com );
	}

	//std::cout<<"breakp1";
	//std::ifstream colors("/home/btyukodi/assembly_openmesh/OpenMesh-9.0/src/OpenMesh/Apps/Assembly/vtk_colors.conf");
	std::ifstream colors(vtk_conf_file);
	std::vector<double> red, green, blue, alpha;
	double tmp;
	//std::cout<<"colors reading";	

	while (!colors.eof()){		
		colors >> tmp;
		red.push_back(tmp);
		colors >> tmp;
		green.push_back(tmp);
		colors >>tmp;
		blue.push_back(tmp);
		colors >>tmp;
		alpha.push_back(tmp);
	}
	//std::cout<<"colors read";
	colors.close();
	os<<"# vtk DataFile Version 3.1"<<std::endl;
	os<<"# created by btyukodi"<<std::endl;	
	os<<"ASCII"<<std::endl;
	os<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

	os<<"POINTS "<< mesh.n_vertices() + 3*mesh.n_faces()<<" FLOAT"<<std::endl;	
	MyMesh::Point p;
	//actual vertices
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		p = mesh.point(*v);
		os<<"\t "<<p[0]<<" \t "<<p[1]<<" \t "<<p[2]<<std::endl;
	}
	//3 extra vertex per face for parallel triangles
	mesh.update_normals();
	MyMesh::Point normal;
	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
		normal = mesh.normal(*face);
		for (MyMesh::FaceVertexIter fv = mesh.fv_iter(*face); fv.is_valid(); ++fv ){
			//mesh.update_normal(*fv);

			p = mesh.point(*fv) - 1e-5*normal;
			//std::cout<<normal[0]<<std::endl;
			os<<"\t "<<p[0]<<" \t "<<p[1]<<" \t "<<p[2]<<std::endl;			
		}
	}


	os<<"CELLS "<<2*mesh.n_faces()<<" "<<2*mesh.n_faces()<<std::endl;

	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
		os<<"3";
		for (MyMesh::FaceVertexIter fv = mesh.fv_iter(*face); fv.is_valid(); ++fv ){
			os<<" "<<(*fv).idx();
		}
		os<<std::endl;
	}


	int ixp=0;
	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
		os<<"3";
		for (MyMesh::FaceVertexIter fv = mesh.fv_iter(*face); fv.is_valid(); ++fv ){
			os<<" "<<mesh.n_vertices()+ixp;
			ixp++;
		}
		os<<std::endl;
	}



	os<<"CELL_TYPES "<<2*mesh.n_faces()<<std::endl;
	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
		os<<"5"<<std::endl;
	}	
	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
		os<<"5"<<std::endl;
	}			
	os<<"CELL_DATA "<<2*mesh.n_faces()<<std::endl;
	os<<"COLOR_SCALARS faceColor 4"<<std::endl;	

	bool stripped_mesh;
	if (OpenMesh::hasProperty<OpenMesh::FaceHandle, FaceProp>(mesh, "face_props")) {
		stripped_mesh=false;		
	}
	else{
		stripped_mesh=true;
	}

	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto face_props_stripped =  OpenMesh::FProp<int>(mesh, "face_type");

	int face_type, nft;
	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
		//get face type
		if (stripped_mesh){
			face_type = face_props_stripped[*face];
		}
		else{
			face_type = face_props[*face].face_type;			
		}
	//std::cout<<"strp? "<<stripped_mesh<<std::endl;	
	//std::cout<<"face type  "<<face_type<<std::endl;	
		if (face_type>=0){	
			os<<red[face_type]<<" "<<green[face_type]<<" "<<blue[face_type]<<" "<<alpha[face_type]<<std::endl;
		}
		else{
			os<<"1.0 1.0 0.0 1.0"<<std::endl;

			/* set it back to color negative face types the same
			nft=-(face_type+2);
			os<<red[nft]<<" "<<green[nft]<<" "<<blue[nft]<<" "<<alpha[nft]<<std::endl;		*/	
		}

	}		


	for (MyMesh::FaceIter face = mesh.faces_sbegin(); face!=mesh.faces_end(); ++face){
		//get face type
		if (stripped_mesh){
			face_type = face_props_stripped[*face];
		}
		else{
			face_type = face_props[*face].face_type;			
		}
	//std::cout<<"strp? "<<stripped_mesh<<std::endl;	
	//std::cout<<"face type  "<<face_type<<std::endl;


///!!!! grey inner face    os<<"0.1 0.1 0.1 1.0"<<std::endl;	
		//---------------same colored inner face----------------------
		if (face_type>=0){	
			os<<red[face_type]<<" "<<green[face_type]<<" "<<blue[face_type]<<" "<<alpha[face_type]<<std::endl;
		}
		else{
			os<<"1.0 1.0 0.0 1.0"<<std::endl;		
		}
		//--------------------------------------

//		os<<"1.0 1.0 0.0 1.0"<<std::endl;	
/////		if (face_type>=0){	
/////
/////			//os<<red[face_type]*1.5<<" "<<green[face_type]*1.5<<" "<<blue[face_type]*1.5<<" "<<alpha[face_type]<<std::endl;
/////			//os<<"1.0 0.0 0.0 1.0"<<std::endl;
/////			os<<"0.1 0.1 0.1 1.0"<<std::endl;
/////		}
/////		else{
/////			os<<"1.0 1.0 0.0 1.0"<<std::endl;
/////		}

	}	


	os.close();
	return 1;
}

//convert om file to VTK grid format
int convert_om_to_VTK_snapshots(std::string om_filename, std::string vtk_base, std::string vtk_conf_file){

	std::ifstream omfile(om_filename, std::ios::out | std::ios::binary);

 	//alernatingly read a timestep (long) and a mesh (.om) objects

	const std::string _ext=".OM";
	OpenMesh::IO::Options ropt, wopt;	
	ropt+=OpenMesh::IO::Options::FaceColor; 	
	long timestep;
	MyMesh mesh;

	auto face_props_stripped =  OpenMesh::FProp<int>(mesh, "face_type");	
	while (omfile >> timestep){
   		
   		OpenMesh::IO::read_mesh(mesh, omfile, _ext, ropt);
   		dump_vtk_snapshot(mesh,  vtk_base + std::to_string(timestep) +".vtk", vtk_conf_file);

	}
	omfile.close();
	return 1;	

}


void print_mesh(MyMesh & mesh){
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	int edge_type;
	for (MyMesh::FaceIter f_it = mesh.faces_sbegin(); f_it != mesh.faces_end(); ++f_it){
		std::cout<<"face "<< *f_it<<" face_type "<<face_props[*f_it].face_type;
		for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(*f_it); fh_it.is_valid(); ++fh_it){	
			edge_type = edge_props[*fh_it].edge_type;	
			//std::cout<<"e_b"<< e_b<<" theta0 "<<theta0;
			std::cout<<" edge_type "<<edge_type;//<<std::endl;
		}	
	std::cout<<std::endl; 
}
return;
}


/*void print_mesh_stripped(MyMesh & mesh){
	auto edge_props = OpenMesh::HProp<int>(mesh, "edge_type");
	auto face_props = OpenMesh::FProp<int>(mesh, "face_type");
	int edge_type;
	for (MyMesh::FaceIter f_it = mesh.faces_sbegin(); f_it != mesh.faces_end(); ++f_it){
		std::cout<<"face "<< *f_it<<" face_type "<<face_props[*f_it];
		for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(*f_it); fh_it.is_valid(); ++fh_it){	
			edge_type = edge_props[*fh_it];	
			//std::cout<<"e_b"<< e_b<<" theta0 "<<theta0;
			std::cout<<" edge_type "<<edge_type;//<<std::endl;
		}	
	std::cout<<std::endl; 
}
return;
}*/

void create_single_subunit_om(std::string om_filename, int face_type, int edge_type1, int edge_type2, int edge_type3){
	MyMesh mesh;
	mesh.request_face_status();
  	mesh.request_edge_status();
  	mesh.request_vertex_status();    
  	mesh.request_halfedge_status();
	mesh.request_face_normals();

	auto face_props_bare = OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_bare = OpenMesh::HProp<int>(mesh, "edge_type");



	std::vector<MyMesh::VertexHandle> vhandles;
	MyMesh::FaceHandle face;
	MyMesh::HalfedgeHandle he;
/*	vhandles.push_back(mesh.add_vertex(MyMesh::Point(0,0,0)));
	vhandles.push_back(mesh.add_vertex(MyMesh::Point(0,1,0)));
	vhandles.push_back(mesh.add_vertex(MyMesh::Point(1,0,0)));	
*/
	vhandles.push_back(mesh.add_vertex(MyMesh::Point(-0.5,-0.289,0)));
	vhandles.push_back(mesh.add_vertex(MyMesh::Point(0.0, 0.577,0)));
	vhandles.push_back(mesh.add_vertex(MyMesh::Point(0.5, -0.289,0)));	

	face = mesh.add_face(vhandles);	

	MyMesh::Point n = mesh.calc_normal(face);
	std::cout<<"face normal nx, ny, nz= "<<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl;

	face_props_bare[face] = face_type;
	he = mesh.halfedge_handle(face);
	edge_props_bare[he] = edge_type1;
	he = mesh.next_halfedge_handle(he);
	edge_props_bare[he] = edge_type2;	
	he = mesh.next_halfedge_handle(he);
	edge_props_bare[he] = edge_type3;		



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


	//read the mesh here from .om
	const std::string _ext=".OM";
	//OpenMesh::IO::Options ropt, wopt;
	//ropt+=OpenMesh::IO::Options::FaceColor; 	
	std::ofstream omfile(om_filename, std::ios::out | std::ios::binary);
	if (!omfile.is_open()) std::cout<<"create_single_subunit_om() could not open om file"<<std::endl;
	long timestep=0;
	omfile << timestep;
	if (!OpenMesh::IO::write_mesh(mesh, omfile, _ext)) 			
	{
		std::cerr << "init config write error\n";
    	exit(1);
	}		

/*	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
		face_props_bare[*fit] = face_props[*fit].face_type;
	}

	for (MyMesh::HalfedgeIter hit = mesh.halfedges_sbegin(); hit!=mesh.halfedges_end(); ++hit){
		edge_props_bare[*hit] = edge_props[*hit].edge_type;
	}	
*/
	std::cout<<"face_type "<<face_type<<std::endl;
	return;	
}

void create_hexamer_subunit_om(std::string om_filename, int face_type, int edge_type1, int edge_type2, int edge_type3){
	MyMesh mesh;
	mesh.request_face_status();
  	mesh.request_edge_status();
  	mesh.request_vertex_status();    
  	mesh.request_halfedge_status();
	mesh.request_face_normals();

	auto face_props_bare = OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_bare = OpenMesh::HProp<int>(mesh, "edge_type");



	std::vector<MyMesh::VertexHandle> all_vhandles, vhandles;
	MyMesh::FaceHandle face;
	MyMesh::HalfedgeHandle he;
/*	vhandles.push_back(mesh.add_vertex(MyMesh::Point(0,0,0)));
	vhandles.push_back(mesh.add_vertex(MyMesh::Point(0,1,0)));
	vhandles.push_back(mesh.add_vertex(MyMesh::Point(1,0,0)));	
*/
	all_vhandles.push_back(mesh.add_vertex(MyMesh::Point(0, 0, 0)));
	all_vhandles.push_back(mesh.add_vertex(MyMesh::Point(1, 0, 0)));
	all_vhandles.push_back(mesh.add_vertex(MyMesh::Point(0.5, 0.866025404, 0)));	
	all_vhandles.push_back(mesh.add_vertex(MyMesh::Point(-0.5, 0.866025404, 0)));
	all_vhandles.push_back(mesh.add_vertex(MyMesh::Point(-1, 0, 0)));
	all_vhandles.push_back(mesh.add_vertex(MyMesh::Point(-0.5, -0.866025404, 0)));		
	all_vhandles.push_back(mesh.add_vertex(MyMesh::Point(0.5, -0.866025404, 0)));	

	MyMesh::Point n; 
	for (unsigned int i=1; i<7; i++){
		vhandles.clear();
		vhandles.push_back(all_vhandles[0]);
		vhandles.push_back(all_vhandles[i%6+1]);
		vhandles.push_back(all_vhandles[i]);

		face = mesh.add_face(vhandles);	
		n = mesh.calc_normal(face);
		std::cout<<"face normal nx, ny, nz= "<<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl;

		face_props_bare[face] = face_type;
		he = mesh.halfedge_handle(face);
		edge_props_bare[he] = edge_type1;
		he = mesh.next_halfedge_handle(he);
		edge_props_bare[he] = edge_type2;	
		he = mesh.next_halfedge_handle(he);
		edge_props_bare[he] = edge_type3;		

	}

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


	//read the mesh here from .om
	const std::string _ext=".OM";
	//OpenMesh::IO::Options ropt, wopt;
	//ropt+=OpenMesh::IO::Options::FaceColor; 	
	std::ofstream omfile(om_filename, std::ios::out | std::ios::binary);
	if (!omfile.is_open()) std::cout<<"create_hexamer_subunit_om() could not open om file"<<std::endl;
	long timestep=0;
	omfile << timestep;
	if (!OpenMesh::IO::write_mesh(mesh, omfile, _ext)) 			
	{
		std::cerr << "init config write error\n";
    	exit(1);
	}		

/*	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
		face_props_bare[*fit] = face_props[*fit].face_type;
	}

	for (MyMesh::HalfedgeIter hit = mesh.halfedges_sbegin(); hit!=mesh.halfedges_end(); ++hit){
		edge_props_bare[*hit] = edge_props[*hit].edge_type;
	}	
*/
	std::cout<<"face_type "<<face_type<<std::endl;
	return;	
}


void convert_dat_to_om(std::string vertices_dat, std::string faces_dat, std::string edges_dat, std::string omfile){
	MyMesh mesh;
  	mesh.request_face_status();
  	mesh.request_edge_status();
  	mesh.request_vertex_status();    
  	mesh.request_halfedge_status();
	mesh.request_face_normals();

	auto face_props_bare = OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_bare = OpenMesh::HProp<int>(mesh, "edge_type");	

	//std::cout<<"run_debug()"<<std::endl;
	//std::ifstream vertices("/home/btyukodi/assembly_openmesh/sandbox/carlos/vertices.dat");
	//std::ifstream faces("/home/btyukodi/assembly_openmesh/sandbox/carlos/faces.dat");
	//std::ifstream edges("/home/btyukodi/assembly_openmesh/sandbox/carlos/edges.dat");	

	std::ifstream vertices(vertices_dat);
	std::ifstream faces(faces_dat);
	std::ifstream edges(edges_dat);		

	double x, y, z;
	long ix, ix1, ix2, ix3;
	int face_type, edge_type;
	long t=0;
	std::vector<MyMesh::VertexHandle> vhandles;
	while(vertices >> ix){
		vertices >> x;
		vertices >> y;
		vertices >> z;

		//std::cout<<ix<<" "<<x<<" "<<y<<" "<<z<<std::endl;
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
		face_props_bare[fh] = face_type-1;
		//std::cout<<ix1<<" "<<ix2<<" "<<ix3<<" "<<face_type<<std::endl;
	}
	faces.close();

	MyMesh::HalfedgeHandle he;
	while (edges >> ix){
		edges >> ix1;
		edges >> ix2;
		edges >> edge_type;
		he = mesh.find_halfedge(vhandles[ix1-1], vhandles[ix2-1]);
		edge_props_bare[he] = edge_type-1;
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

	//std::ofstream os("/home/btyukodi/assembly_openmesh/sandbox/carlos/minimal.om", std::ios::out | std::ios::binary);
	std::ofstream os(omfile, std::ios::out | std::ios::binary);	
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
	return;
}

//only works on stripped mesh
void dump_dat_snapshot(MyMesh mesh, std::string vertices_dat, std::string faces_dat, std::string edges_dat){
	//std::cout<<"dump dat snapshot"<< mesh.n_faces() << std::endl;
   	std::ofstream os_vertices( vertices_dat );
   	std::ofstream os_faces( faces_dat );
   	std::ofstream os_edges( edges_dat );   	
	os_vertices.precision(6);	
	os_vertices.setf(std::ios::fixed);
	os_vertices.setf(std::ios::showpoint);

	os_faces.precision(6);	
	os_faces.setf(std::ios::fixed);
	os_faces.setf(std::ios::showpoint);

	os_edges.precision(6);	
	os_edges.setf(std::ios::fixed);
	os_edges.setf(std::ios::showpoint);		

	MyMesh::Point com = MyMesh::Point(0,0,0);
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		com+=mesh.point(*v) / mesh.n_vertices();
	}
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		mesh.set_point(*v, mesh.point(*v) - com );
	}

	std::vector<MyMesh::VertexHandle> vhandles;
	std::vector<int> idx;
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		vhandles.push_back(*v);
		idx.push_back((*v).idx());
	}

	//now create a map from internal indices to file indices
	std::map<int, int> idx_i;
	for (unsigned int i=0; i<vhandles.size(); i++){
		idx_i[ idx[i] ] = i+1;
	}

	//write the vertex file
	MyMesh::VertexHandle v;
	for (unsigned int i=0; i<vhandles.size(); i++){
		v = vhandles[i];
		os_vertices<<i+1<<" "<<mesh.point(v)[0]<<" "<<mesh.point(v)[1]<<" "<<mesh.point(v)[2]<<std::endl;
		//std::cout<<"vertex";
	}	

	//write the face file
	auto face_props_stripped =  OpenMesh::FProp<int>(mesh, "face_type");
	int fix=1;
	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
		os_faces<<fix;

		for (MyMesh::FaceVertexIter fv = mesh.fv_iter(*fit); fv.is_valid(); ++fv ){
			os_faces<<" "<<idx_i[(*fv).idx()];
		}
		os_faces<<" "<< face_props_stripped[*fit]+1 <<std::endl;
		fix++;
	}

	//write the edge file
	auto edge_props_stripped = OpenMesh::HProp<int>(mesh, "edge_type");
	//skip boundary halfedges
	int hix=1;

	for (MyMesh::HalfedgeIter hit = mesh.halfedges_sbegin(); hit!=mesh.halfedges_end(); ++hit){
		if (!mesh.is_boundary(*hit)){
			os_edges<<hix<<" "<<idx_i[mesh.from_vertex_handle(*hit).idx()]<<" "<<idx_i[mesh.to_vertex_handle(*hit).idx()]<<" "<<edge_props_stripped[*hit]+1<<std::endl;
			hix++;
			//std::cout<<"edge";
		}
		

		
	}
	os_vertices.close();
	os_faces.close();
	os_edges.close();

	return;
}


int convert_om_to_dat_snapshots(std::string om_filename, std::string dat_base){
	//std::cout<<"reading";
	std::ifstream omfile(om_filename, std::ios::out | std::ios::binary);

 	//alernatingly read a timestep (long) and a mesh (.om) objects

	const std::string _ext=".OM";
	OpenMesh::IO::Options ropt, wopt;	
	ropt+=OpenMesh::IO::Options::FaceColor; 	
	long timestep;
	MyMesh mesh;
	auto face_props_stripped =  OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_stripped = OpenMesh::HProp<int>(mesh, "edge_type");	
	while (omfile >> timestep){
   		
   		OpenMesh::IO::read_mesh(mesh, omfile, _ext, ropt);
   		dump_dat_snapshot(mesh,  dat_base + std::to_string(timestep) +"_vertices.dat", dat_base + std::to_string(timestep) +"_faces.dat", dat_base + std::to_string(timestep) +"_edges.dat");

	}
	omfile.close();
	return 1;

}


int convert_om_to_dat_snapshots(std::string om_filename, std::string dat_base, long long om_byte_position){
	//std::cout<<"reading";
	std::ifstream omfile(om_filename, std::ios::out | std::ios::binary);
	omfile.seekg(om_byte_position, omfile.beg); 
 	//alernatingly read a timestep (long) and a mesh (.om) objects

	const std::string _ext=".OM";
	OpenMesh::IO::Options ropt, wopt;	
	ropt+=OpenMesh::IO::Options::FaceColor; 	
	long timestep;
	MyMesh mesh;
	auto face_props_stripped =  OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_stripped = OpenMesh::HProp<int>(mesh, "edge_type");	
	omfile >> timestep;
  		
   	OpenMesh::IO::read_mesh(mesh, omfile, _ext, ropt);
   	dump_dat_snapshot(mesh,  dat_base + std::to_string(timestep) +"_vertices.dat", dat_base + std::to_string(timestep) +"_faces.dat", dat_base + std::to_string(timestep) +"_edges.dat");

	
	omfile.close();
	return 1;

}

//need many inputs:
//mesh, umbrella window, wall, timestep, rng, ...
int create_checkpoint(MyMesh & mesh, UmbrellaWindow & uw, long timestep, std::string ckp_filename){
	std::ofstream checkpoint_file(ckp_filename, std::ios::out | std::ios::binary);

	//checkpoint_file.write((char *)&uw, sizeof(UmbrellaWindow));
	checkpoint_file.write((char *)&(uw.spring_const), sizeof(double) );
	checkpoint_file.write((char *)&(uw.N0), sizeof(double) );

	dump_om(mesh, timestep, checkpoint_file); //checkpoint_file<<t<<mesh			


	checkpoint_file.close();

	return 1;
}

int load_from_checkpoint(MyMesh & mesh, UmbrellaWindow & uw, long & timestep, std::string ckp_filename, std::string input_filename){
	std::ifstream checkpoint_file(ckp_filename, std::ios::in | std::ios::binary);


	checkpoint_file.read((char *)&(uw.spring_const), sizeof(double) );
	checkpoint_file.read((char *)&(uw.N0), sizeof(double) );
	long long current_position = checkpoint_file.tellg();

	checkpoint_file>>timestep;
	init_mesh_from_om(mesh, ckp_filename, current_position, input_filename);						
	
	checkpoint_file.close();	

	return 1;
}


int create_checkpoint_wall(MyMesh & mesh, Wall & wall, long timestep, std::string ckp_filename){
	std::ofstream checkpoint_file(ckp_filename, std::ios::out | std::ios::binary);

	//checkpoint_file.write((char *)&uw, sizeof(UmbrellaWindow));
	checkpoint_file.write((char *)&(wall.hardness), sizeof(double) );
	checkpoint_file.write((char *)&(wall.z_position_1), sizeof(double) );
	checkpoint_file.write((char *)&(wall.z_position_2), sizeof(double) );
	checkpoint_file.write((char *)&(wall.amplitude), sizeof(double) );

	dump_om(mesh, timestep, checkpoint_file); //checkpoint_file<<t<<mesh			


	checkpoint_file.close();

	return 1;
}

int load_from_checkpoint_wall(MyMesh & mesh, Wall & wall, long & timestep, std::string ckp_filename, std::string input_filename){
	std::ifstream checkpoint_file(ckp_filename, std::ios::in | std::ios::binary);


	checkpoint_file.read((char *)&(wall.hardness), sizeof(double) );
	checkpoint_file.read((char *)&(wall.z_position_1), sizeof(double) );
	checkpoint_file.read((char *)&(wall.z_position_2), sizeof(double) );
	checkpoint_file.read((char *)&(wall.amplitude), sizeof(double) );
	long long current_position = checkpoint_file.tellg();

	checkpoint_file>>timestep;
	init_mesh_from_om(mesh, ckp_filename, current_position, input_filename);						
	
	checkpoint_file.close();	

	return 1;
}

bool file_exists (std::string name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


void dump_error(std::string error_msg, std::string filename){
	std::ofstream os(filename );
	os<<error_msg<<std::endl;
	os.close();
	return;

}

//appends the snapshot from omfile at position to omfile_to_append_to
long long append_snapshot(std::string om_filename, long long om_byte_position, std::string omfile_to_append_to){

	std::ofstream os_out(omfile_to_append_to, std::ios::out | std::ios::binary | std::ios::app);

	std::ifstream omfile(om_filename, std::ios::out | std::ios::binary);
	omfile.seekg(om_byte_position, omfile.beg); 
 	//alernatingly read a timestep (long) and a mesh (.om) objects

	const std::string _ext=".OM";
	OpenMesh::IO::Options ropt, wopt;	
	ropt+=OpenMesh::IO::Options::FaceColor; 	
	long timestep;
	MyMesh mesh;
	auto face_props_stripped =  OpenMesh::FProp<int>(mesh, "face_type");
	auto edge_props_stripped = OpenMesh::HProp<int>(mesh, "edge_type");	
	omfile >> timestep;
  		
   	OpenMesh::IO::read_mesh(mesh, omfile, _ext, ropt);
	omfile.close();	



	long long newpos;
	newpos=dump_om_stripped(mesh, om_byte_position, os_out);

	os_out.close();

	return newpos;
}

//passed by copy, not reference here
/*void dump_mesh(MyMesh mesh, std::string om_filename){

	mesh = get_stripped_mesh(mesh);

	mesh.request_face_colors();	
    MyMesh::Color clr;	
    clr[0] = 252;
    clr[1] = 168;
    clr[2] = 3;
	for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++) {    
		mesh.set_color(*f_it, clr);
	}

	OpenMesh::IO::Options ropt, wopt;
	wopt+=OpenMesh::IO::Options::FaceColor; 
	wopt+=OpenMesh::IO::Options::Custom; 

	//std::cout<<"rm has property "<<OpenMesh::hasProperty<OpenMesh::FaceHandle, FaceProp>(mesh, "face_props")<<std::endl;	

//will need .om to restore halfedge types
try
  {
    if ( !OpenMesh::IO::write_mesh(mesh, om_filename) )
    {
      std::cerr << "Cannot write mesh to file 'output.om'" << std::endl;
      return;
    }
  }
  catch( std::exception& x )
  {
    std::cerr << x.what() << std::endl;
    return;
  }


try
  {
    if ( !OpenMesh::IO::write_mesh(mesh, "output.off", wopt) )
    {
      std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
      return;
    }
  }
  catch( std::exception& x )
  {
    std::cerr << x.what() << std::endl;
    return;
  }  


try
  {
    if ( !OpenMesh::IO::write_mesh(mesh, "output.ply", wopt) )
    {
      std::cerr << "Cannot write mesh to file 'output.ply'" << std::endl;
      return;
    }
  }
  catch( std::exception& x )
  {
    std::cerr << x.what() << std::endl;
    return;
  }    


try
  {
    if ( !OpenMesh::IO::write_mesh(mesh, "output.obj", wopt) )
    {
      std::cerr << "Cannot write mesh to file 'output.obj'" << std::endl;
      return;
    }
  }
  catch( std::exception& x )
  {
    std::cerr << x.what() << std::endl;
    return;
  }  

try
  {
    if ( !OpenMesh::IO::write_mesh(mesh, "output.stl", wopt) )
    {
      std::cerr << "Cannot write mesh to file 'output.stl'" << std::endl;
      return;
    }
  }
  catch( std::exception& x )
  {
    std::cerr << x.what() << std::endl;
    return;
  }    

try
  {
    if ( !OpenMesh::IO::write_mesh(mesh, "output.vtk", wopt) )
    {
      std::cerr << "Cannot write mesh to file 'output.vtk'" << std::endl;
      return;
    }
  }
  catch( std::exception& x )
  {
    std::cerr << x.what() << std::endl;
    return;
  }      

  MyMesh readmesh;
  auto face_props_readmesh = OpenMesh::FProp<int>(readmesh, "face_type");
  auto edge_props_readmesh = OpenMesh::HProp<int>(readmesh, "edge_type");
 

  try
  {
    if (OpenMesh::IO::read_mesh( readmesh, "output.ply" ))
      std::cout << "  ok\n";
    else
    {
      std::cout << "  failed!\n";
      return;
    }
    readmesh.property_stats(std::cout);
	std::cout<<"readmesh has property "<<OpenMesh::hasProperty<OpenMesh::FaceHandle, int>(readmesh, "face_type")<<std::endl;	    
  }
  catch( std::exception &x )
  {
    std::cerr << x.what() << std::endl;
    return;
  }    


	int edge_type;
	float e_b,theta0;
	for (MyMesh::FaceIter f_it = readmesh.faces_sbegin(); f_it != readmesh.faces_end(); ++f_it){
		std::cout<<"face "<< *f_it<<" ";
		for (MyMesh::FaceHalfedgeIter fh_it = readmesh.fh_iter(*f_it); fh_it.is_valid(); ++fh_it){	
			edge_type = edge_props_readmesh[*fh_it];	
			//std::cout<<"e_b"<< e_b<<" theta0 "<<theta0;
			std::cout<<" edge_type "<<edge_type;//<<std::endl;


		}	
		std::cout<<std::endl;
 
	}



*/