

/*
struct Excluder
{
	int f;
};
*/

/*double sqr_dist(MyMesh & mesh, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2){
	MyMesh::Point v1p = mesh.point(v1);
	MyMesh::Point v2p = mesh.point(v2);	
	std::cout<<"norm 1, 2 "<<(v1p-v2p).sqrnorm()<<" "<<(v1p[0]-v2p[0])*(v1p[0]-v2p[0]) +  (v1p[1]-v2p[1])*(v1p[1]-v2p[1]) + (v1p[2]-v2p[2])*(v1p[2]-v2p[2])<<std::endl;

	return (v1p[0]-v2p[0])*(v1p[0]-v2p[0]) +  (v1p[1]-v2p[1])*(v1p[1]-v2p[1]) + (v1p[2]-v2p[2])*(v1p[2]-v2p[2]);
}*/



/*
1. Add two more parameters:
- param_filename
- init_config_filename - has to be .om for the edge/face types; optional

2. Watch out for copying from the prototype mesh; get the right subunit with right rotation
*/

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



void init_mesh_pp(MyMesh & mesh){
  	mesh.request_face_status();
  	mesh.request_edge_status();
  	mesh.request_vertex_status();    
  	mesh.request_halfedge_status();
	mesh.request_face_normals();  

	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/single_sub.off")) 	
//	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/two_sub.off")) 	
//	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/hexamer.off")) 			
	{
		std::cerr << "read error\n";
    	exit(1);
	}			

	//read_set_mesh_params(mesh, "/home/btyukodi/assembly_openmesh/sandbox/input_pretty.json");
	read_set_mesh_params(mesh, "/home/btyukodi/assembly_openmesh/sandbox/input_carlos2.json");	

	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto proto_edge_props = OpenMesh::HProp<HalfedgeProp>((*mesh_props).prototypes_mesh, "edge_props");


	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");

	MyMesh::FaceHandle fh;
	MyMesh::HalfedgeHandle anchor_halfedge;
	for (MyMesh::FaceIter face=mesh.faces_begin(); face!=mesh.faces_end(); ++face){
		fh = *face;
		//for now, find any edge as anchor
		anchor_halfedge = mesh.halfedge_handle(fh);
		clone_properties_from_prototype(mesh, 0, fh, anchor_halfedge);

		mesh.calc_face_centroid(fh, face_props[fh].COM );		
	}
	update_full_neighbor_list(mesh);	

	int edge_type;
	float e_b,theta0;
	for (MyMesh::FaceIter f_it = mesh.faces_sbegin(); f_it != mesh.faces_end(); ++f_it){
		std::cout<<"face "<< *f_it<<" ";
		for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(*f_it); fh_it.is_valid(); ++fh_it){
			e_b = edge_props[*fh_it].e_b[0];
			theta0 = edge_props[*fh_it].theta0[0];		
			edge_type = edge_props[*fh_it].edge_type;	
			std::cout<<"e_b"<< e_b<<" theta0 "<<theta0;
			std::cout<<" edge_type "<<edge_type;//<<std::endl;


		}	
		std::cout<<std::endl;
 
	}	
}


//this function will go somewhere else, now it's just used for debugging
void init_mesh_dont_delete(MyMesh & mesh){


  	mesh.request_face_status();
  	mesh.request_edge_status();
  	mesh.request_vertex_status();    
  	mesh.request_halfedge_status();
	mesh.request_face_normals();  		

	read_set_mesh_params(mesh, "/home/btyukodi/assembly_openmesh/sandbox/input_pretty.json");	
	//read a cracked mesh from 
//	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/hex_sheet_large.off")) 
//	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/wedge.off")) 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/single_sub.off")) 	
//	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/two_sub.off")) 	

//	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/hexamer.off")) 			
	{
		std::cerr << "read error\n";
    	exit(1);
	}

	//set all edge types equal for now
	//MyMesh::HalfedgeHandle hedge = mesh.halfedge_handle(0);
	//MyMesh::HalfedgeHandle hedge2 = mesh.halfedge_handle(1);	
	/*
	OpenMesh::HPropHandleT<HalfedgeProp> edge_props;
	mesh.add_property(edge_props,  "edge_props");
	mesh.property(edge_props).set_persistent(true);
	*/

	MyMesh prototypes_mesh;
	std::vector<MyMesh::VertexHandle> vhandles;
	std::vector<MyMesh::FaceHandle> prototype_faces;
	std::vector<MyMesh::HalfedgeHandle> anchor_halfedges;
	MyMesh::FaceHandle proto_face;

	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	(*mesh_props).kT = 1.0;
	(*mesh_props).d_max = 0.2;
	(*mesh_props).R_add = 0.2;	
	(*mesh_props).l_fuse = 0.5;	
	(*mesh_props).L_neighbor = 1.2;
	(*mesh_props).k_fusion = 0.001;//0.01;//0.01;		
	(*mesh_props).k_fusion2 = 0.0001;//0.001;	
	(*mesh_props).k_fusion_edge = 0.01;			

	//add and set prototype faces here
	vhandles.push_back(prototypes_mesh.add_vertex(MyMesh::Point(0,0,0)));
	vhandles.push_back(prototypes_mesh.add_vertex(MyMesh::Point(0,1,0)));
	vhandles.push_back(prototypes_mesh.add_vertex(MyMesh::Point(1,0,0)));	
	proto_face = prototypes_mesh.add_face(vhandles);


	auto proto_edge_props = OpenMesh::HProp<HalfedgeProp>(prototypes_mesh, "edge_props");
	auto proto_face_props = OpenMesh::FProp<FaceProp>(prototypes_mesh, "face_props");

	proto_face_props[proto_face].face_type = 0;
	proto_face_props[proto_face].mu = -3.0;
	proto_face_props[proto_face].R_exc = 0.2;	
	proto_face_props[proto_face].k_insertion = 0.0;		
	proto_face_props[proto_face].number_of_rotational_configs = 1;	



	std::vector<float> e_b, bend_mod, theta0;
	e_b.clear();
	bend_mod.clear();
	theta0.clear();

	e_b.push_back(-8.0);
	bend_mod.push_back(300.0);
	theta0.push_back(0.41);	

	e_b.push_back(100.0);
	bend_mod.push_back(300.0);
	theta0.push_back(0.41);

	MyMesh::HalfedgeHandle halfedge = prototypes_mesh.halfedge_handle(proto_face);
	for (int i=0; i<3; i++){
		proto_edge_props[halfedge].l0 = 1.0;
		proto_edge_props[halfedge].edge_type=0;
		proto_edge_props[halfedge].stretch_mod=200.0;
		proto_edge_props[halfedge].e_b = e_b;
		proto_edge_props[halfedge].bend_mod = bend_mod;
		proto_edge_props[halfedge].theta0 = theta0;

		halfedge = prototypes_mesh.next_halfedge_handle(halfedge);
	}	


//	(*mesh_props).prototypes_mesh = prototypes_mesh;

	prototype_faces.push_back(proto_face);
//	(*mesh_props).prototype_faces = prototype_faces;

	anchor_halfedges.push_back(halfedge);
//	(*mesh_props).anchor_halfedges = anchor_halfedges;

//------------------------

	//add and set prototype faces here
	vhandles.clear();
	vhandles.push_back(prototypes_mesh.add_vertex(MyMesh::Point(0,0,0)));
	vhandles.push_back(prototypes_mesh.add_vertex(MyMesh::Point(0,-1,0)));
	vhandles.push_back(prototypes_mesh.add_vertex(MyMesh::Point(1,0,0)));	
	proto_face = prototypes_mesh.add_face(vhandles);


	proto_face_props[proto_face].face_type = 1;
	proto_face_props[proto_face].mu = -3.0;
	proto_face_props[proto_face].R_exc = 0.2;
	proto_face_props[proto_face].k_insertion = 0.005;//0.005;			
	proto_face_props[proto_face].number_of_rotational_configs = 3;	


	halfedge = prototypes_mesh.halfedge_handle(proto_face);

	e_b.clear();
	bend_mod.clear();
	theta0.clear();

/*	e_b.push_back(-8.0);
	bend_mod.push_back(300.0);
	theta0.push_back(0.41);	

	e_b.push_back(-8.0);
	bend_mod.push_back(300.0);
	theta0.push_back(0.41); *///(0.41);


	//for (int i=0; i<3; i++){
		proto_edge_props[halfedge].l0 = 1.0;
		proto_edge_props[halfedge].edge_type=1;
		proto_edge_props[halfedge].stretch_mod=200.0;


		e_b.clear();
		bend_mod.clear();
		theta0.clear();
		e_b.insert(e_b.end(), 	   		{100.0,  -8.0,  100.0, 100.0} );
		bend_mod.insert(bend_mod.end(), {0.0,    300.0, 0.0,   0.0});
		theta0.insert(theta0.end(),  	{0.0,    0.41,  0.0,   0.0});

		proto_edge_props[halfedge].e_b = e_b;
		proto_edge_props[halfedge].bend_mod = bend_mod;
		proto_edge_props[halfedge].theta0 = theta0;

		halfedge = prototypes_mesh.next_halfedge_handle(halfedge);

		proto_edge_props[halfedge].l0 = 1.0;
		proto_edge_props[halfedge].edge_type=2;
		proto_edge_props[halfedge].stretch_mod=200.0;

		e_b.clear();
		bend_mod.clear();
		theta0.clear();
		e_b.insert(e_b.end(), 	   		{100.0,  100.0, -8.0,  100.0} );
		bend_mod.insert(bend_mod.end(), {0.0,    0.0,   300.0, 0.0});
		theta0.insert(theta0.end(),  	{0.0,    0.0,   0.41,  0.0});

		proto_edge_props[halfedge].e_b = e_b;
		proto_edge_props[halfedge].bend_mod = bend_mod;
		proto_edge_props[halfedge].theta0 = theta0;

		halfedge = prototypes_mesh.next_halfedge_handle(halfedge);

		proto_edge_props[halfedge].l0 = 1.0;
		proto_edge_props[halfedge].edge_type=3;
		proto_edge_props[halfedge].stretch_mod=200.0;

		e_b.clear();
		bend_mod.clear();
		theta0.clear();
		e_b.insert(e_b.end(), 	   		{100.0,  100.0, 100.0, -8.0 } );
		bend_mod.insert(bend_mod.end(), {0.0,    0.0,   0.0,   300.0});
		theta0.insert(theta0.end(),  	{0.0,    0.0,   0.0,   -0.21*1.8});

		proto_edge_props[halfedge].e_b = e_b;
		proto_edge_props[halfedge].bend_mod = bend_mod;
		proto_edge_props[halfedge].theta0 = theta0;

		halfedge = prototypes_mesh.next_halfedge_handle(halfedge);		



	//}	


	(*mesh_props).prototypes_mesh = prototypes_mesh;

	prototype_faces.push_back(proto_face);
	(*mesh_props).prototype_faces = prototype_faces;

	anchor_halfedges.push_back(halfedge);
	(*mesh_props).anchor_halfedges = anchor_halfedges;

//-----------------------	


	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");

	MyMesh::FaceHandle fh;
	MyMesh::HalfedgeHandle anchor_halfedge;
	for (MyMesh::FaceIter face=mesh.faces_begin(); face!=mesh.faces_end(); ++face){
		fh = *face;
		//for now, find any edge as anchor
		anchor_halfedge = mesh.halfedge_handle(fh);
		clone_properties_from_prototype(mesh, 1, fh, anchor_halfedge);

		mesh.calc_face_centroid(fh, face_props[fh].COM );		
	}
	update_full_neighbor_list(mesh);

/*	for (MyMesh::FaceIter face1 = mesh.faces_sbegin(); face1!=mesh.faces_end(); ++face1){	
		std::cout<<"--- face "<<*face1<<std::endl;
		std::cout<<"neighbors ";
		for (auto & neighbor: face_props[ *face1 ].neighbor_list){
	    		std::cout<<" "<<neighbor;
	    }
		std::cout<<std::endl;
	}
*/
/*
	for (MyMesh::HalfedgeIter hedge=mesh.halfedges_begin(); hedge!=mesh.halfedges_end(); ++hedge){

		edge_props[*hedge].l0 = 1.0;
		edge_props[*hedge].edge_type=0;
		edge_props[*hedge].stretch_mod=100.0;
		edge_props[*hedge].e_b = e_b;
		edge_props[*hedge].bend_mod = bend_mod;
		edge_props[*hedge].theta0 = theta0;

	}

	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	std::vector<float> e_b, bend_mod, theta0;

	e_b.push_back(-5.0);
	bend_mod.push_back(300.0);
	theta0.push_back(0.4);





	for (MyMesh::HalfedgeIter hedge=mesh.halfedges_begin(); hedge!=mesh.halfedges_end(); ++hedge){

		edge_props[*hedge].l0 = 1.0;
		edge_props[*hedge].edge_type=0;
		edge_props[*hedge].stretch_mod=100.0;
		edge_props[*hedge].e_b = e_b;
		edge_props[*hedge].bend_mod = bend_mod;
		edge_props[*hedge].theta0 = theta0;

	}

	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	for (MyMesh::FaceIter face=mesh.faces_begin(); face!=mesh.faces_end(); ++face){
		face_props[*face].face_type = 0;
		face_props[*face].mu = -1.3;		
	}	

	//mesh.copy_all_properties(mesh.face_handle(0) ,(*mesh_props).prototype_faces[0]);
	proto_face_props[(*mesh_props).prototype_faces[0]].mu = -1.1;
	//(*mesh_props).prototypes_mesh.copy_all_properties(mesh.face_handle(0) ,(*mesh_props).prototype_faces[0]);	
	proto_face_props[(*mesh_props).prototype_faces[0]] =  face_props[mesh.face_handle(0)];


/*	mesh.property(edge_props, hedge).l0 = 0.5;
	mesh.property(edge_props, hedge).edge_type = 5;	
*/
	//mesh.copy_all_properties(hedge, hedge2);

	//std::cout<<edge_props[hedge2].l0<<std::endl;

	//return mesh;
}





int mainX2(){
	MyMesh mesh1, mesh2;
	std::vector<MyMesh::VertexHandle> vhandles;
	MyMesh::FaceHandle face1, face2;
	mesh1.request_vertex_status();
	mesh2.request_vertex_status();
		
	auto mesh1_props = OpenMesh::FProp<int>(mesh1, "face_props");
	auto mesh2_props = OpenMesh::FProp<int>(mesh2, "face_props");			

	vhandles.push_back(mesh1.add_vertex(MyMesh::Point(0,0,0)));
	vhandles.push_back(mesh1.add_vertex(MyMesh::Point(0,1,0)));
	vhandles.push_back(mesh1.add_vertex(MyMesh::Point(1,0,0)));	
	face1 = mesh1.add_face(vhandles);
	mesh1_props[face1] = 5;

	vhandles.clear();
	vhandles.push_back(mesh2.add_vertex(MyMesh::Point(0,0,0)));
	vhandles.push_back(mesh2.add_vertex(MyMesh::Point(0,1,0)));
	vhandles.push_back(mesh2.add_vertex(MyMesh::Point(1,0,0)));
	face2 = mesh2.add_face(vhandles);
	mesh2.status(vhandles[0]).set_tagged(true);

	mesh1.copy_all_properties(face1, face2);
	std::cout<<"face1 prop "<<mesh1_props[face1]<<std::endl;
	std::cout<<"face2 prop "<<mesh2_props[face2]<<std::endl;

	return 0;
}





bool test_manifold(MyMesh & mesh){
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
			if (!mesh.is_manifold(*v_it)){
				return false;
			}
		}	
	return true;
}


bool test_coord(MyMesh & mesh, int prev_step){
	/*for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
			if (std::abs( mesh.point(*v_it)[2] )>8.0  ){
				std::cout<<"prev step "<<prev_step<<std::endl;
				return false;
			}
	}*/	
	double l;
	for (MyMesh::EdgeIter e_it = mesh.edges_sbegin(); e_it != mesh.edges_end(); ++e_it){
			l = mesh.calc_edge_length(*e_it);
			if ( (l >1.5) || (l<0.5)  ){
				std::cout<<"prev step "<<prev_step<<" l="<<l<<std::endl;
				return false;
			}
	}		
	return true;
}


//independent com calculator
MyMesh::Point get_com(MyMesh & mesh, MyMesh::FaceHandle face){
	MyMesh::Point com = MyMesh::Point(0,0,0);
	for (MyMesh::FaceVertexIter fv = mesh.fv_iter(face); fv.is_valid(); ++fv ){
		com+=mesh.point(*fv);
	}
	com = com/3.0;
	return com;
}

//check if neighbor lists are up to date and symmetric, and com positions
bool test_neighbor_list(MyMesh & mesh){
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	//loop over faces by index, construct 2D neighbor matrix
	int nf = mesh.n_faces();
	int neighbor_matrix[nf][nf], reference_matrix[nf][nf];
	int id1, id2;
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");		
	double L = (*mesh_props).L_neighbor;	
	MyMesh::Point com1, com2, dr;

	for (int i=0; i<nf; i++){
		for (int j=0; j<nf; j++){
			neighbor_matrix[i][j]=0;
			reference_matrix[i][j]=0;
		}
	}



	for (MyMesh::FaceIter face1 = mesh.faces_sbegin(); face1!=mesh.faces_end(); ++face1){
		id1 = (*face1).idx();	
		com1 = get_com(mesh, *face1);
		//std::cout<<"COM "<<com1<<" "<<mesh.calc_face_centroid(*face1)<<std::endl;		
		for (MyMesh::FaceIter face2 = mesh.faces_sbegin(); face2!=mesh.faces_end(); ++face2){
			if ( *face1!= *face2){				
				id2 = (*face2).idx();
				com2 = get_com(mesh, *face2);
				//std::cout<<"ID "<<id1<<" "<<id2<<" "<<com1<<" "<<com2<<std::endl;

				dr = com1 - com2;
				if ( abs(dr[0])<L && abs(dr[1])<L && abs(dr[2])<L ){
					reference_matrix[id1][id2]=1;
					
				}
			}
				
		}
	}


    for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
    	id1 = (*fit).idx();
		for (auto & neighbor: face_props[*fit].neighbor_list){
			id2 = neighbor.idx();
			neighbor_matrix[id1][id2]++;
		}
    }

	for (int i=0; i<nf; i++){
		for (int j=0; j<nf; j++){
			if (neighbor_matrix[i][j]>1){
				std::cout<<" neighbor matrix>1"<<neighbor_matrix[i][j]<<std::endl;
			}
			if (neighbor_matrix[i][j] != neighbor_matrix[j][i] ){
				std::cout<<" neighbor matrix not symmetric"<<neighbor_matrix[i][j]<<std::endl;
			}		
			//this is a bad test; neighbor lists are not up to date but the idea is to still avoid overlap
			//if (neighbor_matrix[i][j] !=reference_matrix[i][j] ){
			//	std::cout<<" reference matrix != neighbor_matrix "<<i<<" "<<j<<" "<<neighbor_matrix[i][j]<<std::endl;
			//}					
		}
	}

//MyMesh::FaceHandle face = mesh.face_handle(0);
//std::cout<<face.idx()<<std::endl;
}

double test_com_updateness(MyMesh & mesh){
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");		
	double com_err = 0.0;
	MyMesh::Point com;
	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
		com = mesh.calc_face_centroid(*fit);
		com_err+=(com - face_props[*fit].COM).norm() ;
	}
	return com_err;
}

double print_angle_potential(MyMesh & mesh){
	MyMesh::VertexHandle v0 = mesh.vertex_handle(0);
	MyMesh::VertexHandle v1 = mesh.vertex_handle(1);
	MyMesh::VertexHandle v2 = mesh.vertex_handle(2);
	MyMesh::VertexHandle v3 = mesh.vertex_handle(3);
	MyMesh::Point com = (mesh.point(v0) + mesh.point(v2))*0.5;
	mesh.set_point(v0, mesh.point(v0)-com);
	mesh.set_point(v1, mesh.point(v1)-com);
	mesh.set_point(v2, mesh.point(v2)-com);
	mesh.set_point(v3, mesh.point(v3)-com);	
	MyMesh::Point ex = (mesh.point(v0) - mesh.point(v2)).normalize();
	MyMesh::Point p1 = mesh.point(v1);
	MyMesh::EdgeHandle edge = find_edge(mesh, v0, v2);
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");


	double theta=-M_PI;
	double dtheta=0.05;
	rotatevec(p1, ex, theta);
	mesh.set_point(v1, p1);		

	for (int i=0; i< 2*abs(theta)/dtheta; i++){
		std::cout<<mesh.calc_dihedral_angle(edge)<<" "<<edge_bending_energy(mesh, edge, edge_props)<<std::endl;	

		//std::cout<<theta+i*dtheta<<" "<<mesh.calc_dihedral_angle(edge)<<std::endl;	
		p1 = mesh.point(v1);
		rotatevec(p1, ex, dtheta);
		mesh.set_point(v1, p1);	
	}


}

/*
//will need to add umbrella window
void propagate(MyMesh & mesh, long t_init, long t_final, RunParameters rp, std::ostream& outfile_om,  std::ostream& outfile_data){

    int nmoves = 10;
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



	for (long t=t_init; t<t_final; t++){
		if (t % rp.dtsave ==0){		
			std::cout<<"t = "<<t<<" n_faces="<<mesh.n_faces()<<std::endl;	
			
			fpos = dump_om(mesh, t, outfile_om);
			dump_data(mesh, t, outfile_data, fpos);
			//dump_om(mesh, "omtest"+std::to_string(t)+".om");
		}

		
		for (int n=0; n<mesh.n_vertices(); n++){
			attempt_move(mesh);
		}


		which_move = rand() % nmoves;
		accepted=false;
		switch (which_move){
			case 0:
				k_fusion = (*mesh_props).k_fusion; 
				fusion_vectors = get_type1_fusion_triplets(mesh, l_fuse);	
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion rate! - type1_fusion"<<std::endl;
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){	
					accepted = attempt_type1_fusion(mesh);
				}
				//if (accepted) std::cout<<"type1 fusion"<<std::endl;
				break;

			case 1:
				k_fusion = (*mesh_props).k_fusion; 
				fission_vectors = get_type1_fission_pairs(mesh);	
    			v1 = std::get<0>(fission_vectors);					
				p_propose = k_fusion * v1.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion rate! - type1_fission"<<std::endl;				
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){	
					accepted = attempt_type1_fission(mesh);	
				}
				//if (accepted) std::cout<<"type1 fission"<<std::endl;
				break;

			case 2:		
				k_fusion = (*mesh_props).k_fusion2; 
				fusion_vectors = get_type2_fusion_triplets(mesh, l_fuse);	
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fusion"<<std::endl;				
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){
					accepted = attempt_type2_fusion(mesh);
				}
				//if (accepted) std::cout<<"type2 fusion "<<p_propose<<std::endl;		
				break;
			case 3:			
				k_fusion = (*mesh_props).k_fusion2; 
				type2_fission_vectors = get_type2_fission_triplets(mesh);	
    			v1 = std::get<0>(type2_fission_vectors);					
				p_propose = k_fusion * v1.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fission"<<std::endl;				
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){
					accepted = attempt_type2_fission(mesh);	
				}
				//if (accepted) std::cout<<"type2 fission "<<p_propose<<std::endl;	
				break;

			case 4:	
				//if (mesh.n_faces()>n_umb_sim ) break;

				//pick subunit type from prototype subunits
				which_type = rand() % (*mesh_props).prototype_faces.size();
				number_of_rotational_configs = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].number_of_rotational_configs;
				which_rotation = rand() % number_of_rotational_configs;
				//----
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				boundary_halfedges = find_boundary_halfedges(mesh);			
				p_propose = k_insertion * boundary_halfedges.size() * number_of_rotational_configs;
				if (p_propose>1) std::cout<<"Warning! Decrease k_insertion rate! - insertion"<<std::endl;				
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){	
					accepted = attempt_insertion(mesh, which_type, which_rotation);	
				}  
				//if (accepted) std::cout<<"insertion type "<< which_type <<" "<<p_propose<<std::endl;					
				break;

			case 5:	
				//if (mesh.n_faces()>n_umb_sim ) break;			
				//pick subunit type from prototype subunits
				which_type = rand() % (*mesh_props).prototype_faces.size();
				number_of_rotational_configs = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].number_of_rotational_configs;
				which_rotation = rand() % number_of_rotational_configs;
				//---
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				wedges = get_open_wedge_triplets(mesh);		
				w1 = std::get<0>(wedges);		
				p_propose = k_insertion * w1.size() * number_of_rotational_configs;
				if (p_propose>1) std::cout<<"Warning! Decrease k_insertion rate! - wedge_insertion"<<std::endl;						
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){
					accepted = attempt_wedge_insertion(mesh, which_type, which_rotation);
				}
				//if (accepted) std::cout<<"wedge insertion type "<< which_type <<" "<<p_propose<<std::endl;				
				break;
			case 6:
				//pick subunit type from prototype subunits
				which_type = rand() % (*mesh_props).prototype_faces.size();
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				removable_faces = get_simply_removable_faces(mesh, which_type);

				p_propose = k_insertion * removable_faces.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_insertion rate! - removal"<<std::endl;						
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){	

					accepted = attempt_removal(mesh, which_type);
				}
				break;
			case 7:
				//pick subunit type from prototype subunits
				which_type = rand() % (*mesh_props).prototype_faces.size();
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				removable_faces = get_wedge_removable_faces(mesh, which_type); 

				p_propose = k_insertion * removable_faces.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_insertion rate! - wedge_removal"<<std::endl;						
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){	
					accepted = attempt_wedge_removal(mesh, which_type);
				}
				break;

			case 8:
				k_fusion = (*mesh_props).k_fusion_edge;
				fusion_halfedges = get_halfedge_fusion_pairs(mesh, l_fuse);
    			h1 = std::get<0>(fusion_halfedges);
    			p_propose = k_fusion * h1.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fusion"<<std::endl;		
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){					
					accepted = attempt_edge_fusion(mesh);
				}
				//if (accepted) std::cout<<"EDGE fusion " <<" "<<p_propose<<std::endl;
				break;

			case 9:
				k_fusion = (*mesh_props).k_fusion_edge;
				fission_edges = get_fission_edges(mesh);
    			p_propose = k_fusion * fission_edges.size();	
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fission"<<std::endl;		
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){    									
					accepted = attempt_edge_fission(mesh);
				}
				//if (accepted) std::cout<<"EDGE fission " <<" "<<p_propose<<std::endl;
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
	}


	return;
}

void run_dynamical(std::string input_file){

	MyMesh mesh;	
	//add parameters to init_mesh; if init_config provided, init_from_om too
    RunParameters rp;	
    read_set_run_params(rp, input_file);   

    //initialize RNG from ensemble seed
    std::srand(rp.ensemble*73 + 17);

	init_mesh(mesh, rp.init_file, input_file, rp.init_file_pos);

    //should I create data_folder here?
    system(("mkdir -p "+rp.data_folder).c_str());    

    std::ofstream outfile_om(rp.data_folder+"snapshots.om", std::ios::out | std::ios::binary);
    std::ofstream outfile_data(rp.data_folder+"data_log.dat");	
    outfile_data<<"t"<<"\t"<<"key" <<"\t" <<"n_f" <<"\t"<<"n_v"<<"\t"<<"n_e"<<"\t"<< "E_el"<<"\t"<<"E_full"<<std::endl;


    long t_init=0;
    long t_final=rp.timesteps;
    propagate(mesh, t_init, t_final, rp, outfile_om,  outfile_data);
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
		convert_om_to_VTK_snapshots(rp.data_folder+"snapshots.om", rp.data_folder+"conversions/vtk_");    

}
*/



int main_old(){
	MyMesh mesh;

	//!!!!init_mesh(mesh);
	std::cout<<"n_faces "<<mesh.n_faces()<<std::endl;
	double l;
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto proto_face_props = OpenMesh::FProp<FaceProp>((*mesh_props).prototypes_mesh, "face_props");	

	OpenMesh::FPropHandleT< FaceProp > fprop;

	/*for (MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it){
			l = mesh.calc_edge_length(*e_it);

				std::cout<<"l= "<<l<<std::endl;

			
	}*/	
	//print_angle_potential(mesh);

	bool accepted, overlap;

	std::vector<MyMesh::VertexHandle*> handle_tracking_v;
   	std::vector<MyMesh::FaceHandle*> handle_tracking_f;
    std::vector<MyMesh::HalfedgeHandle*> handle_tracking_h;
    	MyMesh::FaceHandle ff;




//std::cout<<attempt_type1_fusion_wtf(mesh, 0)<<std::endl;

//std::cout<<attempt_type1_fusion_tmp(mesh)<<std::endl;

	/*for (long t=0; t<100; t++){

    		accepted = attempt_type2_fission_tmp(mesh);		
   	}
   	*/
 /*   std::cout<<"E1 = "<<full_energy(mesh)<<std::endl;
    while (!attempt_wedge_insertion_tmp(mesh, 1, 2));
    std::cout<<"E2 = "<<full_energy(mesh)<<std::endl;  
    std::cout<<"Overlap? "<<check_full_overlap(mesh)<<std::endl;  
*/
/*
    MyMesh::FaceHandle fh1 = mesh.face_handle(0);	
    MyMesh::FaceHandle fh2 = mesh.face_handle(0);	
    std::cout<<"fh1=fh2? "<<(fh1==fh2)<<std::endl;
    mesh.delete_face(fh1);
    std::cout<<"deleted? fh1, fh2 "<<mesh.status(fh1).deleted()<<" "<<mesh.status(fh2).deleted()<<std::endl;
*/

    RunParameters run_params;
    read_set_run_params(run_params, "/home/btyukodi/assembly_openmesh/sandbox/input_carlos2.json");
    std::cout<<run_params.data_folder<<std::endl;
	std::cout.precision(10);
    int which_move, which_type, which_rotation, number_of_rotational_configs;
    int nmoves = 10;
    double k_insertion, k_fusion, p_propose, r;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> , std::vector<int>> wedges;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > fusion_vectors, type2_fission_vectors;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> >	fission_vectors;	
	std::vector <MyMesh::VertexHandle> w1, v1;  
	std::vector<MyMesh::HalfedgeHandle> boundary_halfedges, h1;	  
	std::vector<MyMesh::FaceHandle> removable_faces_geom, removable_faces;	
	std::tuple<std::vector<MyMesh::HalfedgeHandle>, std::vector<MyMesh::HalfedgeHandle> > fusion_halfedges;	
	std::vector<MyMesh::EdgeHandle> fission_edges;	
	std::vector<int> n_faces;

	//auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	float kT = (*mesh_props).kT;
	double l_fuse = (*mesh_props).l_fuse;

	int n_umb_sim = 1;
/*
for (MyMesh::EdgeIter eit = mesh.edges_sbegin(); eit!=mesh.edges_end(); ++eit){

	std::cout<<"dih "<<mesh.calc_dihedral_angle(*eit)<<std::endl;

}

for (MyMesh::VertexIter vit = mesh.vertices_sbegin(); vit!=mesh.vertices_end(); ++vit){
		mesh.set_point(*vit, mesh.point(*vit)+MyMesh::Point(0,0,-2.0) );
	}
*/

	long long byte_position;
/*   std::cout<<"------- TEST restart ----------------"<<std::endl;
	byte_position = 1109;
	init_mesh_from_om(mesh, "outfileom.om", byte_position, "/home/btyukodi/assembly_openmesh/sandbox/input_carlos2.json");
	std::cout<<byte_position <<"\t" <<mesh.n_faces() <<"\t"<<mesh.n_vertices()<<"\t"<<mesh.n_edges()<<"\t"<< full_elastic_energy(mesh)<<"\t"<<full_energy(mesh)<<std::endl;   

	byte_position = 3282;
	init_mesh_from_om(mesh, "outfileom.om", byte_position, "/home/btyukodi/assembly_openmesh/sandbox/input_carlos2.json");
	std::cout<<byte_position <<"\t" <<mesh.n_faces() <<"\t"<<mesh.n_vertices()<<"\t"<<mesh.n_edges()<<"\t"<< full_elastic_energy(mesh)<<"\t"<<full_energy(mesh)<<std::endl;   


	byte_position = 10672;
	init_mesh_from_om(mesh, "outfileom.om", byte_position,"/home/btyukodi/assembly_openmesh/sandbox/input_carlos2.json");
	std::cout<<byte_position <<"\t" <<mesh.n_faces() <<"\t"<<mesh.n_vertices()<<"\t"<<mesh.n_edges()<<"\t"<< full_elastic_energy(mesh)<<"\t"<<full_energy(mesh)<<std::endl;   


	byte_position = 27482;
	init_mesh_from_om(mesh, "outfileom.om", byte_position,"/home/btyukodi/assembly_openmesh/sandbox/input_carlos2.json");
	std::cout<<byte_position <<"\t" <<mesh.n_faces() <<"\t"<<mesh.n_vertices()<<"\t"<<mesh.n_edges()<<"\t"<< full_elastic_energy(mesh)<<"\t"<<full_energy(mesh)<<std::endl;   


 std::cout<<"------- TEST restart end----------------"<<std::endl;
*/
   std::ofstream outfile("outfileom5.om", std::ios::out | std::ios::binary);
   std::ofstream datafile("data_out5.dat");
   long long fpos;

	for (long t=0; t<2000; t++){
		for (int n=0; n<mesh.n_vertices(); n++){
				attempt_move(mesh);
		}
	}	

	for (long t=0; t<7200000*0+2000000; t++){
/*		if (!test_manifold(mesh)){
			std::cout<<"NOT MANIFOLD t="<<t<<std::endl;
			break;
		}
*/
		if (t % 10000 ==0){
			//std::cout<<"------------- mesh copied --------------------------"<<std::endl;				
			std::cout<<"t = "<<t<<" n_faces="<<mesh.n_faces()<<std::endl;	
			//std::cout<<"------------- prop removed --------------------------"<<std::endl;	

			//if (mesh.n_faces()>10){
				//std::cout<<"tellp "<<outfile.tellp()<<std::endl; 			
				fpos = dump_om(mesh, t, outfile);
				dump_data(mesh, t, datafile, fpos);
				fpos=0;
				//outfile<<"-----------------";

			//	break;						
				//std::cout<<"has property "<<OpenMesh::hasProperty<OpenMesh::FaceHandle, FaceProp>(mesh, "face_props")<<std::endl;
			//}

		}

		if (t % 40000 ==0){
			n_umb_sim++;
		}
		
		//std::cout<<attempt_move_tmp(mesh)<<std::endl;
		for (int n=0; n<mesh.n_vertices(); n++){
			attempt_move(mesh);
		}


		which_move = rand() % nmoves;
		accepted=false;
		switch (which_move){
			case 0:
				k_fusion = (*mesh_props).k_fusion; 
				fusion_vectors = get_type1_fusion_triplets(mesh, l_fuse);	
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion rate! - type1_fusion"<<std::endl;
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){	
					accepted = attempt_type1_fusion(mesh);
				}
				//if (accepted) std::cout<<"type1 fusion"<<std::endl;
				break;

			case 1:
				k_fusion = (*mesh_props).k_fusion; 
				fission_vectors = get_type1_fission_pairs(mesh);	
    			v1 = std::get<0>(fission_vectors);					
				p_propose = k_fusion * v1.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion rate! - type1_fission"<<std::endl;				
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){	
					accepted = attempt_type1_fission(mesh);	
				}
				//if (accepted) std::cout<<"type1 fission"<<std::endl;
				break;

			case 2:		
				k_fusion = (*mesh_props).k_fusion2; 
				fusion_vectors = get_type2_fusion_triplets(mesh, l_fuse);	
    			v1 = std::get<0>(fusion_vectors);					
				p_propose = k_fusion * v1.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fusion"<<std::endl;				
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){
					accepted = attempt_type2_fusion(mesh);
				}
				//if (accepted) std::cout<<"type2 fusion "<<p_propose<<std::endl;		
				break;
			case 3:			
				k_fusion = (*mesh_props).k_fusion2; 
				type2_fission_vectors = get_type2_fission_triplets(mesh);	
    			v1 = std::get<0>(type2_fission_vectors);					
				p_propose = k_fusion * v1.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion2 rate! - type2_fission"<<std::endl;				
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){
					accepted = attempt_type2_fission(mesh);	
				}
				//if (accepted) std::cout<<"type2 fission "<<p_propose<<std::endl;	
				break;

			case 4:	
				//if (mesh.n_faces()>n_umb_sim ) break;

				//pick subunit type from prototype subunits
				which_type = rand() % (*mesh_props).prototype_faces.size();
				number_of_rotational_configs = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].number_of_rotational_configs;
				which_rotation = rand() % number_of_rotational_configs;
				//----
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				boundary_halfedges = find_boundary_halfedges(mesh);			
				p_propose = k_insertion * boundary_halfedges.size() * number_of_rotational_configs;
				if (p_propose>1) std::cout<<"Warning! Decrease k_insertion rate! - insertion"<<std::endl;				
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){	
					accepted = attempt_insertion(mesh, which_type, which_rotation);	
				}  
				//if (accepted) std::cout<<"insertion type "<< which_type <<" "<<p_propose<<std::endl;					
				break;

			case 5:	
				//if (mesh.n_faces()>n_umb_sim ) break;			
				//pick subunit type from prototype subunits
				which_type = rand() % (*mesh_props).prototype_faces.size();
				number_of_rotational_configs = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].number_of_rotational_configs;
				which_rotation = rand() % number_of_rotational_configs;
				//---
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				wedges = get_open_wedge_triplets(mesh);		
				w1 = std::get<0>(wedges);		
				p_propose = k_insertion * w1.size() * number_of_rotational_configs;
				if (p_propose>1) std::cout<<"Warning! Decrease k_insertion rate! - wedge_insertion"<<std::endl;						
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){
					accepted = attempt_wedge_insertion(mesh, which_type, which_rotation);
				}
				//if (accepted) std::cout<<"wedge insertion type "<< which_type <<" "<<p_propose<<std::endl;				
				break;
			case 6:
				//pick subunit type from prototype subunits
				which_type = rand() % (*mesh_props).prototype_faces.size();
				//----
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				removable_faces = get_simply_removable_faces(mesh, which_type);

				//FILTER IT FOR THE GIVEN TYPE!!!!!	
				/*removable_faces.clear();
				for (auto & rface : removable_faces_geom){
					if (face_props[rface].face_type == which_type){
						removable_faces.push_back(rface);
					}
				}*/

				p_propose = k_insertion * removable_faces.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_insertion rate! - removal"<<std::endl;						
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){	

					attempt_removal(mesh, which_type);
				}
				break;
			case 7:
				//pick subunit type from prototype subunits
				which_type = rand() % (*mesh_props).prototype_faces.size();
				//----
				k_insertion = proto_face_props[ (*mesh_props).prototype_faces[which_type] ].k_insertion;
				removable_faces = get_wedge_removable_faces(mesh, which_type); //by geometry only

				//FILTER IT FOR THE GIVEN TYPE!!!!!	
				/*removable_faces.clear();
				for (auto & rface : removable_faces_geom){
					if (face_props[rface].face_type == which_type){
						removable_faces.push_back(rface);
					}
				}*/

				p_propose = k_insertion * removable_faces.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_insertion rate! - wedge_removal"<<std::endl;						
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){	
					attempt_wedge_removal(mesh, which_type);
				}
				break;

			case 8:
				//need to add k_fusion_edge!!!
				k_fusion = (*mesh_props).k_fusion_edge;
				fusion_halfedges = get_halfedge_fusion_pairs(mesh, l_fuse);
    			h1 = std::get<0>(fusion_halfedges);
    			p_propose = k_fusion * h1.size();
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fusion"<<std::endl;		
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){					
					accepted = attempt_edge_fusion(mesh);
				}
				if (accepted) std::cout<<"EDGE fusion " <<" "<<p_propose<<std::endl;
				//std::cout<<"COM "<< test_com_updateness(mesh)<<std::endl;
				//if (accepted) goto QUIT;
				break;

			case 9:
				//need to add k_fusion_edge!!!
				k_fusion = (*mesh_props).k_fusion_edge;
				fission_edges = get_fission_edges(mesh);
    			p_propose = k_fusion * fission_edges.size();	
				if (p_propose>1) std::cout<<"Warning! Decrease k_fusion_edge rate! - edge_fission"<<std::endl;		
				r = rand()/(RAND_MAX + 1.0);
				if (r<p_propose){    									
					accepted = attempt_edge_fission(mesh);
				}
				if (accepted) std::cout<<"EDGE fission " <<" "<<p_propose<<std::endl;
				//std::cout<<"COM "<< test_com_updateness(mesh)<<std::endl;
				//if (accepted) goto QUIT;
				break;

		}

    	handle_tracking_f.clear();
    	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
    		for (auto & neighbor: face_props[*fit].neighbor_list){
    			handle_tracking_f.push_back( &neighbor);
    		}


    		/*for (int xx=0; xx<10; xx++){
    			handle_tracking_f.push_back( &ff );
    		}*/
    	}
		mesh.garbage_collection<std::vector<MyMesh::VertexHandle*>, std::vector<MyMesh::HalfedgeHandle*>, std::vector<MyMesh::FaceHandle*> >(handle_tracking_v, handle_tracking_h, handle_tracking_f);


		//test_neighbor_list(mesh);				
		//mesh.garbage_collection();
		//if (accepted){
		//	std::cout<<accepted<<std::endl;
		//}
		//test_neighbor_list(mesh);		
		//std::cout<<"COM "<< test_com_updateness(mesh)<<std::endl;


		if (t % 150 ==0){
			update_full_neighbor_list(mesh);
			//overlap = check_full_overlap(mesh);
			//if (overlap)
			//	std::cout<<"Overlap? "<<overlap<<std::endl;
			//CHECK OVERLAPS FOR FULL CAPSID HERE
		}
		n_faces.push_back( mesh.n_faces() );

		//std::cout<<"GARBAGE COLL length="<<handle_tracking_f.size()<<std::endl;

		//mesh.garbage_collection();

		//!!MAYBE ONLY DO GARBAGE COLLECTION AT THE END OF EACH TIMESTEP??
	//std::cout<<"n_faces "<<mesh.n_faces()<<std::endl;
/*		accepted = attempt_type2_fission_tmp(mesh);
		if (accepted){
			std::cout<<"fission "<<accepted<<std::endl;
		}
*/		
		//std::cout<<attempt_type2_fission_tmp(mesh)<<std::endl;
		//get_type2_fission_triplets(mesh);		
	}

    	//attempt_removal_tmp(mesh);

/*	for (MyMesh::VertexIter vit = mesh.vertices_sbegin(); vit!=mesh.vertices_end(); ++vit){
		std::cout<<"valence "<<mesh.valence(*vit)<<std::endl;
	}
*/

QUIT:std::ofstream outFile("/home/btyukodi/assembly_openmesh/tests/n_faces_mu3.txt");
// the important part
for (const auto &nf : n_faces) outFile << nf << "\n";

MyMesh::FaceHandle fh0 = mesh.face_handle(0);
MyMesh::FaceHandle fh10 = mesh.face_handle(1);
MyMesh::FaceHandle fhd = mesh.face_handle(3);
/*mesh.delete_face(fh0);
mesh.delete_face(fh10);
mesh.delete_face(fhd);*/
double zz = 0.05;
std::cout.precision(20);
std::cout<<"E before "<<full_elastic_energy(mesh)<<" "<<zz<<std::endl;
mesh.garbage_collection();
std::cout<<"E after "<<full_elastic_energy(mesh)<<std::endl;



	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");

	int edge_type;
	float e_b,theta0;
	for (MyMesh::FaceIter f_it = mesh.faces_sbegin(); f_it != mesh.faces_end(); ++f_it){
		std::cout<<"face "<< *f_it<<" ; type: "<<face_props[*f_it].face_type;
		for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(*f_it); fh_it.is_valid(); ++fh_it){
			e_b = edge_props[*fh_it].e_b[0];
			theta0 = edge_props[*fh_it].theta0[0];		
			edge_type = edge_props[*fh_it].edge_type;	
			//std::cout<<"e_b"<< e_b<<" theta0 "<<theta0;
			std::cout<<" edge_type "<<edge_type;//<<std::endl;


		}	
		std::cout<<std::endl;
 
	}
		

/*	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	for (MyMesh::HalfedgeIter h_it=mesh.halfedges_begin(); h_it!=mesh.halfedges_end(); ++h_it){
		std::cout<<"prop "<<*h_it<<" "<<edge_props[*h_it].l0<<std::endl;
	}

	MyMesh::HalfedgeHandle hedge = mesh.halfedge_handle(1);
	double l0 = edge_props[hedge].l0;
		std::cout<<"propf "<<hedge<<" "<< l0<<std::endl;
	double E;
	E = edge_stretch_energy(mesh, mesh.edge_handle(hedge));
//	std::cout<<1.0/(std::numeric_limits<double>::max()*2.0)<<std::endl;
	std::cout<<E<<std::endl;
*/

mesh.garbage_collection();


outfile.close();
datafile.close();

//dump_lammps_snapshot(mesh, "lammps_snapshot.dat");

//-------
std::cout<<"TESTREAD "<<std::endl;
convert_om_to_lammps_trajectory("outfileom4.om");
	//std::cout<<"breakpA";
convert_om_to_lammps_snapshots("outfileom4.om", "dump_snapshots/lammps_");
	//std::cout<<"breakpB";
convert_om_to_VTK_snapshots("outfileom4.om", "dump_snapshots/vtk_");
//dump_vtk_snapshot(mesh, "vtktest.vtk");	
return 0;

  MyMesh readmesh;
  auto face_props_readmesh = OpenMesh::FProp<int>(readmesh, "face_type");
  auto edge_props_readmesh = OpenMesh::HProp<int>(readmesh, "edge_type");
const std::string _ext=".OM";
	OpenMesh::IO::Options ropt, wopt;
	ropt+=OpenMesh::IO::Options::FaceColor; 
	//wopt+=OpenMesh::IO::Options::Custom; 
std::ifstream rf("outfileom.om", std::ios::out | std::ios::binary);
   rf.seekg(1564, rf.beg); 
if (!OpenMesh::IO::read_mesh(readmesh, rf, _ext, ropt)) 	
		
	{
		std::cerr << "read error\n";
    	exit(1);
	}



	for (MyMesh::FaceIter f_it = readmesh.faces_sbegin(); f_it != readmesh.faces_end(); ++f_it){
		std::cout<<"face "<< *f_it<<" ";
		for (MyMesh::FaceHalfedgeIter fh_it = readmesh.fh_iter(*f_it); fh_it.is_valid(); ++fh_it){	
			edge_type = edge_props_readmesh[*fh_it];	
			//std::cout<<"e_b"<< e_b<<" theta0 "<<theta0;
			std::cout<<" edge_type "<<edge_type;//<<std::endl;


		}	
		std::cout<<std::endl;
 
	}


//auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
//auto proto_face_props = OpenMesh::FProp<FaceProp>((*mesh_props).prototypes_mesh, "face_props");
std::cout<<"prototype mesh "<< proto_face_props[(*mesh_props).prototype_faces[0]].mu<<std::endl;


	std::cout<<"n_faces "<<mesh.n_faces()<<std::endl;
try
  {
    if ( !OpenMesh::IO::write_mesh(mesh, "output.off") )
    {
      std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
      return 1;
    }
  }
  catch( std::exception& x )
  {
    std::cerr << x.what() << std::endl;
    return 1;
  }
  return 0;

	return 0;
}