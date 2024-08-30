#include <iostream>
#include <stdlib.h>
#include <random>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Apps/Assembly/topology.hh>
#include <OpenMesh/Apps/Assembly/energy.hh>
#include <OpenMesh/Apps/Assembly/geometry.hh>
#include <OpenMesh/Apps/Assembly/monte_carlo.hh>
#include <OpenMesh/Apps/Assembly/excluders.hh>
#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>
#include <OpenMesh/Apps/Assembly/random.hh>
/*
void dump_mesh(MyMesh & mesh){
	try
	  {
	    if ( !OpenMesh::IO::write_mesh(mesh, "err_output.off") )
	    {
	      std::cerr << "Cannot write mesh to file 'err_output.off'" << std::endl;
	      return;
	    }
	  }
	  catch( std::exception& x )
	  {
	    std::cerr << x.what() << std::endl;
	    return;
	  }

}

double get_mesh_comx(MyMesh & mesh){
	double comx=0.0;
	int nv=0;
	for (MyMesh::VertexIter v_it = mesh.vertices_sbegin();v_it != mesh.vertices_end(); ++v_it ){
		comx+=mesh.point(*v_it)[0];
		nv++;
	}
	return comx;
}


double get_mesh_face_comx(MyMesh & mesh){
	double comx=0.0;
	int nv=0;
	MyMesh::Point COM;
	for (MyMesh::FaceIter f_it = mesh.faces_sbegin();f_it != mesh.faces_end(); ++f_it ){

		mesh.calc_face_centroid(*f_it, COM );
		comx+=COM[0];
		nv++;
	}
	return comx;
}
*/

//returns the equilibrium position of the unbound vertex of face
//if the unbound vertex is set at this position, there should be no bending energy on the bound edge and no stretching energy on the unbound edges
//meant to be used for simple insertion and simple removal
MyMesh::Point get_equilibrium_position(MyMesh & mesh, MyMesh::FaceHandle & face){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	

	MyMesh::HalfedgeHandle hedge1, hedge2, he;
	double l_AC, l_BC, l_AB, xC, yC;

	MyMesh::VertexHandle vC, v_from, v_to;
	for (MyMesh::FaceVertexIter fv = mesh.fv_iter(face); fv.is_valid(); ++fv ){
		//std::cout<<"valence "<<mesh.valence(*fv)<<std::endl;
		if (mesh.valence(*fv)==2){
			vC = *fv;
		}
	}	

	//find the non-boundary edge of face (that's the one that is attached)
	for (MyMesh::FaceHalfedgeIter fh = mesh.fh_iter(face); fh.is_valid(); ++fh ){
		//std::cout<<"valence "<<mesh.valence(*fv)<<std::endl;
		//if parent edge is not boundary, we found it
		if( !mesh.is_boundary(  mesh.edge_handle(*fh) ) ){
			he = *fh;
		}
	
	}		

	v_from = mesh.from_vertex_handle(he);
	v_to = mesh.to_vertex_handle(he);

	//get edge lengths and in-plane position of new vertex
	hedge1 = mesh.find_halfedge(v_to, vC);
	hedge2 = mesh.find_halfedge( vC, v_from);
	if (mesh.is_boundary(hedge1) || mesh.is_boundary(hedge2)){
		std::cout<<"Warning! attempt_insertion; not OK, one of them is boundary"<<std::endl;
	}


	l_AC = edge_props[hedge1].l0;
	l_BC = edge_props[hedge2].l0;	
	l_AB = mesh.calc_edge_length( mesh.edge_handle(he));
	//get in-plane position of new_vertex, assuming the general case of non-equilateral triangles, l01!=l02 necessarily
	std::tie(xC, yC) = calculate_triangle_third_vertex_coordinate(l_AB, l_BC, l_AC);

	//need a Point to std::vector<double> converter
	//or just use Point for r vectors

	MyMesh::VertexHandle v0 = mesh.to_vertex_handle(mesh.next_halfedge_handle( mesh.opposite_halfedge_handle(he) ) );
	MyMesh::Point r0 = mesh.point(v0);
	MyMesh::Point ex = (mesh.point(v_from) - mesh.point(v_to)).normalize();
	MyMesh::Point b = r0 - mesh.point(v_from);
	MyMesh::Point ey = -(b - (b.dot(ex)) * ex).normalize();
	MyMesh::Point rC = xC*ex + yC*ey;

	hedge1 =  mesh.find_halfedge(v_to, v_from);
	hedge2 = mesh.opposite_halfedge_handle(hedge1);
	int edge_type1 = edge_props[hedge1].edge_type;
	int edge_type2 = edge_props[hedge2].edge_type;
	double theta0 = 0.5*(edge_props[hedge1].theta0[edge_type2] + edge_props[hedge2].theta0[edge_type1] );

	//now need to pick a rotation
	rotatevec(rC, ex, theta0);



	//REPLACE WITH PROPER VOLUME RAND
	rC+=mesh.point(v_to);

	return rC;
}

/*
Attempts insertion of a new face of type face_type with rotation rotational_config=0,1,2 compared to the
anchor edge of the corresponding prototype face. 
Position (the boundary edge) where the new subunit is inserted is chosen randomly.
Positions (boundary edges) are currently searched twice, in this function and when computing p_propose.
This could be optimized by searching only outside and passing boundary_halfedges to this function.
Doesn't seem to make much of a speedup now.

- test if randvec() works as expected
*/
bool attempt_insertion_old_bak(MyMesh & mesh, int face_type, int rotational_config, UmbrellaWindow & uw, std::mt19937 & eng){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");		
	//get a prototype face of kind face_type for properties
	MyMesh::FaceHandle proto_face = (*mesh_props).prototype_faces[face_type];
	std::vector<MyMesh::HalfedgeHandle> boundary_halfedges = find_boundary_halfedges(mesh);
	//int which = int(rand() % boundary_halfedges.size());
	int which = randint(boundary_halfedges.size(), eng);	
	MyMesh::VertexHandle v_from, v_to, new_vertex;
	std::vector<MyMesh::VertexHandle> vhandles;
	std::vector<MyMesh::EdgeHandle> affected_edges;	
	MyMesh::FaceHandle new_face;
	MyMesh::Point new_vertex_point;
	double E_before, E_after, dE, p, r, kT, d_max, R_add, v_add, l_AC, l_BC, l_AB, xC, yC, D;
	MyMesh::HalfedgeHandle he = boundary_halfedges[which];
	MyMesh::HalfedgeHandle anchor_halfedge, hedge1, hedge2;
	//should be up to date, after garbage collection!!
	int n_faces_before = mesh.n_faces();
	int n_faces_after = n_faces_before + 1;

	//---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------

	kT = (*mesh_props).kT;
	R_add = (*mesh_props).R_add;	
	v_add = 4.0*M_PI*R_add*R_add*R_add/3.0;	
	D = (*mesh_props).DOS;
	//d_max = (*mesh_props).d_max;		

	//affected edges before: edge(he)
	affected_edges.push_back( mesh.edge_handle(he) );
	E_before = energy(mesh, affected_edges);
	E_before+=uw.bias_potential(n_faces_before);


	v_from = mesh.from_vertex_handle(he);
	v_to = mesh.to_vertex_handle(he);

	//need to compute new vertex position: find eq. position then rotate
	//for now use some crap; needs to be repositioned after rotation is decided, according to theta0 between old hedge and new hedge
	new_vertex_point = 0.5*(mesh.point( v_from ) + mesh.point( v_to ));

	new_vertex = mesh.add_vertex(new_vertex_point);
	vhandles.push_back(v_from);
	vhandles.push_back(v_to);
	vhandles.push_back(new_vertex);
	new_face = mesh.add_face(vhandles);

	//need to set new_face 1) face properties, 2) edge properties from prototype
	//adjust edge rotation/permutation
	anchor_halfedge = mesh.find_halfedge(v_from, v_to);
	for (int ir=0; ir<rotational_config; ir++){
		anchor_halfedge = mesh.next_halfedge_handle(anchor_halfedge);
	}
	clone_properties_from_prototype(mesh, face_type, new_face, anchor_halfedge);

	//get edge lengths and in-plane position of new vertex
	hedge1 = mesh.find_halfedge(v_to, new_vertex);
	hedge2 = mesh.find_halfedge( new_vertex, v_from);
	if (mesh.is_boundary(hedge1) || mesh.is_boundary(hedge2)){
		std::cout<<"Warning! attempt_insertion; not OK, one of them is boundary"<<std::endl;
	}


	l_AC = edge_props[hedge1].l0;
	l_BC = edge_props[hedge2].l0;	
	l_AB = mesh.calc_edge_length( mesh.edge_handle(he));
	//get in-plane position of new_vertex, assuming the general case of non-equilateral triangles, l01!=l02 necessarily
	std::tie(xC, yC) = calculate_triangle_third_vertex_coordinate(l_AB, l_BC, l_AC);

	//need a Point to std::vector<double> converter
	//or just use Point for r vectors

	MyMesh::VertexHandle v0 = mesh.to_vertex_handle(mesh.next_halfedge_handle( mesh.opposite_halfedge_handle(he) ) );
	MyMesh::Point r0 = mesh.point(v0);
	MyMesh::Point ex = (mesh.point(v_from) - mesh.point(v_to)).normalize();
	MyMesh::Point b = r0 - mesh.point(v_from);
	MyMesh::Point ey = -(b - (b.dot(ex)) * ex).normalize();
	MyMesh::Point rC = xC*ex + yC*ey;

	hedge1 =  mesh.find_halfedge(v_to, v_from);
	hedge2 = mesh.opposite_halfedge_handle(hedge1);
	int edge_type1 = edge_props[hedge1].edge_type;
	int edge_type2 = edge_props[hedge2].edge_type;
	double theta0 = 0.5*(edge_props[hedge1].theta0[edge_type2] + edge_props[hedge2].theta0[edge_type1] );

	//now need to pick a rotation
	rotatevec(rC, ex, theta0);



	//REPLACE WITH PROPER VOLUME RAND
	rC+=mesh.point(v_to);

	/*std::cout<<"rC1 "<<rC<<std::endl;
	std::cout<<"rC2 "<<get_equilibrium_position(mesh, new_face)<<std::endl;
	mesh.set_point(new_vertex, rC);
	std::cout<<"E_bend1 "<<edge_bending_energy(mesh, mesh.edge_handle(anchor_halfedge), edge_props)<<std::endl;
	//std::cout<<"E_bend2 "<<edge_bending_energy(mesh, mesh.edge_handle(hedge2), edge_props)<<std::endl;	
	std::cout<<"E_stretch1 "<<edge_stretch_energy(mesh, mesh.edge_handle(mesh.find_halfedge(v_to, new_vertex)), edge_props)<<std::endl;
	std::cout<<"E_stretch2 "<<edge_stretch_energy(mesh, mesh.edge_handle(mesh.find_halfedge( new_vertex, v_from)), edge_props)<<std::endl;		
	*/
	//std::cout<<"rC "<<rC<<" "<<(rC - mesh.point(v_to)).norm()<<" "<< (rC - mesh.point(v_from)).norm() <<std::endl;		
	rC+=randvec(eng)*R_add;


	mesh.set_point(new_vertex, rC);


	//========= excluder =========================================
	//update new face's COM
	mesh.calc_face_centroid(new_face, face_props[new_face].COM );
	//update new face's neighbor list
	update_neighbor_list(mesh, new_face);
	//check new face's overlap with its neighbors
	bool overlap = check_neighbor_overlap(mesh, new_face);
	//if (overlap) std::cout<<"OVERLAP "<<std::endl;
	//============================================================


	//affected_edges_after: of get_affected_edges(new_vertex)
	affected_edges = get_affected_edges(mesh, new_vertex);
	E_after = energy(mesh, affected_edges);
	E_after+=uw.bias_potential(n_faces_after);
	dE = E_after-E_before;
	//std::cout<<"dE "<<dE<<std::endl;	


	float mu=face_props[new_face].mu;	
	p = D*v_add*exp(-(dE-mu)/kT);
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if (  (r<p) && (!overlap) ){
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//std::cout<<"INSERTION acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//------------------------

		//========= excluder =========================================
		//if accepted, add the new face to neighbors' neighbor lists
		for (auto & neighbor: face_props[new_face].neighbor_list){
			face_props[neighbor].neighbor_list.push_back(new_face);
		}
		//============================================================

		return true;
	}
	else{
		mesh.delete_face(new_face);
		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//std::cout<<"INSERTION rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------

		return false;
	}
}


bool attempt_insertion(MyMesh & mesh, int face_type, int rotational_config, UmbrellaWindow & uw, std::mt19937 & eng){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");		
	//get a prototype face of kind face_type for properties
	MyMesh::FaceHandle proto_face = (*mesh_props).prototype_faces[face_type];
	std::vector<MyMesh::HalfedgeHandle> boundary_halfedges = find_boundary_halfedges(mesh);
	//int which = int(rand() % boundary_halfedges.size());
	int which = randint(boundary_halfedges.size(), eng);	
	MyMesh::VertexHandle v_from, v_to, new_vertex;
	std::vector<MyMesh::VertexHandle> vhandles;
	std::vector<MyMesh::EdgeHandle> affected_edges;	
	MyMesh::FaceHandle new_face;
	MyMesh::Point new_vertex_point;
	double E_before, E_after, dE, p, r, kT, d_max, R_add, v_add, D;
	MyMesh::HalfedgeHandle he = boundary_halfedges[which];
	MyMesh::HalfedgeHandle anchor_halfedge;
	//should be up to date, after garbage collection!!
	int n_faces_before = mesh.n_faces();
	int n_faces_after = n_faces_before + 1;

	//---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------

	kT = (*mesh_props).kT;
	R_add = (*mesh_props).R_add;	
	v_add = 4.0*M_PI*R_add*R_add*R_add/3.0;	
	D = (*mesh_props).DOS;
	//d_max = (*mesh_props).d_max;		

	//affected edges before: edge(he)
	affected_edges.push_back( mesh.edge_handle(he) );
	E_before = energy(mesh, affected_edges);
	E_before+=uw.bias_potential(n_faces_before);


	v_from = mesh.from_vertex_handle(he);
	v_to = mesh.to_vertex_handle(he);

	//need to compute new vertex position: find eq. position then rotate
	//for now use some crap; needs to be repositioned after rotation is decided, according to theta0 between old hedge and new hedge
	new_vertex_point = 0.5*(mesh.point( v_from ) + mesh.point( v_to ));

	new_vertex = mesh.add_vertex(new_vertex_point);
	vhandles.push_back(v_from);
	vhandles.push_back(v_to);
	vhandles.push_back(new_vertex);
	new_face = mesh.add_face(vhandles);

	//need to set new_face 1) face properties, 2) edge properties from prototype
	//adjust edge rotation/permutation
	anchor_halfedge = mesh.find_halfedge(v_from, v_to);
	for (int ir=0; ir<rotational_config; ir++){
		anchor_halfedge = mesh.next_halfedge_handle(anchor_halfedge);
	}
	clone_properties_from_prototype(mesh, face_type, new_face, anchor_halfedge);


	MyMesh::Point rC = get_equilibrium_position(mesh, new_face);

	/*std::cout<<"rC1 "<<rC<<std::endl;
	std::cout<<"rC2 "<<get_equilibrium_position(mesh, new_face)<<std::endl;
	mesh.set_point(new_vertex, rC);
	std::cout<<"E_bend1 "<<edge_bending_energy(mesh, mesh.edge_handle(anchor_halfedge), edge_props)<<std::endl;
	//std::cout<<"E_bend2 "<<edge_bending_energy(mesh, mesh.edge_handle(hedge2), edge_props)<<std::endl;	
	std::cout<<"E_stretch1 "<<edge_stretch_energy(mesh, mesh.edge_handle(mesh.find_halfedge(v_to, new_vertex)), edge_props)<<std::endl;
	std::cout<<"E_stretch2 "<<edge_stretch_energy(mesh, mesh.edge_handle(mesh.find_halfedge( new_vertex, v_from)), edge_props)<<std::endl;		
	*/
	//std::cout<<"rC "<<rC<<" "<<(rC - mesh.point(v_to)).norm()<<" "<< (rC - mesh.point(v_from)).norm() <<std::endl;	
	MyMesh::Point dr = 	randvec_normal(eng)*R_add;
	rC+=dr;


	mesh.set_point(new_vertex, rC);


	//========= excluder =========================================
	//update new face's COM
	mesh.calc_face_centroid(new_face, face_props[new_face].COM );
	//update new face's neighbor list
	update_neighbor_list(mesh, new_face);
	//check new face's overlap with its neighbors
	bool overlap = check_neighbor_overlap(mesh, new_face);
	//if (overlap) std::cout<<"OVERLAP "<<std::endl;
	//============================================================


	//affected_edges_after: of get_affected_edges(new_vertex)
	affected_edges = get_affected_edges(mesh, new_vertex);
	E_after = energy(mesh, affected_edges);
	E_after+=uw.bias_potential(n_faces_after);
	dE = E_after-E_before;
	//std::cout<<"dE "<<dE<<std::endl;	


	float mu=face_props[new_face].mu;	
	//p = D*v_add*exp(-(dE-mu)/kT);

	double pdf = normal_pdf(dr[0]/R_add)*normal_pdf(dr[1]/R_add)*normal_pdf(dr[2]/R_add) / (R_add*R_add*R_add);
	p = (1.0/pdf)*exp(-(dE-mu)/kT);

	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if (  (r<p) && (!overlap) ){
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//std::cout<<"INSERTION acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//------------------------

		//========= excluder =========================================
		//if accepted, add the new face to neighbors' neighbor lists
		for (auto & neighbor: face_props[new_face].neighbor_list){
			face_props[neighbor].neighbor_list.push_back(new_face);
		}
		//============================================================

		return true;
	}
	else{
		mesh.delete_face(new_face);
		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//std::cout<<"INSERTION rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------

		return false;
	}
}

/*
Attempts insertion of a new face of type face_type with rotation rotational_config=0,1,2 compared to the
anchor edge of the corresponding prototype face. 
Position (the open  wedge) where the new subunit is inserted is chosen randomly.
Positions (wedges) are currently searched twice, in this function and when computing p_propose.
This could be optimized by searching only outside and passing wedges to this function.
Doesn't seem to make much of a speedup now.
*/

bool attempt_wedge_insertion(MyMesh & mesh, int face_type, int rotational_config, UmbrellaWindow & uw, std::mt19937 & eng){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	//get_wedges
	std::vector <MyMesh::VertexHandle> w1, w2, w3;
    std::vector<MyMesh::EdgeHandle> affected_edges1, affected_edges2;	
    double kT, E_before, E_after, dE, r, p;

	//should be up to date, after garbage collection!!
	int n_faces_before = mesh.n_faces();
	int n_faces_after = n_faces_before + 1;

	auto wedges = get_open_wedge_triplets(mesh);
    w1 = std::get<0>(wedges);
    w2 = std::get<1>(wedges);
    w3 = std::get<2>(wedges); 
    if (w1.size()==0){
    	return false;
    }

    //---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------

	//int which = int(rand() % w1.size());
	int which = randint( w1.size(), eng );

	affected_edges1 = get_affected_edges(mesh, w1[which]);
	affected_edges2 = get_affected_edges(mesh, w3[which]);

	E_before = energy(mesh, affected_edges1) + energy(mesh, affected_edges2);
	E_before+=uw.bias_potential(n_faces_before);	

	std::vector<MyMesh::VertexHandle> vhandles;
	vhandles.push_back(w1[which]);
	vhandles.push_back(w2[which]);
	vhandles.push_back(w3[which]);

	MyMesh::FaceHandle new_face = mesh.add_face(vhandles);

	kT = (*mesh_props).kT;

	//anchor halfedge belongs to newly added face
	MyMesh::HalfedgeHandle anchor_halfedge = mesh.find_halfedge(w1[which], w2[which]);
	for (int ir=0; ir<rotational_config; ir++){
		anchor_halfedge = mesh.next_halfedge_handle(anchor_halfedge);
	}
	clone_properties_from_prototype(mesh, face_type, new_face, anchor_halfedge);

	//========= excluder =========================================
	//update new face's COM
	mesh.calc_face_centroid(new_face, face_props[new_face].COM );
	//update new face's neighbor list
	update_neighbor_list(mesh, new_face);
	//check new face's overlap with its neighbors
	bool overlap = check_neighbor_overlap(mesh, new_face);
	//if (overlap) std::cout<<"OVERLAP "<<std::endl;
	//============================================================	

	//affected_edges1 = get_affected_edges(mesh, w1[which]);
	//affected_edges2 = get_affected_edges(mesh, w3[which]);
	MyMesh::EdgeHandle new_edge = find_edge(mesh, w1[which], w3[which]);
	E_after = energy(mesh, affected_edges1) + energy(mesh, affected_edges2) + edge_stretch_energy(mesh, new_edge, edge_props );
	E_after+=uw.bias_potential(n_faces_after);

	dE = E_after-E_before;
	//std::cout<<"wedge dE "<<dE<<std::endl;	


	float mu=face_props[new_face].mu;	
	p = exp(-(dE-mu)/kT);
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if ( (r<p)  && (!overlap)){
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//std::cout<<"WEDGE INSERTION acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//------------------------

		//========= excluder =========================================
		//if accepted, add the new face to neighbors' neighbor lists
		for (auto & neighbor: face_props[new_face].neighbor_list){
			face_props[neighbor].neighbor_list.push_back(new_face);
		}
		//============================================================

		return true;
	}
	else{
		mesh.delete_face(new_face);
		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//std::cout<<"WEDGE INSERTION rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<" "<<overlap<<std::endl;
		//std::cout<<"OL? "<<overlap<<std::endl;
		//--------------------------------------

		return false;
	}
}


/*
Attempts simple (non-wedge) removal of a new face of type face_type
Removable face is chosen randomly.
Removable faces are currently searched twice, in this function and when computing p_propose.
This could be optimized by searching only outside and passing removable_faces to this function.
Doesn't seem to make much of a speedup now.
*/
bool attempt_removal(MyMesh & mesh, int face_type, UmbrellaWindow & uw, std::mt19937 & eng){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	double E_before, E_after, dE, kT, mu, r, p, R_add, v_add, D;

	//should be up to date, after garbage collection!!
	int n_faces_before = mesh.n_faces();
	int n_faces_after = n_faces_before - 1;	

	//find removable subunits geometrically
/*	std::vector<MyMesh::FaceHandle> removable_faces_geom = get_simply_removable_faces(mesh, face_type);
	//filter them for the current subunit type only
	std::vector<MyMesh::FaceHandle> removable_faces;
	for (auto & rface : removable_faces_geom){
		if (face_props[rface].face_type == face_type){
			removable_faces.push_back(rface);
		}
	}	
*/
	std::vector<MyMesh::FaceHandle> removable_faces = get_simply_removable_faces(mesh, face_type);
	if (removable_faces.size()==0){
		//std::cout<<"nothing to remove"<<std::endl;
		return false;
	}
	kT = (*mesh_props).kT;
	R_add = (*mesh_props).R_add;	
	v_add = 4.0*M_PI*R_add*R_add*R_add/3.0;
	D = (*mesh_props).DOS;
    //---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------	

	//int which = int(rand() % removable_faces.size());
	int which = randint( removable_faces.size(), eng);

	mu=face_props[ removable_faces[which] ].mu;		
	//need to find hanging vertex to find affected edges; or just use the faces' edges
	MyMesh::VertexHandle v0, v1, v2;
	MyMesh::HalfedgeHandle he;
	for (MyMesh::FaceVertexIter fv = mesh.fv_iter(removable_faces[which]); fv.is_valid(); ++fv ){
		//std::cout<<"valence "<<mesh.valence(*fv)<<std::endl;
		if (mesh.valence(*fv)==2){
			v0 = *fv;
		}
	}

	he = mesh.halfedge_handle(v0);
	if (mesh.is_boundary(he)){
		he = mesh.opposite_halfedge_handle(he);
	}
	if (mesh.from_vertex_handle(he)==v0){
		he = mesh.next_halfedge_handle(he);
		v1 = mesh.from_vertex_handle(he);
		v2 = mesh.to_vertex_handle(he);
	}
	if (mesh.to_vertex_handle(he)==v0){
		he = mesh.prev_halfedge_handle(he);
		v1 = mesh.from_vertex_handle(he);
		v2 = mesh.to_vertex_handle(he);
	}	

    std::vector<MyMesh::EdgeHandle> affected_edges = get_affected_edges(mesh, v0);
    E_before = energy(mesh, affected_edges);
	E_before+=uw.bias_potential(n_faces_before);    
   
    E_after = halfedge_stretch_energy(mesh, mesh.opposite_halfedge_handle(he), edge_props); 	
	E_after+=uw.bias_potential(n_faces_after);    

    dE = E_after - E_before;


	MyMesh::Point rC = get_equilibrium_position(mesh, removable_faces[which]);
	MyMesh::Point dr = rC - mesh.point(v0);
    double pdf = normal_pdf(dr[0]/R_add)*normal_pdf(dr[1]/R_add)*normal_pdf(dr[2]/R_add)/ (R_add*R_add*R_add);
	//p = (1.0/ (D*v_add) )*exp(-(dE+mu)/kT);
	p = pdf*exp(-(dE+mu)/kT);
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if (r<p){

    	//============== excluder ============================  	
    	//remove the to-be-removed face from its neighbors neighbor lists
    	//going through neighbors of the removable face
    	std::vector<MyMesh::FaceHandle> * vec;
    	for (auto & neighbor: face_props[ removable_faces[which] ].neighbor_list){
    		vec = &(face_props[neighbor].neighbor_list);
    		//this line supposedly removes from vector
			(*vec).erase(std::remove((*vec).begin(), (*vec).end(), removable_faces[which] ), (*vec).end());    	
		}
    	//====================================================

		//std::cout<<"REMOVED"<<std::endl;
    	mesh.delete_face(removable_faces[which]);		
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//	std::cout<<"REMOVAL acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//------------------------

		return true;
	}    

	else{

		//std::cout<<"NOT REMOVED"<<std::endl;
		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//std::cout<<"REMOVAL rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------		
		return false;
	}

}


bool attempt_removal_old_bak(MyMesh & mesh, int face_type, UmbrellaWindow & uw, std::mt19937 & eng){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	double E_before, E_after, dE, kT, mu, r, p, R_add, v_add, D;

	//should be up to date, after garbage collection!!
	int n_faces_before = mesh.n_faces();
	int n_faces_after = n_faces_before - 1;	

	//find removable subunits geometrically
/*	std::vector<MyMesh::FaceHandle> removable_faces_geom = get_simply_removable_faces(mesh, face_type);
	//filter them for the current subunit type only
	std::vector<MyMesh::FaceHandle> removable_faces;
	for (auto & rface : removable_faces_geom){
		if (face_props[rface].face_type == face_type){
			removable_faces.push_back(rface);
		}
	}	
*/
	std::vector<MyMesh::FaceHandle> removable_faces = get_simply_removable_faces(mesh, face_type);
	if (removable_faces.size()==0){
		//std::cout<<"nothing to remove"<<std::endl;
		return false;
	}
	kT = (*mesh_props).kT;
	R_add = (*mesh_props).R_add;	
	v_add = 4.0*M_PI*R_add*R_add*R_add/3.0;
	D = (*mesh_props).DOS;
    //---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------	

	//int which = int(rand() % removable_faces.size());
	int which = randint( removable_faces.size(), eng);

	mu=face_props[ removable_faces[which] ].mu;		
	//need to find hanging vertex to find affected edges; or just use the faces' edges
	MyMesh::VertexHandle v0, v1, v2;
	MyMesh::HalfedgeHandle he;
	for (MyMesh::FaceVertexIter fv = mesh.fv_iter(removable_faces[which]); fv.is_valid(); ++fv ){
		//std::cout<<"valence "<<mesh.valence(*fv)<<std::endl;
		if (mesh.valence(*fv)==2){
			v0 = *fv;
		}
	}

	he = mesh.halfedge_handle(v0);
	if (mesh.is_boundary(he)){
		he = mesh.opposite_halfedge_handle(he);
	}
	if (mesh.from_vertex_handle(he)==v0){
		he = mesh.next_halfedge_handle(he);
		v1 = mesh.from_vertex_handle(he);
		v2 = mesh.to_vertex_handle(he);
	}
	if (mesh.to_vertex_handle(he)==v0){
		he = mesh.prev_halfedge_handle(he);
		v1 = mesh.from_vertex_handle(he);
		v2 = mesh.to_vertex_handle(he);
	}	

    std::vector<MyMesh::EdgeHandle> affected_edges = get_affected_edges(mesh, v0);
    E_before = energy(mesh, affected_edges);
	E_before+=uw.bias_potential(n_faces_before);    
   
    E_after = halfedge_stretch_energy(mesh, mesh.opposite_halfedge_handle(he), edge_props); 	
	E_after+=uw.bias_potential(n_faces_after);    

    dE = E_after - E_before;
	p = (1.0/ (D*v_add) )*exp(-(dE+mu)/kT);
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if (r<p){

    	//============== excluder ============================  	
    	//remove the to-be-removed face from its neighbors neighbor lists
    	//going through neighbors of the removable face
    	std::vector<MyMesh::FaceHandle> * vec;
    	for (auto & neighbor: face_props[ removable_faces[which] ].neighbor_list){
    		vec = &(face_props[neighbor].neighbor_list);
    		//this line supposedly removes from vector
			(*vec).erase(std::remove((*vec).begin(), (*vec).end(), removable_faces[which] ), (*vec).end());    	
		}
    	//====================================================

		//std::cout<<"REMOVED"<<std::endl;
    	mesh.delete_face(removable_faces[which]);		
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//	std::cout<<"REMOVAL acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//------------------------

		return true;
	}    

	else{

		//std::cout<<"NOT REMOVED"<<std::endl;
		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//std::cout<<"REMOVAL rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------		
		return false;
	}

}

/*
Attempts wedge removal of a new face of type face_type
Removable wedge face is chosen randomly.
Removable faces are currently searched twice, in this function and when computing p_propose.
This could be optimized by searching only outside and passing removable_faces to this function.
Doesn't seem to make much of a speedup now.
*/
bool attempt_wedge_removal(MyMesh & mesh, int face_type, UmbrellaWindow & uw, std::mt19937 & eng){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	double E_before, E_after, dE, kT, mu, r, p;
	std::vector<MyMesh::EdgeHandle> affected_edges_before;
	std::vector<MyMesh::HalfedgeHandle> affected_halfedges_after;
	MyMesh::EdgeHandle edge;

	//should be up to date, after garbage collection!!
	int n_faces_before = mesh.n_faces();
	int n_faces_after = n_faces_before - 1;	

	//find removable subunits geometrically
/*	std::vector<MyMesh::FaceHandle> removable_faces_geom = get_wedge_removable_faces(mesh);
	//filter them for the current subunit type only
	std::vector<MyMesh::FaceHandle> removable_faces;
	for (auto & rface : removable_faces_geom){
		if (face_props[rface].face_type == face_type){
			removable_faces.push_back(rface);
		}
	}
*/
	std::vector<MyMesh::FaceHandle> removable_faces = get_wedge_removable_faces(mesh, face_type);
	if (removable_faces.size()==0){
		//std::cout<<"W nothing to remove"<<std::endl;
		return false;
	}
	kT = (*mesh_props).kT;

    //---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------		

	//int which = int(rand() % removable_faces.size());
	int which = randint( removable_faces.size(), eng);	

	mu=face_props[ removable_faces[which] ].mu;	
	for (MyMesh::FaceHalfedgeIter fh = mesh.fh_iter(removable_faces[which]); fh.is_valid(); ++fh ){
		edge = mesh.edge_handle(*fh);
		affected_edges_before.push_back( edge );
		if (!mesh.is_boundary(edge)){
			affected_halfedges_after.push_back( mesh.opposite_halfedge_handle(*fh) );
		}
	}

	E_before = energy(mesh, affected_edges_before);
	E_before+=uw.bias_potential(n_faces_before);

	E_after = halfedge_stretch_energy(mesh, affected_halfedges_after[0], edge_props) + halfedge_stretch_energy(mesh, affected_halfedges_after[1], edge_props);
	E_after+=uw.bias_potential(n_faces_after);

    dE = E_after - E_before;
	p = exp(-(dE+mu)/kT);
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if (r<p){

    	//============== excluder ============================  	
    	//remove the to-be-removed face from its neighbors neighbor lists
    	//going through neighbors of the removable face
    	std::vector<MyMesh::FaceHandle> * vec;
    	for (auto & neighbor: face_props[ removable_faces[which] ].neighbor_list){
    		vec = &(face_props[neighbor].neighbor_list);
    		//this line supposedly removes from vector
			(*vec).erase(std::remove((*vec).begin(), (*vec).end(), removable_faces[which] ), (*vec).end());    	
		}
    	//====================================================

		//std::cout<<"W REMOVED"<<std::endl;
    	mesh.delete_face(removable_faces[which]);	
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//	std::cout<<"WEDGE REMOVAL acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//------------------------

		return true;
	}    

	else{

		//std::cout<<"W NOT REMOVED"<<std::endl;
		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//std::cout<<"WEDGE REMOVAL rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------			
		return false;
	}

}

/*
Attempts simple move of a vertex. 
*/
bool attempt_move(MyMesh & mesh, std::mt19937 & eng){
	//std::cout<<"attempt_move_tmp.."<<std::endl;
	//int which = int(rand() % mesh.n_vertices());
	int which = randint( mesh.n_vertices(), eng);	
	MyMesh::VertexHandle v = mesh.vertex_handle(which);
	std::vector<MyMesh::EdgeHandle> affected_edges = get_affected_edges(mesh, v);
	double E_before = elastic_energy(mesh, affected_edges);
	//double E_after_ae, dE_ae;
	//double E_before = full_elastic_energy(mesh);
	double E_after, dE, p, r, kT, d_max, dx, dy, dz;
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	kT = (*mesh_props).kT;
	d_max = (*mesh_props).d_max;

	/*dx = d_max*((2.0*rand()/(RAND_MAX + 1.0))-1.0) ;
	dy = d_max*((2.0*rand()/(RAND_MAX + 1.0))-1.0) ;
	dz = d_max*((2.0*rand()/(RAND_MAX + 1.0))-1.0) ;*/
	dx = randdouble(-d_max, d_max, eng);
	dy = randdouble(-d_max, d_max, eng);
	dz = randdouble(-d_max, d_max, eng);
	//std::cout<<"dz "<<dz<<std::endl;

	mesh.set_point(v, mesh.point(v)+MyMesh::Point(dx, dy, dz));

	//E_after = full_elastic_energy(mesh);
	E_after = elastic_energy(mesh, affected_edges);

	dE = E_after - E_before;
	//dE_ae = E_after_ae- E_before_ae;
	//std::cout<<dE<<" "<<dE_ae<<std::endl;


	//========= excluder =========================================
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	bool overlap=false;
	for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v); vf.is_valid(); ++vf ){
		//update affected faces' COMs
		mesh.calc_face_centroid(*vf, face_props[*vf].COM );

		//check affected face's overlap with its neighbors
		overlap = check_neighbor_overlap(mesh, *vf, face_props);
		if (overlap) break;
	}
	//============================================================	

	p = exp(-dE/kT);	
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);
	if ( (r<p) && (!overlap)){
		return true;
	}
	else{
		mesh.set_point(v, mesh.point(v)-MyMesh::Point(dx, dy, dz));
		//========= excluder =========================================		
		//reset affected faces' COMs back		
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		}
		//============================================================
		return false;
	}
	//return;
}


/*
Attempts type1 (wedge) fusion. Fusable pairs are searched here as well, they could be passed as a variable.
*/
bool attempt_type1_fusion(MyMesh & mesh, std::mt19937 & eng){
	//for now, find fusion pairs here every time an attempt is made
	//but this will be changed:
	//1. wedges will be found outside and passed as argument to get_type1_fusion_triplets
	//2. type1_fusion_triplets are also found outside and filtered for distance between pairs
	//3. eligible fusion triplets will be passed as arguments to this function 
	//this will be problematic because of either garbage collection or vertex replacements upon attempts
	//wedges can't really be preserved in between (rejected) attempts unless the involved wedges are updated.
	//in that case, garbage collection is only done upon accepted attempts

    int which;
    MyMesh::VertexHandle vA, vB, vC, vp, vp2;
    MyMesh::Point vAinit, vBinit;
    double p, r, E_before, E_after, dE, v_fuse, D;
    std::vector<MyMesh::EdgeHandle> affected_edges, affected_edges1, affected_edges2;

	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	float kT = (*mesh_props).kT;
	double l_fuse = (*mesh_props).l_fuse;
	//float l_fuse2 = (*mesh_props).l_fuse * (*mesh_props).l_fuse;

	v_fuse = 4.0*M_PI*l_fuse*l_fuse*l_fuse*0.125/3.0;	
	D = (*mesh_props).DOS;

	/*std::vector<MyMesh::VertexHandle> v1, v2, v3, v1f, v2f, v3f;	
    auto fusion_vectors = get_type1_fusion_triplets(mesh);
    v1 = std::get<0>(fusion_vectors);
    v2 = std::get<1>(fusion_vectors);
    v3 = std::get<2>(fusion_vectors); 
*/

	std::vector<MyMesh::VertexHandle> v1f, v2f, v3f;	
    auto fusion_vectors = get_type1_fusion_triplets(mesh, l_fuse);
    v1f = std::get<0>(fusion_vectors);
    v2f = std::get<1>(fusion_vectors);
    v3f = std::get<2>(fusion_vectors);     




/*
	for (int i=0; i<v1.size(); i++){
		if ( (mesh.point(v1[i]) - mesh.point(v3[i])).sqrnorm() < l_fuse2){			
			v1f.push_back(v1[i]);
			v3f.push_back(v3[i]);
			v2f.push_back(v2[i]);			
		}
	}
*/
	if (v1f.size()==0){
		return false;
	}

	//which = int(rand() % v1f.size());
	which = randint( v1f.size(), eng);	
	vC = v2f[which];
	//no need for this, v1f-v2f-v3f are ordered
/*	if (mesh.is_boundary(mesh.find_halfedge(v2f[which], v1f[which]))){	
		vA = v1f[which];
		vB = v3f[which];
	}
	else{
*/		
		vB = v1f[which];
		vA = v3f[which];		
//	}

    //---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------


	//std::cout<<vA<<" "<<vB<<" "<<vC<<std::endl;
	affected_edges1 = get_affected_edges(mesh, vA);
	affected_edges2 = get_affected_edges(mesh, vB);

	E_before = energy(mesh, affected_edges1) + energy(mesh, affected_edges2);

	vAinit = mesh.point(vA);
	vBinit = mesh.point(vB);

	vp = merge_vertices(mesh, vB, vA);

	if ( (vp!=vA) && (vp!=vB)){
		std::cout<<"both vertices lost in fusion"<<std::endl;
	}
	mesh.set_point(vp, 0.5*(vAinit + vBinit ));

	//========= excluder =========================================
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	bool overlap=false;
	for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vp); vf.is_valid(); ++vf ){
		//update affected faces' COMs
		mesh.calc_face_centroid(*vf, face_props[*vf].COM );

		//check affected face's overlap with its neighbors
		overlap = check_neighbor_overlap(mesh, *vf);
		if (overlap) break;
	}
	//============================================================	


	affected_edges = get_affected_edges(mesh, vp);
	E_after = energy(mesh, affected_edges);
	dE = E_after - E_before;

	p = (1.0/ (D*v_fuse) )*exp(-dE/kT);

	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if ( (r<p) && (!overlap)){
		//mesh.garbage_collection();
		//mesh.garbage_collection<std::vector<MyMesh::VertexHandle*>, std::vector<MyMesh::HalfedgeHandle*>, std::vector<MyMesh::FaceHandle*> >(handle_tracking_v, handle_tracking_h, handle_tracking_f);

		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//std::cout<<"TYPE1 FUSION acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//---------------------------

		return true;
	}
	else{

		vp2 = split_vertices(mesh, vp, vC);
	
		if (mesh.is_boundary(mesh.find_halfedge(vC, vp))){
			mesh.set_point(vp, vAinit);
			mesh.set_point(vp2, vBinit);			
		}
		else{		
			mesh.set_point(vp, vBinit);
			mesh.set_point(vp2, vAinit);	
		}

		//========= excluder =========================================		
		//reset affected faces' COMs back		
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vp); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		}

		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vp2); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		}		
		//============================================================

		
		//mesh.garbage_collection();
    	//mesh.garbage_collection<std::vector<MyMesh::VertexHandle*>, std::vector<MyMesh::HalfedgeHandle*>, std::vector<MyMesh::FaceHandle*> >(handle_tracking_v, handle_tracking_h, handle_tracking_f);

		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//std::cout<<"TYPE1 FUSION rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------


		return false;
	}	

}

/*
Attempts type2 (cracks) fusion. Fusable triplets are searched here as well, they could be passed as a variable.
*/
bool attempt_type2_fusion(MyMesh & mesh, std::mt19937 & eng){

/*	std::vector<MyMesh::VertexHandle> v1, v2, v3, v1f, v2f, v3f;
    auto fusion_vectors = get_type2_fusion_triplets(mesh);
    v1 = std::get<0>(fusion_vectors);
    v2 = std::get<1>(fusion_vectors);
    v3 = std::get<2>(fusion_vectors); 
*/
    int which;
    MyMesh::VertexHandle vA, vB, vC, vp, vp2, vD;
    MyMesh::Point vAinit, vBinit;
    double p, r, E_before, E_after, dE, v_fuse, D;
    std::vector<MyMesh::EdgeHandle> affected_edges, affected_edges1, affected_edges2;

	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	float kT = (*mesh_props).kT;
	double l_fuse = (*mesh_props).l_fuse;	
	//float l_fuse2 = (*mesh_props).l_fuse * (*mesh_props).l_fuse;
	v_fuse = 4.0*M_PI*l_fuse*l_fuse*l_fuse*0.125/3.0;
	D = (*mesh_props).DOS;

	std::vector<MyMesh::VertexHandle> v1f, v2f, v3f;
    auto fusion_vectors = get_type2_fusion_triplets(mesh, l_fuse);
    v1f = std::get<0>(fusion_vectors);
    v2f = std::get<1>(fusion_vectors);
    v3f = std::get<2>(fusion_vectors);     

		
/*
	for (int i=0; i<v1.size(); i++){
//		if (sqr_dist(mesh, v1[i], v3[i]) < l_fuse2){
		if ( (mesh.point(v1[i]) - mesh.point(v3[i])).sqrnorm() < l_fuse2){					
			v1f.push_back(v1[i]);
			v3f.push_back(v3[i]);
			v2f.push_back(v2[i]);			
		}
	}
*/
	if (v1f.size()==0){
		return false;
	}

	//which = int(rand() % v1f.size());
	which = randint( v1f.size(), eng);	
	vC = v2f[which];
	//no need for this, v1f-v2f-v3f are ordered
/*	if (mesh.is_boundary(mesh.find_halfedge(v2f[which], v1f[which]))){	
		vA = v1f[which];
		vB = v3f[which];
	}
	else{
*/		

	//---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------
		vB = v1f[which];
		vA = v3f[which];		
//	}

		vD = mesh.to_vertex_handle( mesh.next_halfedge_handle(mesh.find_halfedge(vC, vA)) );

	//std::cout<<vA<<" "<<vB<<" "<<vC<<std::endl;
	affected_edges1 = get_affected_edges(mesh, vA);
	affected_edges2 = get_affected_edges(mesh, vB);

	E_before = energy(mesh, affected_edges1) + energy(mesh, affected_edges2);

	vAinit = mesh.point(vA);
	vBinit = mesh.point(vB);

	vp = merge_vertices(mesh, vB, vA);

	//========= excluder =========================================
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	bool overlap=false;
	for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vp); vf.is_valid(); ++vf ){
		//update affected faces' COMs
		mesh.calc_face_centroid(*vf, face_props[*vf].COM );

		//check affected face's overlap with its neighbors
		overlap = check_neighbor_overlap(mesh, *vf);
		if (overlap) break;
	}
	//============================================================


	if ( (vp!=vA) && (vp!=vB)){
		std::cout<<"both vertices lost in fusion"<<std::endl;
	}
	mesh.set_point(vp, 0.5*(vAinit + vBinit ));
	affected_edges = get_affected_edges(mesh, vp);
	E_after = energy(mesh, affected_edges);
	dE = E_after - E_before;

	p = (1.0/ (D*v_fuse) )*exp(-dE/kT);

	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if ( (r<p) && (!overlap) ){

		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//std::cout<<"TYPE2 FUSION acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//---------------------------		

		return true;
	}
	else{

		vp2 = split_vertices(mesh, vp, vC, vD);
	
		if (mesh.is_boundary(mesh.find_halfedge(vC, vp))){
			mesh.set_point(vp, vAinit);
			mesh.set_point(vp2, vBinit);			
		}
		else{		
			mesh.set_point(vp, vBinit);
			mesh.set_point(vp2, vAinit);	
		}

		//========= excluder =========================================		
		//reset affected faces' COMs back		
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vp); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		}

		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vp2); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		}		
		//============================================================


		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t) > 1e-8)		
		//std::cout<<"TYPE2 FUSION rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------
		
		return false;
	}	

}

/*
Attempts type2 (crack) fission. Fission triplets are searched here as well, they could be passed as a variable.
*/
bool attempt_type2_fission(MyMesh & mesh, std::mt19937 & eng){

	std::vector<MyMesh::VertexHandle> v1, v2, v3;
	MyMesh::VertexHandle vp, vback;	
	MyMesh::Point v2init;
	int which;
    double p, r, E_before, E_after, dE, v_fuse, D;
    std::vector<MyMesh::EdgeHandle> affected_edges, affected_edges1, affected_edges2;

	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	float kT = (*mesh_props).kT;
	float l_fuse = (*mesh_props).l_fuse;
	float l_fuse2 = l_fuse * l_fuse;
	v_fuse = 4.0*M_PI*l_fuse*l_fuse*l_fuse*0.125/3.0;	
	D = (*mesh_props).DOS;	

	auto fission_vectors = get_type2_fission_triplets(mesh);
	v1 = std::get<0>(fission_vectors);
	v2 = std::get<1>(fission_vectors);
	v3 = std::get<2>(fission_vectors);


	if (v1.size()==0){
		return false;
	}
	//which = int(rand() % v1.size());
	which = randint( v1.size(), eng);	
	v2init = mesh.point(v2[which]);	

	//---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------	


	affected_edges = get_affected_edges(mesh, v2[which]);
	E_before = energy(mesh, affected_edges);	

	vp = split_vertices(mesh, v2[which], v1[which], v3[which]);

	//!!! Need to set to random from sphere
	MyMesh::Point dr = randvec(eng)*0.5*l_fuse;
	mesh.set_point(vp, mesh.point(vp)+dr);
	mesh.set_point(v2[which], mesh.point(v2[which])-dr);

	//========= excluder =========================================		
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	bool overlap=false;	
	//update affected faces' COMs	
	for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vp); vf.is_valid(); ++vf ){
		//update affected faces' COMs
		mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		//check affected face's overlap with its neighbors
		overlap = check_neighbor_overlap(mesh, *vf);
		if (overlap) break;

	}
	if (!overlap){
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v2[which]); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
			//check affected face's overlap with its neighbors
			overlap = check_neighbor_overlap(mesh, *vf);
			if (overlap) break;				
		}		
	}
	//============================================================


	affected_edges1 = get_affected_edges(mesh, v2[which]);
	affected_edges2 = get_affected_edges(mesh, vp);
	E_after = energy(mesh, affected_edges1) + energy(mesh, affected_edges2);

	dE = E_after - E_before;
	p = D*v_fuse * exp(-dE/kT);
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if ( (r<p) && (!overlap) ){

		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//std::cout<<"TYPE2 FISSION acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//---------------------------	


		return true;
	}
	else{

		vback = merge_vertices(mesh, v2[which], vp);
		mesh.set_point(vback, v2init);
		//need to put back where v2 was

		//========= excluder =========================================		
		//reset affected faces' COMs back		
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vback); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		}	
		//============================================================


		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t) > 1e-8)		
		//std::cout<<"TYPE2 FISSION rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------		

		return false;
	}


}

/*
Attempts type1 (wedge) fission. Fission pairs are searched here as well, they could be passed as a variable.
*/
bool attempt_type1_fission(MyMesh & mesh, std::mt19937 & eng){

	std::vector<MyMesh::VertexHandle> v1, v2;
	MyMesh::VertexHandle vp, vback;	
	MyMesh::Point v1init;
	int which;
    double p, r, E_before, E_after, dE, v_fuse, D;
    std::vector<MyMesh::EdgeHandle> affected_edges, affected_edges1, affected_edges2;

	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	float kT = (*mesh_props).kT;
	float l_fuse = (*mesh_props).l_fuse;	
	float l_fuse2 = l_fuse * l_fuse;

	v_fuse = 4.0*M_PI*l_fuse*l_fuse*l_fuse*0.125/3.0;	
	D = (*mesh_props).DOS;

	auto fission_vectors = get_type1_fission_pairs(mesh);
	v1 = std::get<0>(fission_vectors);
	v2 = std::get<1>(fission_vectors);

	if (v1.size()==0){
		return false;
	}

	//---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------

	//which = int(rand() % v1.size());
	which = randint( v1.size(), eng);

	v1init = mesh.point(v1[which]);	


	affected_edges = get_affected_edges(mesh, v1[which]);
	E_before = energy(mesh, affected_edges);	

	vp = split_vertices(mesh, v1[which], v2[which]);

	//!!! Need to set to random from sphere
	MyMesh::Point dr = randvec(eng)*0.5*l_fuse;	
	mesh.set_point(vp, mesh.point(vp)+dr);
	mesh.set_point(v1[which], mesh.point(v1[which])-dr);	

	//========= excluder =========================================		
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	bool overlap=false;	
	//update affected faces' COMs	
	for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vp); vf.is_valid(); ++vf ){
		//update affected faces' COMs
		mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		//check affected face's overlap with its neighbors
		overlap = check_neighbor_overlap(mesh, *vf);
		if (overlap) break;

	}
	if (!overlap){
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v1[which]); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
			//check affected face's overlap with its neighbors
			overlap = check_neighbor_overlap(mesh, *vf);
			if (overlap) break;				
		}		
	}
	//============================================================


	affected_edges1 = get_affected_edges(mesh, v1[which]);
	affected_edges2 = get_affected_edges(mesh, vp);
	E_after = energy(mesh, affected_edges1) + energy(mesh, affected_edges2);

	dE = E_after - E_before;
	p = D*v_fuse*exp(-dE/kT);
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);
	//std::cout<<"Ebefore "<<E_before<<" Eafter "<<E_after<<" p "<<p<<" r "<<r<<std::endl;

	if ( (r<p) && (!overlap)){
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//std::cout<<"TYPE1 FISSION acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//---------------------------	

		return true;
	}
	else{	
		vback = merge_vertices(mesh, v1[which], vp);
		mesh.set_point(vback, v1init);
		//need to put back where v2 was

		//========= excluder =========================================		
		//reset affected faces' COMs back		
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vback); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		}	
		//============================================================		

		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t) > 1e-8)		
		//std::cout<<"TYPE1 FISSION rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------			

		return false;
	}
}


bool attempt_edge_fusion(MyMesh & mesh, std::mt19937 & eng){
	std::vector<MyMesh::HalfedgeHandle> h1, h2, h1f, h2f;
	MyMesh::VertexHandle v1, v2, v3, v4, v1p, v2p, v3p, v4p;
	//std::vector<MyMesh::VertexHandle> v1f, v2f, v3f, v4f;	
/*    auto fusion_halfedge = get_halfedge_fusion_pairs(mesh, l_fuse);
    h1 = std::get<0>(fusion_halfedge);
    h2 = std::get<1>(fusion_halfedge);
*/


    int which;
    MyMesh::VertexHandle vn1, vn2, vB1, vB2, vA1, vA2;
    double p, r, E_before, E_after, dE, v_fuse, D;
    std::vector<MyMesh::EdgeHandle> affected_edges, affected_edges1, affected_edges2, affected_edges3, affected_edges4;
    MyMesh::EdgeHandle edge12, edge34, edge;

	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	float kT = (*mesh_props).kT;
	double l_fuse = (*mesh_props).l_fuse;
	//float l_fuse2 = (*mesh_props).l_fuse * (*mesh_props).l_fuse;    

	v_fuse = 4.0*M_PI*l_fuse*l_fuse*l_fuse*0.125/3.0;	
	D = (*mesh_props).DOS;

    auto fusion_halfedge = get_halfedge_fusion_pairs(mesh, l_fuse);
    h1f = std::get<0>(fusion_halfedge);
    h2f = std::get<1>(fusion_halfedge);    	

    //filter pairs for distance
/*    for (int i=0; i<h1.size(); i++){
    	v1 = mesh.from_vertex_handle(h1[i]);
    	v2 = mesh.to_vertex_handle(h1[i]);
    	v3 = mesh.from_vertex_handle(h2[i]);
    	v4 = mesh.to_vertex_handle(h2[i]);
    	if ( ((mesh.point(v1) - mesh.point(v4)).sqrnorm() < l_fuse2) && ((mesh.point(v2) - mesh.point(v3)).sqrnorm() < l_fuse2) ){
    		h1f.push_back(h1[i]);
    		h2f.push_back(h2[i]);
    	}
    }
*/

	if (h1f.size()==0){
		return false;
	}
    //which = int(rand() % h1f.size());
	which = randint( h1f.size(), eng);    
    v1 = mesh.from_vertex_handle(h1f[which]);
    v2 = mesh.to_vertex_handle(h1f[which]);
    v3 = mesh.from_vertex_handle(h2f[which]);
    v4 = mesh.to_vertex_handle(h2f[which]);

	vB1 = mesh.from_vertex_handle(mesh.prev_halfedge_handle(h1f[which]));
	vB2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(h2f[which]));
	vA1 = mesh.to_vertex_handle(mesh.next_halfedge_handle(h1f[which]));
	vA2 = mesh.from_vertex_handle(mesh.prev_halfedge_handle(h2f[which]));   

    //have to be careful with affected edges here because v1--v2 is the same as v2--v1 and will be counted 2 times (for each v1 and v2 too). Same for v3--v4.
    affected_edges1 = get_affected_edges(mesh, v1);
    affected_edges2 = get_affected_edges(mesh, v2);
    affected_edges3 = get_affected_edges(mesh, v3);
    affected_edges4 = get_affected_edges(mesh, v4);     

    edge12 = find_edge(mesh, v1, v2);
    edge34 = find_edge(mesh, v3, v4);


    //gather all affected edges in one vector
	affected_edges.insert( affected_edges.end(), affected_edges1.begin(), affected_edges1.end() );
	affected_edges.insert( affected_edges.end(), affected_edges2.begin(), affected_edges2.end() );
	affected_edges.insert( affected_edges.end(), affected_edges3.begin(), affected_edges3.end() );
	affected_edges.insert( affected_edges.end(), affected_edges4.begin(), affected_edges4.end() );

    //remove duplicates from affected edges; nearby vertices neighborhoods overlap thus edges are counted multiple times
	std::sort( affected_edges.begin(), affected_edges.end() );
	affected_edges.erase( unique( affected_edges.begin(), affected_edges.end() ), affected_edges.end() );
		
	//---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------

	E_before = energy(mesh, affected_edges);

	MyMesh::Point v1init = mesh.point(v1);
	MyMesh::Point v2init = mesh.point(v2);
	MyMesh::Point v3init = mesh.point(v3);
	MyMesh::Point v4init = mesh.point(v4);		

	MyMesh::Point new_point1 = (mesh.point(v1)+mesh.point(v4))*0.5;
	MyMesh::Point new_point2 = (mesh.point(v2)+mesh.point(v3))*0.5;	


	std::tie(vn1, vn2) = merge_halfedges(mesh, h1f[which], h2f[which]);	

	mesh.set_point(vn1, new_point1);
	mesh.set_point(vn2, new_point2);

	//========= excluder =========================================		
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	bool overlap=false;	
	//may double iterate faces (via vn1 or vn2 because they may share common faces) but doesn't matter, 
	//we'll just double check if they overlap with their neighbors; it'll just make it a little slower
	//update affected faces' COMs - around vn1	
	for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vn1); vf.is_valid(); ++vf ){
		//update affected faces' COMs
		mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		//check affected face's overlap with its neighbors
		overlap = check_neighbor_overlap(mesh, *vf);
		if (overlap) break;

	}

	if (!overlap){
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vn2); vf.is_valid(); ++vf ){
			//update affected faces' COMs - around vn2
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
			//check affected face's overlap with its neighbors
			overlap = check_neighbor_overlap(mesh, *vf);
			if (overlap) break;				
		}		
	}
	//============================================================	


    //have to be careful with affected edges here because v1n--vn2 is the same as vn2--vn1 and will be counted 2 times (for vn1 and vn2 too)	  
    //and not only that; other, opposite edges may be double counted as well
	affected_edges1 = get_affected_edges(mesh, vn1);
    affected_edges2 = get_affected_edges(mesh, vn2);
    affected_edges.clear();
    //gather all affected edges in one vector
	affected_edges.insert( affected_edges.end(), affected_edges1.begin(), affected_edges1.end() );
	affected_edges.insert( affected_edges.end(), affected_edges2.begin(), affected_edges2.end() );
    //remove duplicates from affected edges; nearby vertices neighborhoods overlap thus edges are counted multiple times
	std::sort( affected_edges.begin(), affected_edges.end() );
	affected_edges.erase( unique( affected_edges.begin(), affected_edges.end() ), affected_edges.end() );	

	E_after = energy(mesh, affected_edges);	

	dE = E_after - E_before;
	p = ( 1.0/ (D*D*v_fuse*v_fuse) )*exp(-dE/kT);
	//r = rand()/(RAND_MAX + 1.0);	
	r = randdouble(eng);

	if ( (r<p) && (!overlap)){
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//std::cout<<"EDGE FUSION acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<" "<<r<<std::endl;
		//---------------------------	

		return true;
	}

	else{

		//split them and put them back
		edge = find_edge(mesh, vn1, vn2);
		auto vertices = split_edge(mesh, edge);
		v1p = mesh.to_vertex_handle( find_outgoing_boundary_he(mesh, vB1) );
		v2p = mesh.from_vertex_handle( find_incoming_boundary_he(mesh, vA1) );
		v3p = mesh.to_vertex_handle( find_outgoing_boundary_he(mesh, vA2) );	
		v4p = mesh.from_vertex_handle( find_incoming_boundary_he(mesh, vB2) );	

    	mesh.set_point(v1p, v1init);
    	mesh.set_point(v2p, v2init);
    	mesh.set_point(v3p, v3init);
    	mesh.set_point(v4p, v4init);			

		//========= excluder =========================================		
		//reset affected faces' COMs back		
		std::vector<MyMesh::VertexHandle> vps = {v1p, v2p, v3p, v4p};
		for (auto & vp: vps){
			for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vp); vf.is_valid(); ++vf ){
				//update affected faces' COMs
				mesh.calc_face_centroid(*vf, face_props[*vf].COM );
			}	
		}
		//============================================================	    		

		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t) > 1e-8)		
		//std::cout<<"EDGE FUSION rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------	


    	return false;

	}

}


bool attempt_edge_fission(MyMesh & mesh, std::mt19937 & eng){
    MyMesh::VertexHandle v1, v1new, v2, v2new, vn1, vn2, vA, vB;
    MyMesh::Point v1init, v2init;
    MyMesh::HalfedgeHandle h1, h2, he;
	MyMesh::EdgeHandle edge;   
    std::vector<MyMesh::EdgeHandle> affected_edges, affected_edges1, affected_edges2, affected_edges3, affected_edges4;	 
    double p, r, E_before, E_after, dE, v_fuse, D;	

	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	float kT = (*mesh_props).kT;
	double l_fuse = (*mesh_props).l_fuse;
	float l_fuse2 = (*mesh_props).l_fuse * (*mesh_props).l_fuse;    

	v_fuse = 4.0*M_PI*l_fuse*l_fuse*l_fuse*0.125/3.0;	
	D = (*mesh_props).DOS;	

	auto fission_edges = get_fission_edges(mesh);
	if (fission_edges.size()==0){
		return false;
	}

	int which;
    //which = int(rand() % fission_edges.size());
    which = randint( fission_edges.size(), eng);    
    edge = fission_edges[which];
    he = mesh.halfedge_handle(edge, 0);
    vA = mesh.from_vertex_handle(he);
    vB = mesh.to_vertex_handle(he);
    //get affected edges - get 2 vertices first, get before energy
    affected_edges1 = get_affected_edges(mesh, vA);
    affected_edges2 = get_affected_edges(mesh, vB);
    //gather all affected edges in one vector
	affected_edges.insert( affected_edges.end(), affected_edges1.begin(), affected_edges1.end() );    
	affected_edges.insert( affected_edges.end(), affected_edges2.begin(), affected_edges2.end() );	
    //remove duplicates from affected edges; nearby vertices neighborhoods overlap thus edges are counted multiple times
	std::sort( affected_edges.begin(), affected_edges.end() );
	affected_edges.erase( unique( affected_edges.begin(), affected_edges.end() ), affected_edges.end() );


	//---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------
	E_before = energy(mesh, affected_edges);


    auto vertices = split_edge(mesh, edge);

    v1 = std::get<0>(vertices);
    v1new = std::get<1>(vertices);
    v2 = std::get<2>(vertices); 
    v2new = std::get<3>(vertices); 

    h1 = mesh.find_halfedge(v1, v2);
    if (mesh.is_valid_handle(h1)){
	    if (!mesh.is_boundary(h1)){
	    	h1 = mesh.opposite_halfedge_handle(h1);
	    }
	    h2 = mesh.find_halfedge(v1new, v2new);
	    if (!mesh.is_boundary(h2)){
	    	h2 = mesh.opposite_halfedge_handle(h2);
	    }
	}
	else{
		h1 = mesh.find_halfedge(v1, v2new);
	    if (!mesh.is_boundary(h1)){
	    	h1 = mesh.opposite_halfedge_handle(h1);
	    }
	    h2 = mesh.find_halfedge(v1new, v2);
	    if (!mesh.is_boundary(h2)){
	    	h2 = mesh.opposite_halfedge_handle(h2);
	    }		
	}

    v1init = mesh.point( mesh.from_vertex_handle(h1) );
    v2init = mesh.point( mesh.from_vertex_handle(h2) ); 

    MyMesh::Point dr1 = randvec(eng)*0.5*l_fuse;
    MyMesh::Point dr2 = randvec(eng)*0.5*l_fuse;
    mesh.set_point(v1, mesh.point(v1)+dr1);
    mesh.set_point(v1new, mesh.point(v1new)-dr1);
    mesh.set_point(v2, mesh.point(v2)+dr2);
    mesh.set_point(v2new, mesh.point(v2new)-dr2);

	//========= excluder =========================================		
	//update affected faces' COMs back; vertices share common faces, those will be updated multiple times and that's OK	
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	bool overlap=false;			
	std::vector<MyMesh::VertexHandle> vps = {v1, v1new, v2, v2new};
	for (auto & vp: vps){
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vp); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
			overlap = check_neighbor_overlap(mesh, *vf);
			if (overlap) break;			
		}	
		if (overlap) break;
	}
	//============================================================	     

    affected_edges1 = get_affected_edges(mesh, v1);
    affected_edges2 = get_affected_edges(mesh, v1new);
    affected_edges3 = get_affected_edges(mesh, v2);
    affected_edges4 = get_affected_edges(mesh, v2new);  

    affected_edges.clear();
	affected_edges.insert( affected_edges.end(), affected_edges1.begin(), affected_edges1.end() );
	affected_edges.insert( affected_edges.end(), affected_edges2.begin(), affected_edges2.end() );
	affected_edges.insert( affected_edges.end(), affected_edges3.begin(), affected_edges3.end() );
	affected_edges.insert( affected_edges.end(), affected_edges4.begin(), affected_edges4.end() );           

    //remove duplicates from affected edges; nearby vertices neighborhoods overlap thus edges are counted multiple times
	std::sort( affected_edges.begin(), affected_edges.end() );
	affected_edges.erase( unique( affected_edges.begin(), affected_edges.end() ), affected_edges.end() );	

	E_after = energy(mesh, affected_edges);	

	dE = E_after - E_before;
	//p = 1.1+v_fuse*v_fuse*exp(-dE/kT); //wtf is this
	p = D*D*v_fuse*v_fuse*exp(-dE/kT);
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);
	if ( (r<p) && (!overlap)){
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//std::cout<<"EDGE FISSION acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<" "<<r<<std::endl;
		//---------------------------	

		return true;
	}
	else{

		std::tie(vn1, vn2) = merge_halfedges(mesh, h1, h2);	    
		mesh.set_point(vn1, v1init);
		mesh.set_point(vn2, v2init);


		//========= excluder =========================================	
		//reset affected faces (neighbors of vn1 and vn2) COM back
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vn1); vf.is_valid(); ++vf ){
			//reset back affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		}	
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(vn2); vf.is_valid(); ++vf ){
			//reset back affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		}		
		//============================================================

		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t) > 1e-8)		
		//std::cout<<"EDGE FISSION rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------			

		return false;
	}

}


/*
Attempts simple move of a vertex assuming kT=0. 
*/
bool attempt_move_athermal(MyMesh & mesh, std::mt19937 & eng){
	//std::cout<<"attempt_move_tmp.."<<std::endl;
	//int which = int(rand() % mesh.n_vertices());
	int which = randint( mesh.n_vertices(), eng);	
	MyMesh::VertexHandle v = mesh.vertex_handle(which);
	std::vector<MyMesh::EdgeHandle> affected_edges = get_affected_edges(mesh, v);
	double E_before = elastic_energy(mesh, affected_edges);
	//double E_after_ae, dE_ae;
	//double E_before = full_elastic_energy(mesh);
	double E_after, dE, p, r, kT, d_max, dx, dy, dz;
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	//kT = (*mesh_props).kT;
	d_max = (*mesh_props).d_max;

	/*dx = d_max*((2.0*rand()/(RAND_MAX + 1.0))-1.0) ;
	dy = d_max*((2.0*rand()/(RAND_MAX + 1.0))-1.0) ;
	dz = d_max*((2.0*rand()/(RAND_MAX + 1.0))-1.0) ;*/
	dx = randdouble(-d_max, d_max, eng);
	dy = randdouble(-d_max, d_max, eng);
	dz = randdouble(-d_max, d_max, eng);
	//std::cout<<"dz "<<dz<<std::endl;

	mesh.set_point(v, mesh.point(v)+MyMesh::Point(dx, dy, dz));

	//E_after = full_elastic_energy(mesh);
	E_after = elastic_energy(mesh, affected_edges);

	dE = E_after - E_before;
	//dE_ae = E_after_ae- E_before_ae;
	//std::cout<<dE<<" "<<dE_ae<<std::endl;


	//========= excluder =========================================
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	bool overlap=false;
	for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v); vf.is_valid(); ++vf ){
		//update affected faces' COMs
		mesh.calc_face_centroid(*vf, face_props[*vf].COM );

		//check affected face's overlap with its neighbors
		overlap = check_neighbor_overlap(mesh, *vf, face_props);
		if (overlap) break;
	}
	//============================================================	

	//p = exp(-dE/kT);	
	if (dE>0){
		p=0.0;
	}
	else{
		p=1.0;
	}
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);
	if ( (r<p) && (!overlap)){
		return true;
	}
	else{
		mesh.set_point(v, mesh.point(v)-MyMesh::Point(dx, dy, dz));
		//========= excluder =========================================		
		//reset affected faces' COMs back		
		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v); vf.is_valid(); ++vf ){
			//update affected faces' COMs
			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
		}
		//============================================================
		return false;
	}
	//return;
}



/*
Attempts insertion of a new face of type face_type with rotation rotational_config=0,1,2 compared to the
anchor edge of the corresponding prototype face. 
Position (the open  hole) where the new subunit is inserted is chosen randomly.
Positions (holes) are currently searched twice, in this function and when computing p_propose.
This could be optimized by searching only outside and passing holes to this function.
Doesn't seem to make much of a speedup now.
*/

bool attempt_hole_insertion(MyMesh & mesh, int face_type, int rotational_config, UmbrellaWindow & uw, std::mt19937 & eng){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	//get_wedges
	std::vector <MyMesh::VertexHandle> w1, w2, w3;
    //std::vector<MyMesh::EdgeHandle> affected_edges1, affected_edges2, affected_edges3;	
    std::vector<MyMesh::EdgeHandle> affected_edges;
    MyMesh::EdgeHandle edge1, edge2, edge3;
    double kT, E_before, E_after, dE, r, p;

	//should be up to date, after garbage collection!!
	int n_faces_before = mesh.n_faces();
	int n_faces_after = n_faces_before + 1;

	auto holes = get_hole_triplets(mesh);
    w1 = std::get<0>(holes);
    w2 = std::get<1>(holes);
    w3 = std::get<2>(holes); 
    if (w1.size()==0){
    	return false;
    }

    //---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------

	//int which = int(rand() % w1.size());
	int which = randint( w1.size(), eng );

	/*affected_edges1 = get_affected_edges(mesh, w1[which]);
	affected_edges2 = get_affected_edges(mesh, w3[which]);
	affected_edges3 = get_affected_edges(mesh, w2[which]);*/

	edge1 = find_edge(mesh, w1[which], w2[which]);
	edge2 = find_edge(mesh, w2[which], w3[which]);
	edge3 = find_edge(mesh, w3[which], w1[which]);
	affected_edges.push_back(edge1);
	affected_edges.push_back(edge2);
	affected_edges.push_back(edge3);		



	//E_before = energy(mesh, affected_edges1) + energy(mesh, affected_edges2)+energy(mesh, affected_edges3);
	E_before = energy(mesh, affected_edges);
	E_before+=uw.bias_potential(n_faces_before);	

	std::vector<MyMesh::VertexHandle> vhandles;
	vhandles.push_back(w1[which]);
	vhandles.push_back(w2[which]);
	vhandles.push_back(w3[which]);

	MyMesh::FaceHandle new_face = mesh.add_face(vhandles);

	kT = (*mesh_props).kT;

	//anchor halfedge belongs to newly added face
	MyMesh::HalfedgeHandle anchor_halfedge = mesh.find_halfedge(w1[which], w2[which]);
	for (int ir=0; ir<rotational_config; ir++){
		anchor_halfedge = mesh.next_halfedge_handle(anchor_halfedge);
	}
	clone_properties_from_prototype(mesh, face_type, new_face, anchor_halfedge);

	//========= excluder =========================================
	//update new face's COM
	mesh.calc_face_centroid(new_face, face_props[new_face].COM );
	//update new face's neighbor list
	update_neighbor_list(mesh, new_face);
	//check new face's overlap with its neighbors
	bool overlap = check_neighbor_overlap(mesh, new_face);
	//if (overlap) std::cout<<"OVERLAP "<<std::endl;
	//============================================================	

	//affected_edges1 = get_affected_edges(mesh, w1[which]);
	//affected_edges2 = get_affected_edges(mesh, w3[which]);
	//MyMesh::EdgeHandle new_edge = find_edge(mesh, w1[which], w3[which]);
	//E_after = energy(mesh, affected_edges1) + energy(mesh, affected_edges2) + edge_stretch_energy(mesh, new_edge, edge_props );
	E_after = energy(mesh, affected_edges);
	E_after+=uw.bias_potential(n_faces_after);

	dE = E_after-E_before;
	//std::cout<<"wedge dE "<<dE<<std::endl;	


	float mu=face_props[new_face].mu;	
	p = exp(-(dE-mu)/kT);
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if ( (r<p)  && (!overlap)){
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//std::cout<<"WEDGE INSERTION acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//------------------------

		//========= excluder =========================================
		//if accepted, add the new face to neighbors' neighbor lists
		for (auto & neighbor: face_props[new_face].neighbor_list){
			face_props[neighbor].neighbor_list.push_back(new_face);
		}
		//============================================================

		return true;
	}
	else{
		mesh.delete_face(new_face);
		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//std::cout<<"WEDGE INSERTION rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<" "<<overlap<<std::endl;
		//std::cout<<"OL? "<<overlap<<std::endl;
		//--------------------------------------

		return false;
	}
}


/*
Attempts hole removal of a face of type face_type
Removable hole face is chosen randomly.
Removable faces are currently searched twice, in this function and when computing p_propose.
This could be optimized by searching only outside and passing removable_faces to this function.
Doesn't seem to make much of a speedup now.
*/
bool attempt_hole_removal(MyMesh & mesh, int face_type, UmbrellaWindow & uw, std::mt19937 & eng){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	double E_before, E_after, dE, kT, mu, r, p;
	std::vector<MyMesh::EdgeHandle> affected_edges_before;
	std::vector<MyMesh::HalfedgeHandle> affected_halfedges_after;
	MyMesh::EdgeHandle edge;

	//should be up to date, after garbage collection!!
	int n_faces_before = mesh.n_faces();
	int n_faces_after = n_faces_before - 1;	

	//find removable subunits geometrically
/*	std::vector<MyMesh::FaceHandle> removable_faces_geom = get_wedge_removable_faces(mesh);
	//filter them for the current subunit type only
	std::vector<MyMesh::FaceHandle> removable_faces;
	for (auto & rface : removable_faces_geom){
		if (face_props[rface].face_type == face_type){
			removable_faces.push_back(rface);
		}
	}
*/
	std::vector<MyMesh::FaceHandle> removable_faces = get_hole_removable_faces(mesh, face_type);
	if (removable_faces.size()==0){
		//std::cout<<"W nothing to remove"<<std::endl;
		return false;
	}
	kT = (*mesh_props).kT;

    //---------------------
	//double Eb_t, Ea_t, dE_t;
	//Eb_t = full_energy(mesh);
	//----------------------		

	//int which = int(rand() % removable_faces.size());
	int which = randint( removable_faces.size(), eng);	

	mu=face_props[ removable_faces[which] ].mu;	
	for (MyMesh::FaceHalfedgeIter fh = mesh.fh_iter(removable_faces[which]); fh.is_valid(); ++fh ){
		edge = mesh.edge_handle(*fh);
		affected_edges_before.push_back( edge );
		if (!mesh.is_boundary(edge)){ //<--- condition preserved from wedge insertion, irrelevant here as none of the edges should be boundary
			affected_halfedges_after.push_back( mesh.opposite_halfedge_handle(*fh) );
		}
	}

	E_before = energy(mesh, affected_edges_before);
	E_before+=uw.bias_potential(n_faces_before);

	E_after = halfedge_stretch_energy(mesh, affected_halfedges_after[0], edge_props) + halfedge_stretch_energy(mesh, affected_halfedges_after[1], edge_props)+ halfedge_stretch_energy(mesh, affected_halfedges_after[2], edge_props);
	E_after+=uw.bias_potential(n_faces_after);

    dE = E_after - E_before;
	p = exp(-(dE+mu)/kT);
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);

	if (r<p){

    	//============== excluder ============================  	
    	//remove the to-be-removed face from its neighbors neighbor lists
    	//going through neighbors of the removable face
    	std::vector<MyMesh::FaceHandle> * vec;
    	for (auto & neighbor: face_props[ removable_faces[which] ].neighbor_list){
    		vec = &(face_props[neighbor].neighbor_list);
    		//this line supposedly removes from vector
			(*vec).erase(std::remove((*vec).begin(), (*vec).end(), removable_faces[which] ), (*vec).end());    	
		}
    	//====================================================

		//std::cout<<"W REMOVED"<<std::endl;
    	mesh.delete_face(removable_faces[which]);	
		//--------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//if (std::abs(dE_t -dE) > 1e-8)
		//	std::cout<<"WEDGE REMOVAL acc, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//------------------------

		return true;
	}    

	else{

		//std::cout<<"W NOT REMOVED"<<std::endl;
		//------------------------------------
		//Ea_t = full_energy(mesh);
		//dE_t = Ea_t - Eb_t;
		//std::cout<<"WEDGE REMOVAL rej, dE, dE_t, dE-dE_t "<<dE<<" "<<dE_t<<" "<<dE-dE_t <<std::endl;
		//--------------------------------------			
		return false;
	}

}