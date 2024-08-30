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
#include <OpenMesh/Apps/Assembly/monte_carlo_external_potential.hh>
#include <OpenMesh/Apps/Assembly/thermodynamic_integration_move.hh>
#include <OpenMesh/Apps/Assembly/excluders.hh>
#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>
#include <OpenMesh/Apps/Assembly/random.hh>




/*
Attempts simple move of a vertex. 

 - 3 vertex handles to the restricted vertices
 - vertex properties for the Einstein spring center coordinates
 - Einstein spring constant
 - lambda

 - need function to compute U0 and U1

*/
/*
e1 : unit vector connecting v_fix1 -- v_fix2
e2: unit vector perpendicular to e1, in the v_fix1, v_fix2, v_fix3 plane


Excluder overlap removed for speed. 
*/

bool attempt_move_thermodynamic_integration(MyMesh & mesh, std::mt19937 & eng, MyMesh::VertexHandle v_fix1, MyMesh::VertexHandle v_fix2, MyMesh::VertexHandle v_fix3, 
	MyMesh::Point e1, MyMesh::Point e2, double k_einstein, double lambda_ti){
	//std::cout<<"attempt_move_tmp.."<<std::endl;
	//int which = int(rand() % mesh.n_vertices());
	int which = randint( mesh.n_vertices(), eng);	
	MyMesh::VertexHandle v = mesh.vertex_handle(which);


	//the fixed vertex does not move
	if (v==v_fix1){
		return false;
	}

	std::vector<MyMesh::EdgeHandle> affected_edges = get_affected_edges(mesh, v);

	auto vertex_props = OpenMesh::VProp<VertexProp>(mesh, "vertex_props");
	MyMesh::Point r0=vertex_props[v].r0;
	MyMesh::Point rp = mesh.point(v);
	//MyMesh::Point e1 = vertex_props[v].e1; //unit vector connecting v_fix1 -- v_fix2
	//MyMesh::Point e2 = vertex_props[v].e2; //unit vector perpendicular to e1
	double dl2 = (rp-r0).sqrnorm();
	double dru, drv;
	
	double E_before = lambda_ti * elastic_energy(mesh, affected_edges) + (1.0-lambda_ti)*k_einstein*dl2;


	//double E_after_ae, dE_ae;
	//double E_before = full_elastic_energy(mesh);
	double E_after, dE, p, r, kT, d_max, dx, dy, dz;
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	kT = (*mesh_props).kT;
	d_max = (*mesh_props).d_max;

	/*dx = d_max*((2.0*rand()/(RAND_MAX + 1.0))-1.0) ;
	dy = d_max*((2.0*rand()/(RAND_MAX + 1.0))-1.0) ;
	dz = d_max*((2.0*rand()/(RAND_MAX + 1.0))-1.0) ;*/

	if (v==v_fix2){
		dru = randdouble(-d_max, d_max, eng);
		dx = dru*e1[0];
		dy = dru*e1[1];
		dz = dru*e1[2];

	}
	else if (v==v_fix3){
		dru = randdouble(-d_max, d_max, eng);
		drv= randdouble(-d_max, d_max, eng);		
	
		dx = dru*e1[0] + drv*e2[0];
		dy = dru*e1[1] + drv*e2[1];
		dz = dru*e1[2] + drv*e2[2];		

	}

	else{
		dx = randdouble(-d_max, d_max, eng);
		dy = randdouble(-d_max, d_max, eng);
		dz = randdouble(-d_max, d_max, eng);
	}
	//std::cout<<"dz "<<dz<<std::endl;

	//======= wall =====================
	//std::vector<MyMesh::FaceHandle> affected_faces = get_affected_faces(mesh, v);
	//auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	//@@ E_before+=externalpotential.energy(mesh, v);
	//===================================

	mesh.set_point(v, mesh.point(v)+MyMesh::Point(dx, dy, dz));


	//dE_ae = E_after_ae- E_before_ae;
	//std::cout<<dE<<" "<<dE_ae<<std::endl;


	//========= excluder =========================================
//====	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
		bool overlap=false;
//====	for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v); vf.is_valid(); ++vf ){
//====		//update affected faces' COMs
//====		mesh.calc_face_centroid(*vf, face_props[*vf].COM );
//====
//====		//check affected face's overlap with its neighbors
//====		overlap = check_neighbor_overlap(mesh, *vf, face_props);
//====		if (overlap) break;
//====	}
	//============================================================	

	r0=vertex_props[v].r0;
	rp = mesh.point(v);
	dl2 = (rp-r0).sqrnorm();

	//E_after = full_elastic_energy(mesh);
	E_after = lambda_ti * elastic_energy(mesh, affected_edges)+ (1.0-lambda_ti)*k_einstein*dl2;

	//========= wall ==============
	//#@E_after+=externalpotential.energy(mesh, v);
	//=============================

	dE = E_after - E_before;	

	//std::cout<<E_after<<" "<<E_before<<" "<<dE<<std::endl;

	p = exp(-dE/kT);	
	//r = rand()/(RAND_MAX + 1.0);
	r = randdouble(eng);
	if ( (r<p) && (!overlap)){
		return true;
	}
	else{
		mesh.set_point(v, mesh.point(v)-MyMesh::Point(dx, dy, dz));
		//========= excluder =========================================		
//===		//reset affected faces' COMs back		
//===		for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v); vf.is_valid(); ++vf ){
//===			//update affected faces' COMs
//===			mesh.calc_face_centroid(*vf, face_props[*vf].COM );
//===		}
		//============================================================
		return false;
	}
	//return;
}
