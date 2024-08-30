#include <iostream>
#include <stdlib.h>
#include <math.h> 
#include<limits>
// -------------------- OpenMesh
//#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
//#include <OpenMesh/Apps/Assembly/topology.hh>
//#include <OpenMesh/Apps/Assembly/geometry.hh>
#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>
#include <OpenMesh/Apps/Assembly/energy.hh>



//this and edge_stretch energy go together, they are implemented TWICE for performance
//passing edge_props by reference to avoid new edge_props instance in here
double halfedge_stretch_energy(MyMesh & mesh, MyMesh::HalfedgeHandle hedge, OpenMesh::HProp<HalfedgeProp> & edge_props){
	double l0, stretch_mod;	
	//auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	double l = mesh.calc_edge_length(hedge);
	//if changed here, change it also in edge_stretch_energy() !!!!			
	const double hc = 0.5;	
	double E=0.0;

	//boundary halfedges are undefined and bear no energy
	if (!mesh.is_boundary(hedge)){
		l0 = edge_props[hedge].l0;	
		if ( (l>l0*(1.0+hc)) || (l<l0*(1.0-hc)) ){
			return 100000.0;//not returning std::numeric_limits<double>::infinity() because other values may be added to it
		}	
		stretch_mod = edge_props[hedge].stretch_mod;
		E+=0.5*stretch_mod*(l-l0)*(l-l0);
	}
	return E;
}

//this and halfedge_stretch energy go together, they are implemented TWICE for performance
//should this be edge or halfedge stretch energy? TBD, maybe should implement both?
//maybe just the full edge energy; that way edge length only has to be computed once
//passing edge_props by reference to avoid new edge_props instance in here
double edge_stretch_energy(MyMesh & mesh, MyMesh::EdgeHandle edge, OpenMesh::HProp<HalfedgeProp> & edge_props){
	//
	MyMesh::HalfedgeHandle hedge1 = mesh.halfedge_handle(edge, 0);
	MyMesh::HalfedgeHandle hedge2 = mesh.halfedge_handle(edge, 1);

	//USE THIS SECTION FOR A SINGLE LENGTH CALCULATION FOR EACH EDGE
	double l = mesh.calc_edge_length(edge);	
	double l0, stretch_mod;
	double E=0.0;
	//hard constraint max/min length ratio
	//if changed here, change it also in halfedge_stretch_energy() !!!!
	const double hc = 0.5;
	//auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	
	//boundary halfedges are undefined and bear no energy
	if (!mesh.is_boundary(hedge1)){
		l0 = edge_props[hedge1].l0;	
		if ( (l>l0*(1.0+hc)) || (l<l0*(1.0-hc)) ){
			return 100000.0;
		}	
		stretch_mod = edge_props[hedge1].stretch_mod;
		E+=0.5*stretch_mod*(l-l0)*(l-l0);
	}
	if (!mesh.is_boundary(hedge2)){
		l0 = edge_props[hedge2].l0;	
		if ( (l>l0*(1.0+hc)) || (l<l0*(1.0-hc)) ){
			return 100000.0;
		}		
		stretch_mod = edge_props[hedge2].stretch_mod;			
		E+=0.5*stretch_mod*(l-l0)*(l-l0);
	}	
	//std::cout<<E-(halfedge_stretch_energy(mesh, hedge1) + halfedge_stretch_energy(mesh, hedge2))<<std::endl;	
	return E;

	//This is a compact version but I spelled it out above because only requires one edge length calculation
	//return halfedge_stretch_energy(mesh, hedge1) + halfedge_stretch_energy(mesh, hedge2);	
}

//passing edge_props by reference to avoid new edge_props instance in here
double edge_bending_energy(MyMesh & mesh, MyMesh::EdgeHandle edge, OpenMesh::HProp<HalfedgeProp> & edge_props){
	//no faces meet at a surface edge
	if (mesh.is_boundary(edge)){
		return 0.0;
	}
	MyMesh::HalfedgeHandle hedge1 = mesh.halfedge_handle(edge, 0);
	MyMesh::HalfedgeHandle hedge2 = mesh.halfedge_handle(edge, 1);
	//auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	int edge_type1 = edge_props[hedge1].edge_type;
	int edge_type2 = edge_props[hedge2].edge_type;

	double theta0 = 0.5*(edge_props[hedge1].theta0[edge_type2] + edge_props[hedge2].theta0[edge_type1] );
	double bend_mod = 0.5* ( edge_props[hedge1].bend_mod[edge_type2] + edge_props[hedge2].bend_mod[edge_type1] );
	double theta = mesh.calc_dihedral_angle(edge);
	//double thetap = mesh.calc_dihedral_angle(edge);
	//std::cout<<"theta "<<theta<<" "<<std::endl;
	return 0.5*bend_mod*(theta - theta0)*(theta - theta0);//bend_mod*(1.0-cos(theta - theta0));

}

//passing edge_props by reference to avoid new edge_props instance in here
double edge_binding_energy(MyMesh & mesh, MyMesh::EdgeHandle edge, OpenMesh::HProp<HalfedgeProp> & edge_props){
	//no faces meet at a surface edge
	if (mesh.is_boundary(edge)){
		return 0.0;
	}
	MyMesh::HalfedgeHandle hedge1 = mesh.halfedge_handle(edge, 0);
	MyMesh::HalfedgeHandle hedge2 = mesh.halfedge_handle(edge, 1);
	//auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	int edge_type1 = edge_props[hedge1].edge_type;
	int edge_type2 = edge_props[hedge2].edge_type;
	double e_b = 0.5*(edge_props[hedge1].e_b[edge_type2] + edge_props[hedge2].e_b[edge_type1] );
	return e_b;
}

double full_elastic_energy(MyMesh & mesh){
	double E=0.0;
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	for (MyMesh::EdgeIter e_it=mesh.edges_sbegin(); e_it!=mesh.edges_end(); ++e_it){
		E+=edge_stretch_energy(mesh, *e_it, edge_props) + edge_bending_energy(mesh, *e_it, edge_props);
	}
	
	return E;
}

double full_binding_energy(MyMesh & mesh){
	double E=0.0;
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	for (MyMesh::EdgeIter e_it=mesh.edges_sbegin(); e_it!=mesh.edges_end(); ++e_it){
		E+=edge_binding_energy(mesh, *e_it, edge_props);
	}	
	return E;
}

double full_mu_N(MyMesh & mesh){
	double mu_N = 0.0;
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	for (MyMesh::FaceIter f_it=mesh.faces_sbegin(); f_it!=mesh.faces_end(); ++f_it){
		mu_N+=face_props[*f_it].mu;
	}
	return mu_N;	
}

double full_energy(MyMesh & mesh){
	double E=0.0;
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	for (MyMesh::EdgeIter e_it=mesh.edges_sbegin(); e_it!=mesh.edges_end(); ++e_it){
		E+=edge_stretch_energy(mesh, *e_it, edge_props) + edge_bending_energy(mesh, *e_it, edge_props) + edge_binding_energy(mesh, *e_it, edge_props);
	}
	
	return E;
}

double elastic_energy(MyMesh & mesh, std::vector<MyMesh::EdgeHandle> affected_edges){
	double E=0.0;
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	for (auto & edge: affected_edges){
		E+=edge_stretch_energy(mesh, edge, edge_props) + edge_bending_energy(mesh, edge, edge_props);
	}
	return E;	
}

double energy(MyMesh & mesh, std::vector<MyMesh::EdgeHandle> affected_edges){
	double E=0.0;
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	for (auto & edge: affected_edges){
		E+=edge_stretch_energy(mesh, edge, edge_props) + edge_bending_energy(mesh, edge, edge_props) + edge_binding_energy(mesh, edge, edge_props);
	}
	return E;	
}

double energy(MyMesh & mesh, MyMesh::EdgeHandle edge){
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	return edge_stretch_energy(mesh, edge, edge_props) + edge_bending_energy(mesh, edge, edge_props) + edge_binding_energy(mesh, edge, edge_props);
}



double full_stretch_energy(MyMesh & mesh){
	double E=0.0;
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	for (MyMesh::EdgeIter e_it=mesh.edges_sbegin(); e_it!=mesh.edges_end(); ++e_it){
		E+=edge_stretch_energy(mesh, *e_it, edge_props);
	}
	
	return E;
}


double full_bending_energy(MyMesh & mesh){
	double E=0.0;
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");	
	for (MyMesh::EdgeIter e_it=mesh.edges_sbegin(); e_it!=mesh.edges_end(); ++e_it){
		E+= edge_bending_energy(mesh, *e_it, edge_props);
	}
	
	return E;
}


double einstein_solid_energy(MyMesh & mesh, double k_einstein){
	double E=0.0;
	MyMesh::Point r0, r;
	auto vertex_props = OpenMesh::VProp<VertexProp>(mesh, "vertex_props");
    for(MyMesh::VertexIter it = mesh.vertices_sbegin(); it != mesh.vertices_end(); ++it) {
    	r0 = vertex_props[*it].r0;
    	r  = mesh.point(*it);

    	E+=k_einstein*(r-r0).sqrnorm(); 
    }
	
	return E;
}