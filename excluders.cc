
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Apps/Assembly/excluders.hh>
#include <stdlib.h>
#include <iostream>


//for now only for FaceProps, etc.
//#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>


//checks overlap between two faces' excluders
//this can be extended to multiple excluder types,
//looping over all "overlappable" excluders;
//If extended, need to change mesh.calc_face_centroid(*vf, face_props[*vf].COM ) updates to "update_excluders_position() - could be just a wrapper for now"
//!!! assumes that f1 and f2 COM are up to date
bool check_overlap(MyMesh & mesh, MyMesh::FaceHandle face1, MyMesh::FaceHandle face2, OpenMesh::FProp<FaceProp> & face_props ){
	if (face1==face2){
		return false;
	}
	//auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	MyMesh::Point com1, com2;
	double R_exc1, R_exc2;
	com1 = face_props[face1].COM;
	com2 = face_props[face2].COM;
	R_exc1 = face_props[face1].R_exc;
	R_exc2 = face_props[face2].R_exc;
	//std::cout<<"OL "<<com1<<" "<<com2<<" "<<R_exc1<<" "<<R_exc2<<std::endl;
	if ( (com1-com2).sqrnorm()< (R_exc1 + R_exc2 )*(R_exc1 + R_exc2 ) ){
		return true;
	}
	return false;
}


//assumes all COM are up to date
//use for testing only
bool check_full_overlap(MyMesh & mesh){
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	for (MyMesh::FaceIter face1 = mesh.faces_sbegin(); face1!=mesh.faces_end(); ++face1){
		for (MyMesh::FaceIter face2 = mesh.faces_sbegin(); face2!=mesh.faces_end(); ++face2){
			if (check_overlap(mesh, *face1, *face2, face_props)){
				return true;
			}
		}
	}
	return false;
}

//checks overlap between a face and its neighbors
//!!! assumes COM of face and neighbors are all up to date
bool check_neighbor_overlap(MyMesh & mesh, MyMesh::FaceHandle face){
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	std::vector<MyMesh::FaceHandle> neighbor_list = face_props[face].neighbor_list;		
	for (auto & neighbor : neighbor_list){
		if (check_overlap(mesh, face, neighbor, face_props)){
			return true;
		}
	}
	return false;
}




//checks overlap between a face and its neighbors
//!!! assumes COM of face and neighbors are all up to date
//same as check_neighbor_overlap(MyMesh & mesh, MyMesh::FaceHandle face) but takes face_props as parameter by reference for some (small) performance improvement
bool check_neighbor_overlap(MyMesh & mesh, MyMesh::FaceHandle face, OpenMesh::FProp<FaceProp> & face_props){
	//auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	//std::vector<MyMesh::FaceHandle> neighbor_list = face_props[face].neighbor_list;		
	for (auto & neighbor : face_props[face].neighbor_list){
		if (check_overlap(mesh, face, neighbor, face_props)){
			return true;
		}
	}
	return false;
}



//updates neighbor list of face but NOT its neighbors?
//SHOULD ENSURE ELSEWHERE THAT ALL NEIGHBOR LISTS ARE SYMMETRIC AT ALL TIMES
void update_neighbor_list(MyMesh & mesh, MyMesh::FaceHandle face){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");		
	double L = (*mesh_props).L_neighbor; //<-- should be a mesh parameter instead of hardcoding
	MyMesh::Point com1, com2, dr;
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");

	com2 = face_props[face].COM;
	face_props[face].neighbor_list.clear();		
	for (MyMesh::FaceIter face1 = mesh.faces_sbegin(); face1!=mesh.faces_end(); ++face1){	
		if ( face!= *face1){				
			com1 = face_props[*face1].COM;
			dr = com1 - com2;
			if ( abs(dr[0])<L && abs(dr[1])<L && abs(dr[2])<L ){
				//should be done symmetrically for a 1/2 speedup
				//need to be careful not to add twice
				face_props[face].neighbor_list.push_back(*face1);
			
			}			
		}
	}	

}

//remove face from all of its neighbors' neighbor lists
void remove_from_neighbor_lists(MyMesh & mesh, MyMesh::FaceHandle face){

}

//add the new face to all necessary neighbor lists and also populate face's neighbor list
void add_to_neighbor_lists(MyMesh & mesh, MyMesh::FaceHandle face){

}

void update_all_com(MyMesh & mesh){

}



//assumes all COM are up to date
void update_full_neighbor_list(MyMesh & mesh){
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");		
	double L = (*mesh_props).L_neighbor; //<-- should be a mesh parameter instead of hardcoding
	MyMesh::Point com1, com2, dr;
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");

	for (MyMesh::FaceIter face1 = mesh.faces_sbegin(); face1!=mesh.faces_end(); ++face1){	
		face_props[*face1].neighbor_list.clear();				
		com1 = face_props[*face1].COM;
		for (MyMesh::FaceIter face2 = mesh.faces_sbegin(); face2!=mesh.faces_end(); ++face2){
			if ( *face1!= *face2){	
				com2 = face_props[*face2].COM;
				dr = com1 - com2;
				if ( abs(dr[0])<L && abs(dr[1])<L && abs(dr[2])<L ){
					//should be done symmetrically for a 1/2 speedup
					//need to be careful not to add twice
					face_props[*face1].neighbor_list.push_back(*face2);
					//std::cout<<"FACE "<<*face2<<std::endl;
				
				}
			}	
		}
	}

/*	for (MyMesh::FaceIter face1 = mesh.faces_sbegin(); face1!=mesh.faces_end(); ++face1){	
		std::cout<<"face "<<*face1<<std::endl;
		std::cout<<"neighbors ";
		for (auto & neighbor: face_props[ *face1 ].neighbor_list){
	    		std::cout<<" "<<neighbor;
	    }
		std::cout<<std::endl;
	}
*/	

}


