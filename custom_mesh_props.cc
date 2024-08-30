#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>
#include <iostream>

//ideally I would have used copy_all_properties but that one doesn't seem to work in between meshes
//proto_face_type: integer, type of the face to be taken from the proto mesh (stored in prototype_faces[proto_face_type])
//face: the face to copy the proto_face properties
//anchor_halfedge: halfedge of face that should match the corresponding anchor halfedge of the prototype. This is to keep track of rotations
void clone_properties_from_prototype(MyMesh & mesh, int proto_face_type, MyMesh::FaceHandle & face, MyMesh::HalfedgeHandle & anchor_halfedge){

	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");	
	MyMesh proto_mesh = (*mesh_props).prototypes_mesh;	
	auto proto_edge_props = OpenMesh::HProp<HalfedgeProp>(proto_mesh, "edge_props");
	auto proto_face_props = OpenMesh::FProp<FaceProp>(proto_mesh, "face_props");
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	MyMesh::FaceHandle proto_face = (*mesh_props).prototype_faces[proto_face_type];
	face_props[face] = proto_face_props[proto_face];

	MyMesh::HalfedgeHandle halfedge = anchor_halfedge;//mesh.halfedge_handle(face);
	MyMesh::HalfedgeHandle proto_halfedge = (*mesh_props).anchor_halfedges[proto_face_type];//proto_mesh.halfedge_handle(proto_face);
	for (int i=0; i<3; i++){
		edge_props[halfedge] = proto_edge_props[proto_halfedge];

		halfedge = mesh.next_halfedge_handle(halfedge);
		proto_halfedge = proto_mesh.next_halfedge_handle(proto_halfedge);
	}

	//proto_face_props[(*mesh_props).prototype_faces[0]] =  face_props[mesh.face_handle(0)];		
	return;
}


/*
Update all mesh element properties from its prototypes, except neighbor lists, based on edge and face types.
*/
void update_properties_from_prototypes(MyMesh & mesh){

	auto edge_props = OpenMesh::HProp<HalfedgeProp>(mesh, "edge_props");
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	auto mesh_props = OpenMesh::MProp<MeshProp>(mesh, "mesh_props");
	MyMesh proto_mesh = (*mesh_props).prototypes_mesh;			
	auto proto_edge_props = OpenMesh::HProp<HalfedgeProp>(proto_mesh, "edge_props");
	auto proto_face_props = OpenMesh::FProp<FaceProp>(proto_mesh, "face_props");	
	


	MyMesh::FaceHandle proto_face;
	MyMesh::HalfedgeHandle anchor_halfedge, anchor_halfedge_proto;
	int face_type, edge_type;
	for (MyMesh::FaceIter face=mesh.faces_begin(); face!=mesh.faces_end(); ++face){
		face_type = face_props[*face].face_type;

		anchor_halfedge_proto = (*mesh_props).anchor_halfedges[face_type];
		//for now, find any edge as anchor
		//NEEDS to match anchor he
		anchor_halfedge = mesh.halfedge_handle(*face);
		for (int ir=0; ir<3; ir++){
			edge_type = edge_props[anchor_halfedge].edge_type;
			if ( proto_edge_props[anchor_halfedge_proto].edge_type == edge_type){
				//std::cout<<"ROT POS FOUND facetype "<<face_type <<std::endl;
				break;
			}

			anchor_halfedge = mesh.next_halfedge_handle(anchor_halfedge);
		}

		proto_face = (*mesh_props).prototype_faces[face_type];

		face_props[*face].mu = proto_face_props[proto_face].mu;		
		face_props[*face].k_insertion = proto_face_props[proto_face].k_insertion;	
		face_props[*face].number_of_rotational_configs = proto_face_props[proto_face].number_of_rotational_configs;
		face_props[*face].R_exc = proto_face_props[proto_face].R_exc;	

		for (int i=0; i<3; i++){
			edge_props[anchor_halfedge] = proto_edge_props[anchor_halfedge_proto];

			anchor_halfedge = mesh.next_halfedge_handle(anchor_halfedge);
			anchor_halfedge_proto = proto_mesh.next_halfedge_handle(anchor_halfedge_proto);
		}	


/*		for (int i=0; i<3; i++){
					//anchor_halfedge = mesh.halfedge_handle(fh);
			std::cout<<edge_props[anchor_halfedge].e_b[0]<<std::endl;

			anchor_halfedge = mesh.next_halfedge_handle(anchor_halfedge);

		}			
		std::cout<<"----"<<std::endl;									
*/		//std::cout<<edge_props[anchor_halfedge].e_b[0]<<std::endl;
		//clone_properties_from_prototype(mesh, face_type, fh, anchor_halfedge);

		
	}


	return;
}

//little helper to switch values of two double variables
void swap_floats(float & a, float & b){
	float tmp;
	tmp = a;
	a = b;
	b = tmp;
	return;
}

/*
Need to swap properties of the prototype meshes to keep handles, etc. valid
proto_mesh1 = proto_mesh2 does a deep copy and invalidates handles.
*/

void swap_prototype_mesh_properties(MyMesh & mesh1, MyMesh & mesh2){
	auto mesh_props1 = OpenMesh::MProp<MeshProp>(mesh1, "mesh_props");	
	auto mesh_props2 = OpenMesh::MProp<MeshProp>(mesh2, "mesh_props");	
	auto edge_props1 = OpenMesh::HProp<HalfedgeProp>( (*mesh_props1).prototypes_mesh, "edge_props");
	auto face_props1 = OpenMesh::FProp<FaceProp>((*mesh_props1).prototypes_mesh, "face_props");		
	auto edge_props2 = OpenMesh::HProp<HalfedgeProp>( (*mesh_props2).prototypes_mesh, "edge_props");
	auto face_props2 = OpenMesh::FProp<FaceProp>((*mesh_props2).prototypes_mesh, "face_props");	


	//for face in prototype_faces
	//	for halfedge in neighborhood; 
	//	simply copy all properties for all (neighbor list irrelevant for proto faces)
	MyMesh::FaceHandle proto_face1, proto_face2;
	MyMesh::HalfedgeHandle he1, he2;
	FaceProp tmp_faceprop;
	HalfedgeProp tmp_edgeprop;

	for (unsigned int i=0; i<(*mesh_props1).prototype_faces.size(); i++){
		//swap prototype face properties
		proto_face1 = (*mesh_props1).prototype_faces[i];
		proto_face2 = (*mesh_props2).prototype_faces[i];
		tmp_faceprop = face_props1[ proto_face1 ];
		face_props1[ proto_face1 ] = face_props2[ proto_face2 ];
		face_props2[ proto_face2 ] = tmp_faceprop;

		//swap prototype edge properties
		he1 = (*mesh_props1).prototypes_mesh.halfedge_handle(proto_face1);
		he2 = (*mesh_props2).prototypes_mesh.halfedge_handle(proto_face2);
		//std::cout<<edge_props1[he1].e_b[0]<<" "<<edge_props2[he2].e_b[0]<<std::endl;		
		for (int ir=0; ir<3; ir++){
			if (edge_props1[he1].edge_type == edge_props2[he2].edge_type) break;
			he2 = (*mesh_props2).prototypes_mesh.next_halfedge_handle(he2);
		}

		for (int ir=0; ir<3; ir++){
			tmp_edgeprop = edge_props1[he1];
			edge_props1[he1] = edge_props2[he2];
			edge_props2[he2] = tmp_edgeprop;
			he2 = (*mesh_props2).prototypes_mesh.next_halfedge_handle(he2);
			he1 = (*mesh_props1).prototypes_mesh.next_halfedge_handle(he1);			

		}

	}

/*	for (int ir=0; ir<3; ir++){
		std::cout<<edge_props1[he1].e_b[0]<<" "<<edge_props2[he2].e_b[0]<<std::endl;
		he1 = (*mesh_props1).prototypes_mesh.next_halfedge_handle(he1);	
		he2 = (*mesh_props2).prototypes_mesh.next_halfedge_handle(he2);					
	}	
*/
}

//switch configurations between two meshes.
//assumes face_types and edge_types are the same, but parameters are different
//to be used for parallel tempering
//use each others prototype meshes to cross-copy properties
//a little tricky because need to preserve neighbor lists, COM, etc. and only copy parameters
//also need to cross-copy prototype meshes
void swap_configurations(MyMesh & mesh1, MyMesh & mesh2){
	auto mesh_props1 = OpenMesh::MProp<MeshProp>(mesh1, "mesh_props");	
	auto mesh_props2 = OpenMesh::MProp<MeshProp>(mesh2, "mesh_props");			

	//swap mesh properties
	swap_floats( (*mesh_props1).kT, (*mesh_props2).kT );
	swap_floats( (*mesh_props1).d_max, (*mesh_props2).d_max );	
	swap_floats( (*mesh_props1).l_fuse, (*mesh_props2).l_fuse );
	swap_floats( (*mesh_props1).R_add, (*mesh_props2).R_add );
	swap_floats( (*mesh_props1).L_neighbor, (*mesh_props2).L_neighbor );
	swap_floats( (*mesh_props1).k_fusion, (*mesh_props2).k_fusion );
	swap_floats( (*mesh_props1).k_fusion2, (*mesh_props2).k_fusion2 );
	swap_floats( (*mesh_props1).k_fusion_edge, (*mesh_props2).k_fusion_edge );
	swap_floats( (*mesh_props1).k_vertex_move, (*mesh_props2).k_vertex_move );	


	swap_prototype_mesh_properties(mesh1, mesh2);
/*
	(*mesh_props1).prototypes_mesh;
	(*mesh_props2).prototypes_mesh;
	std::cout<<"mesh addr "<< proto1<<" "<< proto2<<std::endl;
	MyMesh tmp_proto_mesh;
	tmp_proto_mesh = (*mesh_props1).prototypes_mesh;
	(*mesh_props1).prototypes_mesh = (*mesh_props2).prototypes_mesh;
	(*mesh_props2).prototypes_mesh = tmp_proto_mesh;
	std::cout<<"mesh addr "<< proto1<<" "<< proto2<<std::endl;


	std::vector<MyMesh::FaceHandle> prototype_faces_tmp;
	prototype_faces_tmp = (*mesh_props1).prototype_faces;
	(*mesh_props1).prototype_faces	= (*mesh_props2).prototype_faces;
	(*mesh_props2).prototype_faces	= prototype_faces_tmp;

	std::vector<MyMesh::HalfedgeHandle> prototype_anchor_halfedges_tmp;
	 prototype_anchor_halfedges_tmp = (*mesh_props1).anchor_halfedges;
	(*mesh_props1).anchor_halfedges	= (*mesh_props2).anchor_halfedges;
	(*mesh_props2).anchor_halfedges	=  prototype_anchor_halfedges_tmp;
*/
	//update face and edge properties, except neighbor lists
	update_properties_from_prototypes(mesh1);
	update_properties_from_prototypes(mesh2);

	//switch face properties
	/* Simply clone from prototype faces? No, want to keep COM and neighbor_list 
	So copy 
	face_type
	mu
	k_insertion
	R_exc 
	from prototypes

	*/

	//switch halfedge properties
	/*
	l0
	stretch_mod
	e_b
	theta0
	bend_mod
	
	*/
}