#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Apps/Assembly/topology.hh>
#include <stdlib.h>
#include <OpenMesh/Tools/SmartTagger/SmartTaggerT.hh>
/*
make && ./Build/bin/Assembly && python -c "from convert_off_to_lammps import *; convert_off_to_lammps('output.off', 'lammps_0.dat')" && ovito lammps_0.dat
*/


/*
Helper function to set the half edge type of halfedges connecting v01 -- v1
To be used for boundary edges. 
The boundary halfedge edge type is set to -1
The non-boundary halfedge (of the boundary edge) is set to value
*/
void set_halfedge_type(MyMesh & mesh, MyMesh::VertexHandle & v01, MyMesh::VertexHandle & v1, OpenMesh::HPropHandleT<int> & edge_type, int value){
  //mark halfedges
  MyMesh::HalfedgeHandle he1;
  he1 = mesh.find_halfedge(v01, v1);
  if (mesh.is_boundary(he1)){
  	mesh.property(edge_type, he1)=-1;
  	//edge_type[he1] = -1;
  	he1 = mesh.opposite_halfedge_handle(he1);
  	mesh.property(edge_type, he1) = value;
  }
  else{
  	mesh.property(edge_type, he1) = value;
  	he1 = mesh.opposite_halfedge_handle(he1);
  	mesh.property(edge_type, he1) = -1;  	
  }
	return;
}


OpenMesh::HPropHandleT<int> init_edge_types(MyMesh & mesh){

	OpenMesh::HPropHandleT<int> edge_type;
	mesh.add_property(edge_type,  "edge_type");
	mesh.property(edge_type).set_persistent(true);

	std::vector<MyMesh::FaceHandle> faces;
	std::vector<std::vector<int>> halfedges;
	std::vector<std::vector<int>> halfedges_new;	

	int ntype=0;
	//iterate over all faces
  	for(MyMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit) {
  		faces.push_back(*fit);
  		halfedges.push_back(std::vector<int>());
  		//iterate over adjacent halfedges (3) of *fit
  		for (MyMesh::FaceHalfedgeIter hit = mesh.fh_iter(*fit); hit.is_valid(); ++hit ){
  			mesh.property(edge_type, *hit)=ntype;
  			halfedges.back().push_back(ntype);
  			ntype++;
  		}
	}

	return edge_type;
}

void print_edge_types(MyMesh & mesh, OpenMesh::HPropHandleT<int> edge_type){

	std::vector<MyMesh::FaceHandle> faces;
	std::vector<std::vector<int>> halfedges;
	std::vector<std::vector<MyMesh::VertexHandle>> v_from;
	std::vector<std::vector<MyMesh::VertexHandle>> v_to;		
  	for(MyMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit) {
  		faces.push_back(*fit);
  		halfedges.push_back(std::vector<int>());
  		v_from.push_back(std::vector<MyMesh::VertexHandle>());
  		v_to.push_back(std::vector<MyMesh::VertexHandle>());  		  		
  		//iterate over adjacent halfedges (3) of *fit
  		for (MyMesh::FaceHalfedgeIter hit = mesh.fh_iter(*fit); hit.is_valid(); ++hit ){
  			halfedges.back().push_back( mesh.property(edge_type, *hit) );
  			v_from.back().push_back( mesh.from_vertex_handle(*hit) );  	
  			v_to.back().push_back( mesh.to_vertex_handle(*hit) );  		  					
  		}
	}		

	//std::cout<<"-----Final edge types-------------"<<std::endl;
	for (std::size_t i=0; i<faces.size(); ++i){
		std::cout<<"face "<<faces[i]<<":";
		for (int j=0; j<3; j++){
			std::cout<<" "<<halfedges[i][j]<<" ("<<v_from[i][j]<<","<<v_to[i][j]<<")";
		}
		std::cout<<std::endl;

	}	
	return;
}

int test_vector_properties(MyMesh & mesh){

	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}	


	MyMesh::HalfedgeHandle hedge = mesh.halfedge_handle(1);
/*	OpenMesh::HPropHandleT<int> edge_type;
	mesh.add_property(edge_type,  "edge_type");
	mesh.property(edge_type).set_persistent(true);
	mesh.property(edge_type, hedge)= 2;
*/

/*	OpenMesh::HPropHandleT<std::vector<float>> preferred_angles_prop;
	std::vector<float> preferred_angles;
	mesh.add_property(preferred_angles_prop,  "preferred_angles_prop");
	mesh.property(preferred_angles_prop).set_persistent(true);	
	mesh.property(preferred_angles_prop, hedge)= preferred_angles;

	mesh.property(preferred_angles_prop, hedge).push_back(10.0);
	mesh.property(preferred_angles_prop, hedge).push_back(15.0);	


	hedge = mesh.halfedge_handle(1);
	std::cout<<mesh.property(preferred_angles_prop, hedge)[1]<<std::endl;
*/


	OpenMesh::HPropHandleT<float*> preferred_angles_prop;
	float preferred_angles[10];
	mesh.add_property(preferred_angles_prop,  "preferred_angles_prop");
	mesh.property(preferred_angles_prop).set_persistent(true);	
	mesh.property(preferred_angles_prop, hedge)= preferred_angles;

	mesh.property(preferred_angles_prop, hedge)[0]=10.0;
	mesh.property(preferred_angles_prop, hedge)[1]=15.0;	


	hedge = mesh.halfedge_handle(1);
	std::cout<<mesh.property(preferred_angles_prop, hedge)[1]<<std::endl;


  try
  {
    if ( !OpenMesh::IO::write_mesh(mesh, "persistence-check.om") )
    {
      std::cerr << "Cannot write mesh to file 'persistence-check.om'" << std::endl;
      //return 1;
    }
  }
  catch( std::exception& x )
  {
    std::cerr << x.what() << std::endl;
    //return 1;
  }

  mesh.clear();
  std::cout << "Read back mesh..";
  try
  {
    if (OpenMesh::IO::read_mesh( mesh, "persistence-check.om" ))
      std::cout << "  ok\n";
    else
    {
      std::cout << "  failed!\n";
      //return 1;
    }
    //mesh_stats(mesh, "  ");
  }
  catch( std::exception &x )
  {
    std::cerr << x.what() << std::endl;
    //return 1;
  }  

	std::cout<<"read back "<<mesh.property(preferred_angles_prop, hedge)[1]<<std::endl;

return 0;
}

/*
Test function to test crack fusion
*/
void test_crack_fusion(MyMesh & mesh){

	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}

	//define the edge_type property, stored on halfedges
	//set it to 99 for all halfedges
	OpenMesh::HPropHandleT<int> edge_type = init_edge_types(mesh);
	std::cout<<"-----Initial edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);


	//get the 4 vertices of the crack
	//MyMesh::VertexHandle v01 = mesh.vertex_handle(13);
	//MyMesh::VertexHandle v02 = mesh.vertex_handle(12);  
	MyMesh::VertexHandle v1 = mesh.vertex_handle(6);
	MyMesh::VertexHandle v2 = mesh.vertex_handle(20); 

	merge_vertices(mesh, v1, v2);

	std::cout<<"-----Final edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);

	return;
}

/*
Test function to test wedge fusion
*/
void test_wedge_fusion(MyMesh & mesh){

	//read a cracked mesh from 
	//if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/wedge.off")) 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/hex_sheet_large.off")) 		
	{
		std::cerr << "read error\n";
    	exit(1);
	}

	MyMesh::Point v1init, v2init;
	OpenMesh::HPropHandleT<int> edge_type = init_edge_types(mesh);
	std::cout<<"-----Initial edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);

	//get the 3 vertices of the crack
	MyMesh::VertexHandle v01 = mesh.vertex_handle(24); //(0);
	MyMesh::VertexHandle v1 = mesh.vertex_handle(230); //(1);
	MyMesh::VertexHandle v2 = mesh.vertex_handle(18); //(6); 
	MyMesh::VertexHandle vC;

	v1init = mesh.point(v1);
	v2init = mesh.point(v2);


	vC = merge_vertices(mesh, v1, v2);
	mesh.set_point(vC, 0.5*(v1init + v2init ));
	std::cout<<"-----Return vertex-------------"<<std::endl;	
	std::cout<<"vC, v01 "<<vC<<" "<<v01<<std::endl;
	std::cout<<"HE handles "<<mesh.halfedge_handle(vC)<<" "<<mesh.halfedge_handle(v01)<<std::endl;

	std::cout<<"-----Final edge types-------------"<<std::endl;
	//print_edge_types(mesh, edge_type);

	return;
}


void test_delete_circulators(MyMesh & mesh){
	//read a cracked mesh from 
	//if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/wedge.off")) 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/hex_sheet_large.off")) 		
	{
		std::cerr << "read error\n";
    	exit(1);
	}

	MyMesh::VertexHandle v0 = mesh.vertex_handle(110); //(0);
	for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(v0); vv_it.is_valid(); ++vv_it){
		std::cout<<*vv_it<<std::endl;

	}

std::cout<<"---faces-----"<<std::endl;
	for (MyMesh::VertexFaceIter vf_it = mesh.vf_iter(v0); vf_it.is_valid(); ++vf_it){
		std::cout<<*vf_it<<std::endl;

	}

	MyMesh::VertexHandle v1 = mesh.vertex_handle(111);	
	//mesh.delete_vertex(v1);
	MyMesh::EdgeHandle eh;
	eh = find_edge(mesh, v0, v1);
	mesh.delete_edge(eh);

std::cout<<"---after deleteion-----"<<std::endl;
	for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(v0); vv_it.is_valid(); ++vv_it){
		std::cout<<*vv_it<<std::endl;

	}
std::cout<<"---faces-----"<<std::endl;
	for (MyMesh::VertexFaceIter vf_it = mesh.vf_iter(v0); vf_it.is_valid(); ++vf_it){
		std::cout<<*vf_it<<std::endl;

	}	

std::vector<MyMesh::VertexHandle> vh;
MyMesh::VertexHandle v2 = mesh.vertex_handle(178);	
vh.push_back(v0);
vh.push_back(v2);
vh.push_back(v1);
mesh.add_face(vh);	

mesh.garbage_collection();	
	return;
}

/*
Test function to test crack fission
*/
void test_crack_fission(MyMesh & mesh){
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}


	//get the 4 vertices of the crack
	MyMesh::VertexHandle v01 = mesh.vertex_handle(13);
	MyMesh::VertexHandle v02 = mesh.vertex_handle(12);  
	MyMesh::VertexHandle v1 = mesh.vertex_handle(6);
	MyMesh::VertexHandle v2 = mesh.vertex_handle(20); 

	merge_vertices(mesh, v1, v2);

	OpenMesh::HPropHandleT<int> edge_type = init_edge_types(mesh);
	std::cout<<"-----Initial edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);


	MyMesh::VertexHandle v2f = mesh.vertex_handle(18);
	MyMesh::VertexHandle v3f = mesh.vertex_handle(12);  
	MyMesh::VertexHandle v1f = mesh.vertex_handle(6);
	MyMesh::VertexHandle newv = split_vertices(mesh, v1f, v2f, v3f);
	mesh.set_point(v1f, mesh.point(v1f)+MyMesh::Point(0.1, 0.1,  0.3));

	std::cout<<"-----Final edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);
	return;
}

/*
Test function to test wedge fission
*/
void test_wedge_fission(MyMesh & mesh){
	//read a cracked mesh from 
	//if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/wedge.off")) 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 	
	{
		std::cerr << "read error\n";
    	exit(1);
	}
	//get the 3 vertices of the crack
	MyMesh::VertexHandle v01 = mesh.vertex_handle(0);
	MyMesh::VertexHandle v1 = mesh.vertex_handle(1);
	MyMesh::VertexHandle v2 = mesh.vertex_handle(6); 


	//merge_vertices(mesh, v2, v1);


	OpenMesh::HPropHandleT<int> edge_type = init_edge_types(mesh);
	std::cout<<"-----Initial edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);

	//MyMesh::VertexHandle v2f = mesh.vertex_handle(0);
	MyMesh::VertexHandle v3f = mesh.vertex_handle(4);  
	//MyMesh::VertexHandle v1f = mesh.vertex_handle(5);

	MyMesh::VertexHandle v1f = mesh.vertex_handle(0);
	MyMesh::VertexHandle v2f = mesh.vertex_handle(19); 

	mesh.set_point(v2f, mesh.point(v2f)+MyMesh::Point(0.1, 0.1,  -0.3));	
	MyMesh::VertexHandle newv = split_vertices(mesh, v1f, v2f);
	mesh.set_point(v1f, mesh.point(v1f)+MyMesh::Point(0.1, 0.1,  0.3));	

	std::cout<<"-----Final edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);


	return;
}


/*
Test function to test finding open wedges
*/
void test_open_wedge_triplets(MyMesh & mesh){
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}
    std::vector<int> number_of_common_neighbor_vertices;
    std::vector< MyMesh::VertexHandle > wedge1, wedge2, wedge3;
    auto wedge = get_open_wedge_triplets(mesh);
    wedge1 = std::get<0>(wedge);
    wedge2 = std::get<1>(wedge);
    wedge3 = std::get<2>(wedge);    
    number_of_common_neighbor_vertices = std::get<3>(wedge);  
    for (unsigned int i=0; i<wedge1.size(); ++i){
    	std::cout<<"wedge "<<wedge1[i]<<" "<<wedge2[i]<<" "<<wedge3[i]<<" | "<<number_of_common_neighbor_vertices[i]<<std::endl;
    }
    return;

}

/*
Test function to test finding wedge fission pairs
*/
void test_type1_fission_pairs(MyMesh & mesh){
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}

	std::vector<MyMesh::VertexHandle> v1, v2;
	std::tie(v1, v2) = get_type1_fission_pairs(mesh);
	for (std::size_t i=0; i<v1.size(); ++i){
		std::cout<<"Type1 fission "<<v1[i]<<" "<<v2[i]<<std::endl;
	}
	return;
}


/*
Test function to test finding crack fission triplets
*/
void test_type2_fission_triplets(MyMesh & mesh){
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}
	std::vector<MyMesh::VertexHandle> v1, v2, v3;
    auto fusion_vectors = get_type2_fission_triplets(mesh);
    v1 = std::get<0>(fusion_vectors);
    v2 = std::get<1>(fusion_vectors);
    v3 = std::get<2>(fusion_vectors); 
    std::cout<<"----- Pre fusion -----------"<<std::endl;	
	for (std::size_t i=0; i<v1.size(); ++i){
		std::cout<<"Type2 fission "<<v1[i]<<" "<<v2[i]<<" "<<v3[i]<<std::endl;
	}



	MyMesh::VertexHandle vf1 = mesh.vertex_handle(6);
	MyMesh::VertexHandle vf2 = mesh.vertex_handle(20);   
	merge_vertices(mesh, vf1, vf2);

	fusion_vectors = get_type2_fission_triplets(mesh);
    v1 = std::get<0>(fusion_vectors);
    v2 = std::get<1>(fusion_vectors);
    v3 = std::get<2>(fusion_vectors); 

    std::cout<<"----- Post fusion -----------"<<std::endl;	
	for (std::size_t i=0; i<v1.size(); ++i){
		std::cout<<"Type2 fission "<<v1[i]<<" "<<v2[i]<<" "<<v3[i]<<std::endl;
	}

	return;
}

/*
Test finding wedge fusion triplets
*/
void test_type1_fusion_triplets(MyMesh & mesh){
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}

	std::vector<MyMesh::VertexHandle> v1, v2, v3;
    auto fusion_vectors = get_type1_fusion_triplets(mesh);
    v1 = std::get<0>(fusion_vectors);
    v2 = std::get<1>(fusion_vectors);
    v3 = std::get<2>(fusion_vectors); 

	for (std::size_t i=0; i<v1.size(); ++i){
		std::cout<<"Type1 fusion "<<v1[i]<<" "<<v2[i]<<" "<<v3[i]<<std::endl;
	}
	return;
}

/*
Test finding crack fusion triplets
*/
void test_type2_fusion_triplets(MyMesh & mesh){
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}

	std::vector<MyMesh::VertexHandle> v1, v2, v3;
    auto fusion_vectors = get_type2_fusion_triplets(mesh);
    v1 = std::get<0>(fusion_vectors);
    v2 = std::get<1>(fusion_vectors);
    v3 = std::get<2>(fusion_vectors); 

	for (std::size_t i=0; i<v1.size(); ++i){
		std::cout<<"Type2 fusion "<<v1[i]<<" "<<v2[i]<<" "<<v3[i]<<std::endl;
	}
	return;
}

/*
Test finding crack fusion triplets
*/
void test_type2_fusion_triplets_2(MyMesh & mesh){
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/open_crack.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}


	MyMesh::VertexHandle v1f = mesh.vertex_handle(19);
	MyMesh::VertexHandle v2f = mesh.vertex_handle(13);  
	MyMesh::VertexHandle v3f = mesh.vertex_handle(6);
	MyMesh::VertexHandle newv = split_vertices(mesh, v1f, v3f, v2f);
	mesh.set_point(v1f, mesh.point(v1f)+MyMesh::Point(0.1, 0.1,  0.5));


	std::vector<MyMesh::VertexHandle> v1, v2, v3;
    auto fusion_vectors = get_type2_fusion_triplets(mesh);
    v1 = std::get<0>(fusion_vectors);
    v2 = std::get<1>(fusion_vectors);
    v3 = std::get<2>(fusion_vectors); 

	for (std::size_t i=0; i<v1.size(); ++i){
		std::cout<<"Type2 fusion "<<v1[i]<<" "<<v2[i]<<" "<<v3[i]<<std::endl;
	}
	return;
}

/*
Test if mesh can break to disconnected faces via face removal
looks like it can
*/
void test_disconnected(MyMesh & mesh){
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/hex_sheet_long.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}

	MyMesh::VertexHandle v1 =  mesh.vertex_handle(17);
	MyMesh::VertexHandle v2 =  mesh.vertex_handle(25);	

	for(MyMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit) {
		for (MyMesh::FaceVertexIter vit = mesh.fv_iter(*fit); vit.is_valid(); ++vit){
			if (*vit==v1 || *vit==v2){
				mesh.delete_face(*fit);
				break;
			}
		}
	}
	mesh.garbage_collection();

	return;
}

/*
Test edge fusion
*/
/*
void test_merge_edges(MyMesh & mesh){
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/disjoint.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}	

	OpenMesh::HPropHandleT<int> edge_type;
	mesh.add_property(edge_type,  "edge_type");
	mesh.property(edge_type).set_persistent(true);

	std::vector<MyMesh::FaceHandle> faces;
	std::vector<std::vector<int>> halfedges;
	std::vector<std::vector<int>> halfedges_new;	

	int ntype=0;
	//iterate over all faces
  	for(MyMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit) {
  		faces.push_back(*fit);
  		halfedges.push_back(std::vector<int>());
  		//iterate over adjacent halfedges (3) of *fit
  		for (MyMesh::FaceHalfedgeIter hit = mesh.fh_iter(*fit); hit.is_valid(); ++hit ){
  			mesh.property(edge_type, *hit)=ntype;
  			halfedges.back().push_back(ntype);
  			ntype++;
  		}
	}

	std::cout<<"-----Initial edge types-------------"<<std::endl;
	for (std::size_t i=0; i<faces.size(); ++i){
		std::cout<<"face "<<faces[i]<<":";
		for (int j=0; j<3; j++){
			std::cout<<" "<<halfedges[i][j];
		}
		std::cout<<std::endl;

	}

	MyMesh::VertexHandle v1 =  mesh.vertex_handle(5);
	MyMesh::VertexHandle v2 =  mesh.vertex_handle(37);
	MyMesh::VertexHandle v3 =  mesh.vertex_handle(29); //(7 2);
	MyMesh::VertexHandle v4 =  mesh.vertex_handle(38); //(4 14);	

	MyMesh::EdgeHandle e1 = find_edge(mesh, v2, v1);
	MyMesh::EdgeHandle e2 = find_edge(mesh, v3, v4);
	merge_edges(mesh, e1, e2);

	std::cout<<"-----Final edge types-------------"<<std::endl;
	for (std::size_t i=0; i<faces.size(); ++i){
		std::cout<<"face "<<faces[i]<<":";
		for (int j=0; j<3; j++){
			std::cout<<" "<<halfedges[i][j];
		}
		std::cout<<std::endl;

	}	


	return;
}
*/
/*
Test edge fusion
*/
void test_merge_halfedges(MyMesh & mesh){
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/disjoint.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}	

	OpenMesh::HPropHandleT<int> edge_type = init_edge_types(mesh);
	std::cout<<"-----Initial edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);

	MyMesh::VertexHandle v1 =  mesh.vertex_handle(37);
	MyMesh::VertexHandle v2 =  mesh.vertex_handle(5);
	MyMesh::VertexHandle v3 =  mesh.vertex_handle(29); //(7 2);
	MyMesh::VertexHandle v4 =  mesh.vertex_handle(38); //(4 14);	

	MyMesh::Point new_point1 = (mesh.point(v1)+mesh.point(v4))*0.5;
	MyMesh::Point new_point2 = (mesh.point(v2)+mesh.point(v3))*0.5;	

	MyMesh::HalfedgeHandle e1 = mesh.find_halfedge(v1, v2);
	MyMesh::HalfedgeHandle e2 = mesh.find_halfedge(v3, v4);
	MyMesh::VertexHandle vn1, vn2;
	std::tie(vn1, vn2) = merge_halfedges(mesh, e1, e2);



	mesh.set_point(vn1, new_point1);	
	mesh.set_point(vn2, new_point2);	
	std::cout<<"vn1 "<<vn1<<" vn2 "<<vn2<<std::endl;

	std::cout<<"-----Final edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);


	return;
}

/*
Test edge fission
*/
void test_split_edge(MyMesh & mesh){
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/hex_sheet_long.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}	

	OpenMesh::HPropHandleT<int> edge_type = init_edge_types(mesh);
	std::cout<<"-----Initial edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);

	MyMesh::VertexHandle v1 =  mesh.vertex_handle(11);
	MyMesh::VertexHandle v2 =  mesh.vertex_handle(20);
	MyMesh::VertexHandle v1new; 
	MyMesh::VertexHandle v2new; 

	

	MyMesh::EdgeHandle edge = find_edge(mesh, v1, v2);

	//std::tie(vn1, vn2) = merge_halfedges(mesh, e1, e2);

    auto vertices = split_edge(mesh, edge);
    v1 = std::get<0>(vertices);
    v1new = std::get<1>(vertices);
    v2 = std::get<2>(vertices); 
    v2new = std::get<3>(vertices);     


	MyMesh::Point new_point1 = mesh.point(v1)+MyMesh::Point(0, 0.0,  0.5);
	MyMesh::Point new_point2 = mesh.point(v2)+MyMesh::Point(0.0, 0.0,  -0.5);


	mesh.set_point(v1, new_point1);	
	mesh.set_point(v2, new_point2);	
	//std::cout<<"vn1 "<<vn1<<" vn2 "<<vn2<<std::endl;

	std::cout<<"-----Final edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);


	return;
}

/*
To each halfedge of a face, assign each a different edge types.
Do a bunch of fusion/fission moves and see if they are conserved.
This is to verify that edge properties are properly copied in 
merge_vertices() and split_vertices()
*/
void test_edge_type_conservation(MyMesh & mesh){

	//read a cracked mesh from 
	//if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/hex_sheet_large.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}


	OpenMesh::HPropHandleT<int> edge_type = init_edge_types(mesh);
	std::cout<<"-----Initial edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);

	std::vector<MyMesh::VertexHandle> v1, v2, v1f, v2f, v3f, surface_neighbors;
	MyMesh::VertexHandle v3, v1p, v2p;
	int len, which;
	MyMesh::VertexHandle newv;
	std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > fusion_vectors;

	std::tuple<std::vector<MyMesh::HalfedgeHandle>, std::vector<MyMesh::HalfedgeHandle> > edge_fusion_vectors;
	std::vector<MyMesh::HalfedgeHandle> h1, h2;
	MyMesh::HalfedgeHandle hp;
	std::vector<MyMesh::EdgeHandle> edge_fission;

	MyMesh::VertexHandle vA1, vA2, vB1, vB2, vm;
	MyMesh::FaceHandle dummy_face;
	std::vector<MyMesh::VertexHandle> dummy_v;



	//do a bunch of topology operations
	for (int t=0; t<100000*4; t++){


		for (int m=0; m<mesh.n_vertices(); m++ ){

			vm = mesh.vertex_handle(m);
			mesh.set_point(vm, mesh.point(vm)+MyMesh::Point(0, 0.0,  0.5));

		}
		//std::cout<<"t "<<t<<std::endl;
		std::tie(v1, v2) = get_type1_fission_pairs(mesh);
		len = v1.size();
		if (len>0){

			which = int(rand() % len);
			std::cout<<"type1 fission... "<<t<<" v1"<<v1[which]<<" v2"<<v2[which]<<std::endl;				
	
			newv = split_vertices(mesh, v1[which], v2[which]);
			mesh.set_point(v1[which], mesh.point(v1[which])+MyMesh::Point(0.1, 0.1,  0.3));		
		
					
		}		
							
		if (t%12==0){
		    fusion_vectors = get_type2_fusion_triplets(mesh);
		    v1f = std::get<0>(fusion_vectors);
		    v2f = std::get<1>(fusion_vectors);
		    v3f = std::get<2>(fusion_vectors); 
			len = v1f.size();
			if (len>0){
				std::cout<<"type2 fusion... "<<t<<std::endl;	
				which = int(rand() % len);
				merge_vertices(mesh, v1f[which], v3f[which]);
				
			}	
		}   

		if (t%17==0){
		    fusion_vectors = get_type1_fusion_triplets(mesh);
		    v1f = std::get<0>(fusion_vectors);
		    v2f = std::get<1>(fusion_vectors);
		    v3f = std::get<2>(fusion_vectors); 
			len = v1f.size();
			if (len>0){
				std::cout<<"type1 fusion... "<<t<<std::endl;					
				which = int(rand() % len);
				merge_vertices(mesh, v1f[which], v3f[which]);
				
			}	
		}  	

		if (t%7==0){
		    fusion_vectors = get_type2_fission_triplets(mesh);
		    v1f = std::get<0>(fusion_vectors);
		    v2f = std::get<1>(fusion_vectors);
		    v3f = std::get<2>(fusion_vectors); 
			len = v1f.size();
			if (len>0){
				std::cout<<"type2 fission... "<<t<<std::endl;					
				which = int(rand() % len);
				split_vertices(mesh, v2f[which], v1f[which], v3f[which]);
				
			}	
		} 


		if (t%27==0){
		    edge_fusion_vectors = get_halfedge_fusion_pairs(mesh);
    		h1 = std::get<0>(edge_fusion_vectors);
    		h2 = std::get<1>(edge_fusion_vectors);	
    		len = h1.size();
    		if (len>0){
				std::cout<<"edge fusion... "<<t<<std::endl;	    			
 				which = int(rand() % len);
 				merge_halfedges(mesh, h1[which], h2[which]);   	
 						
    		}		 

    	}

		for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
			if (!mesh.is_manifold(*v_it)){
    			std::cout<<"is manifold  "<<mesh.is_manifold(*v_it)<<std::endl;
    		}
    	}

		if (t%5==0){
		    edge_fission = get_fission_edges(mesh);
    		len = edge_fission.size();
    		if (len>0){
				std::cout<<"edge fission... "<<t<<std::endl;	    			
 				which = int(rand() % len);
 				//if (t<85){
 					split_edge(mesh, edge_fission[which]);   
 				//}
 				/*else{
					std::cout<<"problem edge... "<<which<<" "<<edge_fission[which] <<std::endl;	 	
					hp = mesh.halfedge_handle(edge_fission[which],0);
					v1p = mesh.from_vertex_handle(hp);
					v2p = mesh.to_vertex_handle(hp);
					std::cout<<"problem edge vertices... "<<v1p<<" "<<v2p <<" "<<mesh.is_boundary(v1p)<<" "<<mesh.is_boundary(v2p)<<std::endl;

					surface_neighbors = find_surface_neighbors(mesh, v1p);
					for (auto &sn: surface_neighbors){
						std::cout<<"v1 surface surface_neighbors: "<<sn<<std::endl;
					}
					surface_neighbors = find_surface_neighbors(mesh, v2p);
					for (auto &sn: surface_neighbors){
						std::cout<<"v2 surface surface_neighbors: "<<sn<<std::endl;
					}			

					vA1 =  mesh.vertex_handle(16);	
					vA2 =  mesh.vertex_handle(7);			
					vB1 =  mesh.vertex_handle(14);	
					vB2 =  mesh.vertex_handle(0);	

					//std::cout<<"is manifold  "<<mesh.is_manifold(vA1)<<std::endl;							

					dummy_v.push_back(vA1);	
					dummy_v.push_back(vA2);	
					dummy_v.push_back(v1p);				
					//dummy_face = mesh.add_face(dummy_v);		
					//std::cout<<"vA1--vA2 : "<<mesh.find_halfedge(vA2, vA1)<<std::endl;	
					std::cout<<"vA1--vA2 : "<<mesh.find_halfedge(vB2, vB1)<<std::endl;	


					 //mesh.property(edge_type, *hit)						

 				}	*/
 						
    		}		 

    	}  
    	  	
	}


	std::cout<<"-----Final edge types-------------"<<std::endl;
	print_edge_types(mesh, edge_type);

	return;
}


void test_halfedge_fusion_pairs(MyMesh & mesh){
	//(ParticleIdentifier ==5) || (ParticleIdentifier ==8) || ParticleIdentifier ==2 || ParticleIdentifier==11
	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}
	MyMesh::VertexHandle vf1 = mesh.vertex_handle(6);
	MyMesh::VertexHandle vf2 = mesh.vertex_handle(20);   
	merge_vertices(mesh, vf1, vf2);

	std::vector<MyMesh::HalfedgeHandle> h1, h2;
	MyMesh::VertexHandle v1, v2, v3, v4;
    auto fusion_halfedge = get_halfedge_fusion_pairs(mesh);
    h1 = std::get<0>(fusion_halfedge);
    h2 = std::get<1>(fusion_halfedge);

    for (int i; i<h1.size(); i++){
    	v1 = mesh.from_vertex_handle(h1[i]);
    	v2 = mesh.to_vertex_handle(h1[i]);
    	v3 = mesh.from_vertex_handle(h2[i]);
    	v4 = mesh.to_vertex_handle(h2[i]);

    	std::cout<<"edge fusion "<<v1<<" "<<v2<<" "<<v3<<" "<<v4<<std::endl;
    }

}

void test_fission_edges(MyMesh & mesh){
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/cracked.off")) 
	//if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/hex_sheet_long.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}
	std::vector<MyMesh::EdgeHandle> fission_edges = get_fission_edges(mesh);
	MyMesh::HalfedgeHandle hedge;

	for (auto& edge: fission_edges){
		hedge = mesh.halfedge_handle(edge, 0);
		std::cout<<"fission edge "<<mesh.from_vertex_handle(hedge)<<" "<<mesh.to_vertex_handle(hedge)<<std::endl;
	}

}

struct Excluder
{
	int f;
};


int mainX()
{
  std::cout << "Running..." << std::endl;
  MyMesh mesh;
  // the request has to be called before a vertex/face/edge can be deleted. it grants access to the status attribute
  mesh.request_face_status();
  mesh.request_edge_status();
  mesh.request_vertex_status();    
  OpenMesh::SmartTaggerFT< MyMesh > tagger(mesh,50);


//test_delete_circulators(mesh);
// test_crack_fusion(mesh);
// test_wedge_fusion(mesh);
// test_crack_fission(mesh);
// test_wedge_fission(mesh);
// test_open_wedge_triplets(mesh);
// test_type1_fission_pairs(mesh);  
// test_type2_fission_triplets(mesh);
// test_type1_fusion_triplets(mesh);
// test_type2_fusion_triplets(mesh);
// test_type2_fusion_triplets_2(mesh);   
// test_edge_type_conservation(mesh);

//  test_merge_halfedges(mesh);  
//  test_split_edge(mesh);
//  test_disconnected(mesh);
//  test_halfedge_fusion_pairs(mesh);
//  test_fission_edges(mesh);


	//read a cracked mesh from 
	if (!OpenMesh::IO::read_mesh(mesh, "/home/btyukodi/assembly_openmesh/test_meshes/wedge.off")) 
	{
		std::cerr << "read error\n";
    	exit(1);
	}

	auto face_props = OpenMesh::FProp<Excluder>(mesh, "face_props");

// ok, this doesn't work with vectors  test_vector_properties(mesh);
MyMesh::FaceHandle fh0 = mesh.face_handle(0);
MyMesh::FaceHandle fh10 = mesh.face_handle(1);
MyMesh::FaceHandle fhd = mesh.face_handle(3);
//tagger.set_tag(fh10);
//tagger.set_tag(fh0, 4);
//mesh.delete_face(fhd);
//mesh.garbage_collection();
	int nf=0;
  for (MyMesh::FaceIter fit = mesh.faces_begin(); fit!=mesh.faces_end(); ++fit){
  	tagger.set_tag(*fit, nf);
  	face_props[*fit].f = nf;
  	nf++;
  	//std::cout<<"tagged? "<<*fit<<" " <<tagger.get_tag(*fit)<<std::endl;
  	/*if (!tagger.get_tag(*fit)){
  		mesh.delete_face(*fit);
  	}*/
  } 

  for (MyMesh::FaceIter fit = mesh.faces_begin(); fit!=mesh.faces_end(); ++fit){
  	//std::cout<<"tagged? "<<*fit<<" " <<tagger.get_tag(*fit)<<" fprop "<< face_props[*fit].f<<" "<<&face_props[*fit]<<std::endl;
  	std::cout<<"tagged? "<<*fit<<" " <<tagger.get_tag(*fit)<<" fprop "<< face_props[*fit].f<<" "<<&mesh.face(*fit)<<std::endl;  	

  }   
    	std::cout<<"----------------------"<<std::endl;
    	std::cout<<"nfaces "<<mesh.n_faces()<<std::endl;
  mesh.delete_face(fhd);
  mesh.delete_face(fh10);  
    	std::cout<<"nfaces "<<mesh.n_faces()<<std::endl;

  for (MyMesh::FaceIter fit = mesh.faces_begin(); fit!=mesh.faces_end(); ++fit){
//  	std::cout<<"tagged? "<<*fit<<" " <<tagger.get_tag(*fit)<<" fprop "<< face_props[*fit].f<<" "<<&face_props[*fit]<<std::endl;
  	std::cout<<"tagged? "<<*fit<<" " <<tagger.get_tag(*fit)<<" fprop "<< face_props[*fit].f<<" "<<&mesh.face(*fit)<<std::endl;    	

  }   

  mesh.garbage_collection();
    	std::cout<<"nfaces "<<mesh.n_faces()<<std::endl;
    	std::cout<<"----------------------"<<std::endl;  
  for (MyMesh::FaceIter fit = mesh.faces_begin(); fit!=mesh.faces_end(); ++fit){
  	//std::cout<<"tagged? "<<*fit<<" " <<tagger.get_tag(*fit)<<" fprop "<< face_props[*fit].f<<" "<<&face_props[*fit]<<std::endl;
  	std::cout<<"tagged? "<<*fit<<" " <<tagger.get_tag(*fit)<<" fprop "<< face_props[*fit].f<<" "<<&mesh.face(*fit)<<std::endl;    	

  }   

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
 