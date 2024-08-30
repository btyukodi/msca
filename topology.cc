#include <iostream>
#include <stdlib.h>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
//#include <OpenMesh/Tools/SmartTagger/SmartTaggerT.hh>
//#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>
#include <OpenMesh/Core/Mesh/Status.hh>

#include <OpenMesh/Apps/Assembly/topology.hh>
//need custom_mesh_props for face type filtering
#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>

//typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

/*
Find edge connecting v1 and v2
*/
MyMesh::EdgeHandle find_edge(MyMesh & mesh, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v2 ){

	MyMesh::HalfedgeHandle heh = mesh.find_halfedge(v1, v2);
    if (heh.is_valid()) {
        return mesh.edge_handle(heh);
    }
    else {
        return MyMesh::InvalidEdgeHandle;
    } 
}

/*
Find an outgoing boundary halfedge of v. v should only have maximum one outgoing boundary halfedge.
*/
MyMesh::HalfedgeHandle find_outgoing_boundary_he(MyMesh & mesh, MyMesh::VertexHandle & v){
	MyMesh::HalfedgeHandle heh;
	for(MyMesh::VertexOHalfedgeIter oh = mesh.voh_iter(v); oh.is_valid(); ++oh) {
		if (mesh.is_boundary(*oh)){
			return *oh;
		}
	}
    return MyMesh::InvalidHalfedgeHandle;
}

/*
Find an incoming boundary halfedge of v. v should only have maximum one incoming boundary halfedge.
*/
MyMesh::HalfedgeHandle find_incoming_boundary_he(MyMesh & mesh, MyMesh::VertexHandle & v){
	MyMesh::HalfedgeHandle heh;
	for(MyMesh::VertexIHalfedgeIter ih = mesh.vih_iter(v); ih.is_valid(); ++ih) {
		if (mesh.is_boundary(*ih)){
			return *ih;
		}
	}
    return MyMesh::InvalidHalfedgeHandle;
}


//find all (undeleted) boundary halfedges
std::vector<MyMesh::HalfedgeHandle> find_boundary_halfedges(MyMesh & mesh){
	std::vector<MyMesh::HalfedgeHandle> surface_halfedges;
	for (MyMesh::HalfedgeIter h_it = mesh.halfedges_sbegin();h_it != mesh.halfedges_end(); ++h_it ){
		if (mesh.is_boundary(*h_it)){
			surface_halfedges.push_back(*h_it);
		}
	}
	return surface_halfedges;
}

/*
Returns a vector with the common neighbors of v1 and v2, i.e. vertices that are directly connected both to v1 and v2
*/
std::vector<MyMesh::VertexHandle> find_common_neighbors(MyMesh & mesh, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v2){
  std::vector<MyMesh::VertexHandle> common_neighbors;
  //find common neighbors of v1 and v2
  //iterate over neighbors of v1
  for(MyMesh::VertexVertexIter it1 = mesh.vv_iter(v1); it1.is_valid(); ++it1) {
    //neighbors of v2
    for(MyMesh::VertexVertexIter it2 = mesh.vv_iter(v2); it2.is_valid(); ++it2) {
      if ((*it1)==(*it2)){
        common_neighbors.push_back(*it1);     
      }
    }
  }
  return common_neighbors;
}



/*
Returns nearest surface neighbors of v1
*/
std::vector<MyMesh::VertexHandle> find_surface_neighbors(MyMesh & mesh, MyMesh::VertexHandle & v1){
	std::vector<MyMesh::VertexHandle> surface_neighbors;
	MyMesh::HalfedgeHandle he;
	MyMesh::VertexHandle v_from, v_to;

	//find the surface edges out of all edges of v1
/*    for(MyMesh::VertexEdgeIter edge = mesh.ve_iter(v1); edge.is_valid(); ++edge) {
    	if (mesh.is_boundary(*edge)){
    		he = (*edge).h0();
    		v_from = mesh.from_vertex_handle(he);
    		v_to = mesh.to_vertex_handle(he);
    		if (v_from != v1){
    			surface_neighbors.push_back(v_from);    			
    		}
    		if (v_to != v1){
    			surface_neighbors.push_back(v_to);    			
    		}    		
    	}
    }
*/
	he = find_incoming_boundary_he(mesh, v1);
	surface_neighbors.push_back(mesh.from_vertex_handle(he));
	he = find_outgoing_boundary_he(mesh, v1);
	surface_neighbors.push_back(mesh.to_vertex_handle(he));	    

    if (surface_neighbors.size()!=2){
    	std::cout<<"something's wrong, v1 should have 2 surface neighbor vertices; "<<surface_neighbors.size()<<std::endl;
    }
    return surface_neighbors;
}

/*Does the merge by inserting a dummy face (or two dummy faces if crack) between (v1, v2, v0) and collapses the half edge connecting v1 and v2
returns the common neighbor of v1, v2 for wedge, and the first common neighbor for crack
*/
MyMesh::VertexHandle merge_vertices(MyMesh & mesh, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v2){
    std::vector<MyMesh::VertexHandle> common_neighbors = find_common_neighbors(mesh, v1, v2);
  
    if (common_neighbors.size()==0){
      std::cout<<"no common neighbors, vertices can't be merged"<<std::endl;
    }

    if (common_neighbors.size()>2){
      std::cout<<"more than two common neighbors, vertices can't be merged"<<std::endl;
    }  

    if (mesh.is_valid_handle (mesh.find_halfedge(v1, v2))){
      std::cout<<"v1 and v2 are connected!!"<<std::endl;    	
    }
  
    std::vector<MyMesh::VertexHandle>  dummy_face_vhandles;
   // std::vector<MyMesh::VertexHandle*> handle_tracking_v;
   // std::vector<MyMesh::FaceHandle*> handle_tracking_f;
   // std::vector<MyMesh::HalfedgeHandle*> handle_tracking_h;
    MyMesh::FaceHandle dummy_face;
    MyMesh::VertexHandle remaining_vertex;
    MyMesh::HalfedgeHandle he, he10_surface, he10_nonsurface, he20_surface, he20_nonsurface, he_tmp;

    //need to copy edge properties
    for (auto& v0: common_neighbors){
        he_tmp = mesh.find_halfedge(v0, v1);
        if (mesh.is_boundary(he_tmp)){
        	he10_surface = he_tmp;
        	he10_nonsurface = mesh.opposite_halfedge_handle(he_tmp);
        }
        else{
        	he10_nonsurface = he_tmp;
        	he10_surface = mesh.opposite_halfedge_handle(he_tmp);    	
        }
    
        he_tmp = mesh.find_halfedge(v0, v2);
        if (mesh.is_boundary(he_tmp)){
        	he20_surface = he_tmp;
        	he20_nonsurface = mesh.opposite_halfedge_handle(he_tmp);
        }
        else{
        	he20_nonsurface = he_tmp;
        	he20_surface = mesh.opposite_halfedge_handle(he_tmp);    	
        }  
  	    mesh.copy_all_properties(he10_nonsurface, he20_surface);
	    mesh.copy_all_properties(he20_nonsurface, he10_surface);      
    }  

    //for all common neighbors v0 of v1 and v2, insert dummy faces (1 for wedge, 2 for crack);
    //auto v0 = common_neighbors[0];
    for (auto& v0: common_neighbors){
    	he = mesh.find_halfedge(v0, v1);
    	//need to get the orientation of dummy face right
    	if (mesh.is_boundary(he)){
    		he = mesh.opposite_halfedge_handle(he);
    	}
        dummy_face_vhandles.clear();  
        dummy_face_vhandles.push_back( mesh.to_vertex_handle(he) );
        dummy_face_vhandles.push_back( mesh.from_vertex_handle(he) );
        dummy_face_vhandles.push_back(v2);
        dummy_face = mesh.add_face(dummy_face_vhandles);  	
  
    }

    he = mesh.find_halfedge(v1, v2);   
    remaining_vertex = mesh.to_vertex_handle(he); 

    mesh.collapse(he);
  
 	//std::cout<<"@@@ "<<mesh.deleted(v1)<<" "<<mesh.deleted(v2)<<" "<<mesh.deleted(remaining_vertex)<<" he"<<he<<std::endl; 
    //handle_tracking_v.push_back(& (common_neighbors[0]) );
    //handle_tracking_v.push_back(&v1);
    //handle_tracking_v.push_back(&v2);    
    //mesh.garbage_collection<std::vector<MyMesh::VertexHandle> >(handle_tracking);
/*    mesh.garbage_collection<std::vector<MyMesh::VertexHandle*>, std::vector<MyMesh::HalfedgeHandle*>, std::vector<MyMesh::FaceHandle*> >(handle_tracking_v, handle_tracking_h, handle_tracking_f);


 	//mesh.garbage_collection();

    if (mesh.is_valid_handle(v1)){
    	remaining_vertex = v1;
    }
    if (mesh.is_valid_handle(v2)){
    	remaining_vertex = v2;
    }  
 	//std::cout<<"**** "<<mesh.halfedge_handle(v1)<<" "<<mesh.halfedge_handle(v2)<<" "<<mesh.halfedge_handle(remaining_vertex)<<" he"<<he<<std::endl;
 	//std::cout<<"### "<<mesh.is_valid_handle(v1)<<" "<<mesh.is_valid_handle(v2)<<" "<<mesh.is_valid_handle(remaining_vertex)<<" he"<<he<<std::endl; 
*/		
    return remaining_vertex;//common_neighbors[0];

}


/*Does the wedge/crack closure by inserting a dummy face (or two dummy faces if crack) between (w1, w3, v0) and collapses the half edge connecting w1 and w3
returns the vertices w13p = w1 or w3 (whichever is left) and w2 after the merge
Note: w1-w2-w3 must be ordered wedge triplets!!
*/
/*std::tuple<MyMesh::VertexHandle, MyMesh::VertexHandle> close_wedge(MyMesh & mesh, MyMesh::VertexHandle & w1, MyMesh::VertexHandle & w2, MyMesh::VertexHandle & w3){
	//one of the common neighbors should be w2; if crack there are two common neighbors
    std::vector<MyMesh::VertexHandle> common_neighbors = find_common_neighbors(mesh, w1, w3);
  
  
    std::vector<MyMesh::VertexHandle>  dummy_face_vhandles;
    std::vector<MyMesh::VertexHandle*> handle_tracking_v;
    std::vector<MyMesh::FaceHandle*> handle_tracking_f;
    std::vector<MyMesh::HalfedgeHandle*> handle_tracking_h;
    MyMesh::FaceHandle dummy_face;
    MyMesh::VertexHandle w13p;
    MyMesh::HalfedgeHandle he, he10_surface, he10_nonsurface, he20_surface, he20_nonsurface, he_tmp;

    //need to copy edge properties
    for (auto& v0: common_neighbors){
        he_tmp = mesh.find_halfedge(v0, w1);
        if (mesh.is_boundary(he_tmp)){
        	he10_surface = he_tmp;
        	he10_nonsurface = mesh.opposite_halfedge_handle(he_tmp);
        }
        else{
        	he10_nonsurface = he_tmp;
        	he10_surface = mesh.opposite_halfedge_handle(he_tmp);    	
        }
    
        he_tmp = mesh.find_halfedge(v0, w3);
        if (mesh.is_boundary(he_tmp)){
        	he20_surface = he_tmp;
        	he20_nonsurface = mesh.opposite_halfedge_handle(he_tmp);
        }
        else{
        	he20_nonsurface = he_tmp;
        	he20_surface = mesh.opposite_halfedge_handle(he_tmp);    	
        }  
  	    mesh.copy_all_properties(he10_nonsurface, he20_surface);
	    mesh.copy_all_properties(he20_nonsurface, he10_surface);      
    }  

    //for all common neighbors v0 of v1 and v2, insert dummy faces (1 for wedge, 2 for crack);
    //auto v0 = common_neighbors[0];
    for (auto& v0: common_neighbors){
    	he = mesh.find_halfedge(v0, w1);
    	//need to get the orientation of dummy face right
    	if (mesh.is_boundary(he)){
    		he = mesh.opposite_halfedge_handle(he);
    	}
        dummy_face_vhandles.clear();  
        dummy_face_vhandles.push_back( mesh.to_vertex_handle(he) );
        dummy_face_vhandles.push_back( mesh.from_vertex_handle(he) );
        dummy_face_vhandles.push_back(w3);
        dummy_face = mesh.add_face(dummy_face_vhandles);  	
  
    }

    he = mesh.find_halfedge(w1, w3);   
    //remaining_vertex = mesh.to_vertex_handle(he); 

    mesh.collapse(he);
  
 	//std::cout<<"@@@ "<<mesh.deleted(v1)<<" "<<mesh.deleted(v2)<<" "<<mesh.deleted(remaining_vertex)<<" he"<<he<<std::endl; 
    handle_tracking_v.push_back(&w1);
    handle_tracking_v.push_back(&w2);
    handle_tracking_v.push_back(&w3);    
    //mesh.garbage_collection<std::vector<MyMesh::VertexHandle> >(handle_tracking);
    mesh.garbage_collection<std::vector<MyMesh::VertexHandle*>, std::vector<MyMesh::HalfedgeHandle*>, std::vector<MyMesh::FaceHandle*> >(handle_tracking_v, handle_tracking_h, handle_tracking_f);

 	//mesh.garbage_collection();

    if (mesh.is_valid_handle(w1)){
    	w13p = w1;
    }
    if (mesh.is_valid_handle(w3)){
    	w13p = w3;
    }  
 	//std::cout<<"**** "<<mesh.halfedge_handle(v1)<<" "<<mesh.halfedge_handle(v2)<<" "<<mesh.halfedge_handle(remaining_vertex)<<" he"<<he<<std::endl;
 	std::cout<<"### "<<mesh.is_valid_handle(w13p)<<" "<<mesh.is_valid_handle(w2)<<" "<<mesh.is_valid_handle(w1)<<" he"<<he<<std::endl; 	
    return std::make_tuple(w13p, w2);

}
*/

/*
Creates a crack along v2 - v1 - v3. Can only be used for cracks, i.e. v1 has to be non-boundary
Returns the two new vertices v1--> (v0, v1')
*/
MyMesh::VertexHandle split_vertices(MyMesh & mesh, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v2, MyMesh::VertexHandle & v3){

	if (mesh.is_boundary(v1)){
		std::cout<<"This is not a crack! Try split_vertices(mesh, v1, v2) for wedge openings"<<std::endl;
	}

	//std::vector<MyMesh::VertexHandle*> handle_tracking_v;
   //std::vector<MyMesh::FaceHandle*> handle_tracking_f;
    //std::vector<MyMesh::HalfedgeHandle*> handle_tracking_h;

    //vertex_split returns a halfedge; should be able to use to get the newly created vertex
    MyMesh::HalfedgeHandle new_halfedge = mesh.vertex_split(MyMesh::Point(mesh.point(v1)), v1,  v3, v2); 

    //v0 will be the newly inserted vertex
    MyMesh::VertexHandle v0 = mesh.from_vertex_handle(new_halfedge);

	mesh.delete_edge( find_edge(mesh, v0, v1) );    


    //remaining edge's surface halfedge (hold_surface) properties should be copied to the
    //newly created edge's non-surface halfedge (hnew_nonsurface)
    MyMesh::HalfedgeHandle hold_surface;
    MyMesh::HalfedgeHandle hnew_nonsurface;    
    MyMesh::HalfedgeHandle he;          

    //assuming halfedge v2-v1 stays
    he = mesh.find_halfedge(v2, v1);

    if (mesh.is_boundary(he)){
    	hold_surface = he;
    }
    else{
    	hold_surface = mesh.opposite_halfedge_handle(he); 	
    }

    he = mesh.find_halfedge(v2, v0);
    if (mesh.is_boundary(he)){
    	hnew_nonsurface = mesh.opposite_halfedge_handle(he);
    }
    else{
    	hnew_nonsurface = he;    	
    }


    mesh.copy_all_properties(hold_surface, hnew_nonsurface);

    //assuming halfedge v3-v1 stays
    he = mesh.find_halfedge(v3, v1);

    if (mesh.is_boundary(he)){
    	hold_surface = he;
    }
    else{
    	hold_surface = mesh.opposite_halfedge_handle(he); 	
    }

    he = mesh.find_halfedge(v3, v0);
    if (mesh.is_boundary(he)){
    	hnew_nonsurface = mesh.opposite_halfedge_handle(he);
    }
    else{
    	hnew_nonsurface = he;    	
    }

    mesh.copy_all_properties(hold_surface, hnew_nonsurface);    
		
    //mesh.garbage_collection();   
    //mesh.garbage_collection<std::vector<MyMesh::VertexHandle*>, std::vector<MyMesh::HalfedgeHandle*>, std::vector<MyMesh::FaceHandle*> >(handle_tracking_v, handle_tracking_h, handle_tracking_f);

    //std::cout<<"OH3 "<<mesh.halfedge_handle(v0)<<std::endl;     
    return v0;
}

//---------
/*
MyMesh::VertexHandle open_crack(MyMesh & mesh, MyMesh::VertexHandle & v2, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v3){

	if (mesh.is_boundary(v1)){
		std::cout<<"This is not a crack! Try split_vertices(mesh, v1, v2) for wedge openings"<<std::endl;
	}

    //vertex_split returns a halfedge; should be able to use to get the newly created vertex
    MyMesh::HalfedgeHandle new_halfedge = mesh.vertex_split(MyMesh::Point(mesh.point(v1)), v1,  v3, v2); 

    //v0 will be the newly inserted vertex
    MyMesh::VertexHandle v0 = mesh.from_vertex_handle(new_halfedge);

	mesh.delete_edge( find_edge(mesh, v0, v1) );    


    //remaining edge's surface halfedge (hold_surface) properties should be copied to the
    //newly created edge's non-surface halfedge (hnew_nonsurface)
    MyMesh::HalfedgeHandle hold_surface;
    MyMesh::HalfedgeHandle hnew_nonsurface;    
    MyMesh::HalfedgeHandle he;          

    //assuming halfedge v2-v1 stays
    he = mesh.find_halfedge(v2, v1);

    if (mesh.is_boundary(he)){
    	hold_surface = he;
    }
    else{
    	hold_surface = mesh.opposite_halfedge_handle(he); 	
    }

    he = mesh.find_halfedge(v2, v0);
    if (mesh.is_boundary(he)){
    	hnew_nonsurface = mesh.opposite_halfedge_handle(he);
    }
    else{
    	hnew_nonsurface = he;    	
    }


    mesh.copy_all_properties(hold_surface, hnew_nonsurface);

    //assuming halfedge v3-v1 stays
    he = mesh.find_halfedge(v3, v1);

    if (mesh.is_boundary(he)){
    	hold_surface = he;
    }
    else{
    	hold_surface = mesh.opposite_halfedge_handle(he); 	
    }

    he = mesh.find_halfedge(v3, v0);
    if (mesh.is_boundary(he)){
    	hnew_nonsurface = mesh.opposite_halfedge_handle(he);
    }
    else{
    	hnew_nonsurface = he;    	
    }

    mesh.copy_all_properties(hold_surface, hnew_nonsurface);    
		
    mesh.garbage_collection();   
    std::cout<<"OH3 "<<mesh.halfedge_handle(v0)<<std::endl;     
    return v0;
}
*/

/*
Creates a wedge crack along v2 - v1. Can only be used for wedges, i.e. v1 has to be boundary
Returns the newly created vertex
*/
MyMesh::VertexHandle split_vertices(MyMesh & mesh, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v2){
	/*if (!mesh.is_boundary(v1) || mesh.is_boundary(v2)){
		std::cout<<"v1 has to be boundary, v2 non-boundary"<<std::endl;
	}*/

	//insert dummy face
    std::vector<MyMesh::VertexHandle> surface_neighbors = find_surface_neighbors(mesh, v1);
    MyMesh::VertexHandle vA, vB, vC, v0, vD;
    MyMesh::HalfedgeHandle he;
    std::vector<MyMesh::VertexHandle> dummy_face_vhandles;
    MyMesh::FaceHandle dummy_face1, dummy_face2, dummy_face3, dummy_face4;    
  

    vA = surface_neighbors[0];
    vB = surface_neighbors[1];
    vC = mesh.add_vertex(0.5*(mesh.point(vA) + mesh.point(vB))  +MyMesh::Point(0.0,0.0,1.5) );
    //vD = mesh.add_vertex( mesh.point(vC) - 3*(mesh.point(vC) - mesh.point(v1) )  +MyMesh::Point(0.0,0.0,1.5));
    //std::cout<<"v1 "<<v1<<" v2 "<<v2<<" vA "<<vA<<" vB "<<vB<<" vC "<<vC<<std::endl;


    dummy_face_vhandles.clear(); 
    he = mesh.find_halfedge(vA, v1);
    //need to get the orientation of dummy face right
    //this could be done probably by the order of surface neighbors returned by find_surface_neighbors()
    if (mesh.is_boundary(he)){   	 
	    dummy_face_vhandles.push_back( vA);
	    dummy_face_vhandles.push_back( v1);
	    dummy_face_vhandles.push_back(vC);
	    dummy_face1 = mesh.add_face(dummy_face_vhandles);  	  
    }
    else{
	    dummy_face_vhandles.push_back( v1);
	    dummy_face_vhandles.push_back( vA);
	    dummy_face_vhandles.push_back(vC);
	    dummy_face1 = mesh.add_face(dummy_face_vhandles);  	      
    }

    dummy_face_vhandles.clear(); 
    he = mesh.find_halfedge(v1, vB);
    //need to get the orientation of dummy face right
    if (mesh.is_boundary(he)){   	 
	    dummy_face_vhandles.push_back( vC);
	    dummy_face_vhandles.push_back( v1);
	    dummy_face_vhandles.push_back(vB);
	    dummy_face2 = mesh.add_face(dummy_face_vhandles);  	  
    }
    else{
	    dummy_face_vhandles.push_back( v1);
	    dummy_face_vhandles.push_back( vC);
	    dummy_face_vhandles.push_back(vB);
	    dummy_face2 = mesh.add_face(dummy_face_vhandles);  	      
    } 



/*
   dummy_face_vhandles.clear(); 
    he = mesh.find_halfedge(vC, vB);
    //need to get the orientation of dummy face right
    if (mesh.is_boundary(he)){   	 
	    dummy_face_vhandles.push_back( vD);
	    dummy_face_vhandles.push_back( vC);
	    dummy_face_vhandles.push_back(vB);
	    dummy_face3 = mesh.add_face(dummy_face_vhandles);  	  
    }
    else{
	    dummy_face_vhandles.push_back( vC);
	    dummy_face_vhandles.push_back( vD);
	    dummy_face_vhandles.push_back(vB);
	    dummy_face3 = mesh.add_face(dummy_face_vhandles);  	      
    }   


    dummy_face_vhandles.clear(); 
    he = mesh.find_halfedge(vD, vA);
    //need to get the orientation of dummy face right
    if (mesh.is_boundary(he)){   	 
	    dummy_face_vhandles.push_back( vA);
	    dummy_face_vhandles.push_back( vC);
	    dummy_face_vhandles.push_back(vD);
	    dummy_face4 = mesh.add_face(dummy_face_vhandles);  	  
    }
    else{
	    dummy_face_vhandles.push_back( vC);
	    dummy_face_vhandles.push_back( vA);
	    dummy_face_vhandles.push_back(vD);
	    dummy_face4 = mesh.add_face(dummy_face_vhandles);  	      
    }   
*/

	//crack it
	v0 = split_vertices(mesh, v1, v2, vC);
    //std::cout<<"OH2 "<<mesh.halfedge_handle(v0)<<std::endl;
	//remove dummy face
	mesh.delete_face(dummy_face1);
    //mesh.garbage_collection(); 	
	mesh.delete_face(dummy_face2);
	/*mesh.delete_edge( find_edge(mesh, v0, vC), false );
    mesh.garbage_collection(); 		
	mesh.delete_edge( find_edge(mesh, v1, vC), false );	
	*/
	//mesh.delete_vertex(vC);
   // mesh.garbage_collection(); 	
   // std::cout<<"OH "<<mesh.halfedge_handle(v0)<<std::endl;
    return v0;
/*    surface_neighbors = find_surface_neighbors(mesh, v2);
    if (surface_neighbors[0]==v1){
    	return surface_neighbors[1];
    }
    else{
    	return surface_neighbors[0];
    }
   */ 
}


/*Returns wedge crack fission pair vertices
  returns a tuple of two vectors
  first vector contains surface vertices
  second vector contains corresponding inner vertices
  wedge understood to be open along v_pair_surface[i] -- v_pair_not_surface[i] edge
*/
std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > get_type1_fission_pairs(MyMesh & mesh){

    std::vector<MyMesh::VertexHandle> v_pair_surface;
    std::vector<MyMesh::VertexHandle> v_pair_not_surface;  
    std::vector<MyMesh::VertexHandle>  surface_vertices;

    //find boundary vertices
    for(MyMesh::VertexIter it = mesh.vertices_sbegin(); it != mesh.vertices_end(); ++it) {
      if(mesh.is_boundary(*it)){
        surface_vertices.push_back(*it);
      }
    }

    //iterate over all surface vertices
    for (auto& v1: surface_vertices){
      //iterate over all neighbors v2 of surface vertex v1
      for(MyMesh::VertexVertexIter v2 =  mesh.vv_iter(v1); v2.is_valid(); ++v2){
        //if v2 is not a boundary, v2 - v1 can be wedge-open
        if (!mesh.is_boundary(*v2)){
          v_pair_surface.push_back(v1);
          v_pair_not_surface.push_back(*v2);
        }
      }
    }
    return std::make_pair(v_pair_surface,v_pair_not_surface);
}


/*	type2 fission is crack opening
	all 3 vertices of the prospective crack have to be non-boundary
	returns each triplet TWICE
	could have been combined with type1_fission_pairs, but here all triplets appear twice; if duplicates are removed, it could be combined with type1 fission
*/
std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> >  get_type2_fission_triplets(MyMesh & mesh){

    std::vector <MyMesh::VertexHandle> fission1;
    std::vector <MyMesh::VertexHandle> fission2;
    std::vector <MyMesh::VertexHandle> fission3;

    //iterate over all vertices
    for (MyMesh::VertexIter v1 = mesh.vertices_sbegin(); v1 != mesh.vertices_end(); ++v1){
      if (mesh.is_boundary(*v1)){
        continue;
      }

      //iterate over v2 neighbors of v1
      for (MyMesh::VertexVertexIter v2 = mesh.vv_iter(*v1); v2.is_valid(); ++v2){
        if (mesh.is_boundary(*v2)){
          continue;
        }

        //iterate over v3 neighbors of v2
        for (MyMesh::VertexVertexIter v3 = mesh.vv_iter(*v2); v3.is_valid(); ++v3){
          if (mesh.is_boundary(*v3) || (*v3==*v1)){
            continue;
          }
          else{
            fission1.push_back(*v1);
            fission2.push_back(*v2);
            fission3.push_back(*v3);                      
          }
        }
      }
    }
	return std::make_tuple(fission1, fission2, fission3);
}


/*	returns wedges v1--v2--v3
	three vectors containing v1[i], v2[i], v3[i]
	also returns the number of common neighbor vertices of v1 and v2 to be used for wedge closure test
	!! includes crack wedges as well, separately
*/

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> , std::vector<int>>  get_open_wedge_triplets(MyMesh & mesh){
    std::vector <MyMesh::VertexHandle> wedge1, wedge2, wedge3, surface_vertices, common_neighbors;
    std::vector <int> number_of_common_neighbor_vertices; 
    MyMesh::VertexHandle w2, w3;
    MyMesh::HalfedgeHandle he;    
  
  	//find surface vertices
    for (MyMesh::VertexIter it = mesh.vertices_sbegin(); it != mesh.vertices_end(); ++it){
      if (mesh.is_boundary(*it)){
        surface_vertices.push_back(*it);
      }
    }    
  	//loop over surface vertices
    for (auto &w1 : surface_vertices){

    	he = find_outgoing_boundary_he(mesh, w1);
    	w2 = mesh.to_vertex_handle(he);

    	he = mesh.next_halfedge_handle(he);//find_outgoing_boundary_he(mesh, w2);
    	w3 = mesh.to_vertex_handle(he);

    	//only if w1 and w3 are not connected; this solves the triangle hole problem as well (which is not a wedge)
    	if (!mesh.is_valid_handle( mesh.find_halfedge(w1, w3) )){
    		wedge1.push_back(w1);
    		wedge2.push_back(w2);
    		wedge3.push_back(w3);  
			common_neighbors = find_common_neighbors(mesh, w1, w3);
			number_of_common_neighbor_vertices.push_back(common_neighbors.size());


    	}
    }
      
    return std::make_tuple(wedge1, wedge2, wedge3, number_of_common_neighbor_vertices);
}


std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> , std::vector<int>>  old_get_open_wedge_triplets(MyMesh & mesh){
    std::vector <MyMesh::VertexHandle> wedge1, wedge2, wedge3;
    std::vector <int> number_of_common_neighbor_vertices;   
    std::vector <MyMesh::VertexHandle> surface_vertices;

    bool qualifies;
    std::vector<MyMesh::VertexHandle> common_neighbors;
    MyMesh::VertexHandle vertex1, vertex2;
    MyMesh::HalfedgeHandle he1, he2;    
  
  	//find surface vertices
    for (MyMesh::VertexIter it = mesh.vertices_sbegin(); it != mesh.vertices_end(); ++it){
      if (mesh.is_boundary(*it)){
        surface_vertices.push_back(*it);
      }
    }    
  	//loop over surface vertex pairs, check if they have common neighbors attached by surface edges
    for (unsigned int ix1 =0; ix1<surface_vertices.size(); ++ix1){
      //this could possibly be made faster by iterating over the neighbors of surface_vertices[ix1]
      //but then wedges may appear multiple times
        for (unsigned int ix2 =ix1+1; ix2<surface_vertices.size(); ++ix2){
            vertex1 = surface_vertices[ix1];
            vertex2 = surface_vertices[ix2];
            //only if vertex1 and vertex2 are not connected
            if (!mesh.is_valid_handle( mesh.find_halfedge(vertex1, vertex2) )) {
            	common_neighbors = find_common_neighbors(mesh, vertex1, vertex2);
            	//have to check if vertex1--cn and vertex2--cn are both surface edges
            	for (auto &cn : common_neighbors){
            		qualifies=true;
            		he1 = mesh.find_halfedge(cn, vertex1);
            		he2 = mesh.opposite_halfedge_handle(he1);
            		//if neither halfedges of vertex1--cn are boundary
            		if ( (!mesh.is_boundary(he1)) && !(mesh.is_boundary(he2))){
            			qualifies=false;
            		}
            		he1 = mesh.find_halfedge(cn, vertex2);
            		he2 = mesh.opposite_halfedge_handle(he1);
            		//if neither halfedges of vertex2--cn are boundary
            		if ( (!mesh.is_boundary(he1)) && !(mesh.is_boundary(he2))){
            			qualifies=false;
            		}    
            		if (qualifies){
		                wedge1.push_back(vertex1);
	                    wedge2.push_back(cn);
	                    wedge3.push_back(vertex2);
	                    number_of_common_neighbor_vertices.push_back(common_neighbors.size());
                	}        		
            	}
            }
        }
    }
    return std::make_tuple(wedge1, wedge2, wedge3, number_of_common_neighbor_vertices);
}


/*	type1 fusion is wedge closure; open_wedge_triplets could be passed as parameters; might be faster if they're reused for more moves
	returns fusion1, fusion2, fusion3 vertices where fusion2 is the common neighbor and fusion1 and fusion2 can be fused
	returning fusion2 is merely convenience, fusion1 and fusion3 only have one common neighbor
*/
std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > get_type1_fusion_triplets(MyMesh & mesh){
    //if two triangles only, should return
    //get open wedges
    std::vector <MyMesh::VertexHandle> wedge1;
    std::vector <MyMesh::VertexHandle> wedge2;
    std::vector <MyMesh::VertexHandle> wedge3; 
  
    std::vector <MyMesh::VertexHandle> fusion1;
    std::vector <MyMesh::VertexHandle> fusion2;
    std::vector <MyMesh::VertexHandle> fusion3;
     
    bool qualifies;
    std::vector<int> number_of_common_neighbor_vertices;
    auto wedge = get_open_wedge_triplets(mesh);
    wedge1 = std::get<0>(wedge);
    wedge2 = std::get<1>(wedge);
    wedge3 = std::get<2>(wedge);    
    number_of_common_neighbor_vertices = std::get<3>(wedge);  
  
    MyMesh::HalfedgeHandle crack_halfedge;
    MyMesh::VertexHandle next_vertex, start_vertex;
    //not all open wedges qualify for type1 fusion; need to filter them:
    //1. neighboring faces can't fold onto each other
    //2. has to be wedge closure, not crack closure
    //3. single triangle with dangling others should not be able to close, i.e. w1 - w4 can't be connected
    // looks like all these are actually satisfied if the number of common neighbor vertices of w1 and w3 is precisely 1
  
    //loop over all open wedge pairs
    for (unsigned int i=0; i<wedge1.size(); ++i){
      qualifies=true;
  
      //3. avoid single triangle holes; the problem is that large structures can be attached to such holes via a single edge only;
      //those still identified as open wedges, but fusion should not be allowed
      if (number_of_common_neighbor_vertices[i] !=1 ){
        qualifies=false;
      }
  
      if (qualifies){
        fusion1.push_back(wedge1[i]);
        fusion2.push_back(wedge2[i]);
        fusion3.push_back(wedge3[i]);            
      }
    }
    return std::make_tuple(fusion1, fusion2, fusion3);
}

/*
Overloads get_type1_fusion_triplets(mesh); Filters for fusion pairs within l_fuse distance
*/

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > get_type1_fusion_triplets(MyMesh & mesh, double l_fuse){
	double l_fuse2 = l_fuse*l_fuse;
	std::vector<MyMesh::VertexHandle> v1, v2, v3, v1f, v2f, v3f;
	//get topologically fusable triplets
    auto fusion_vectors = get_type1_fusion_triplets(mesh);
    v1 = std::get<0>(fusion_vectors);
    v2 = std::get<1>(fusion_vectors);
    v3 = std::get<2>(fusion_vectors); 

    //filter them for within l_fuse distance
	for (unsigned int i=0; i<v1.size(); i++){
		if ( (mesh.point(v1[i]) - mesh.point(v3[i])).sqrnorm() < l_fuse2){			
			v1f.push_back(v1[i]);
			v3f.push_back(v3[i]);
			v2f.push_back(v2[i]);			
		}
	}
    return std::make_tuple(v1f, v2f, v3f);
}


/*	type2 fusion is crack closure
	this is done by filtering open_wedge_triplets for cracks
	returns each fusion triplet TWICE - the reason is that if v1-v2-v3-v4 is a crack, v1-v2-v3 and v3-v4-v1 are the same fusion move but they are detected twice
*/
std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> >  get_type2_fusion_triplets(MyMesh & mesh){

    //need to test if mesh is composed of 2 triangles only
    //then there is no type2 fusion
    std::vector <MyMesh::VertexHandle> wedge1, wedge2, wedge3, fusion1, fusion2, fusion3;
    std::vector<int> number_of_common_neighbor_vertices;
    auto wedge = get_open_wedge_triplets(mesh);
    wedge1 = std::get<0>(wedge);
    wedge2 = std::get<1>(wedge);
    wedge3 = std::get<2>(wedge);    
    number_of_common_neighbor_vertices = std::get<3>(wedge);  
  
    bool qualifies;
    MyMesh::HalfedgeHandle crack_halfedge;
    MyMesh::VertexHandle next_vertex, start_vertex;
  
    //loop over all open wedge pairs
    for (unsigned int i=0; i<wedge1.size(); ++i){
        qualifies=true;
        //crack fusion requires two common neighbors; but not all with 2 common neighbors are crack fusion
        if (number_of_common_neighbor_vertices[i] !=2 ){
          qualifies=false;
        }
        //if 2 common neighbors, have to do the crack loop 
        if (qualifies){
          //find a boundary halfedge incident to wedge2
            for (MyMesh::VertexEdgeIter edge = mesh.ve_iter(wedge2[i]); edge.is_valid(); ++edge ){
              if (mesh.is_boundary((*edge).h0())){
                crack_halfedge = (*edge).h0();
              }
              if (mesh.is_boundary((*edge).h1())){
                crack_halfedge = (*edge).h1();
              }      
            }
            //do 4 steps and check if can get back to surface_vertex
            start_vertex = mesh.from_vertex_handle(crack_halfedge);   
            //loop_vertices.push_back( start_vertex );        
            for (int j=0; j<4; ++j){
                crack_halfedge = mesh.next_halfedge_handle(crack_halfedge);
                next_vertex = mesh.from_vertex_handle(crack_halfedge);   
         
            }         
            //if back to starting vertex after 4 steps
            if (next_vertex == start_vertex){
                //will need to filter for duplicates; each crack appears 2 times meaning 2 different fusions
                //each fusion will be counted twice;
                fusion1.push_back(wedge1[i]);
                fusion2.push_back(wedge2[i]);
                fusion3.push_back(wedge3[i]);        
          }
        }
    }
    return std::make_tuple(fusion1, fusion2, fusion3);
}


std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> >  get_type2_fusion_triplets(MyMesh & mesh, double l_fuse){
	double l_fuse2 = l_fuse*l_fuse;	
	std::vector<MyMesh::VertexHandle> v1, v2, v3, v1f, v2f, v3f;
    auto fusion_vectors = get_type2_fusion_triplets(mesh);
    v1 = std::get<0>(fusion_vectors);
    v2 = std::get<1>(fusion_vectors);
    v3 = std::get<2>(fusion_vectors); 

	for (unsigned int i=0; i<v1.size(); i++){
		if ( (mesh.point(v1[i]) - mesh.point(v3[i])).sqrnorm() < l_fuse2){					
			v1f.push_back(v1[i]);
			v3f.push_back(v3[i]);
			v2f.push_back(v2[i]);			
		}
	}    
    return std::make_tuple(v1f, v2f, v3f);
}



std::vector<MyMesh::FaceHandle> get_simply_removable_faces(MyMesh & mesh){
	std::vector<MyMesh::FaceHandle> removable_faces;
	int n_surface_edges;
	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
		n_surface_edges = 0;
		for (MyMesh::FaceEdgeIter fe = mesh.fe_iter(*fit); fe.is_valid(); ++fe ){
			if (mesh.is_boundary(*fe)){
				n_surface_edges++;
			}

		}

		if (n_surface_edges==2){
			removable_faces.push_back(*fit);
		}
    }
	return removable_faces;

}


std::vector<MyMesh::FaceHandle> get_simply_removable_faces(MyMesh & mesh, int face_type){
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
	//find removable subunits topologically
	std::vector<MyMesh::FaceHandle> removable_faces_geom = get_simply_removable_faces(mesh);
	//filter them for the current subunit type only
	std::vector<MyMesh::FaceHandle> removable_faces;
	for (auto & rface : removable_faces_geom){
		if (face_props[rface].face_type == face_type){
			removable_faces.push_back(rface);
		}
	}
	return removable_faces;
}



std::vector<MyMesh::FaceHandle> get_wedge_removable_faces(MyMesh & mesh){
	std::vector<MyMesh::FaceHandle> removable_faces;
	int n_surface_edges, n_surface_vertices;
	for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
		n_surface_edges = 0;
		n_surface_vertices = 0;
		for (MyMesh::FaceEdgeIter fe = mesh.fe_iter(*fit); fe.is_valid(); ++fe ){
			if (mesh.is_boundary(*fe)){
				n_surface_edges++;
			}
		}

		for (MyMesh::FaceVertexIter fv = mesh.fv_iter(*fit); fv.is_valid(); ++fv ){
			if (mesh.is_boundary(*fv)){
				n_surface_vertices++;
			}
		}		
		//have to test to make sure structure doesn't fall apart by removing the wedge, eg. /\/\/\/\ zigzag stripe
        //looks like this is only possible if wedge hub is a non-surface hub
		if ( (n_surface_edges==1) && (n_surface_vertices==2)){
			removable_faces.push_back(*fit);
		}
    }
	return removable_faces;

}


std::vector<MyMesh::FaceHandle> get_wedge_removable_faces(MyMesh & mesh, int face_type){
	auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");
	//find removable subunits geometrically
	std::vector<MyMesh::FaceHandle> removable_faces_geom = get_wedge_removable_faces(mesh);
	//filter them for the current subunit type only
	std::vector<MyMesh::FaceHandle> removable_faces;
	for (auto & rface : removable_faces_geom){
		if (face_props[rface].face_type == face_type){
			removable_faces.push_back(rface);
		}
	}
	return removable_faces;
}


/*	(Half)Edge fusion.	
	Barely any consistency check, will fuse any two boundary halfedges;
	Will fuse hedge1.from() to hedge2.to() and hedge1.to() fo hedge2.from() 
	Returns tuple of two vertices, 
	1st vertex: hedge1 from_vertex
	2nd vertex: hedge2 from_vertex
*/
std::tuple< MyMesh::VertexHandle, MyMesh::VertexHandle > merge_halfedges(MyMesh & mesh, MyMesh::HalfedgeHandle & hedge1, MyMesh::HalfedgeHandle & hedge2){
  //both edge1 and edge2 have to be surface edges
    if (!(mesh.is_boundary(hedge1) && mesh.is_boundary(hedge2))){
        std::cout<<"one of the edges is not boundary, cannot merge"<<std::endl;
    }
    //they can't be neighboring edges either; then it would be a crack/wedge closure
    //may test it here
  
    std::vector<MyMesh::VertexHandle>  dummy_face_vhandles;
    MyMesh::FaceHandle dummy_face1, dummy_face2;
    MyMesh::HalfedgeHandle nonboundary_hedge1, nonboundary_hedge2, dummy_edge1, dummy_edge2;
    MyMesh::VertexHandle vA1, vA2, vB1, vB2, vremain1, vremain2;
    std::vector<MyMesh::VertexHandle> vh;
    //identify the boundary and non-boundary halfedges
    nonboundary_hedge1 = mesh.opposite_halfedge_handle(hedge1);
    nonboundary_hedge2 = mesh.opposite_halfedge_handle(hedge2);

    //need to copy edge properties
    mesh.copy_all_properties(nonboundary_hedge1, hedge2);
    mesh.copy_all_properties(nonboundary_hedge2, hedge1);
  
    vA1 = mesh.from_vertex_handle(hedge1);
    vA2 = mesh.to_vertex_handle(hedge1);  
    vB1 = mesh.from_vertex_handle(hedge2);
    vB2 = mesh.to_vertex_handle(hedge2);  
  
    vh.push_back(vB2);
    vh.push_back(vA1);
    vh.push_back(vB1);   
  
    dummy_face1 = mesh.add_face(vh); 
  
    vh.clear();
    vh.push_back(vA1);
    vh.push_back(vA2);
    vh.push_back(vB1);     
  
    dummy_face2 = mesh.add_face(vh); 
  
    dummy_edge1 = mesh.find_halfedge(vA1, vB2); //mesh.halfedge_handle( find_edge(mesh, vA1, vB2), 0 );
    if (!mesh.is_boundary(dummy_edge1)){
    	dummy_edge1 = mesh.opposite_halfedge_handle(dummy_edge1);
    }
    dummy_edge2 = mesh.find_halfedge(vA2, vB1); //mesh.halfedge_handle( find_edge(mesh, vA2, vB1), 0 );
    if (!mesh.is_boundary(dummy_edge2)){
    	dummy_edge2 = mesh.opposite_halfedge_handle(dummy_edge2);
    }   
    vremain1 = mesh.to_vertex_handle(dummy_edge1);
    vremain2 = mesh.to_vertex_handle(dummy_edge2);
    mesh.collapse(dummy_edge1);
    mesh.collapse(dummy_edge2);
    //mesh.garbage_collection();  
    return std::make_tuple(vremain1, vremain2);//find_edge(mesh, vremain1, vremain2);
}


/*	
Edge fission.
Edge is between v1 and v2.
Returns v1, v1new, v2, v2new
*/
std::tuple<MyMesh::VertexHandle, MyMesh::VertexHandle, MyMesh::VertexHandle, MyMesh::VertexHandle > split_edge(MyMesh & mesh, MyMesh::EdgeHandle & edge){

	MyMesh::HalfedgeHandle he = mesh.halfedge_handle(edge, 0);
	MyMesh::VertexHandle v1, v2, vA1, vA2, vB1, vB2, v1new, v2new, vAC, vBC;
    std::vector<MyMesh::VertexHandle> dummy_face_vhandles;
    MyMesh::FaceHandle dummy_face1, dummy_face2, dummy_face3, dummy_face4;  


	v2 = mesh.from_vertex_handle(he);
	v1 = mesh.to_vertex_handle(he);
    std::vector<MyMesh::VertexHandle> surface_neighborsA = find_surface_neighbors(mesh, v1);
    std::vector<MyMesh::VertexHandle> surface_neighborsB = find_surface_neighbors(mesh, v2);  
    vA1 = surface_neighborsA[0];
    vA2 = surface_neighborsA[1];    
    vB1 = surface_neighborsB[0];
    vB2 = surface_neighborsB[1];    

    vAC = mesh.add_vertex(0.5*(mesh.point(vA1) + mesh.point(vA2)));
    vBC = mesh.add_vertex(0.5*(mesh.point(vB1) + mesh.point(vB2)));    

    dummy_face_vhandles.clear(); 
    he = mesh.find_halfedge(vA1, v1);
    //need to get the orientation of dummy face right
    //this could be done probably by the order of surface neighbors returned by find_surface_neighbors()
    if (mesh.is_boundary(he)){   	 
	    dummy_face_vhandles.push_back( vA1);
	    dummy_face_vhandles.push_back( v1);
	    dummy_face_vhandles.push_back(vAC);
	    dummy_face1 = mesh.add_face(dummy_face_vhandles);  	 

		dummy_face_vhandles.clear(); 
	    dummy_face_vhandles.push_back( vAC);
	    dummy_face_vhandles.push_back( v1);
	    dummy_face_vhandles.push_back(vA2);
	    dummy_face3 = mesh.add_face(dummy_face_vhandles);  	     
    }
    else{
	    dummy_face_vhandles.push_back( v1);
	    dummy_face_vhandles.push_back( vA1);
	    dummy_face_vhandles.push_back(vAC);
	    dummy_face1 = mesh.add_face(dummy_face_vhandles);  	 

	    dummy_face_vhandles.clear(); 
	    dummy_face_vhandles.push_back( v1);
	    dummy_face_vhandles.push_back( vAC);
	    dummy_face_vhandles.push_back(vA2);
	    dummy_face3 = mesh.add_face(dummy_face_vhandles);  		         
    }

    dummy_face_vhandles.clear(); 
    he = mesh.find_halfedge(vB1, v2);
    //need to get the orientation of dummy face right
    //this could be done probably by the order of surface neighbors returned by find_surface_neighbors()
    if (mesh.is_boundary(he)){   	 
	    dummy_face_vhandles.push_back( vB1);
	    dummy_face_vhandles.push_back( v2);
	    dummy_face_vhandles.push_back(vBC);
	    dummy_face2 = mesh.add_face(dummy_face_vhandles); 

    	dummy_face_vhandles.clear(); 
	    dummy_face_vhandles.push_back( vBC);
	    dummy_face_vhandles.push_back( v2);
	    dummy_face_vhandles.push_back(vB2);
	    dummy_face4 = mesh.add_face(dummy_face_vhandles); 	     	  
    }
    else{
	    dummy_face_vhandles.push_back( v2);
	    dummy_face_vhandles.push_back( vB1);
	    dummy_face_vhandles.push_back(vBC);
	    dummy_face2 = mesh.add_face(dummy_face_vhandles);  	

	    dummy_face_vhandles.clear();  
	    dummy_face_vhandles.push_back( v2);
	    dummy_face_vhandles.push_back( vBC);
	    dummy_face_vhandles.push_back(vB2);
	    dummy_face4 = mesh.add_face(dummy_face_vhandles); 	         
    }    

    v1new = split_vertices(mesh, v1, vA1, v2);
    v2new = split_vertices(mesh, v2, vB2);


    mesh.delete_face(dummy_face1);
    mesh.delete_face(dummy_face2);
    mesh.delete_face(dummy_face3);
    mesh.delete_face(dummy_face4);  

//************************
//    mesh.delete_vertex(vAC);
//    mesh.delete_vertex(vBC);  
    //mesh.garbage_collection();
//*****************************    

    return std::make_tuple(v1, v1new, v2, v2new);
}

/*
Returns two parallel vectors with fusable halfedges.
Each pair appears TWICE (direct and swapped between the parallel vectors)
*/
std::tuple<std::vector<MyMesh::HalfedgeHandle>, std::vector<MyMesh::HalfedgeHandle> > get_halfedge_fusion_pairs(MyMesh & mesh){
	//find surface halfedges
	std::vector<MyMesh::HalfedgeHandle> surface_halfedges, halfedge_fusion_pairsA, halfedge_fusion_pairsB;
	MyMesh::VertexHandle v1A, v2A, v1B, v2B;
	std::vector<MyMesh::VertexHandle> common_neighbors1, common_neighbors2;
	for (MyMesh::HalfedgeIter h_it = mesh.halfedges_sbegin();h_it != mesh.halfedges_end(); ++h_it ){
		if (mesh.is_boundary(*h_it)){
			surface_halfedges.push_back(*h_it);
		}
	}

	//std::cout<<"number of surf halfedges "<<surface_halfedges.size()<<std::endl;
	for (auto & hedgeA: surface_halfedges) {
		v1A = mesh.from_vertex_handle(hedgeA);
		v2A = mesh.to_vertex_handle(hedgeA); 			
		for (auto & hedgeB: surface_halfedges) {   
			v1B = mesh.to_vertex_handle(hedgeB);
			v2B = mesh.from_vertex_handle(hedgeB); 	
		
			if ( (v1A==v1B) || (v2A==v2B) ){
				continue;
			}

			common_neighbors1 = find_common_neighbors(mesh, v1A, v1B);
			common_neighbors2 = find_common_neighbors(mesh, v2A, v2B);

			if ( (common_neighbors1.size()>0) || (common_neighbors2.size()>0)){
				continue;
			}

			else{
				halfedge_fusion_pairsA.push_back(hedgeA);
				halfedge_fusion_pairsB.push_back(hedgeB);			
			}
		}
	}

	return std::make_tuple(halfedge_fusion_pairsA, halfedge_fusion_pairsB);
}


std::tuple<std::vector<MyMesh::HalfedgeHandle>, std::vector<MyMesh::HalfedgeHandle> > get_halfedge_fusion_pairs(MyMesh & mesh, double l_fuse){
	double l_fuse2 = l_fuse*l_fuse;
	std::vector<MyMesh::HalfedgeHandle> h1, h2, h1f, h2f;
	MyMesh::VertexHandle v1, v2, v3, v4;	
    auto fusion_halfedge = get_halfedge_fusion_pairs(mesh);
    h1 = std::get<0>(fusion_halfedge);
    h2 = std::get<1>(fusion_halfedge);	

    for (unsigned int i=0; i<h1.size(); i++){
    	v1 = mesh.from_vertex_handle(h1[i]);
    	v2 = mesh.to_vertex_handle(h1[i]);
    	v3 = mesh.from_vertex_handle(h2[i]);
    	v4 = mesh.to_vertex_handle(h2[i]);
    	if ( ((mesh.point(v1) - mesh.point(v4)).sqrnorm() < l_fuse2) && ((mesh.point(v2) - mesh.point(v3)).sqrnorm() < l_fuse2) ){
    		h1f.push_back(h1[i]);
    		h2f.push_back(h2[i]);
    	}
    }    
    return std::make_tuple(h1f, h2f);
}


/*
Returns vector with edges eligible for edge fission
*/
std::vector<MyMesh::EdgeHandle> get_fission_edges(MyMesh & mesh){
	//find edges that are non-boundary, but both of its vertices are boundary
	//Those are candidates. They have to be filtered so that the structure doesn't fall apart 
	//upon fusion

	std::vector<MyMesh::EdgeHandle> fission_edges;
	MyMesh::VertexHandle v_from, v_to;
	MyMesh::HalfedgeHandle hedge, hedge1;
	MyMesh::EdgeHandle edge;

	//go over all halfedges
	for (MyMesh::EdgeIter e_it = mesh.edges_sbegin();e_it != mesh.edges_end(); ++e_it ){
		//find non-boundary halfedges
		if (!mesh.is_boundary(*e_it)){
			hedge1 = mesh.halfedge_handle(*e_it, 0);
			v_from = mesh.from_vertex_handle(hedge1);
			v_to = mesh.to_vertex_handle(hedge1);
			//edge = find_edge(mesh, v_to, v_from);
			//find boundary vertices of non-boundary halfedge
			if (mesh.is_boundary(v_from) && mesh.is_boundary(v_to) ){
				//check if edge fission would lead to disjoint structures
				//1. find a boundary halfedge of either v_from or v_to
				for (MyMesh::VertexOHalfedgeIter oh = mesh.voh_iter(v_from); oh.is_valid(); ++oh ){
					if (mesh.is_boundary(*oh)){
						hedge = *oh; 
					}
				}

				//fission_edges.push_back( *e_it );
				//2. start from the boundary halfedge and step through the boundary 
				for (int i=0; i<10000; i++){           //should be a while loop but maybe safer this way
					hedge = mesh.next_halfedge_handle(hedge);

					//if looped to the other vertex, structure will break; no need to continue iteration
					if (mesh.to_vertex_handle(hedge)==v_to){
						break;
					}

					//if looped back to where we started from, structure won't break
					if (mesh.to_vertex_handle(hedge)==v_from){
						fission_edges.push_back( *e_it );
						break;
					}

					
				}
				

			}
		}
	}

	return fission_edges;
}


//get all affected edges if vertex v is moved
//these are all incident edges to v and face-opposite edges of v (because the dihedrals on those will change too if v is moved)
std::vector<MyMesh::EdgeHandle> get_affected_edges(MyMesh & mesh, MyMesh::VertexHandle v){
	std::vector<MyMesh::EdgeHandle> affected_edges;
	//get all incident edges to v
  	for (MyMesh::VertexEdgeIter ve = mesh.ve_iter(v); ve.is_valid(); ++ve ){	
  		affected_edges.push_back(*ve);
  	}
  	//also get the opposite edges to v
  	//all faces of v
  	for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v); vf.is_valid(); ++vf ){
  		//all halfedges of all faces of v	
  		for (MyMesh::FaceHalfedgeIter fh = mesh.fh_iter(*vf); fh.is_valid(); ++fh){
  			if (mesh.opposite_vh(*fh) == v){
  				affected_edges.push_back( mesh.edge_handle(*fh));
  			}
  		}
  	}  	
  	return affected_edges;
}


//get all affected faces if vertex v is moved
//these are all faces adjacent to v; faces that move when v moves
std::vector<MyMesh::FaceHandle> get_affected_faces(MyMesh & mesh, MyMesh::VertexHandle v){
    std::vector<MyMesh::FaceHandle> affected_faces;
    for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v); vf.is_valid(); ++vf ){
        affected_faces.push_back(*vf);
    }   
    return affected_faces;
}


//get faces that can be removed while creating a hole
//these are triangles attached by all their 3 edges, basically all non-boundary faces
std::vector<MyMesh::FaceHandle> get_hole_removable_faces(MyMesh & mesh){
    std::vector<MyMesh::FaceHandle> removable_faces;
    for (MyMesh::FaceIter fit = mesh.faces_sbegin(); fit!=mesh.faces_end(); ++fit){
        if (!mesh.is_boundary(*fit, true)){ //true to check vertices too
            removable_faces.push_back(*fit);
        }
    }
    return removable_faces;

}


std::vector<MyMesh::FaceHandle> get_hole_removable_faces(MyMesh & mesh, int face_type){
    auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");    
    //find removable subunits topologically
    std::vector<MyMesh::FaceHandle> removable_faces_geom = get_hole_removable_faces(mesh);
    //filter them for the current subunit type only
    std::vector<MyMesh::FaceHandle> removable_faces;
    for (auto & rface : removable_faces_geom){
        if (face_props[rface].face_type == face_type){
            removable_faces.push_back(rface);
        }
    }
    return removable_faces;
}



/*  returns holes v1--v2--v3
    three vectors containing v1[i], v2[i], v3[i]
    RETURNS ALL HOLES 3 TIMES becase if v1-v2-v3 is a hole, then v2-v3-v1 and v3-v1-v2 are also counted but they are the same holes!!!!
    need to avoid 3 fold counting in rates!!!
*/

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> >  get_hole_triplets(MyMesh & mesh){
    std::vector <MyMesh::VertexHandle> vhole1, vhole2, vhole3;
    std::vector <MyMesh::HalfedgeHandle> surface_halfedges;

    MyMesh::VertexHandle v1, v2, v3, v4;
    MyMesh::HalfedgeHandle he;    
  
    //find surface vertices
    for (MyMesh::HalfedgeIter hit = mesh.halfedges_sbegin(); hit != mesh.halfedges_end(); ++hit){
      if (mesh.is_boundary(*hit)){
         surface_halfedges.push_back(*hit);
       }
    }    
    //loop over surface halfedges; check if can get back in 3 steps
    for (auto &he_start : surface_halfedges){


        he = he_start;
        v1 = mesh.to_vertex_handle(he);

        he = mesh.next_halfedge_handle(he);
        v2 = mesh.to_vertex_handle(he);        

        he = mesh.next_halfedge_handle(he);
        v3 = mesh.to_vertex_handle(he);

        he = mesh.next_halfedge_handle(he);                
          

        if (he==he_start){
            vhole1.push_back(v1);
            vhole2.push_back(v2);
            vhole3.push_back(v3);

        }


    }
      
    return std::make_tuple(vhole1, vhole2, vhole3);
}