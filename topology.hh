#include <iostream>
#include <stdlib.h>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

//#ifndef TOPOLOGYHH
//#define TOPOLOGYHH

/*
Everything is passed by reference. It would probably work even if handles weren't passed by reference.
Handles seem to already do that, they refer to the same elements even if passed by copy.
Example: set two handles to same face, delete one. Both will show deleted:

    MyMesh::FaceHandle fh1 = mesh.face_handle(0);	
    MyMesh::FaceHandle fh2 = mesh.face_handle(0);	
    std::cout<<"fh1=fh2? "<<(fh1==fh2)<<std::endl;
    mesh.delete_face(fh1);
    std::cout<<"deleted? fh1, fh2 "<<mesh.status(fh1).deleted()<<" "<<mesh.status(fh2).deleted()<<std::endl;

Probably is nonsense to pass by reference, correct these once everything else works.
*/

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

MyMesh::EdgeHandle find_edge(MyMesh & mesh, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v2 );

MyMesh::VertexHandle merge_vertices(MyMesh & mesh, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v2);

std::vector<MyMesh::VertexHandle> find_common_neighbors(MyMesh & mesh, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v2);

std::vector<MyMesh::VertexHandle> find_surface_neighbors(MyMesh & mesh, MyMesh::VertexHandle & v1);

std::vector<MyMesh::HalfedgeHandle> find_boundary_halfedges(MyMesh & mesh);

MyMesh::HalfedgeHandle find_outgoing_boundary_he(MyMesh & mesh, MyMesh::VertexHandle & v);

MyMesh::HalfedgeHandle find_incoming_boundary_he(MyMesh & mesh, MyMesh::VertexHandle & v);

MyMesh::VertexHandle split_vertices(MyMesh & mesh, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v2, MyMesh::VertexHandle & v3);

MyMesh::VertexHandle split_vertices(MyMesh & mesh, MyMesh::VertexHandle & v1, MyMesh::VertexHandle & v2);

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> , std::vector<int>>  get_open_wedge_triplets(MyMesh & mesh);

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > get_type1_fission_pairs(MyMesh & mesh);

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > get_type2_fission_triplets(MyMesh & mesh);

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > get_type1_fusion_triplets(MyMesh & mesh);

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > get_type1_fusion_triplets(MyMesh & mesh, double l_fuse);

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > get_type2_fusion_triplets(MyMesh & mesh);

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> > get_type2_fusion_triplets(MyMesh & mesh, double l_fuse);

std::vector<MyMesh::FaceHandle> get_simply_removable_faces(MyMesh & mesh);

std::vector<MyMesh::FaceHandle> get_simply_removable_faces(MyMesh & mesh, int face_type);

std::vector<MyMesh::FaceHandle> get_wedge_removable_faces(MyMesh & mesh);

std::vector<MyMesh::FaceHandle> get_wedge_removable_faces(MyMesh & mesh, int face_type);

std::tuple< MyMesh::VertexHandle, MyMesh::VertexHandle > merge_halfedges(MyMesh & mesh, MyMesh::HalfedgeHandle & hedge1, MyMesh::HalfedgeHandle & hedge2);

std::tuple<MyMesh::VertexHandle, MyMesh::VertexHandle, MyMesh::VertexHandle, MyMesh::VertexHandle > split_edge(MyMesh & mesh, MyMesh::EdgeHandle & edge);

std::tuple<std::vector<MyMesh::HalfedgeHandle>, std::vector<MyMesh::HalfedgeHandle> > get_halfedge_fusion_pairs(MyMesh & mesh);

std::tuple<std::vector<MyMesh::HalfedgeHandle>, std::vector<MyMesh::HalfedgeHandle> > get_halfedge_fusion_pairs(MyMesh & mesh, double l_fuse);

std::vector<MyMesh::EdgeHandle> get_fission_edges(MyMesh & mesh);

std::vector<MyMesh::EdgeHandle> get_affected_edges(MyMesh & mesh, MyMesh::VertexHandle v);

std::vector<MyMesh::FaceHandle> get_affected_faces(MyMesh & mesh, MyMesh::VertexHandle v);

std::tuple< std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle>, std::vector<MyMesh::VertexHandle> >  get_hole_triplets(MyMesh & mesh);

std::vector<MyMesh::FaceHandle> get_hole_removable_faces(MyMesh & mesh);

std::vector<MyMesh::FaceHandle> get_hole_removable_faces(MyMesh & mesh, int face_type);

//#endif