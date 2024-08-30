
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
//for now only for FaceProps, etc.
#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;


bool check_overlap(MyMesh & mesh, MyMesh::FaceHandle face1, MyMesh::FaceHandle face2, OpenMesh::FProp<FaceProp> & face_props );

bool check_full_overlap(MyMesh & mesh);

bool check_neighbor_overlap(MyMesh & mesh, MyMesh::FaceHandle face);

bool check_neighbor_overlap(MyMesh & mesh, MyMesh::FaceHandle face, OpenMesh::FProp<FaceProp> & face_props);

void update_neighbor_list(MyMesh & mesh, MyMesh::FaceHandle face);

void update_full_neighbor_list(MyMesh & mesh);


