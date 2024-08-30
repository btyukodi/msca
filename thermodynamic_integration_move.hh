#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <random>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;


bool attempt_move_thermodynamic_integration(MyMesh & mesh, std::mt19937 & eng, MyMesh::VertexHandle v_fix1, MyMesh::VertexHandle v_fix2, 
	MyMesh::VertexHandle v_fix3, MyMesh::Point e1, MyMesh::Point e2, double k_einstein, double lambda_ti);

