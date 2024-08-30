
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <random>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

/*
bool attempt_insertion(MyMesh & mesh, int face_type, int rotational_config, UmbrellaWindow & uw, std::mt19937 & eng);

bool attempt_wedge_insertion(MyMesh & mesh, int face_type, int rotational_config, UmbrellaWindow & uw, std::mt19937 & eng);

bool attempt_removal(MyMesh & mesh, int face_type, UmbrellaWindow & uw, std::mt19937 & eng);

bool attempt_wedge_removal(MyMesh & mesh, int face_type, UmbrellaWindow & uw, std::mt19937 & eng);
*/
bool attempt_insertion_external_potential(MyMesh & mesh, int face_type, int rotational_config, UmbrellaWindow & uw, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_wedge_insertion_external_potential(MyMesh & mesh, int face_type, int rotational_config, UmbrellaWindow & uw, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_removal_external_potential(MyMesh & mesh, int face_type, UmbrellaWindow & uw, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_wedge_removal_external_potential(MyMesh & mesh, int face_type, UmbrellaWindow & uw, std::mt19937 & eng, ExternalPotential & externalpotential);


bool attempt_move_external_potential(MyMesh & mesh, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_type1_fusion_external_potential(MyMesh & mesh, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_type2_fusion_external_potential(MyMesh & mesh, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_type2_fission_external_potential(MyMesh & mesh, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_type1_fission_external_potential(MyMesh & mesh, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_edge_fusion_external_potential(MyMesh & mesh, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_edge_fission_external_potential(MyMesh & mesh, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_hole_insertion_external_potential(MyMesh & mesh, int face_type, int rotational_config, UmbrellaWindow & uw, std::mt19937 & eng, ExternalPotential & externalpotential);

bool attempt_hole_removal_external_potential(MyMesh & mesh, int face_type, UmbrellaWindow & uw, std::mt19937 & eng, ExternalPotential & externalpotential);