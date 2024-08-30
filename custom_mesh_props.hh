//properties and functions referring to triangular subunits
//these customize the Mesh to be usable for self assembly
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#ifndef EXCLUDERSHH
#define EXCLUDERSHH

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

//data to be stored on halfedges
//this will go somewhere else too
struct HalfedgeProp
{
	//equilibrium length of edge
	float l0;
	//stretching modulus of edge
	float stretch_mod;
	//type of edge
	int edge_type;
	//binding energy with edge type i is e_b[i]
	std::vector<float> e_b;
	//preferred angle with edge type i is theta0[i]
	std::vector<float> theta0;
	//bending modulus with edge type i is bend_mod[i]
	std::vector<float> bend_mod;	
};



//subunit properties
struct FaceProp
{
	//type of the subunit
	int face_type;

	//subunit chemical potential
	float mu;

	//insertion rate
	float k_insertion;

	//rotational symmetry number of subunit
	int number_of_rotational_configs;

	//neighbors of current face
	std::vector<MyMesh::FaceHandle> neighbor_list;

	//center of mass coordinates of the face
	MyMesh::Point COM;	

	//for now, only center-of-mass excluder
	//more sophisticated data structure will be needed for multiple excluders 
	//and a list of all excluders here
	//also an overlap matrix - which type overlaps which

	//excluder radius
	double R_exc;


};


//mesh properties
struct MeshProp
{
	//mesh temperature
	float kT;
	//d_max
	float d_max;
	//maximum fusion distance; also, fusion sphere DIAMETER
	float l_fuse;

	//radius of sphere to add new vertex, sphere centered on the equilibrium position of the new vertex
	float R_add;

	//neighbor list box length
	float L_neighbor;

	//fusion/fission rate for type1 fusions/fissions
	float k_fusion;

	//fusion/fission rate for type2 fusions/fissions
	float k_fusion2;	

	//fusion/fission rate for edge fusion/fission
	float k_fusion_edge;

	//hold a prototype of each subunit type for property storage
	//std::vector<SubunitProp> subunit_prototypes;
	MyMesh prototypes_mesh;

	//hold a handle to each subunit type, so that prototype_faces[subunit_type] returns a sample subunit(face) of the given type
	//wouldn't need this if there was guarantee that prototypes_mesh.face_handle(ix) is ordered in the order of face additions
	std::vector<MyMesh::FaceHandle> prototype_faces;

	//holds a handle to one halfedge of each prototype face to keep track of orientation when copying from it
	std::vector<MyMesh::HalfedgeHandle> anchor_halfedges;

	//density of states to use upon placements; for instance, (DOS*v_fuse) is the number of states (for a vertex) within volume v_fuse
	float DOS;

	//rate for vertex moves
	float k_vertex_move;

};

//used for thermodynamic integration	
struct VertexProp
{

	//center of Einstein springs
	MyMesh::Point r0;

	int vertex_type;
};


void clone_properties_from_prototype(MyMesh & mesh, int proto_face_type, MyMesh::FaceHandle & face, MyMesh::HalfedgeHandle & anchor_halfedge);

void swap_configurations(MyMesh & mesh1, MyMesh & mesh2);

#endif