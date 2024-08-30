#include <OpenMesh/Apps/Assembly/topology.hh>
#include <OpenMesh/Apps/Assembly/monte_carlo.hh>
#include <stdlib.h>


int main_(){

	std::cout << "Running..." << std::endl;
	MyMesh mesh;
	// the request has to be called before a vertex/face/edge can be deleted. it grants access to the status attribute
	mesh.request_face_status();
	mesh.request_edge_status();
	mesh.request_vertex_status();  

	//attempt_move(mesh, 0);

	
}