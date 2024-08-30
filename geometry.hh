
#include <iostream>
#include <stdlib.h>
#include <random>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Apps/Assembly/topology.hh>


/*
double norm(MyMesh::Point v);

double dot(std::vector<double> v1, std::vector<double> v2);

std::vector<double> cross(std::vector<double> v1, std::vector<double> v2);
*/
MyMesh::Point randvec(std::mt19937 & eng);

MyMesh::Point randvec_normal(std::mt19937 & eng);

void rotatevec(MyMesh::Point & vec, MyMesh::Point axis, double angle);

std::tuple<double, double> calculate_triangle_third_vertex_coordinate(double AB, double BC, double AC);

void zero_com_shift_mesh(MyMesh & mesh);

void random_rotate_mesh(MyMesh & mesh, std::mt19937 & eng);

//void normalize(std::vector<double> & v);