#include <stdlib.h>
#include <math.h> 
#include <random>
#include <OpenMesh/Apps/Assembly/random.hh>
#include <OpenMesh/Apps/Assembly/geometry.hh>


MyMesh::Point randvec(std::mt19937 & eng)
{
	double psi1 = randdouble(eng);//rand()/(RAND_MAX + 1.0);
	double psi2 = randdouble(eng);//rand()/(RAND_MAX + 1.0);
	double theta = 2 * M_PI * psi2;
	double phi = acos(1 - 2 * psi1);
	MyMesh::Point v;	
	double u = randdouble(eng);//rand()/(RAND_MAX + 1.0);
	double up = pow(u, 1.0/3.0);
	v[0] = up*( sin(phi) * cos(theta) );
	v[1] = up*( sin(phi) * sin(theta) );
	v[2] = up* cos(phi);

	return v;
}


MyMesh::Point randvec_normal(std::mt19937 & eng)
{
	MyMesh::Point v;	
	v[0] = randnormal(eng);
	v[1] = randnormal(eng);
	v[2] = randnormal(eng);

	return v;
}


void rotatevec(MyMesh::Point & vec, MyMesh::Point axis, double angle)
{
	double vx, vy, vz, ex, ey, ez, nm; //,vxr,vyr,vzr;
	vx = vec[0];
	vy = vec[1];
	vz = vec[2];
	//finding normalized vector
	ex = axis[0];
	ey = axis[1];
	ez = axis[2];

	//cout << "hevec is" << he[heindex0].hevec[0] << " " << he[heindex0].hevec[1] << " " << he[heindex0].hevec[2] <<endl;
	nm = norm(axis);
	ex /= nm;
	ey /= nm;
	ez /= nm;
	double ct = cos(angle);
	double st = sin(angle);
	double mct = 1 - cos(angle);
	//rotating

	vec[0] = vx * (ct + ex * ex * mct) + vy * (ex * ey * mct - ez * st) + vz * (ex * ez * mct + ey * st) ;
	vec[1] = vx * (ey * ex * mct + ez * st) + vy * (ct + ey * ey * mct) + vz * (ey * ez * mct - ex * st) ;
	vec[2] = vx * (ez * ex * mct - ey * st) + vy * (ez * ey * mct + ex * st) + vz * (ct + ez * ez * mct) ;
	return;
}

/*
Given a triangle ABC and its edge lengths, calculate the in-plane coordinate of vertex C
assuming that A = (0,0) and B = (AB, 0) (Ox is directed towards AB)
https://math.stackexchange.com/questions/2156851/calculate-the-coordinates-of-the-third-vertex-of-triangle-given-the-other-two-an
*/
std::tuple<double, double> calculate_triangle_third_vertex_coordinate(double AB, double BC, double AC){
	double xC, yC;
	xC = (AB*AB + AC*AC - BC*BC) / (2.0 * AB);
	yC = sqrt( (AB+AC+BC)*(AB+AC-BC)*(AB-AC+BC)*(-AB+AC+BC) ) / (2.0 * AB);
	return std::make_tuple(xC, yC);
}



void zero_com_shift_mesh(MyMesh & mesh){
	MyMesh::Point com = MyMesh::Point(0,0,0);
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		com+=mesh.point(*v) / mesh.n_vertices();
	}
	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		mesh.set_point(*v, mesh.point(*v) - com );
	}	
	return;
}

//rotates a mesh randomly so that all orientations are equally probable
void random_rotate_mesh(MyMesh & mesh, std::mt19937 & eng){
	//pick a point from the mesh
	MyMesh::Point p_from = mesh.point( mesh.vertex_handle(0) );
	MyMesh::Point axis;
	MyMesh::Point p_to, p;	
	double angle;

	p_to = randvec(eng);
	p_to = p_to.normalize();
	p_to*=p_from.norm();

	axis = p_from.cross(p_to);
	axis = axis.normalize();
	
	angle = OpenMesh::angle( p_from.dot(p_to), (p_from.cross(p_to)).norm() );


	for (MyMesh::VertexIter v = mesh.vertices_sbegin(); v!=mesh.vertices_end(); ++v){
		p = mesh.point(*v);
		rotatevec(p, axis, angle);
		mesh.set_point(*v, p);
	}	

	return;
}