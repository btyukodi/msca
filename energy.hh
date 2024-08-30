
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Apps/Assembly/custom_mesh_props.hh>
#include <math.h> 
#ifndef ENERGYHH
#define ENERGYHH

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;





double halfedge_stretch_energy(MyMesh & mesh, MyMesh::HalfedgeHandle hedge, OpenMesh::HProp<HalfedgeProp> & edge_props);

//should this be edge or halfedge stretch energy? TBD, maybe should implement both?
//maybe just the full edge energy; that way edge length only has to be computed once
double edge_stretch_energy(MyMesh & mesh, MyMesh::EdgeHandle edge, OpenMesh::HProp<HalfedgeProp> & edge_props);


double edge_bending_energy(MyMesh & mesh, MyMesh::EdgeHandle edge, OpenMesh::HProp<HalfedgeProp> & edge_props);

double edge_binding_energy(MyMesh & mesh, MyMesh::EdgeHandle edge, OpenMesh::HProp<HalfedgeProp> & edge_props);

double full_elastic_energy(MyMesh & mesh);

double full_binding_energy(MyMesh & mesh);

double full_mu_N(MyMesh & mesh);

double full_energy(MyMesh & mesh);

double elastic_energy(MyMesh & mesh, std::vector<MyMesh::EdgeHandle> affected_edges);

double energy(MyMesh & mesh, std::vector<MyMesh::EdgeHandle> affected_edges);

double energy(MyMesh & mesh, MyMesh::EdgeHandle edge);

double full_stretch_energy(MyMesh & mesh);

double full_bending_energy(MyMesh & mesh);

double einstein_solid_energy(MyMesh & mesh, double k_einstein);

struct UmbrellaWindow{
	double spring_const;
	double N0;
	std::vector<int> values;
	UmbrellaWindow(double spring_const, double N0){
		this->spring_const = spring_const;
		this->N0 = N0;
	}
	double bias_potential(int N){
		return 0.5*this->spring_const*(N - this->N0)*(N - this->N0);
	}
	void add_value(int v){
		this->values.push_back(v);
	}
};

//this should (have been) a derived class from a more general, external potential class
struct Wall{
	double hardness;
	double z_position_1;
	double z_position_2;
	double amplitude;
	Wall(double hardness, double z_position_1, double z_position_2, double amplitude){
		this->hardness = hardness;
		this->z_position_1 = z_position_1;
		this->z_position_2 = z_position_2;
		this->amplitude = amplitude;
	}

	double potential_energy_1(double z){
		return this->amplitude * exp(this->hardness*(z - z_position_1));
	}

	double potential_energy_2(double z){
		return this->amplitude * exp(this->hardness*(z_position_2 - z)) ;
	}

	double potential_energy(double z){
		//return this->amplitude*( exp(this->hardness*(z - z_position_1) ) + exp(this->hardness*(z_position_2 - z)) );
		return potential_energy_1(z)+potential_energy_2(z);
	}

	//needs to be changed if the potential is changed!!!
	double force_1(double z){
		return -this->hardness * this->amplitude * exp(this->hardness*(z - this->z_position_1));
	}
	double force_2(double z){
		return this->hardness * this->amplitude * exp(this->hardness*(this->z_position_2 - z)) ;		
	}

	double wall_force_1(MyMesh & mesh, MyMesh::VertexHandle v){
		double z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1;
		return n_v*this->force_1(z);
	}
	double wall_force_2(MyMesh & mesh, MyMesh::VertexHandle v){
		double z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1;
		return n_v*this->force_2(z);
	}	
	double wall_force_1(MyMesh & mesh){
		double F=0.0;
	    for(MyMesh::VertexIter v = mesh.vertices_sbegin(); v != mesh.vertices_end(); ++v) {
	    	F+=wall_force_1(mesh, *v);
	    }
	    return F;
	}
	double wall_force_2(MyMesh & mesh){
		double F=0.0;
	    for(MyMesh::VertexIter v = mesh.vertices_sbegin(); v != mesh.vertices_end(); ++v) {
	    	F+=wall_force_2(mesh, *v);
	    }
	    return F;
	}			


	double wall_energy_1(MyMesh & mesh, MyMesh::VertexHandle v){
		double z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1;
		return n_v*this->potential_energy_1(z);
	}

	double wall_energy_2(MyMesh & mesh, MyMesh::VertexHandle v){
		double z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1;
		return n_v*this->potential_energy_2(z);
	}	

	double wall_energy(MyMesh & mesh, MyMesh::VertexHandle v){

		//should be identical to this:
		//return wall_energy_1(mesh, v)+wall_energy_2(mesh, v);

		double z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1;
		//multiply by "number of vertices in the hub" to make sure simple bond breaks don't change the wall energy
		//valence should be verified
		return n_v*this->potential_energy(z);
		/*MyMesh::Point com;
		double E=0.0;
		double z;
		for (auto & face: affected_faces){
			com = face_props[face].COM;
			z = com[2];
			E+=this->potential_energy(z);
		}
		return E;
		*/
	}

	//overloading for all vertices
	double wall_energy_1(MyMesh & mesh){
		double E=0.0;
	    for(MyMesh::VertexIter v = mesh.vertices_sbegin(); v != mesh.vertices_end(); ++v) {
	    	E+=wall_energy_1(mesh, *v);
	    }
	    return E;
	}

	//overloading for all vertices
	double wall_energy_2(MyMesh & mesh){
		double E=0.0;
	    for(MyMesh::VertexIter v = mesh.vertices_sbegin(); v != mesh.vertices_end(); ++v) {
	    	E+=wall_energy_2(mesh, *v);
	    }
	    return E;
	}



	//moves the upper wall by dz
	void move_1(double dz){
		this->z_position_1+=dz;
	}
	//moves the lower wall by dz
	void move_2(double dz){
		this->z_position_2+=dz;
	}	

};

//static AcceptanceMonitor default_AM = AcceptanceMonitor();
static MyMesh::FaceHandle default_facehandle;

//a generic external potential; all external potentials should be 
//inherited from this
struct ExternalPotential{
	virtual double energy(MyMesh & mesh, MyMesh::VertexHandle v, int valence_shift=0, MyMesh::FaceHandle face_to_remove=default_facehandle)=0;
	//full energy, overloading per-vertex energy
	//should be virtual and implemented in all subsequent potentials
	double energy(MyMesh & mesh){
		double E=0.0;
	    for(MyMesh::VertexIter v = mesh.vertices_sbegin(); v != mesh.vertices_end(); ++v) {
	    	E+=energy(mesh, *v);
	    }
	    return E;

	}
};

struct CorePotential : ExternalPotential{
	double radius;
	double hardness;

	CorePotential(double radius, double hardness){
		this->radius = radius;
		this->hardness = hardness;
	}

	//potential energy function
	double potential_energy(double x, double y, double z){
		//choose whatever potential with minimum at x^2+y^2+z^2=radius
		double r = sqrt(x*x+y*y+z*z);
		//return -this->hardness*exp(-(r-this->radius)*(r-this->radius));
		//std::cout<<"CorePotential "<<this->radius<<" "<<r<<std::endl;
		return this->hardness*(r-this->radius)*(r-this->radius);
	}

	//valence shift used in case the energy of removal is needed but the face is not yet removed
	double energy(MyMesh & mesh, MyMesh::VertexHandle v, int valence_shift=0, MyMesh::FaceHandle face_to_remove=default_facehandle){
		double x, y, z;
		x = mesh.point(v)[0];
		y = mesh.point(v)[1];		
		z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1+valence_shift;
		return n_v*this->potential_energy(x, y, z);
	}	



};


struct MorsePotential : ExternalPotential{
	double De;
	double a;
	double re;
	double f_restrain;
	double a_restrain;
	double L;

	using ExternalPotential::energy;

	MorsePotential(double De, double a, double re, double f_restrain, double a_restrain, double L=0){
		this->De = De;
		this->a = a;
		this->re =re;
		this->f_restrain=f_restrain;
		this->a_restrain=a_restrain;
		this->L=L;
	}

	//potential energy function
	double potential_energy(double x, double y, double z, double nx, double ny, double nz){
		//choose whatever potential with minimum at x^2+y^2+z^2=radius
		double r = sqrt(x*x+y*y+z*z);
		//return -this->hardness*exp(-(r-this->radius)*(r-this->radius));
		//std::cout<<"CorePotential "<<this->radius<<" "<<r<<std::endl;
		double sq = ( 1 - exp(-this->a*(r-this->re) ));
		//incorrect: return this->De*(1 - exp(-a*(r-this->re)*(r-this->re))  );

		double erx, ery, erz;
		erx=x/r;
		ery=y/r;
		erz=z/r;

		//std::cout<<"dotprod " <<(erx*nx + ery*ny + erz*nz)<<std::endl;
		//return this->De*(sq*sq - 1);

		//not OK, falls into the core return this->De*(sq*sq - 1) * (erx*nx + ery*ny + erz*nz);	
		return this->De*(sq*sq - 1) * (2.0 + erx*nx + ery*ny + erz*nz);	
	}

	//valence shift used in case the energy of removal is needed but the face is not yet removed
	/*double energy(MyMesh & mesh, MyMesh::VertexHandle v, int valence_shift=0){
		double x, y, z;
		x = mesh.point(v)[0];
		y = mesh.point(v)[1];		
		z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1+valence_shift;
		return n_v*this->potential_energy(x, y, z);
	}*/

	//can use the same Morse Potential as restraining, but make it deeper and possibly stiffer
	double restraining_potential(double x, double y, double z){
		//choose whatever potential with minimum at x^2+y^2+z^2=radius
		double r = sqrt(x*x+y*y+z*z);
		//return -this->hardness*exp(-(r-this->radius)*(r-this->radius));
		//std::cout<<"CorePotential "<<this->radius<<" "<<r<<std::endl;
		//double sq = ( 1 - exp(-this->a_restrain*(r-this->re) ));
		//incorrect: return this->De*(1 - exp(-a*(r-this->re)*(r-this->re))  );
		return this->De*exp(this->a_restrain* (r-this->f_restrain*this->re) );//this->D_restrain*(sq*sq - 1);

	}


	//double torque_energy(double x, double y, double z, double nx, double ny, double nz){
	//	return 0.0;
	//}

	//new version ensures that negative (i.e. non-removable) face types have an extra wall
	//to prevent them drifting away from the core
	double energy(MyMesh & mesh, MyMesh::VertexHandle v, int valence_shift=0, MyMesh::FaceHandle face_to_remove=default_facehandle){
		auto face_props = OpenMesh::FProp<FaceProp>(mesh, "face_props");	
		double x, y, z, nx, ny, nz;
		double E_core=0.0;
		double E_restrain=0.0;
		x = mesh.point(v)[0];
		y = mesh.point(v)[1];		
		z = mesh.point(v)[2];

		if (this->L > 0.0){
			//if on the cylinder portion
			if ((z<=0.5*this->L) && (z>=-0.5*this->L)){
				z=0.0; //work with the in-plane radius r only
			}
			//if above the top cap
			if (z>=0.5*this->L){
				z-=0.5*this->L; //center it to the top
			}
			if (z<=-0.5*this->L){
				z-=-0.5*this->L; //center it to the bottom
			}
		}

		int n_v = mesh.valence(v)-1+valence_shift;
		MyMesh::Point n;

		//E_core = n_v*this->potential_energy(x, y, z);


		//add an extra restraining potential to the non-removable subunit vertex
		//iterate over adjacent faces, test if any of them is non-removable; 
	    for (MyMesh::VertexFaceIter vf = mesh.vf_iter(v); vf.is_valid(); ++vf ){
			//not very nice, but it is assumed that whenever valence_shift<0, the removable handle is also passed to this function
			if ( (valence_shift<0) && (*vf==face_to_remove)){
				//do nothing if current face is the removable one
			}	
			else{
				n = mesh.calc_normal(*vf);
				//std::cout<<"n "<<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl;
				E_core+=this->potential_energy(x, y, z, n[0], n[1], n[2]);
				//mesh.update_normal(*vf);
				//
				//E_core+=this->torque_energy(x, y, z, n[0], n[1], n[2]);
			}    	


	    	if (face_props[*vf].face_type<0 ){
	    		//std::cout<<"restraining "<<std::endl;
	    		E_restrain+=restraining_potential(x, y, z) ;

	    	}
	    }   

		return E_core+E_restrain;
	}


};


/*struct SpheroCylinderPotential : MorsePotential{
	//length of spherocylinder
	double L;
	SpheroCylinderPotential( double De, double a, double re, double f_restrain, double a_restrain, double L)
	:MorsePotential(De, a, re, f_restrain, a_restrain),
	L(L)
	 {
	}


};*/


struct DummyPotential : ExternalPotential{
	DummyPotential(){}

	double energy(MyMesh & mesh, MyMesh::VertexHandle v, int valence_shift=0, MyMesh::FaceHandle face_to_remove=default_facehandle){
		return 0.0;
	}

};



struct WallPotential : ExternalPotential{
	double hardness;
	double z_position_1;
	double z_position_2;
	double amplitude;
	double move_rate;
	double dz_max;
	double weight;
	WallPotential(double hardness, double z_position_1, double z_position_2, double amplitude, double move_rate, double dz_max, double weight){
		this->hardness = hardness;
		this->z_position_1 = z_position_1;
		this->z_position_2 = z_position_2;
		this->amplitude = amplitude;
		this->move_rate = move_rate;
		this->dz_max = dz_max;
		this->weight = weight;
	}

	using ExternalPotential::energy;


	virtual double potential_energy_1(double z){
		return this->amplitude * exp(this->hardness*(z - z_position_1));
	}

	virtual double potential_energy_2(double z){
		return this->amplitude * exp(this->hardness*(z_position_2 - z)) ;
	}

	double potential_energy(double z){
		//return this->amplitude*( exp(this->hardness*(z - z_position_1) ) + exp(this->hardness*(z_position_2 - z)) );
		return potential_energy_1(z)+potential_energy_2(z);
	}

	//needs to be changed if the potential is changed!!!
	double force_1(double z){
		return -this->hardness * this->amplitude * exp(this->hardness*(z - this->z_position_1));
	}
	double force_2(double z){
		return this->hardness * this->amplitude * exp(this->hardness*(this->z_position_2 - z)) ;		
	}

	double wall_force_1(MyMesh & mesh, MyMesh::VertexHandle v){
		double z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1;
		return n_v*this->force_1(z);
	}
	double wall_force_2(MyMesh & mesh, MyMesh::VertexHandle v){
		double z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1;
		return n_v*this->force_2(z);
	}	
	double wall_force_1(MyMesh & mesh){
		double F=0.0;
	    for(MyMesh::VertexIter v = mesh.vertices_sbegin(); v != mesh.vertices_end(); ++v) {
	    	F+=wall_force_1(mesh, *v);
	    }
	    return F;
	}
	double wall_force_2(MyMesh & mesh){
		double F=0.0;
	    for(MyMesh::VertexIter v = mesh.vertices_sbegin(); v != mesh.vertices_end(); ++v) {
	    	F+=wall_force_2(mesh, *v);
	    }
	    return F;
	}			


	double wall_energy_1(MyMesh & mesh, MyMesh::VertexHandle v){
		double z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1;
		return n_v*this->potential_energy_1(z);
	}

	double wall_energy_2(MyMesh & mesh, MyMesh::VertexHandle v){
		double z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1;
		return n_v*this->potential_energy_2(z);
	}	

	double energy(MyMesh & mesh, MyMesh::VertexHandle v, int valence_shift=0, MyMesh::FaceHandle face_to_remove=default_facehandle){

		//should be identical to this:
		//return wall_energy_1(mesh, v)+wall_energy_2(mesh, v);

		double z = mesh.point(v)[2];
		int n_v = mesh.valence(v)-1+valence_shift;;
		//multiply by "number of vertices in the hub" to make sure simple bond breaks don't change the wall energy
		//valence should be verified
		return n_v*this->potential_energy(z);
		/*MyMesh::Point com;
		double E=0.0;
		double z;
		for (auto & face: affected_faces){
			com = face_props[face].COM;
			z = com[2];
			E+=this->potential_energy(z);
		}
		return E;
		*/
	}

	//full wall energy, overloading per-vertex energy
	/*double energy(MyMesh & mesh){
		double E=0.0;
	    for(MyMesh::VertexIter v = mesh.vertices_sbegin(); v != mesh.vertices_end(); ++v) {
	    	E+=energy(mesh, *v);
	    }
	    return E;

	}*/

	//overloading for all vertices
	double wall_energy_1(MyMesh & mesh){
		double E=0.0;
	    for(MyMesh::VertexIter v = mesh.vertices_sbegin(); v != mesh.vertices_end(); ++v) {
	    	E+=wall_energy_1(mesh, *v);
	    }
	    return E;
	}

	//overloading for all vertices
	double wall_energy_2(MyMesh & mesh){
		double E=0.0;
	    for(MyMesh::VertexIter v = mesh.vertices_sbegin(); v != mesh.vertices_end(); ++v) {
	    	E+=wall_energy_2(mesh, *v);
	    }
	    return E;
	}




	//moves the upper wall by dz
	void move_1(double dz){
		this->z_position_1+=dz;
	}
	//moves the lower wall by dz
	void move_2(double dz){
		this->z_position_2+=dz;
	}	

};


//single wall adhesion for testing. The potential energy then should be added to either (or both) walls in WallPotential
struct WallPotentialWithAdhesion : WallPotential{
	
	using WallPotential::WallPotential;

	//overriding the lower wall potential energy to make it Morse
	//might need a different potential with 3 parameters to keep the exponential close to the upper wall but tune the adhesion strength separately.
	//or just a Morse with independent parameters from the upper wall potential parameters. Then needs a constructor
	double potential_energy_2(double z){

		double a=0.5*this->hardness;
		double De=this->amplitude;
		double z0=this->z_position_2;
		double sq = ( 1 - exp(-a*( z - z0) ));

		return De*sq*sq;
	}


/*	double potential_energy_1(double z){
		return 0.0; //no upper wall for now
	}
*/



};

#endif