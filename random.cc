#include <random>
#include <OpenMesh/Apps/Assembly/random.hh>

int randint(int max, std::mt19937 & engine){
	std::uniform_int_distribution<int> distribution(0,max-1);
	return distribution(engine);
}

double randdouble(double min, double max, std::mt19937 & engine){
	std::uniform_real_distribution<double> distribution(min,max);
	return distribution(engine);
}

double randdouble(std::mt19937 & engine){
	static std::uniform_real_distribution<double> distribution(0.0, 1.0);
	return distribution(engine);
}

double randnormal(std::mt19937 & engine){
	static std::normal_distribution<double> distribution(0.0, 1.0);
	return distribution(engine);
}

double normal_pdf(double x){
	const double C = 1.0/sqrt(2.0*3.1415926535);
	return C*exp(-0.5*x*x);
}

