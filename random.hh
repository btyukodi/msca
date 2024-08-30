#include <random>

int randint(int max, std::mt19937 & engine);

double randdouble(double min, double max, std::mt19937 & engine);

double randdouble(std::mt19937 & engine);

double randnormal(std::mt19937 & engine);

double normal_pdf(double x);