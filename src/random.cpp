#include"random.h"

random_device rd;  // Will be used to obtain a seed for the random number engine
mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> dis(0.0, 1.0);



// output dir 
std::string outdir = "output/article/1/3.03_15_3272";