#pragma once
//随机数发生器
#include<random>
using namespace std;
extern random_device rd;
extern mt19937 gen;
extern uniform_real_distribution<> dis;

