#pragma once
//随机数发生器
#include<random>
using namespace std;
extern random_device rd;
extern mt19937 gen;
extern uniform_real_distribution<> dis;
//add head file
#ifndef CONFIG_H
#define CONFIG_H

#include <string>

// 声明（不是定义）
extern std::string outdir;

#endif

