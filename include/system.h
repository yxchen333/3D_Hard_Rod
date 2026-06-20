#pragma once
#include "hyperboloid.h"
#include "particles.h"
#include <vector>

struct BestFitData
{
	double rmse = 1e100;
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double packing_fraction = 0.0;
	int frame = -1;
	std::vector<RigidRod> rods;
};

class System
{
public:
	System();
	System(ifstream& initial);
	void mc_npt(double P_end, int annealing_step);
	void output(ofstream& structure, int timestep);
	void mc_nvt(int annealing_step);

private:
	RigidRod* particles;
	int N;// 棒的数量
	double R_in;// 系统内半径
	double R_out;// 系统外半径
	double L_R = 0.0;// 棒长与容器外直径之比：L/(2*R_out)
	double P = 0.0;// 体系压强，供NPT使用
	vector<double> pf_history;
	bool PFStable();

	BestFitData best_fit;
	HyperboloidFit fitter;

	// 模拟相关参数
	double packing_fraction;
	double delta_x, delta_theta, delta_phi;// 蒙卡移动范围

	bool CheckBoundary(const RigidRod& rod);// 边界重叠判定
	bool if_overlap_unit(int pointer);
	int if_overlap_unit_answer(int pointer);
	bool if_overlap_sys();
	bool attempt_move_x(int pointer, double delta_x);// 尝试平移棒
	bool attempt_move_orientation(int pointer, double delta_theta, double delta_phi);// 尝试转动棒
};
