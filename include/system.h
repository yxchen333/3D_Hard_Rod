
#include"particles.h"

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
	int N;//方块的数量
	double R_in;//系统的半径
	double R_out;
	double P = 0.1;//体系压强,npt相关

	//模拟相关
	double packing_fraction;
	double delta_x, delta_theta,delta_phi;//mc移动范围

	bool CheckBoundary(const RigidRod& rod);//边界重叠判定
	bool if_overlap_unit(int pointer);
	int if_overlap_unit_answer(int pointer);
	bool if_overlap_sys();
	bool attempt_move_x(int pointer, double delta_x);//尝试平移rod
	bool attempt_move_orientation(int pointer, double delta_theta, double delta_phi);
};



