#include"system.h"

System::System()
{
	;
}

System::System(ifstream& initial)
{
	double r = 0;
	double pp1;
	double rl;//rod_length
	double rr = 0.5;//rod_radius
	double Routl;
	initial >> Routl >> pp1 >> rl;
	//three-dimension particle number
	//N = int((0.02* 4*Routl* Routl * Routl * rl * rl* rl* constants::pi/3)/(constants::pi/6+ rl* constants::pi/4));
	N = 14;
	double pf_init = 0.002;
	//cbrt :Cube root
	R_out = cbrt(N * (rl * rr* rr + 4*rr * rr* rr/3) / (pf_init) / (4 * (1 - pp1 * pp1)/3));
	R_in = R_out * pp1;
	particles = new RigidRod[N];
	for (int i = 0; i < N; i++)
	{
		double c_t1 = constants::pi * 2 * dis(gen);
		double c_t2 = constants::pi * dis(gen);
		double c_r = R_in + dis(gen) * (R_out - R_in);
		double theta_ = dis(gen) * constants::pi;
		// 随机取向：各向同性
		double phi = 2.0 * constants::pi * dis(gen);                // 方位角 [0, 2π)
		double theta = acos(2.0 * dis(gen) - 1.0);                    // 极角 [0, π]，保证各向同性
		particles[i] = RigidRod(c_r *sin(c_t2) *cos(c_t1), c_r * sin(c_t2)*sin(c_t1),c_r*cos(c_t2), theta_,phi, rl, rr);//center.x = x_;center.y = y_;	theta = theta_;length = l_;radius = r_;

	}
	//for (int i = 0; i < N; i++)
	//{
	//	bool ok = false;

	//	while (!ok)
	//	{
	//		double c_t1 = constants::pi * 2 * dis(gen);
	//		double c_t2 = constants::pi * dis(gen);
	//		double c_r = R_in + dis(gen) * (R_out - R_in);
	//		double theta_ = dis(gen) * constants::pi;
	//		// 随机取向：各向同性
	//		double phi = 2.0 * constants::pi * dis(gen);                // 方位角 [0, 2π)
	//		double theta = acos(2.0 * dis(gen) - 1.0); 
	//		particles[i] = RigidRod(c_r * sin(c_t2) * cos(c_t1), c_r * sin(c_t2) * sin(c_t1), c_r * cos(c_t2), theta_, phi, rl, rr);//center.x = x_;center.y = y_;	theta = theta_;length = l_;radius = r_;

	//		ok = true;

	//		for (int j = 0; j < i; j++)
	//		{
	//			if (CheckCollision(
	//				particles[i],
	//				particles[j]))
	//			{
	//				ok = false;
	//				break;
	//			}
	//		}
	//	}
	//}
	if (rl > 6)P = 0.003;
	else if(rl>3)P = 0.005;
	else P=0.008;
	//输出本次运行参数信息
	cout << "粒子数：" << N << endl;
	std::cout << "圆环内外径比: " << pp1 << endl;
	cout << "棒长：" << rl << endl;

	double packing_fraction_now = N * (rl * rr * rr + rr * rr * rr * 4/3) / (4*(R_out * R_out * R_out - R_in * R_in * R_in))/3 ;
	//cout << "初始化密度为：" << packing_fraction_now<<endl;

}

////produce module
//System::System(ifstream& initial_test)
//{
//	double r = 0;
//	N = 1;
//	double pp1 = 0;
//	R_out = 7.5;
//
//	double b = 4;
//	double radius = 0.5;
//	double length = b * radius * 2;//15*0.4*2=12;
//	std::cout << "长径比：L/D:" << b << endl;
//
//
//	R_in = R_out * pp1;
//	particles = new RigidRod[N];
//	for (int i = 0; i < N; i++)
//	{
//		initial_test >> particles[i].center.x >> particles[i].center.y >> particles[i].theta;
//
//		particles[i] = RigidRod(particles[i].center.x, particles[i].center.y, particles[i].theta, length, radius);//center.x = x_;center.y = y_;	theta = theta_;length = l_;radius = r_;
//
//	}
//}

//prewrite finish initial configuration
//System::System(ifstream& initial_test)
//{
//	double r = 0;
//	N = 1;
//	double pp1 = 0;
//	R_out = 68;
//
//	double b = 5;
//	double radius = 0.5;
//	double length = b * radius * 2;//15*0.4*2=12;
//	std::cout << "长径比：L/D:" << b << endl;
//
//
//	R_in = R_out * pp1;
//	particles = new RigidRod[N];
//	for (int i = 0; i < N; i++)
//	{
//		initial_test >> particles[i].center.x >> particles[i].center.y >> particles[i].theta;
//
//		particles[i] = RigidRod(particles[i].center.x, particles[i].center.y, particles[i].theta, length, radius);//center.x = x_;center.y = y_;	theta = theta_;length = l_;radius = r_;
//
//	}
//}
void System::output(ofstream& outfile_snap, int timestep)
{
	//时间步
	outfile_snap << "ITEM: TIMESTEP\n" << timestep << "\nITEM: NUMBER OF ATOMS\n" << N  << "\nITEM: BOX BOUNDS pp pp pp\n";
	//box尺寸
	outfile_snap << -R_out << " " << R_out << endl << -R_out << " " << R_out << endl << -R_out << " " << R_out << endl;
	//parameter
	outfile_snap << "ITEM: ATOMS id type x y z theta phi c_shape[1] c_shape[2] c_shape[3] c_orient[1] c_orient[2] c_orient[3] c_orient[4] lap\n";
	for (int i = 0; i < N; i++)
	{
		//double phi_snap = 1.5 * constants::pi - particles[i].theta;

		//粒子id
		outfile_snap << i + 1 << ' ';

		//粒子种类
		if (if_overlap_unit(i)) outfile_snap << 2;
		else outfile_snap << 1;

		//坐标
		outfile_snap << ' ' << particles[i].center.x << " " << particles[i].center.y << " " << particles[i].center.z << " " << particles[i].getTheta() << " " << particles[i].getPhi()<< " ";

		//形状
		outfile_snap << particles[i].radius << ' ' << particles[i].length << ' ' << particles[i].length << ' ';

		//四元数
		double* q = particles[i].compute_quat();
		outfile_snap << q[1] << ' ' << q[2] << ' ' << q[3] << " " << q[0] << " ";
		//outfile_snap<< 0 << ' ' << sin(0.5 * phi_snap) << ' ' << 0 << ' ' << cos(0.5 * phi_snap) << endl;

		//测试信息
		outfile_snap << if_overlap_unit_answer(i) << endl;
	}
	//输出中心小圆
	//outfile_snap << N + 1 << " " << 3 << " " << 0 << ' ' << 0 << " " << 0 << " " << 0 << " " << R_in << " " << R_in << " " << R_in << " " << 0 << " " << 0 << ' ' << 0 << " " << 0 << " " << 0 << endl;
	//这段代码的功能是向文件outfile_snap中写入一个中心位于原点、半径为R_in的圆形或球体的参数

}
void System::mc_nvt(int annealing_step)
{
	//mc相关参数
	int pointer = 0;//定义当前抽样对象
	int accept_x = 0;//位移接受数
	int accept_orien = 0;//旋转接受数
	int change_step = 100 * N;//调整delta的间隔次数,100个MC步
	double accept_min = 0.08 * change_step;//最小接受率
	double accept_max = 0.12 * change_step;//最大接受率

	//初始化，后续会调整
	delta_x = sqrt(N);
	delta_theta = constants::pi / 8;
	delta_phi = constants::pi / 12;
	for (int i = 0; i < annealing_step * N; i++)
	{
		pointer = floor(dis(gen) * N);
		//尝试移动(平动转动抽取一种)7.16尝试，抽取顺序
		if (dis(gen) < 0.5)//先平移
		{
			if (attempt_move_x(pointer, delta_x)) accept_x++;
			if (attempt_move_orientation(pointer, delta_theta, delta_phi)) accept_orien++;
		}
		else//先转动
		{
			if (attempt_move_orientation(pointer, delta_theta, delta_phi)) accept_orien++;
			if (attempt_move_x(pointer, delta_x)) accept_x++;
		}


		//调整delta范围
		if (i % change_step == 0)//隔change_step步根据接受率调整一次delta_x,并且输出一些信息来观测模拟
		{
			if (accept_x < accept_min)//过小/2
			{
				delta_x *= 0.95;
			}
			else if (accept_x > accept_max)//过大*2
			{
				delta_x *= 1.15;
				if (delta_x > R_out)delta_x = R_out;
			}
			//acc_xi = accept_x;
			//cout << "此时位置接受率是：" << 100 * accept_x / change_step  << "%，" << "delta_x :" << delta_x << endl;
			accept_x = 0;//重置接受数
		}
		if (i % change_step == 0)//隔一百步根据接受率调整一次delta_theta,并且输出一些信息来观测模拟
		{
			if (accept_orien < accept_min)//过小/2
			{
				delta_theta *= 0.95;
				delta_phi *= 0.95;
			}
			else if (accept_orien > accept_max)//过大*2
			{
				delta_theta *= 1.15;
				delta_phi *= 1.15;
				if (delta_theta > constants::pi) delta_theta = constants::pi;
				if (delta_phi > constants::pi) delta_phi = constants::pi;  // 与 theta 对称


			}
			//acc_thetai = accept_theta;
			//cout <<"此时角度接受率是：" << 100 * accept_theta / change_step <<"%"<< endl;
			accept_orien = 0;//重置接受数
		}

	}
	
}


void System::mc_npt(double P_end, int annealing_step)
{
	//初始化
	// 
	//npt相关参数
	
	double annealing_interval = 0.01;
	double V = (R_out * R_out - R_in * R_in) * constants::pi;//初始面积
	double beta = 1;//KT^-1,暂定为1
	int step_v = 0;//体积变化的次数
	int accept_v = 0;//体积变化接受数
	double delta_v = 0.001;//体积移动步长
	//
	//mc相关参数
	int pointer = 0;//定义当前抽样对象
	int accept_x = 0;//位移接受数
	int accept_orien = 0;//旋转接受数
	int change_step = 100 * N;//调整delta的间隔次数,100个MC步
	double accept_min = 0.08 * change_step;//最小接受率
	double accept_max = 0.12 * change_step;//最大接受率

	delta_x = sqrt(N);
	delta_theta = constants::pi / 8;
	delta_phi = constants::pi / 12;
	//其他参数
	long int current_squeez = 0;//记录上一次packing_fraction增大0.01时的步数，用来判断是否已经达到平衡，可以输出体系信息
	bool npt = true;//控制npt扩张的开关
	double pkf_point = 0.25;//控制packingfraction输出的点
	ofstream structure("structure.txt");
	ofstream p_observe("p.txt");
	//使用lammps计算自由能，需要这两个参量
	


	std::cout << "初始压强： " << P << " pkf_point:" << pkf_point << endl;
	cout << "最终压强为: " << P_end << " annealing_step: " << annealing_step << endl;
	//std::cout << "N: " << N << " R_out: " << R_out << " R_in： " << R_in << endl;
	


	//MC主循环
	int ci = 0;
	long int i = 0;//循环计步器
	long int P_i = 0;//退火计步器
	int timestep = 0;
	int p_step = 0;
	//while (P <= P_end)
	int cycle = 10;//尝试在压强为0.1下跑平衡，循环参数
	int cycle_i = 0;//统计

	double acc_xi = 0;//统计接受率
	double acc_thetai = 0;
	P = 0.09;//pressure test2026/6/13

	while (P <= P_end)
	{

		i++;
		//随机选点
		pointer = floor(dis(gen) * N);
		//尝试移动(平动转动抽取一种)7.16尝试，抽取顺序
		for (int i = 0; i < annealing_step * N; i++)
		{
			pointer = floor(dis(gen) * N);
			//尝试移动(平动转动抽取一种)7.16尝试，抽取顺序
			if (dis(gen) < 0.5)//先平移
			{
				if (attempt_move_x(pointer, delta_x)) accept_x++;
				if (attempt_move_orientation(pointer, delta_theta, delta_phi)) accept_orien++;
			}
			else//先转动
			{
				if (attempt_move_orientation(pointer, delta_theta, delta_phi)) accept_orien++;
				if (attempt_move_x(pointer, delta_x)) accept_x++;
			}


			//调整delta范围
			if (i % change_step == 0)//隔change_step步根据接受率调整一次delta_x,并且输出一些信息来观测模拟
			{
				if (accept_x < accept_min)//过小/2
				{
					delta_x *= 0.95;
				}
				else if (accept_x > accept_max)//过大*2
				{
					delta_x *= 1.15;
					if (delta_x > R_out)delta_x = R_out;
				}
				//acc_xi = accept_x;
				//cout << "此时位置接受率是：" << 100 * accept_x / change_step  << "%，" << "delta_x :" << delta_x << endl;
				accept_x = 0;//重置接受数
			}
			if (i % change_step == 0)//隔一百步根据接受率调整一次delta_theta,并且输出一些信息来观测模拟
			{
				if (accept_orien < accept_min)//过小/2
				{
					delta_theta *= 0.95;
					delta_phi *= 0.95;
				}
				else if (accept_orien > accept_max)//过大*2
				{
					delta_theta *= 1.15;
					delta_phi *= 1.15;
					if (delta_theta > constants::pi) delta_theta = constants::pi;
					if (delta_phi > constants::pi) delta_phi = constants::pi;  // 与 theta 对称


				}
				//acc_thetai = accept_theta;
				//cout <<"此时角度接受率是：" << 100 * accept_theta / change_step <<"%"<< endl;
				accept_orien = 0;//重置接受数
			}

		}


		if (npt == true)//npt进行的情况
		{
			P_i++;
			//——————————————压强退火——————————————————————————————————//
			if (P_i % (annealing_step * N) == 0)
			{
				double r = particles[pointer].radius;
				double packing_fraction_now = N * (particles[pointer].length * r * r + r * r *r* 4/3) / (4*(R_out * R_out * R_out - R_in * R_in * R_in)/3);
				cout << "进度是：" << cycle_i * 100 / cycle << "%" << endl;
				cout << "Current Packing Fraction is :" << packing_fraction_now << ", P = " << P << endl;
				//cout << "此时位置接受率是：" << 100 * acc_xi / (change_step) << "%" << "位置移动步长：" << delta_x << endl;
				//cout << "此时角度接受率是：" << 100 * acc_thetai / (change_step) << "%" << "角度移动步长：" << delta_theta << endl;
				acc_xi = 0;
				acc_thetai = 0;

				/*if (packing_fraction_now > 0.29)cycle = 100;
				if (packing_fraction_now > 0.29 && packing_fraction_now < 0.291)cout << "pf=0.29,此后压缩慢十倍" << endl;*/

				output(p_observe, p_step++);
				cycle_i++;
			}
			if (P_i % (annealing_step * N * cycle) == 0)
			{

				P += annealing_interval;
				cycle_i = 0;
			}


			//———————————npt蒙卡中，每一步都判断是否进行扩大，扩大后隔1万步平衡输出一次体系信息——————//
			if (dis(gen) < 1.0 / double(N + 1) && npt)
			{

				double save_l = particles[pointer].length;//保存扩大之前的边长，如果扩大失败就恢复
				double save_r = particles[pointer].radius;
				double delta_ln = delta_v * (dis(gen) - 0.5);//体积V在ln空间的位移

				step_v++;//体积变化步数+1

				//扩大边长,和体积的变化是相反的
				double new_l = save_l * exp(-delta_ln / 2);
				double new_r = save_r * exp(-delta_ln / 2);
				for (int k = 0; k < N; k++)
				{
					particles[k].modify_size(new_l, new_r);
				}

				//判断体积变化之后是否会导致重叠，注意只有在delta_ln小于零，即体积减小，粒子变大的时候才判断，粒子变小的时候百分百不重叠
				//定义一个判断是否重叠的信号,默认不重叠。
				bool flag = false;
				if (delta_ln < 0)//体积减小，判断一下会不会发生重叠
				{
					flag = if_overlap_sys();//判断整个体系是否存在重叠
				}
				if (flag)//如果重叠了，恢复原状态
				{
					for (int k = 0; k < N; k++)
					{
						particles[k].modify_size(save_l, save_r);
					}
				}
				else//没重叠，计算体积变化接受率并决定最终是否变化体积
				{
					double V_n = V * exp(delta_ln);//体积变化后的实际大小
					double acc = min(1.0, exp(-1 * beta * (P * (V_n - V) - (N + 1) * delta_ln / beta)));//计算接受概率
					//如果抽到接受，正式改变体积
					if (dis(gen) < acc)
					{
						accept_v++;
						double r = particles[pointer].radius;

						packing_fraction = N * (particles[pointer].length * r * r + r * r * r * 4 / 3) / (4 * (R_out * R_out * R_out - R_in * R_in * R_in) / 3);
						/*std::cout << "npt调整粒子" << ", P = " << P<<endl;*/
						//cout << "npt调整粒子，半径：" << particles[pointer].radius << "；棒长：" << particles[pointer].length << endl;
						 // 检查 packing_fraction 是否超过0.95
						if (packing_fraction > 0.95)
						{
							break; // 跳出while循环
						}
						if (packing_fraction >= pkf_point)//每1%的packing_fraction增加输出一次,每输出依次pkf增大0.01
						{
							current_squeez = i;
							pkf_point += 0.01;//增大1%再输出
							npt = false;//大于0.5的时候，停止压缩准备输出平衡信息
						}
					}
					//如果没抽到，就恢复原状
					else
					{
						for (int k = 0; k < N; k++)
						{
							particles[k].modify_size(save_l, save_r);
						}
					}
				}
				//动态设置体积变化步长
				if ((step_v) % 100 == 0)//隔一百步根据接受率调整一次delta_v,并且输出一些信息来观测模拟
				{
					if (accept_v < 20)//过小/2
					{
						delta_v /= 2;
					}
					else if (accept_v > 50)//过大*2
					{
						delta_v *= 2;
					}
					accept_v = 0;//重置接受数
				}
			}
		}
		else//npt关闭，需要输出平衡构象信息的时候
		{
			//——————————————平衡信息输出————————————//
			if ((i - current_squeez) % (annealing_step * N / 5) == 0)//平衡后输出，平衡时间为annealing step个MC步
			{
				//std::cout << "Current Packing Fraction is :" << packing_fraction << ", P = " << P << "，npt调整粒子后，nvt跑平衡" << endl;
				output(structure, timestep++);

				if (i - current_squeez == annealing_step * N)
				{
					npt = true;
					std::cout << "Current Packing Fraction is :" << packing_fraction << ", P = " << P << "，npt调整粒子后，nvt跑平衡" << endl;
				}
			}
			////关了npt后的nvt
			//output(structure, timestep++);
			//npt = true;

		}


	}
	//最后输出一次
	std::cout << "finally Packing Fraction is :" << packing_fraction << endl;
	//std::cout << "棒长" << particles[pointer].length << endl;
	std::cout << "外径棒长比：" << R_out / particles[pointer].length << endl;
	output(structure, timestep++);
}

bool System::CheckBoundary(const RigidRod& rod) {
	double effective_outer = R_out - rod.radius;
	double effective_inner = R_in + rod.radius;
	auto p1 = rod.end_point[0];
	auto p2 = rod.end_point[1];
	const double r_sq = effective_inner * effective_inner;
	const double R_sq = effective_outer * effective_outer;

	// overlap return true
	return (p1.x * p1.x + p1.y * p1.y+ p1.z * p1.z < r_sq) ||
		(p1.x * p1.x + p1.y * p1.y + p1.z * p1.z > R_sq) ||
		(p2.x * p2.x + p2.y * p2.y + p2.z * p2.z < r_sq) ||
		(p2.x * p2.x + p2.y * p2.y + p2.z * p2.z > R_sq) ||
		(rod.center.x * rod.center.x + rod.center.y * rod.center.y + rod.center.z * rod.center.z < r_sq) ||
		(rod.center.x * rod.center.x + rod.center.y * rod.center.y + rod.center.z * rod.center.z > R_sq);
	/*return 	(p1.x * p1.x + p1.y * p1.y > R_sq) || (p2.x * p2.x + p2.y * p2.y > R_sq) ||
		(rod.center.x * rod.center.x + rod.center.y * rod.center.y > R_sq);*/
}

bool System::if_overlap_unit(int pointer)
{
	//粒子之间重叠判断
	for (int j = 0; j < N; j++)
	{
		if (j == pointer) continue;
		if (CheckCollision(particles[j], particles[pointer])) return true;
	}
	//边界重叠判断
	if (CheckBoundary(particles[pointer])) return true;
	//无重叠返回false
	return false;
}

int System::if_overlap_unit_answer(int pointer)
{
	//粒子之间重叠判断
	for (int j = 0; j < N; j++)
	{
		if (j == pointer) continue;
		if (CheckCollision(particles[j], particles[pointer])) return j;
	}
	//边界重叠判断
	if (CheckBoundary(particles[pointer])) return 2;//重叠返回1
	//无重叠返回false
	return -1;
}

bool System::if_overlap_sys()
{
	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			if (CheckCollision(particles[i], particles[j])) return true;
		}
		//overlap return true
		if (CheckBoundary(particles[i])) return true;
	}
	//没重叠返回false
	return false;
}

bool System::attempt_move_x(int pointer, double delta_x)
{
	//保留数据准备重叠后恢复
	double save_x = particles[pointer].center.x;
	double save_y = particles[pointer].center.y;
	double save_z = particles[pointer].center.z;

	//先尝试平动，且计算移动后所在分区
	double theta_x = dis(gen) * constants::pi * 2;
	double phi_x = acos(2 * dis(gen) - 1);
	double move_r = dis(gen) * delta_x;
	particles[pointer].move_center(save_x + move_r * cos(theta_x)*sin(phi_x), save_y + move_r * sin(theta_x)*sin(phi_x),save_z+move_r*cos(phi_x));

	//判断是否有重叠
	bool flag_x = if_overlap_unit(pointer);

	if (flag_x == false)//说明没重叠
	{
		return true;
	}
	else//重叠了就恢复原状态
	{
		particles[pointer].move_center(save_x, save_y,save_z);
		//if (if_overlap_unit(pointer)) std::cout << "复原之后还是重叠" << endl;
		return false;
	}
}

bool System::attempt_move_orientation(int pointer, double delta_theta, double delta_phi)
{
	// 保存原始角度
	double old_theta = particles[pointer].getTheta();
	double old_phi = particles[pointer].getPhi();

	// 随机生成新角度（在旧角度附近均匀扰动）
	double new_theta = old_theta + (dis(gen) - 0.5) * delta_theta;
	double new_phi = old_phi + (dis(gen) - 0.5) * delta_phi;

	// 处理周期性边界：theta 范围 [0, π]，phi 范围 [-π, π] 或 [0, 2π)
	// 这里将 phi 映射到 [-π, π] 区间
	new_phi = fmod(new_phi + constants::pi, 2.0 * constants::pi);
	if (new_phi < 0) new_phi += 2.0 * constants::pi;
	new_phi -= constants::pi;  // 现在范围 [-π, π]

	// theta 限制在 [0, π]
	new_theta = fmod(new_theta, constants::pi);
	if (new_theta < 0) new_theta += constants::pi;

	// 应用新取向
	particles[pointer].move_theta(new_theta, new_phi);

	// 检查重叠
	bool overlap = if_overlap_unit(pointer);

	if (!overlap) {
		return true;   // 移动成功
	}
	else {
		// 恢复原取向
		particles[pointer].move_theta(old_theta, old_phi);
		return false;
	}
}