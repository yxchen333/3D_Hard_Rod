#include"system.h"
#include <filesystem>


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
	initial >> Routl >> pp1 >> rl>>N;
	//three-dimension particle number
	//N = int((0.02* 4*Routl* Routl * Routl * rl * rl* rl* constants::pi/3)/(constants::pi/6+ rl* constants::pi/4));
	
	double pf_init = 0.002;
	//cbrt :Cube root
	R_out = cbrt(N * (rl * rr * rr + 4 * rr * rr * rr / 3) /
		pf_init / (4 * (1 - pp1 * pp1 * pp1) / 3));
	R_in = R_out * pp1;
	particles = new RigidRod[N];
	for (int i = 0; i < N; i++)
	{
		double c_t1 = constants::pi * 2 * dis(gen);
		double c_t2 = acos(2.0 * dis(gen) - 1.0);
		double c_r = cbrt(
			R_in * R_in * R_in
			+ dis(gen) * (R_out * R_out * R_out - R_in * R_in * R_in));
		// 
		double phi = 2.0 * constants::pi * dis(gen); 
		double theta = acos(2.0 * dis(gen) - 1.0); 
		particles[i] = RigidRod(c_r *sin(c_t2) *cos(c_t1), c_r * sin(c_t2)*sin(c_t1),c_r*cos(c_t2), theta,phi, rl, rr);//center.x = x_;center.y = y_;	theta = theta_;length = l_;radius = r_;

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
	// 
	// 
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
	//output fitter test
	/*HyperboloidFit fitter;

	fitter.TestSamplePoints(particles,N,R_out);*/

	
	cout << "粒子数：" << N << endl;
	std::cout << "内外半径比：" << pp1 << endl;
	cout << "棒长：" << rl << endl;

	double packing_fraction_now =
		N * (rl * rr * rr + 4.0 * rr * rr * rr / 3.0) /
		(4.0 * (R_out * R_out * R_out - R_in * R_in * R_in) / 3.0);
	packing_fraction = packing_fraction_now;
	L_R = rl / R_out / 2;
	cout << "初始化密度为：" << packing_fraction_now << endl;
	cout << "初始化L_R为：" << L_R << endl;



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
// 
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
// 
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
	// 
	outfile_snap << "ITEM: TIMESTEP\n" << timestep << "\nITEM: NUMBER OF ATOMS\n" << N  << "\nITEM: BOX BOUNDS pp pp pp\n";
	// 
	outfile_snap << -R_out << " " << R_out << endl << -R_out << " " << R_out << endl << -R_out << " " << R_out << endl;
	//parameter
	outfile_snap << "ITEM: ATOMS id type x y z theta phi c_shape[1] c_shape[2] c_shape[3] c_orient[1] c_orient[2] c_orient[3] c_orient[4] lap\n";
	for (int i = 0; i < N; i++)
	{
		//double phi_snap = 1.5 * constants::pi - particles[i].theta;

		// 
		outfile_snap << i + 1 << ' ';

		// 
		if (if_overlap_unit(i)) outfile_snap << 2;
		else outfile_snap << 1;

		// 
		outfile_snap << ' ' << particles[i].center.x << " " << particles[i].center.y << " " << particles[i].center.z << " " << particles[i].getTheta() << " " << particles[i].getPhi()<< " ";

		// 
		outfile_snap << particles[i].radius << ' ' << particles[i].length << ' ' << particles[i].length << ' ';

		// 
		double* q = particles[i].compute_quat();
		outfile_snap << q[1] << ' ' << q[2] << ' ' << q[3] << " " << q[0] << " ";
		//outfile_snap<< 0 << ' ' << sin(0.5 * phi_snap) << ' ' << 0 << ' ' << cos(0.5 * phi_snap) << endl;

		// 
		outfile_snap << if_overlap_unit_answer(i) << endl;
	}
	// 
	//outfile_snap << N + 1 << " " << 3 << " " << 0 << ' ' << 0 << " " << 0 << " " << 0 << " " << R_in << " " << R_in << " " << R_in << " " << 0 << " " << 0 << ' ' << 0 << " " << 0 << " " << 0 << endl;
	// 

}
void System::mc_nvt(int annealing_step)
{
	// 
	int pointer = 0; 
	int accept_x = 0; 
	int accept_orien = 0; 
	int change_step = 100 * N; 
	double accept_min = 0.08 * change_step; 
	double accept_max = 0.12 * change_step; 

	// 
	delta_x = sqrt(N);
	delta_theta = constants::pi / 8;
	delta_phi = constants::pi / 12;
	for (int i = 0; i < annealing_step * N; i++)
	{
		pointer = floor(dis(gen) * N);
		// 
		if (dis(gen) < 0.5) 
		{
			if (attempt_move_x(pointer, delta_x)) accept_x++;
			if (attempt_move_orientation(pointer, delta_theta, delta_phi)) accept_orien++;
		}
		else 
		{
			if (attempt_move_orientation(pointer, delta_theta, delta_phi)) accept_orien++;
			if (attempt_move_x(pointer, delta_x)) accept_x++;
		}


		// 
		if (i % change_step == 0) 
		{
			if (accept_x < accept_min) 
			{
				delta_x *= 0.95;
			}
			else if (accept_x > accept_max) 
			{
				delta_x *= 1.15;
				if (delta_x > R_out)delta_x = R_out;
			}
			//acc_xi = accept_x;
			// 
			accept_x = 0; 
		}
		if (i % change_step == 0) 
		{
			if (accept_orien < accept_min) 
			{
				delta_theta *= 0.95;
				delta_phi *= 0.95;
			}
			else if (accept_orien > accept_max) 
			{
				delta_theta *= 1.15;
				delta_phi *= 1.15;
				if (delta_theta > constants::pi) delta_theta = constants::pi;
				if (delta_phi > constants::pi) delta_phi = constants::pi; 


			}
			//acc_thetai = accept_theta;
			// 
			accept_orien = 0; 
		}

	}
	
}


void System::mc_npt(double P_end, int annealing_step)
{
	
	//initial pressure setting
	double P_annealing_interval = 0.02;
	P = 0.5;

	double V = (R_out * R_out - R_in * R_in) * constants::pi; 
	double beta = 1; 
	int step_v = 0; 
	int accept_v = 0; 
	double delta_v = 0.001; 
	//
	// 
	int pointer = 0; 
	int accept_x = 0; 
	int accept_orien = 0; 
	int change_step = 100 * N; 
	double accept_min = 0.08 * change_step;
	double accept_max = 0.12 * change_step; 

	delta_x = sqrt(N);
	delta_theta = constants::pi / 8;
	delta_phi = constants::pi / 12;
	// 
	long int current_squeez = 0; 
	bool npt = true; 
	double pkf_point = 0.25; 
	double L_R_point = 0.7;
	std::filesystem::create_directories(outdir);
	ofstream p_observe(outdir + "p.txt");

	ofstream rmse_out(
		outdir + "rmse.txt");
	rmse_out
		<< "timestep "
		<< "L_R "
		<< "a "
		<< "b "
		<< "c "
		<< "rmse\n";


	bool fit_started = false;
	int fit_count = 0;
	int nvt_sample_count = 0;
	//int fit_count_100 = 0;
	std::string vtk_dir = outdir + "vtk/";
	std::filesystem::create_directories(vtk_dir);
	char surface_name[40960];
	char rod_name[40960];
	
	

	
	std::cout << "初始压强：" << P << "，pkf_point：" << pkf_point << endl;
	cout << "最终压强：" << P_end << "，annealing_step：" << annealing_step << endl;
	
	


	// pressure control1
	int ci = 0;
	long int i = 0; 
	long int P_i = 0; 
	int timestep = 0;
	int p_step = 0;
	int cycle = 10; 
	int cycle_i = 0; 

	double acc_xi = 0; 
	double acc_thetai = 0;
	
	bool first_dump = true;

	while (P <= P_end)
	{

		i++;
		// 
		//pointer = floor(dis(gen) * N);
		// 
		for (int i = 0; i < annealing_step * N; i++)
		{
			pointer = floor(dis(gen) * N);
			// 
			if (dis(gen) < 0.5) 
			{
				if (attempt_move_x(pointer, delta_x)) accept_x++;
				if (attempt_move_orientation(pointer, delta_theta, delta_phi)) accept_orien++;
			}
			else 
			{
				if (attempt_move_orientation(pointer, delta_theta, delta_phi)) accept_orien++;
				if (attempt_move_x(pointer, delta_x)) accept_x++;
			}


			// 
			if (i % change_step == 0) 
			{
				if (accept_x < accept_min) 
				{
					delta_x *= 0.95;
				}
				else if (accept_x > accept_max) 
				{
					delta_x *= 1.15;
					if (delta_x > R_out)delta_x = R_out;
				}
				//acc_xi = accept_x;
				// 
				accept_x = 0; 
			}
			if (i % change_step == 0) 
			{
				if (accept_orien < accept_min) 
				{
					delta_theta *= 0.95;
					delta_phi *= 0.95;
				}
				else if (accept_orien > accept_max) 
				{
					delta_theta *= 1.15;
					delta_phi *= 1.15;
					if (delta_theta > constants::pi) delta_theta = constants::pi;
					if (delta_phi > constants::pi) delta_phi = constants::pi; 


				}
				//acc_thetai = accept_theta;
				// 
				accept_orien = 0; 
			}

		}


		if (npt == true) 
		{
			P_i++;
			
			if (P_i % (annealing_step * N) == 0)
			{
				double r = particles[pointer].radius;

				double packing_fraction_now =
					N *
					(particles[pointer].length * r * r
						+ 4.0 * r * r * r / 3.0)
					/
					(4.0 *
						(R_out * R_out * R_out
							- R_in * R_in * R_in)
						/ 3.0);
				double L_R_now = particles[pointer].length/R_out/2;

				std::cout<< "进度："<< cycle_i * 100 / cycle<< "%"<< std::endl;

				

				std::cout<< "Current Packing Fraction is :"<< packing_fraction_now<< ", P = "<< P<< std::endl;
				std::cout<< "Current L_R_now is :"<< L_R_now<< ", P = "<< P<< std::endl;

				cycle_i++;
				output(p_observe, p_step++);

				
			}
			if (P_i % (annealing_step * N * cycle) == 0)
			{

				P += P_annealing_interval;
				cycle_i = 0;
			}


			// 
			if (dis(gen) < 1.0 / double(N + 1) && npt)
			{

				double save_l = particles[pointer].length; 
				double save_r = particles[pointer].radius;
				double delta_ln = delta_v * (dis(gen) - 0.5); 

				step_v++; 

				// 
				double new_l = save_l * exp(-delta_ln / 2);
				double new_r = save_r * exp(-delta_ln / 2);
				for (int k = 0; k < N; k++)
				{
					particles[k].modify_size(new_l, new_r);
				}

				// 
				// 
				bool flag = false;
				if (delta_ln < 0) 
				{
					flag = if_overlap_sys(); 
				}
				if (flag) 
				{
					for (int k = 0; k < N; k++)
					{
						particles[k].modify_size(save_l, save_r);
					}
				}
				else 
				{
					double V_n = V * exp(delta_ln); 
					double acc = min(1.0, exp(-1 * beta * (P * (V_n - V) - (N + 1) * delta_ln / beta))); 
					// 
					if (dis(gen) < acc)
					{
						accept_v++;
						double r = particles[pointer].radius;


						//packing_fraction = N * (particles[pointer].length * r * r + r * r * r * 4 / 3) / (4 * (R_out * R_out * R_out - R_in * R_in * R_in) / 3);
						L_R = particles[pointer].length / R_out/2;
						
						
						if (L_R >= L_R_point) 
						{
							current_squeez = i;
							// 
							L_R_point += 0.01;
							npt = false; 
							nvt_sample_count = 0;
						}
					}
					// 
					else
					{
						for (int k = 0; k < N; k++)
						{
							particles[k].modify_size(save_l, save_r);
						}
					}
				}
				// 
				if ((step_v) % 100 == 0) 
				{
					if (accept_v < 20) 
					{
						delta_v /= 2;
					}
					else if (accept_v > 50) 
					{
						delta_v *= 2;
					}
					accept_v = 0; 
				}
			}
		}
		else 
		{
			if ((i - current_squeez) % (annealing_step * N ) == 0) {
				
					/*fitter.Fit(particles, N);

					double L_R_now = particles[pointer].length / R_out / 2;
					first_dump = false;*/

					//cout configuration frequence is 1/10

					/*if ((i - current_squeez) % (annealing_step * N*10 ) == 0) {
						sprintf_s(
							surface_name,
							sizeof(surface_name),
							"frame_surface_%06d.vtp",
							fit_count);

						sprintf_s(
							rod_name,
							sizeof(rod_name),
							"frame_rods_%06d.vtp",
							fit_count);
						fitter.OutputSurfaceVTP(
							vtk_dir + std::string(surface_name));
						fitter.OutputRodVTP(
							vtk_dir + std::string(rod_name));


						fitter.OutputPVD(
							vtk_dir + "surface.pvd",
							"frame_surface",
							fit_count);

						fitter.OutputPVD(
							vtk_dir + "rod.pvd",
							"frame_rods",
							fit_count);

					}
					*/

					

					//double current_rmse = fitter.GetRMSE();

					//rmse_out
					//	<< fit_count << " "
					//	<< L_R_now << " "
					//	<< fitter.GetA() << " "
					//	<< fitter.GetB() << " "
					//	<< fitter.GetC() << " "
					//	<< current_rmse
					//	<< std::endl;

					//if (current_rmse < best_fit.rmse)
					//{
					//	std::cout << "\nFind Better Hyperboloid Fit!\n";

					//	best_fit.rmse = current_rmse;

					//	best_fit.a = fitter.GetA();
					//	best_fit.b = fitter.GetB();
					//	best_fit.c = fitter.GetC();

					//	best_fit.frame = fit_count;

					//	//best_fit.packing_fraction = packing_fraction_now;

					//	// 
					//	best_fit.rods.assign(particles, particles + N);

					//	std::cout
					//		<< "Best Frame = "
					//		<< best_fit.frame
					//		<< "  RMSE = "
					//		<< best_fit.rmse
					//		<< std::endl;
					//}

					//fit_count++;
					nvt_sample_count++;




					if (nvt_sample_count == 10)
					{
						npt = true;
						nvt_sample_count = 0;
						//std::cout << "NVT10times完成，L_R = " << L_R << "，P = " << P << endl;
						double packing_fraction_reduced = N * particles[pointer].length * particles[pointer].length * 3 / 4 * constants::pi * R_out * R_out * R_out;
						std::cout << "now  reduced Packing Fraction is :" << packing_fraction_reduced << endl;
					}
				

			}
			
			
			
		



		}


	}
	
	double packing_fraction_reduced = N * particles[pointer].length * particles[pointer].length *3 / 4*constants::pi*R_out * R_out * R_out;
	std::cout << "finally Packing Fraction is :" << packing_fraction_reduced << endl;
	// 
	std::cout << "棒长与直径之比：" << particles[pointer].length / R_out /2<< endl;
	//----------------------------------------------------
// 
//----------------------------------------------------
	//if (!best_fit.rods.empty())
	//{
	//	std::cout << "\n====================================\n";
	//	std::cout << "Output Best Hyperboloid Fit...\n";
	//	std::cout << "Best Frame : " << best_fit.frame << std::endl;
	//	std::cout << "Best RMSE  : " << best_fit.rmse << std::endl;
	//	std::cout << "====================================\n";

	//	// 
	//	fitter.Fit(best_fit.rods.data(), N);

	//	fitter.OutputSurfaceVTP(outdir + "hyperboloid_Surface_best.vtp");
	//	fitter.OutputRodVTP(outdir + "hyperboloid_rod_best.vtp");

	//	
	//}
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
	// 
	for (int j = 0; j < N; j++)
	{
		if (j == pointer) continue;
		if (CheckCollision(particles[j], particles[pointer])) return true;
	}
	// 
	if (CheckBoundary(particles[pointer])) return true;
	// 
	return false;
}

int System::if_overlap_unit_answer(int pointer)
{
	// 
	for (int j = 0; j < N; j++)
	{
		if (j == pointer) continue;
		if (CheckCollision(particles[j], particles[pointer])) return j;
	}
	// 
	if (CheckBoundary(particles[pointer])) return 2; 
	// 
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
	// 
	return false;
}

bool System::attempt_move_x(int pointer, double delta_x)
{
	// 
	double save_x = particles[pointer].center.x;
	double save_y = particles[pointer].center.y;
	double save_z = particles[pointer].center.z;

	// 
	double theta_x = dis(gen) * constants::pi * 2;
	double phi_x = acos(2 * dis(gen) - 1);
	double move_r = dis(gen) * delta_x;
	particles[pointer].move_center(save_x + move_r * cos(theta_x)*sin(phi_x), save_y + move_r * sin(theta_x)*sin(phi_x),save_z+move_r*cos(phi_x));

	// 
	bool flag_x = if_overlap_unit(pointer);

	// 
	if (!flag_x) {
		return true;
	}
	else 
	{
		particles[pointer].move_center(save_x, save_y,save_z);
		// 
		return false;
	}
}

bool System::attempt_move_orientation(int pointer, double delta_theta, double delta_phi)
{
	// 
	double old_theta = particles[pointer].getTheta();
	double old_phi = particles[pointer].getPhi();

	// 
	double new_theta = old_theta + (dis(gen) - 0.5) * delta_theta;
	double new_phi = old_phi + (dis(gen) - 0.5) * delta_phi;

	// 
	// 
	new_phi = fmod(new_phi + constants::pi, 2.0 * constants::pi);
	if (new_phi < 0) new_phi += 2.0 * constants::pi;
	new_phi -= constants::pi; 

	// 
	new_theta = fmod(new_theta, constants::pi);
	if (new_theta < 0) new_theta += constants::pi;

	// 
	particles[pointer].move_theta(new_theta, new_phi);

	// 
	bool overlap = if_overlap_unit(pointer);

	if (!overlap) {
		return true; 
	}
	else {
		// 
		particles[pointer].move_theta(old_theta, old_phi);
		return false;
	}
}

bool System::PFStable()
{
	if (pf_history.size() < 20)
		return false;

	//----------------------------------
	// 
	//----------------------------------

	double mean = 0.0;

	for (double pf : pf_history)
		mean += pf;

	mean /= pf_history.size();

	//----------------------------------
	// 
	//----------------------------------

	double var = 0.0;

	for (double pf : pf_history)
	{
		double d = pf - mean;
		var += d * d;
	}

	var /= pf_history.size();

	double stddev = sqrt(var);

	std::cout
		<< "PF std = "
		<< stddev
		<< std::endl;

	//----------------------------------
	// 
	//----------------------------------

	return stddev < 1e-3;
}
