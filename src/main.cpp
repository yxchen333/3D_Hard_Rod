#include"system.h"
#include<chrono>

using namespace std;



int main()
{
	//记录开始时间
	auto start = chrono::steady_clock::now();
	//-----------------------------------------------------------初始化----------------------------------------------------------------------------------------------
	
	double P_end = 0.6;
	int annealing_step = 200;//默认值2w
	//打开输入输出文件
	ifstream initial("initial.txt");//包含构造函数需要的一切信息，N，r，最后压强和压强退火间隔
	ofstream output("snap_1.dump");
	//构造函数初始化
	System test(initial);
	
	initial >> annealing_step;

	//initial >> P_end >> annealing_step;

	
	//测试直接输出
	test.output(output, 0);
	//-----------------------------------------------------------主程序------------------------------------------------------------------------------------------------
	////压缩蒙卡
	test.mc_npt(P_end, annealing_step);
	//test.mc_nvt( annealing_step);
	//test.output(output, 1);
	//记录结束时间
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed = end - start;
	cout << "Total time: " << elapsed.count() << " seconds." << endl;

	//统计用时
	//cout << "It has been running for " << double(clock()) / 60000 << " minites." << endl;
	return 0;
}

int main1()
{
	//记录开始时间
	auto start = chrono::steady_clock::now();
	ifstream initial_test("initial_test.txt");//包含构造函数需要的一切信息，N，r，最后压强和压强退火间隔
	
	System test(initial_test);
	ofstream output_test("test1.txt");
	test.output(output_test, 0);
	//先测试一下
	double P_end = 2;
	int annealing_step = 30;
	//initial >> P_end >> annealing_step;
	cout << "P_end: " << P_end << " annealing_step: " << annealing_step << endl;;
	

	//压缩蒙卡
	test.mc_npt(P_end, annealing_step);

	//记录结束时间
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed = end - start;
	cout << "Total time: " << elapsed.count() << " seconds." << endl;

	//统计用时
	//cout << "It has been running for " << double(clock()) / 60000 << " minites." << endl;
	return 0;


	/*RigidRod R1(64.6935, 31.0696, 1.51961, 20.0, 1.0);
	RigidRod R2(62.1946, 21.4901, 2.45069, 20.0, 1.0);
	if (CheckCollision(R2, R1)) cout << "重叠" << endl;
	else cout << "不重叠" << endl;
	return 0;*/
}