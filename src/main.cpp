#include "system.h"
#include <chrono>

using namespace std;

int main()
{
	// 记录开始时间
	auto start = chrono::steady_clock::now();

	//-----------------------------------------------------------初始化----------------------------------------------------------------------------------------------
	double P_end = 2.5;
	int annealing_step = 200;// 默认值

	// 打开输入文件
	ifstream initial("input/initial.txt");// 包含构造函数需要的信息：内外径比、棒长和退火间隔

	// 构造函数初始化
	System test(initial);
	initial >> annealing_step;

	//-----------------------------------------------------------主程序------------------------------------------------------------------------------------------------
	// 压缩蒙卡
	test.mc_npt(P_end, annealing_step);

	// 记录结束时间
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed = end - start;
	cout << "Total time: " << elapsed.count() << " seconds." << endl;

	return 0;
}

// 保留的独立测试入口，不参与正常程序执行
int main1()
{
	// 记录开始时间
	auto start = chrono::steady_clock::now();
	ifstream initial_test("initial_test.txt");

	System test(initial_test);
	ofstream output_test("test1.txt");
	test.output(output_test, 0);

	// 测试参数
	double P_end = 2;
	int annealing_step = 30;
	cout << "P_end: " << P_end << " annealing_step: " << annealing_step << endl;

	// 压缩蒙卡
	test.mc_npt(P_end, annealing_step);

	// 记录结束时间
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed = end - start;
	cout << "Total time: " << elapsed.count() << " seconds." << endl;

	return 0;
}
