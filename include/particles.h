
#pragma once
#include<iostream>
#include<fstream>
#include <sstream>
#include <math.h>
#include<cmath>
#include<string>
#include <sys/stat.h>//用来创建文件夹
#include"random.h"
using namespace std;

namespace constants//存放所有全局常量
{
	const double pi = 3.1415926535898;
}

// 三维向量结构（支持球柱体计算）
struct Vec3
{
	double x, y,z;  // 当前位置
	Vec3();//无参数初始化为全0
	Vec3(double dx, double dy, double dz);//初始化成三维向量
	Vec3 operator+(const Vec3& o) const { return { x + o.x, y + o.y ,z+o.z}; }
	Vec3 operator-(const Vec3& o) const { return { x - o.x, y - o.y,z-o.z }; }
	double dot(const Vec3& o) const { return x * o.x + y * o.y+z*o.z; } // 点积
	Vec3 operator*(const double k)const//数乘
	{
		Vec3 c;
		c.x = x * k;
		c.y = y * k;
		c.z = z * k;
		return c;
	};
	//叉乘
	Vec3 operator*(const Vec3& b)const
	{
		Vec3 c;
		c.x = y * b.z - z * b.y;
		c.y = z * b.x - x * b.z;
		c.z = x * b.y - y * b.x;
		return c;
	}
	Vec3 cross(const Vec3& o) const { return {y*o.z-z*o.y,z*o.x-x*o.z,x * o.y - y * o.x }; } // 叉积
	Vec3 rotate(const Vec3& axis, double theta) const { // 向量旋转
		Vec3 k = axis.normalized();

		double c = std::cos(theta);
		double s = std::sin(theta);

		return (*this) * c
			+ (k.cross(*this)) * s
			+ k * (k.dot(*this)) * (1.0 - c);
	};
	//投影到向量r对应的平面上
	Vec3 pro_plane(Vec3& r)const
	{
		Vec3 c;
		Vec3 n = r.normalized();
		c = *this - n * this->dot(n);
		return c;
	}
	// 在 Vec2 结构体内添加
	//Vec3 operator*(double scalar) const {
	//	return { x * scalar, y * scalar ,z*scalar};//位置伸缩
	//}
	double norm() const { return sqrt(x * x + y * y+z*z); }//模长
	Vec3 normalized() const { double n = norm(); return { x / n, y / n ,z/n}; }
};


class RigidRod
{
public:
	RigidRod();//无参数的默认构造函数
	RigidRod(double x_, double y_,double z_,double theta_,double phi_, double l_, double r_);//传入数据初始化
	void move_center(double x_, double y_,double z_);//
	void move_theta(double theta_,double phi_);//
	void modify_size(double length_, double r_);
	double* compute_quat();

	void GetEndpoints();//更新端点坐标
	// 新增：获取极角 theta
	double getTheta() const;

	// 可选：获取方位角 phi
	double getPhi() const;
	static double clamp(double x, double xmin, double xmax);
	//test collision submodule
	static double SegmentSegmentDistance(const RigidRod& a, const RigidRod& b);

	friend class System;
	friend bool CheckCollision(const RigidRod& a, const RigidRod& b);
	friend bool CheckCollision_unit(const RigidRod& a, const RigidRod& b);
	

private:
	Vec3 center;//质心坐标,xyz
	Vec3* end_point;//端点的两个坐标，也就是圆心
	//double theta;//直角边与y轴的夹角，从0到pi/2
	//double phi;
	Vec3 direction;
	double length;
	double radius;
	// 存放链表相关信息。第一个元素代表前一个粒子的编号，第二个元素代表后一个粒子的编号，-1代表前面或后面没有粒子
	// 第三第四个元素代表该链表所对应的cell分区编号，对应到System类中的cell_list数组，但一般只有表头粒子的数据才有意义。
};

bool CheckCollision(const RigidRod& a, const RigidRod& b);

bool CheckCollision_unit(const RigidRod& a, const RigidRod& b);

