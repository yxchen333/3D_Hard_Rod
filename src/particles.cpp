#include"particles.h"

Vec3::Vec3()
{
    x = 0;
    y = 0;
    z = 0;
}


Vec3::Vec3(double dx, double dy, double dz)
{
    x = dx;
    y = dy;
    z = dz;
}


RigidRod::RigidRod()
{
    center = {0,0,0};
    direction = { 1,0,0 };
    radius = 0;
}

RigidRod::RigidRod(double x_, double y_, double z_,double theta_,double phi_, double l_, double r_)
{
   
    
    
    
    center = { x_, y_, z_ };

    

    // 球坐标 -> 单位方向向量
    direction.x = std::sin(theta_) * std::cos(phi_);

    direction.y = std::sin(theta_) * std::sin(phi_);

    direction.z = std::cos(theta_);

    length = l_;

    radius = r_;

    end_point = new Vec3[2];

    GetEndpoints();
}

void RigidRod::GetEndpoints()
{
    Vec3 dir = direction.normalized() * (length / 2);
    end_point[0] = center - dir;
    end_point[1] = center + dir;
}

// 在 .cpp 文件中实现
double RigidRod::getTheta() const {
    double z_clamped = std::max(-1.0, std::min(1.0, direction.z));
    return std::acos(z_clamped);
}

double RigidRod::getPhi() const {
    return std::atan2(direction.y, direction.x);
}
double RigidRod::clamp(double x, double xmin, double xmax)
{
    if(x<xmin)
    {
        return xmin;
    }
    else if (x > xmax)
    {
        return xmax;
    }
    else
    {
        return x;
    }

}
double* RigidRod::compute_quat()const
{
    double* q = new double[4];//四元数数组，分别代表w,x,y,z
    Vec3 zaxis=direction.normalized();//z轴
    Vec3 ref;
    if (fabs(zaxis.z) < 0.9)
        ref = { 0,0,1 };
    else
        ref = { 1,0,0 };
    Vec3 yaxis= ref.cross(zaxis).normalized();
    Vec3 xaxis = yaxis.cross(zaxis).normalized();
    double M[3][3] = { {xaxis.x,yaxis.x,zaxis.x},{ xaxis.y,yaxis.y,zaxis.y },{xaxis.z,yaxis.z,zaxis.z} };//旋转矩阵
    q[0] = 0.5 * sqrt(M[0][0] + M[1][1] + M[2][2] + 1);
    q[1] = 0.5 * sqrt(M[0][0] - M[1][1] - M[2][2] + 1);
    q[2] = 0.5 * sqrt(-1 * M[0][0] + M[1][1] - M[2][2] + 1);
    q[3] = 0.5 * sqrt(-1 * M[0][0] - M[1][1] + M[2][2] + 1);
    double q_max = 0;
    for (int i = 0; i < 4; i++)//寻找最大的q值用于下一步计算，防止数值不稳
    {
        if (q[i] > q_max) q_max = q[i];
    }
    if (q_max == q[0])//
    {
        q[1] = (M[2][1] - M[1][2]) / (4 * q[0]);
        q[2] = (M[0][2] - M[2][0]) / (4 * q[0]);
        q[3] = (M[1][0] - M[0][1]) / (4 * q[0]);
    }
    else if (q_max == q[1])
    {
        q[0] = (M[2][1] - M[1][2]) / (4 * q[1]);
        q[2] = (M[1][0] + M[0][1]) / (4 * q[1]);
        q[3] = (M[0][2] + M[2][0]) / (4 * q[1]);
    }
    else if (q_max == q[2])
    {
        q[0] = (M[0][2] - M[2][0]) / (4 * q[2]);
        q[1] = (M[0][1] + M[1][0]) / (4 * q[2]);
        q[3] = (M[1][2] + M[2][1]) / (4 * q[2]);
    }
    else
    {
        q[0] = (M[1][0] - M[0][1]) / (4 * q[3]);
        q[1] = (M[0][2] + M[2][0]) / (4 * q[3]);
        q[2] = (M[1][2] + M[2][1]) / (4 * q[3]);
    }
    return q;
}

void RigidRod::move_center(double x_, double y_,double z_)
{
    center = { x_,y_,z_ };
    GetEndpoints();
}

void RigidRod::move_theta(
    double theta_,
    double phi_)
{
    

    direction.x =
        std::sin(theta_) * std::cos(phi_);

    direction.y =
        std::sin(theta_) * std::sin(phi_);

    direction.z =
        std::cos(theta_);

    direction = direction.normalized();

    GetEndpoints();
}
void RigidRod::modify_size(double length_, double r_)
{
    radius = r_;
    length = length_;
    GetEndpoints();
}

// 辅助函数：点到线段的最短距离
double pointSegmentDistance(const Vec3& p, const Vec3& a, const Vec3& b) {
    Vec3 ab = b - a;
    Vec3 ap = p - a;
    double t = ap.dot(ab) / ab.dot(ab);
    t = std::max(0.0, std::min(1.0, t));
    Vec3 closest = a + ab * t;
    return (p - closest).norm();
}

//double RigidRod::SegmentSegmentDistance( const RigidRod& a,const RigidRod& b)
//{
//
//    Vec3 u = a.end_point[1] - a.end_point[0];
//    Vec3 v = b.end_point[1] - b.end_point[0];
//    Vec3 w = a.end_point[0] - b.end_point[0];
//
//    double A = u.dot(u);
//    double B = u.dot(v);
//    double C = v.dot(v);
//    double D = u.dot(w);
//    double E = v.dot(w);
//
//    double det = A * C - B * B;
//
//    double s, t;
//    const double eps = 1e-12;
//
//    //if parallel
//    if (fabs(det) < eps)//here need to reconsider,ok
//    {
//        double dist = 1e20;
//        for (int i = 0; i < 2; ++i) {
//            Vec3 p = (i == 0) ? a.end_point[0] : a.end_point[1];
//            dist = std::min(dist, pointSegmentDistance(p, b.end_point[0], b.end_point[1]));
//        }
//        for (int i = 0; i < 2; ++i) {
//            Vec3 p = (i == 0) ? b.end_point[0] : b.end_point[1];
//            dist = std::min(dist, pointSegmentDistance(p, a.end_point[0], a.end_point[1]));
//        }
//        return dist;
//
//
//    }
//    else
//    {
//        s = (B * E - C * D) / det;
//        t = (A * E - B * D) / det;
//        s = clamp(s, 0.0, 1.0);
//        t = clamp(t, 0.0, 1.0);
//
//        Vec3 p = a.end_point[0] + u * s;
//        Vec3 q = b.end_point[0] + v * t;
//        double dist = (p - q).norm();
//        return dist;
//    }
//   
//
//}

double RigidRod::SegmentSegmentDistance(
    const RigidRod& a,
    const RigidRod& b)
{
    Vec3 u = a.end_point[1] - a.end_point[0];
    Vec3 v = b.end_point[1] - b.end_point[0];
    Vec3 w = a.end_point[0] - b.end_point[0];

    double A = u.dot(u);
    double B = u.dot(v);
    double C = v.dot(v);
    double D = u.dot(w);
    double E = v.dot(w);

    double D0 = A * C - B * B;

    double sc, sN, sD = D0;
    double tc, tN, tD = D0;

    const double EPS = 1e-12;

    if (D0 < EPS)
    {
        sN = 0.0;
        sD = 1.0;
        tN = E;
        tD = C;
    }
    else
    {
        sN = (B * E - C * D);
        tN = (A * E - B * D);

        if (sN < 0.0)
        {
            sN = 0.0;
            tN = E;
            tD = C;
        }
        else if (sN > sD)
        {
            sN = sD;
            tN = E + B;
            tD = C;
        }
    }

    if (tN < 0.0)
    {
        tN = 0.0;

        if (-D < 0.0)
            sN = 0.0;
        else if (-D > A)
            sN = sD;
        else
        {
            sN = -D;
            sD = A;
        }
    }
    else if (tN > tD)
    {
        tN = tD;

        if (-D + B < 0.0)
            sN = 0;
        else if (-D + B > A)
            sN = sD;
        else
        {
            sN = -D + B;
            sD = A;
        }
    }

    sc = (fabs(sN) < EPS ? 0.0 : sN / sD);
    tc = (fabs(tN) < EPS ? 0.0 : tN / tD);

    Vec3 dP =
        w +
        u * sc -
        v * tc;

    return dP.norm();
}

//bool CheckCollision_unit(
//    const RigidRod& a,
//    const RigidRod& b)
//{
//    return
//        RigidRod::SegmentSegmentDistance(a, b) <(a.radius + b.radius);
//}
//
//bool CheckCollision(const RigidRod& a, const RigidRod& b)
//{
//    //overlap return true
//    if (CheckCollision_unit(a, b) || CheckCollision_unit(b, a)) return true;
//    else return false;
//}
bool CheckCollision(
    const RigidRod& a,
    const RigidRod& b)
{
    return
        RigidRod::SegmentSegmentDistance(a, b)
        <
        (a.radius + b.radius);
}
const Vec3& RigidRod::GetCenter() const
{
    return center;
}

const Vec3& RigidRod::GetDirection() const
{
    return direction;
}

double RigidRod::GetLength() const
{
    return length;
}
double RigidRod::GetRadius() const
{
    return radius;
}