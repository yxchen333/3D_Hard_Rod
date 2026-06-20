#include "hyperboloid.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

HyperboloidFit::HyperboloidFit()
{
    a = 1.0;
    b = 1.0;
    c = 1.0;

    rmse = 0.0;
}

void HyperboloidFit::GenerateSamplePoints(
    const RigidRod* rods,
    int N)
{
    samples.clear();
    rods_.assign(rods, rods + N);

    const int SAMPLE_NUM = 21;

    for (size_t i = 0; i < N; i++)
    {
        const RigidRod& rod = rods[i];

        //--------------------------------------------------
        // 单位方向
        //--------------------------------------------------

        Vec3 dir = rod.GetDirection().normalized();

        //--------------------------------------------------
        // 从一端采样到另一端
        //--------------------------------------------------

        for (int j = 0; j < SAMPLE_NUM; j++)
        {
            // t : [0,1]
            double t = double(j) / (SAMPLE_NUM - 1);

            // s : [-0.5 , 0.5]
            double s = t - 0.5;

            SamplePoint sp;

            sp.rod_id = (int)i;

            sp.s = s;

            sp.p =
                rod.GetCenter()
                +
                dir * (s * rod.GetLength());

            samples.push_back(sp);
        }
    }

   /* std::cout
        << "Generate "
        << samples.size()
        << " sample points."
        << std::endl;*/

}

//inspect function
void HyperboloidFit::PrintSampleInfo()
{
    std::cout
        << "Sample Number = "
        << samples.size()
        << std::endl;

    for (size_t i = 0; i < std::min(size_t(10), samples.size()); i++)
    {
        std::cout

            << samples[i].rod_id

            << " "

            << samples[i].p.x

            << " "

            << samples[i].p.y

            << " "

            << samples[i].p.z

            << std::endl;
    }
}

//output function
void HyperboloidFit::OutputSamplePoints(const std::string& filename)
{
    std::ofstream out(filename);

    if (!out)
    {
        std::cout << "Cannot open " << filename << std::endl;
        return;
    }

    //-------------------------------------------------
    // 自动计算包围盒（OVITO显示更舒服）
    //-------------------------------------------------

    double xmin = 1e30, xmax = -1e30;
    double ymin = 1e30, ymax = -1e30;
    double zmin = 1e30, zmax = -1e30;

    for (const auto& sp : samples)
    {
        xmin = std::min(xmin, sp.p.x);
        xmax = std::max(xmax, sp.p.x);

        ymin = std::min(ymin, sp.p.y);
        ymax = std::max(ymax, sp.p.y);

        zmin = std::min(zmin, sp.p.z);
        zmax = std::max(zmax, sp.p.z);
    }

    // 给边界留一点余量
    double margin = 2.0;

    xmin -= margin;
    xmax += margin;

    ymin -= margin;
    ymax += margin;

    zmin -= margin;
    zmax += margin;

    //-------------------------------------------------
    // 标准LAMMPS dump格式
    //-------------------------------------------------

    out << "ITEM: TIMESTEP\n";
    out << "0\n";

    out << "ITEM: NUMBER OF ATOMS\n";
    out << samples.size() << "\n";

    out << "ITEM: BOX BOUNDS pp pp pp\n";
    out << xmin << " " << xmax << "\n";
    out << ymin << " " << ymax << "\n";
    out << zmin << " " << zmax << "\n";

    out << "ITEM: ATOMS id type x y z\n";

    for (size_t i = 0; i < samples.size(); i++)
    {
        out
            << i + 1 << " "
            << samples[i].rod_id + 1 << " "
            << samples[i].p.x << " "
            << samples[i].p.y << " "
            << samples[i].p.z << "\n";
    }

    out.close();

    std::cout << "Sample points saved to " << filename << std::endl;
}

void HyperboloidFit::TestSamplePoints(
    const RigidRod* rods,
    int N,
    double R_out)
{
    GenerateSamplePoints(rods, N);
    PrintSampleInfo();
    OutputDebugDump("debug.dump", rods, N, R_out, 0, true, true, true);
    // Optional visualization output can be written by the caller.
    // 测试双曲面方程
    //TestSurfaceValue();
    //test distance()
    a = 8.0;
    b = 8.0;
    c = 15.0;

    std::cout << "Distance = "
        << Distance(samples[0].p)
        << std::endl;


}

void HyperboloidFit::OutputDebugDump(
    const std::string& filename,
    const RigidRod* rods,
    int N,
    double R_out,
    int timestep,
    bool append,
    bool output_sample,
    bool output_surface)
{
    //-------------------------------------------------------
    // 打开文件
    //-------------------------------------------------------

    std::ofstream out;

    if (append)
        out.open(filename, std::ios::app);
    else
        out.open(filename);

    if (!out)
    {
        std::cout << "Cannot open " << filename << std::endl;
        return;
    }

    //-------------------------------------------------------
    // 计算原子数
    //-------------------------------------------------------

    int natoms = N;

    if (output_sample)
        natoms += static_cast<int>(samples.size());

    if (output_surface)
        natoms += static_cast<int>(surface_points.size());

    //-------------------------------------------------------
    // Header
    //-------------------------------------------------------

    out << "ITEM: TIMESTEP\n";
    out << timestep << "\n";

    out << "ITEM: NUMBER OF ATOMS\n";
    out << natoms << "\n";

    out << "ITEM: BOX BOUNDS pp pp pp\n";

    out << -R_out << " " << R_out << "\n";
    out << -R_out << " " << R_out << "\n";
    out << -R_out << " " << R_out << "\n";

    out << "ITEM: ATOMS ";
    out << "id type ";
    out << "x y z ";
    out << "theta phi ";
    out << "c_shape[1] c_shape[2] c_shape[3] ";
    out << "c_orient[1] c_orient[2] c_orient[3] c_orient[4] ";
    out << "lap\n";

    //-------------------------------------------------------
    // 第一部分：Rod
    //-------------------------------------------------------

    int id = 1;

    for (int i = 0; i < N; i++)
    {
        out << id++ << " ";

        out << 1 << " ";

        out
            << rods[i].GetCenter().x << " "
            << rods[i].GetCenter().y << " "
            << rods[i].GetCenter().z << " ";

        out
            << rods[i].getTheta() << " "
            << rods[i].getPhi() << " ";

        out
            << rods[i].GetRadius() << " "
            << rods[i].GetLength() << " "
            << rods[i].GetLength() << " ";

        double* q = rods[i].compute_quat();

        out
            << q[1] << " "
            << q[2] << " "
            << q[3] << " "
            << q[0] << " ";

        out << 0 << "\n";

        // delete[] q;
    }

    //-------------------------------------------------------
    // 第二部分：Sample
    //-------------------------------------------------------

    if (output_sample)
    {
        const double r = 0.15;

        for (const auto& s : samples)
        {
            out << id++ << " ";

            out << 2 << " ";

            out
                << s.p.x << " "
                << s.p.y << " "
                << s.p.z << " ";

            out << "0 0 ";

            out
                << r << " "
                << r << " "
                << r << " ";

            out
                << "0 0 0 1 ";

            out << "0\n";
        }
    }

    //-------------------------------------------------------
    // 第三部分：Surface
    //-------------------------------------------------------

    if (output_surface)
    {
        const double rr = 0.10;

        for (const auto& p : surface_points)
        {
            out << id++ << " ";

            out << 3 << " ";

            out
                << p.x << " "
                << p.y << " "
                << p.z << " ";

            out << "0 0 ";

            out
                << rr << " "
                << rr << " "
                << rr << " ";

            out
                << "0 0 0 1 ";

            out << "0\n";
        }
    }

    out.close();

    std::cout
        << "Debug dump saved (step "
        << timestep
        << ") : "
        << filename
        << std::endl;
}


HyperboloidFitResult HyperboloidFit::FitRods(const RigidRod* rods, int count)
{
    if (rods == nullptr || count <= 0)
    {
        samples.clear(); rods_.clear(); surface_points.clear(); rmse = 0.0;
        return {};
    }
    GenerateSamplePoints(rods, count);
    const bool success = FitPreparedSamples();
    return { success, a, b, c, rmse, success ? static_cast<std::size_t>(count) : 0 };
}

HyperboloidFitResult HyperboloidFit::FitSelectedRods(
    const RigidRod* rods, int count, const std::vector<int>& selected_indices)
{
    if (rods == nullptr || count <= 0 || selected_indices.empty() ||
        !GenerateSelectedSamplePoints(rods, count, selected_indices))
    {
        samples.clear(); rods_.clear(); surface_points.clear(); rmse = 0.0;
        return {};
    }
    const bool success = FitPreparedSamples();
    return { success, a, b, c, rmse, success ? rods_.size() : 0 };
}

bool HyperboloidFit::Fit(const RigidRod* rods, int N)
{
    return FitRods(rods, N).success;
}

bool HyperboloidFit::FitSelected(
    const RigidRod* rods, int N, const std::vector<int>& selected_indices)
{
    return FitSelectedRods(rods, N, selected_indices).success;
}

bool HyperboloidFit::GenerateSelectedSamplePoints(
    const RigidRod* rods, int N, const std::vector<int>& selected_indices)
{
    samples.clear();
    rods_.clear();
    const int sample_num = 21;
    for (int index : selected_indices)
    {
        if (index < 0 || index >= N)
        {
            samples.clear(); rods_.clear();
            return false;
        }
        const RigidRod& rod = rods[index];
        rods_.push_back(rod);
        const Vec3 direction = rod.GetDirection().normalized();
        for (int j = 0; j < sample_num; ++j)
        {
            const double s = static_cast<double>(j) / (sample_num - 1) - 0.5;
            SamplePoint point{};
            point.rod_id = index;
            point.s = s;
            point.p = rod.GetCenter() + direction * (s * rod.GetLength());
            samples.push_back(point);
        }
    }
    return !samples.empty();
}

bool HyperboloidFit::FitPreparedSamples()
{
    if (samples.size() < 6)
    {
        rmse = 0.0; surface_points.clear();
        return false;
    }
    BuildLocalCoordinateSystem();
    SearchABC();
    GenerateSurfaceMesh();
    rmse = ComputeRMSE();
    return std::isfinite(rmse) && std::isfinite(a) &&
        std::isfinite(b) && std::isfinite(c) && a > 0 && b > 0 && c > 0;
}
double HyperboloidFit::ComputeRMSE()
{
    if (samples.empty())
        return 0.0;

    double sum = 0.0;

    for (size_t i = 0; i < samples.size(); i++)
    {
        double d = Distance(samples[i].local);

        sum += d * d;
    }

    return std::sqrt(
        sum /
        samples.size());
}

double HyperboloidFit::GetRMSE() const
{
    return rmse;
}

double HyperboloidFit::SurfaceValue(const Vec3& p)
{

    // 防止参数非法
    if (a <= 0.0 || b <= 0.0 || c <= 0.0)
    {
        return 1e30;
    }

    double xx = p.x * p.x / (a * a);

    double yy = p.y * p.y / (b * b);

    double zz = p.z * p.z / (c * c);

    return xx + yy - zz - 1.0;

}

//test surfacevalue function
void HyperboloidFit::TestSurfaceValue()
{
    Vec3 p;

    p.x = a;
    p.y = 0;
    p.z = 0;

    std::cout
        << "F(a,0,0) = "
        << SurfaceValue(p)
        << std::endl;

    p.x = 0;
    p.y = b;
    p.z = 0;

    std::cout
        << "F(0,b,0) = "
        << SurfaceValue(p)
        << std::endl;

    p.x = 0;
    p.y = 0;
    p.z = c;

    std::cout
        << "F(0,0,c) = "
        << SurfaceValue(p)
        << std::endl;
}

double HyperboloidFit::Distance(const Vec3& p)
{
    //--------------------------------------------------
    // 曲面函数值
    //--------------------------------------------------

    double F = SurfaceValue(p);

    //--------------------------------------------------
    // 梯度
    //--------------------------------------------------

    double Fx = 2.0 * p.x / (a * a);

    double Fy = 2.0 * p.y / (b * b);

    double Fz = -2.0 * p.z / (c * c);

    //--------------------------------------------------
    // 梯度长度
    //--------------------------------------------------

    double grad =
        std::sqrt(
            Fx * Fx +
            Fy * Fy +
            Fz * Fz);

    //--------------------------------------------------
    // 防止除零
    //--------------------------------------------------

    if (grad < 1e-12)
        return 1e30;

    //--------------------------------------------------
    // 一阶近似几何距离
    //--------------------------------------------------

    return std::fabs(F) / grad;
}

void HyperboloidFit::SearchABC()
{
    //----------------------------------------------------
    // 根据局部坐标自动估计搜索范围
    //----------------------------------------------------

    double rmax = 0.0;
    double zmax = 0.0;

    for (const auto& sp : samples)
    {
        double r = std::sqrt(
            sp.local.x * sp.local.x +
            sp.local.y * sp.local.y);

        rmax = std::max(rmax, r);

        zmax = std::max(zmax, std::fabs(sp.local.z));
    }

    //----------------------------------------------------
    // 搜索范围
    //----------------------------------------------------

    double a_min = std::max(0.1, rmax * 0.30);
    double a_max = rmax * 1.20;

    double c_min = std::max(0.1, zmax * 0.30);
    double c_max = zmax * 1.20;

    //----------------------------------------------------
    // Grid Search
    //----------------------------------------------------

    const double step = 0.2;

    double best_rmse = 1e30;

    double best_a = a;
    double best_b = b;
    double best_c = c;

    for (double aa = a_min; aa <= a_max; aa += step)
    {
        for (double cc = c_min; cc <= c_max; cc += step)
        {
            a = aa;
            b = aa;      // 当前仍假设旋转对称
            c = cc;

            double err = ComputeRMSE();

            if (err < best_rmse)
            {
                best_rmse = err;

                best_a = aa;
                best_b = aa;
                best_c = cc;
            }
        }
    }

    //----------------------------------------------------
    // 保存最佳参数
    //----------------------------------------------------

    a = best_a;
    b = best_b;
    c = best_c;

    rmse = best_rmse;

    //----------------------------------------------------
    // 输出
    //----------------------------------------------------

    //std::cout << "\n========== Hyperboloid Fit ==========\n";
    //std::cout << "Center = ("
    //    << center_.x << ", "
    //    << center_.y << ", "
    //    << center_.z << ")\n";

    //std::cout << "Axis Z = ("
    //    << axis_z_.x << ", "
    //    << axis_z_.y << ", "
    //    << axis_z_.z << ")\n";

    //std::cout << "a = " << a << std::endl;
    //std::cout << "b = " << b << std::endl;
    //std::cout << "c = " << c << std::endl;
    //std::cout << "RMSE = " << rmse << std::endl;
    //std::cout << "=====================================\n";
}

//generate  hyperboloid curve //// 供 debug.dump 使用
void HyperboloidFit::GenerateSurfaceMesh() {
    surface_points.clear();
    //modify number of point


    for (int iz = -Nz; iz <= Nz; iz++)
    {
        double z = surface_z_range * c * iz / Nz;

        double r =
            std::sqrt(
                1.0 +
                z * z / (c * c));

        double rx = a * r;
        double ry = b * r;

        for (int it = 0; it < Ntheta; it++)
        {
            double theta =
                2.0 * constants::pi
                * it
                / Ntheta;

            Vec3 local;

            local.x = rx * cos(theta);
            local.y = ry * sin(theta);
            local.z = z;

            // 拟合在 PCA 局部坐标系中完成。输出前转换回世界坐标，
            // 使可视化曲面跟随每一帧重新拟合得到的 Z 轴。
            surface_points.push_back(LocalToWorld(local));
        }
    }


}

void HyperboloidFit::OutputSurfaceVTP(const std::string& filename)
{
    std::ofstream out(filename);

    if (!out)
    {
        std::cerr << "Cannot open " << filename << std::endl;
        return;
    }

    //----------------------------------------------------------
    // Resolution
    //----------------------------------------------------------

    const int Ntheta = 72;
    const int Nz = 80;

    double zmin = -surface_z_range * c;
    double zmax = surface_z_range * c;

    //----------------------------------------------------------
    // Generate Points
    //----------------------------------------------------------

    std::vector<Vec3> points;

    points.reserve(Ntheta * Nz);

    for (int itheta = 0; itheta < Ntheta; itheta++)
    {
        double theta =
            2.0 * constants::pi *
            itheta / Ntheta;

        for (int iz = 0; iz < Nz; iz++)
        {
            double t =
                double(iz) /
                (Nz - 1);

            double z =
                zmin +
                t * (zmax - zmin);

            double r =
                sqrt(
                    1.0 +
                    z * z /
                    (c * c));

            Vec3 local;

            local.x = a * r * cos(theta);
            local.y = b * r * sin(theta);
            local.z = z;

            points.push_back(LocalToWorld(local));
        }
    }

    //----------------------------------------------------------
    // XML Header
    //----------------------------------------------------------

    out << "<?xml version=\"1.0\"?>\n";

    out << "<VTKFile "
        << "type=\"PolyData\" "
        << "version=\"0.1\" "
        << "byte_order=\"LittleEndian\">\n";

    out << "<PolyData>\n";

    //----------------------------------------------------------
    // Number
    //----------------------------------------------------------

    const int NumberOfPoints =
        (int)points.size();

    const int NumberOfPolys =
        Ntheta * (Nz - 1) * 2;

    out
        << "<Piece "
        << "NumberOfPoints=\""
        << NumberOfPoints
        << "\" "
        << "NumberOfPolys=\""
        << NumberOfPolys
        << "\">\n";

    //----------------------------------------------------------
    // Points
    //----------------------------------------------------------

    out << "<Points>\n";

    out << "<DataArray "
        << "type=\"Float64\" "
        << "NumberOfComponents=\"3\" "
        << "format=\"ascii\">\n";

    for (const auto& p : points)
    {
        out
            << p.x << " "
            << p.y << " "
            << p.z << "\n";
    }

    out << "</DataArray>\n";
    out << "</Points>\n";

    //----------------------------------------------------------
    // Polys
    //----------------------------------------------------------

    out << "<Polys>\n";

    //---------------- connectivity ----------------

    out
        << "<DataArray "
        << "type=\"Int32\" "
        << "Name=\"connectivity\" "
        << "format=\"ascii\">\n";

    for (int itheta = 0; itheta < Ntheta; itheta++)
    {
        int nextTheta =
            (itheta + 1) % Ntheta;

        for (int iz = 0; iz < Nz - 1; iz++)
        {
            int A =
                itheta * Nz + iz;

            int B =
                nextTheta * Nz + iz;

            int C =
                nextTheta * Nz + iz + 1;

            int D =
                itheta * Nz + iz + 1;

            out
                << A << " "
                << B << " "
                << C << "\n";

            out
                << A << " "
                << C << " "
                << D << "\n";
        }
    }

    out << "</DataArray>\n";

    //---------------- offsets ----------------

    out
        << "<DataArray "
        << "type=\"Int32\" "
        << "Name=\"offsets\" "
        << "format=\"ascii\">\n";

    int offset = 0;

    for (int i = 0; i < NumberOfPolys; i++)
    {
        offset += 3;
        out << offset << "\n";
    }

    out << "</DataArray>\n";

    out << "</Polys>\n";

    //----------------------------------------------------------
    // CellData
    //----------------------------------------------------------

    out
        << "<CellData Scalars=\"RMSE\">\n";

    out
        << "<DataArray "
        << "type=\"Float64\" "
        << "Name=\"RMSE\" "
        << "format=\"ascii\">\n";

    for (int i = 0; i < NumberOfPolys; i++)
    {
        out << rmse << "\n";
    }

    out << "</DataArray>\n";

    out << "</CellData>\n";

    //----------------------------------------------------------
    // Finish
    //----------------------------------------------------------

    out << "</Piece>\n";
    out << "</PolyData>\n";
    out << "</VTKFile>\n";

    out.close();

   /* std::cout
        << "Surface VTP saved : "
        << filename
        << std::endl;*/
}

void HyperboloidFit::OutputRodVTP(const std::string& filename)
{
    std::ofstream out(filename);

    if (!out) return;

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "<PolyData>\n";

    out << "<Piece NumberOfPoints=\""
        << rods_.size() * 2
        << "\" NumberOfLines=\""
        << rods_.size()
        << "\">\n";

    //---------------- Points ----------------

    out << "<Points>\n";

    out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (auto& rod : rods_)
    {
        Vec3 dir =
            rod.GetDirection().normalized();
        Vec3 p1 =
            rod.GetCenter()
            - rod.GetDirection()
            * rod.GetLength() * 0.5;

        Vec3 p2 =
            rod.GetCenter()
            + rod.GetDirection()
            * rod.GetLength() * 0.5;

        out << p1.x << " "
            << p1.y << " "
            << p1.z << "\n";

        out << p2.x << " "
            << p2.y << " "
            << p2.z << "\n";
    }

    out << "</DataArray>\n";
    out << "</Points>\n";

    //---------------- connectivity ----------------

    out << "<Lines>\n";

    out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    for (int i = 0; i < (int)rods_.size(); i++)
    {
        out << 2 * i
            << " "
            << 2 * i + 1
            << "\n";
    }

    out << "</DataArray>\n";

    //---------------- offsets ----------------

    out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

    int offset = 0;

    for (int i = 0; i < (int)rods_.size(); i++)
    {
        offset += 2;
        out << offset << "\n";
    }

    out << "</DataArray>\n";

    out << "</Lines>\n";

    //---------------- Radius ----------------

    out << "<CellData Scalars=\"Radius\">\n";

    out << "<DataArray "
        << "type=\"Float64\" "
        << "Name=\"Radius\" "
        << "format=\"ascii\">\n";

    for (auto& rod : rods_)
        out << rod.GetRadius() << "\n";

    out << "</DataArray>\n";

    out << "</CellData>\n";

    out << "</Piece>\n";
    out << "</PolyData>\n";
    out << "</VTKFile>\n";
}

void HyperboloidFit::OutputPVD(
    const std::string& filename,
    const std::string& prefix,
    int frame_num)
{
    std::ofstream out(filename);

    if (!out)
    {
        std::cerr
            << "Cannot open "
            << filename
            << std::endl;
        return;
    }

    out << "<?xml version=\"1.0\"?>\n";

    out << "<VTKFile "
        << "type=\"Collection\" "
        << "version=\"0.1\" "
        << "byte_order=\"LittleEndian\">\n";

    out << "  <Collection>\n";

    for (int i = 0; i <= frame_num; i++)
    {
        char vtpname[128];

        sprintf_s(
            vtpname,
            sizeof(vtpname),
            "%s_%06d.vtp",
            prefix.c_str(),
            i);

        out
            << "    <DataSet "
            << "timestep=\"" << i << "\" "
            << "part=\"0\" "
            << "file=\"" << vtpname << "\"/>\n";
    }

    out << "  </Collection>\n";

    out << "</VTKFile>\n";
}

double HyperboloidFit::GetA() const
{
    return a;
}

double HyperboloidFit::GetB() const
{
    return b;
}

double HyperboloidFit::GetC() const
{
    return c;
}

void HyperboloidFit::BuildLocalCoordinateSystem()
{
    //--------------------------------------------------
    // Center
    //--------------------------------------------------

    center_ = Vec3(0, 0, 0);

    for (const auto& s : samples)
        center_ = center_ + s.p;

    center_ = center_ * (static_cast<double>(1) / samples.size());

    //--------------------------------------------------
    // Covariance
    //--------------------------------------------------

    Eigen::Matrix3d cov;

    cov.setZero();

    for (const auto& s : samples)
    {
        Vec3 d = s.p - center_;

        Eigen::Vector3d v(
            d.x,
            d.y,
            d.z);

        cov += v * v.transpose();
    }

    cov /= samples.size();

    //--------------------------------------------------
    // PCA
    //--------------------------------------------------

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);

    Eigen::Matrix3d R =
        solver.eigenvectors();

    axis_x_ = Vec3(
        R(0, 0),
        R(1, 0),
        R(2, 0));

    axis_y_ = Vec3(
        R(0, 1),
        R(1, 1),
        R(2, 1));

    axis_z_ = Vec3(
        R(0, 2),
        R(1, 2),
        R(2, 2));

    //--------------------------------------------------
    // World -> Local
    //--------------------------------------------------

    for (auto& s : samples)
    {
        Vec3 d =
            s.p - center_;

        s.local.x =
            d.dot(axis_x_);

        s.local.y =
            d.dot(axis_y_);

        s.local.z =
            d.dot(axis_z_);
    }
}

Vec3 HyperboloidFit::LocalToWorld(const Vec3& p) const
{
    return center_
        + axis_x_ * p.x
        + axis_y_ * p.y
        + axis_z_ * p.z;
}
