#pragma once
#include <vector>
#include "particles.h"   // Legacy comment normalized to UTF-8.

struct HyperboloidFitResult
{
    bool success = false;
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double rmse = 0.0;
    std::size_t rods_used = 0;
};

class HyperboloidFit
{
public:
    HyperboloidFit();

    // Fit all rods, or only the rods named by indices into the input array.
    // The returned value owns the scalar results and is independent of System.
    HyperboloidFitResult FitRods(const RigidRod* rods, int count);
    HyperboloidFitResult FitSelectedRods(
        const RigidRod* rods,
        int count,
        const std::vector<int>& selected_indices);

    // Legacy stateful entry points retained for existing callers.
    bool Fit(const RigidRod* rods, int N);

    bool FitSelected(
        const RigidRod* rods,
        int N,
        const std::vector<int>& selected_indices);

    double ComputeRMSE();

    void OutputSurfaceVTP(const std::string& filename);
    void OutputRodVTP(const std::string& filename);

    void OutputPVD(
        const std::string& filename,
        const std::string& prefix,
        int frame_num);
    double GetRMSE() const;
    //inspect function
    void PrintSampleInfo();
    //output function
    void OutputSamplePoints(const std::string& filename);

    void TestSamplePoints(
        const RigidRod* rods,
        int N,
        double R_out);
    void OutputDebugDump(
        const std::string& filename,
        const RigidRod* rods,
        int N,
        double R_out,
        int timestep,
        bool append = true,
        bool output_sample = true,
        bool output_surface = true);
    void TestSurfaceValue();

    double GetA() const;
    double GetB() const;
    double GetC() const;
private:
    // surface axis
    Vec3 center_;

    Vec3 axis_x_;
    Vec3 axis_y_;
    Vec3 axis_z_;

private:


    struct SamplePoint
    {
        //------------------------------
        // Legacy comment normalized to UTF-8.
        //------------------------------
        Vec3 p;

        //------------------------------
        // Legacy comment normalized to UTF-8.
        //------------------------------
        Vec3 local;

        //------------------------------
        // Legacy comment normalized to UTF-8.
        //------------------------------
        int rod_id;

        //------------------------------
        // Legacy comment normalized to UTF-8.
        //------------------------------
        double s;

        //------------------------------
        // Legacy comment normalized to UTF-8.
        //------------------------------
        double weight;
    };

    std::vector<SamplePoint> samples;
private:

    std::vector<Vec3> surface_points;
    std::vector<RigidRod> rods_;


    double a, b, c;


    double rmse;

private:

    //debug
    const int Nz = 30;
    const int Ntheta = 30;
    //VTK
    int surface_line_num = 24;
    int surface_point_per_line = 40;
    double surface_z_range = 2;

    void GenerateSamplePoints(
        const RigidRod* rods,
        int N);

    bool GenerateSelectedSamplePoints(
        const RigidRod* rods,
        int N,
        const std::vector<int>& selected_indices);

    bool FitPreparedSamples();

    void BuildLocalCoordinateSystem();
    Vec3 LocalToWorld(const Vec3& p) const;
    double SurfaceValue(const Vec3& p);

    double Distance(const Vec3& p);

    void SearchABC();
    void GenerateSurfaceMesh();
};
