#include <iostream>
// #include "Problem.h"
// #include "Geometry.h"
// #include "Trajectory.h"
#include <variant>
#include <optional>
#include <vector>
#include <memory>
#include <fstream>
#include "mlinterp/mlinterp.hpp"
#include <algorithm>
#include <cmath>
#include <valarray>
#include <iostream>
#include <iomanip>

#include "otap.h"

template <int t>
struct test
{
    int a = t;
};

enum class type
{
    a,
    b,
    c
};

struct data
{
    double d;
    double e;
    data(double a, double b) : d(a), e(b) {}
    ~data()
    {
        d = 0;
        e = 0;
        //  std::cout << "destroyed{" << d << " " << e << "}\n";
    }
};

auto process(std::valarray<double> v)
{
    v = v.apply(std::log10);

    return std::vector(std::begin(v), std::end(v));
}
int main()
{
    using namespace OTAP;

    Geometry geom;
    geom.Push(GeometryPrimitiveType::Arc, 0.0, 20.0, 0.2, 0.0);

    auto traj = VelocityTrajectory(AmbientType::ISA, FluidType::Hansen_Air,
                                   "C:/MINIVER_Aayushi_2022_latest_25Nov/VSCODE_OTAP/otap_pp/docs/trajectory.dat");

    TimeSeries P0(traj.GetTimePoints().size());
    TimeSeries T0(traj.GetTimePoints().size());
    TimeSeries Pe(traj.GetTimePoints().size());
    TimeSeries Te(traj.GetTimePoints().size());
    TimeSeries Me(traj.GetTimePoints().size());

    TimeSeries CPress(traj.GetTimePoints().size());
    TimeSeries Minf(traj.GetTimePoints().size());
    TimeSeries Pinf(traj.GetTimePoints().size());
    TimeSeries Tinf(traj.GetTimePoints().size());
    TimeSeries a1(traj.GetTimePoints().size());

    TimePoints t_result = traj.GetTimePoints();

    auto cpdata = CpData("C:/MINIVER_Aayushi_2022_latest_25Nov/VSCODE_OTAP/otap_pp/docs/test_cp.txt");

    double Position = 7.1;

    auto fluid = make_fluid(FluidType::Hansen_Air);
    double CoefPress;

    for (size_t i = 0; i < t_result.size(); i++)
    {
        Upstream us;
        us.M = traj.GetMinf(t_result[i]);
        us.P = traj.GetPinf(t_result[i]);
        us.T = traj.GetTinf(t_result[i]);

        auto ds = fluid->GetShockDownStream(us);
        CoefPress = cpdata.GetCp(us.M, Position);
        auto edge = fluid->GetEdgeProperties(us, CoefPress);

        P0[i] = ds.P0;
        T0[i] = ds.T0;
        Pe[i] = edge.P;
        Me[i] = edge.M;
        Te[i] = edge.T;
        CPress[i] = CoefPress;
        Minf[i] = us.M;
        Pinf[i] = us.P;
        Tinf[i] = us.T;
        a1[i] = fluid->a(Pinf[i], Tinf[i]);
    }

    // std::ofstream outcheck("C:/MINIVER_Aayushi_2022_latest_25Nov/VSCODE_OTAP/otap_pp/docs/check.txt");

    //  for (size_t i = 0; i < t_result.size(); i++)
    //      {
    //         //outcheck << Minf[i] << " " << CPress[i] << " " << P0[i] << " " << T0[i] << " " << Pe[i] << " " << Pinf[i] << " " << a1[i]*Minf[i] << "\n";

    //          outcheck << Me[i] << " " << "\n";
    //      }

    // std::ofstream out("C:/MINIVER_Aayushi_2022_latest_25Nov/VSCODE_OTAP/otap_pp/docs/props.dat");
    // auto Pr = fluid->Pr(P0, T0);
    // auto mu = fluid->mu(P0, T0);
    // auto Rho = fluid->Rho(P0, T0);
    // auto Z = fluid->Z(P0, T0);
    // auto H = fluid->H(P0, T0);
    // auto Cp = fluid->Cp(P0, T0);

    // for (size_t i = 0; i < t_result.size(); i++)
    // {
    //     out << t_result[i] << " " << P0[i] << " " << T0[i]<< " " << Pr[i] << " " << H[i] << " " << mu[i] << " " << Rho[i]<< '\n';
    // }

    // FayRiddellSolver solver(
    //     CoordinateType::Axisym,
    //     FlowType::Default,
    //     fluid,
    //     geom.GetComponents()[1],
    //     traj.GetPinf(t_result),
    //     traj.GetVinf(t_result),
    //     P0,
    //     T0,
    //     300);



    VanDriestSolver solver(
        CoordinateType::Axisym,
        FlowType::Turbulent,
        fluid,
        geom.GetComponents()[1],
        7.6612,
        traj.GetPinf(t_result),
        traj.GetVinf(t_result),
        P0,
        T0,
        Me,
        Te,
        Pe,
        300,
        1000000,
        0);



// auto res = make_HFSolver( HFSolverType::Eckert,
//         CoordinateType::Axisym,
//         FlowType::Turbulent,
//         fluid,
//         geom.GetComponents()[1],
//         7.6612,
//         traj.GetPinf(t_result),
//         traj.GetVinf(t_result),
//         P0,
//         T0,
//         Me,
//         Te,
//         Pe,
//         300.0,
//         0000000.0,
//         0);



//  res->Solve();

 //   auto res = solver.Solve();

    // std::cout << t_result.size() << "\n";

    // for (size_t i = 0; i <  t_result.size(); i++)
    // {
    // std::cout << Me[i] << "\n" ;
    // }

    // std::ofstream out_res("C:/MINIVER_Aayushi_2022_latest_25Nov/VSCODE_OTAP/otap_pp/docs/trajectory_q.dat");

    // for (size_t i = 0; i < res.h.size(); i++)
    // {
    //     out_res << t_result[i] << " " << res.q[i] << " " << res.h[i] << " " << res.Tg[i] << '\n';
    // }

OTAP::Material mat1(121,2800,960,0.3,5000,5000);

Layer l1(mat1,1,10);
Layer l2(mat1,1,10);


TimeSeries initT(3);

initT[0]=300;
initT[1]=300;
initT[2]=300;

LayerStack LH;
LH.Push(l1);
LH.Push(l2);


LH.InitTemperature(initT);
LH.CreateNodes();

//std::cout << LH.Nodes[1].k() << "\n";

ResponseSolverParams params;

        params.tInit = 0.0;
        params.tFinal;
        params.timestep={{0,0}};
       params.innerRadius;
       params.massRateChar = 0.0;
        params.massRatePyro = 0.0;
        params.propellantMass={{0,0}};

        // Frontwall radiation
        params.qRad={{0,0}};
        params.TRad = 300.0;
        // Table qGen;

        // Backwall radiation and convection
        params.hBackwall = 0.0; // FIXME: Give Options for h and characteristic length for natural convection
        params.TgBackwall = 300.0;
        params.TRadBackwall = 300.0;

 ResponseSolverOptions options_1;

 options_1.coordinateType = CoordinateType::Axisym;
 options_1.Ttol = 1.0;

//  DefaultResponseSolver Solver
//  (
//     &LH, 
//     make_HFSolver(HFSolverType::VanDriest,
//         CoordinateType::Axisym,
//         FlowType::Turbulent,
//         fluid,
//         geom.GetComponents()[1],
//         7.6612,
//         traj.GetPinf(t_result),
//         traj.GetVinf(t_result),
//         P0,
//         T0,
//         Me,
//         Te,
//         Pe,
//         300.0,
//         0000000.0,
//         0),
//        make_HFSolver(HFSolverType::VanDriest,
//         CoordinateType::Axisym,
//         FlowType::Turbulent,
//         fluid,
//         geom.GetComponents()[1],
//         7.6612,
//         traj.GetPinf(t_result),
//         traj.GetVinf(t_result),
//         P0,
//         T0,
//         Me,
//         Te,
//         Pe,
//         300.0,
//         0000000.0,
//         0),
//     initT, 
//     params, 
//     options_1
// );






// CreateNodes();

    return 0;
}
