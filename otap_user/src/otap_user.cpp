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
#include <chrono>

#include "otap.h"

int main()
{
    using namespace OTAP;

    // Geometry geom;
    // geom.Push(GeometryPrimitiveType::Arc, 0.0, 20.0, 0.2, 0.0);

    // auto traj = VelocityTrajectory(AmbientType::ISA, FluidType::Hansen_Air,
    //                                "c:/Users/Arnab Mahanti/source/repos/otap_pp/docs/trajectory.dat");
    // auto ptraj = make_trajectory(TrajectoryType::Velocity,AmbientType::ISA, FluidType::Hansen_Air,
    //                                "c:/Users/Arnab Mahanti/source/repos/otap_pp/docs/trajectory.dat");

    // TimeSeries P0(traj.GetTimePoints().size());
    // TimeSeries T0(traj.GetTimePoints().size());
    // TimeSeries Pe(traj.GetTimePoints().size());
    // TimeSeries Te(traj.GetTimePoints().size());
    // TimeSeries Me(traj.GetTimePoints().size());

    // TimeSeries CPress(traj.GetTimePoints().size());
    // TimeSeries Minf(traj.GetTimePoints().size());
    // TimeSeries Pinf(traj.GetTimePoints().size());
    // TimeSeries Tinf(traj.GetTimePoints().size());
    // TimeSeries a1(traj.GetTimePoints().size());

    // TimePoints t_result = traj.GetTimePoints();

    // // auto cpdata = CpData("C:/MINIVER_Aayushi_2022_latest_25Nov/VSCODE_OTAP/otap_pp/docs/test_cp.txt");

    // // double Position = 7.1;

    // auto fluid = make_fluid(FluidType::Hansen_Air);
    // // double CoefPress;

    // // for (size_t i = 0; i < t_result.size(); i++)
    // // {
    // //     Upstream us;
    // //     us.M = traj.GetMinf(t_result[i]);
    // //     us.P = traj.GetPinf(t_result[i]);
    // //     us.T = traj.GetTinf(t_result[i]);

    // //     auto ds = fluid->GetShockDownStream(us);
    // //     CoefPress = cpdata.GetCp(us.M, Position);
    // //     auto edge = fluid->GetEdgeProperties(us, CoefPress);

    // //     P0[i] = ds.P0;
    // //     T0[i] = ds.T0;
    // //     Pe[i] = edge.P;
    // //     Me[i] = edge.M;
    // //     Te[i] = edge.T;
    // //     CPress[i] = CoefPress;
    // //     Minf[i] = us.M;
    // //     Pinf[i] = us.P;
    // //     Tinf[i] = us.T;
    // //     a1[i] = fluid->a(Pinf[i], Tinf[i]);
    // // }

    // // std::ofstream outcheck("C:/MINIVER_Aayushi_2022_latest_25Nov/VSCODE_OTAP/otap_pp/docs/check.txt");

    // //  for (size_t i = 0; i < t_result.size(); i++)
    // //      {
    // //         //outcheck << Minf[i] << " " << CPress[i] << " " << P0[i] << " " << T0[i] << " " << Pe[i] << " " << Pinf[i] << " " << a1[i]*Minf[i] << "\n";

    // //          outcheck << Me[i] << " " << "\n";
    // //      }

    // // std::ofstream out("C:/MINIVER_Aayushi_2022_latest_25Nov/VSCODE_OTAP/otap_pp/docs/props.dat");
    // // auto Pr = fluid->Pr(P0, T0);
    // // auto mu = fluid->mu(P0, T0);
    // // auto Rho = fluid->Rho(P0, T0);
    // // auto Z = fluid->Z(P0, T0);
    // // auto H = fluid->H(P0, T0);
    // // auto Cp = fluid->Cp(P0, T0);

    // // for (size_t i = 0; i < t_result.size(); i++)
    // // {
    // //     out << t_result[i] << " " << P0[i] << " " << T0[i]<< " " << Pr[i] << " " << H[i] << " " << mu[i] << " " << Rho[i]<< '\n';
    // // }

    // // FayRiddellSolver solver(
    // //     CoordinateType::Axisym,
    // //     FlowType::Default,
    // //     fluid,
    // //     geom.GetComponents()[1],
    // //     traj.GetPinf(t_result),
    // //     traj.GetVinf(t_result),
    // //     P0,
    // //     T0,
    // //     300);

    // // VanDriestSolver solver(
    // //     CoordinateType::Axisym,
    // //     FlowType::Turbulent,
    // //     fluid,
    // //     geom.GetComponents()[1],
    // //     7.6612,
    // //     traj.GetPinf(t_result),
    // //     traj.GetVinf(t_result),
    // //     P0,
    // //     T0,
    // //     Me,
    // //     Te,
    // //     Pe,
    // //     300,
    // //     1000000,
    // //     0);

    // // auto res = make_HFSolver(HFSolverType::Eckert,
    // //                          CoordinateType::Axisym,
    // //                          FlowType::Turbulent,
    // //                          fluid,
    // //                          geom.GetComponents()[1],
    // //                          7.6612,
    // //                          traj.GetPinf(t_result),
    // //                          traj.GetVinf(t_result),
    // //                          P0,
    // //                          T0,
    // //                          Me,
    // //                          Te,
    // //                          Pe,
    // //                          300.0,
    // //                          0000000.0,
    // //                          0.0);

    // //  res->Solve();

    // //   auto res = solver.Solve();

    // // std::cout << t_result.size() << "\n";

    // // for (size_t i = 0; i <  t_result.size(); i++)
    // // {
    // // std::cout << Me[i] << "\n" ;
    // // }

    // // std::ofstream out_res("C:/MINIVER_Aayushi_2022_latest_25Nov/VSCODE_OTAP/otap_pp/docs/trajectory_q.dat");

    // // for (size_t i = 0; i < res.h.size(); i++)
    // // {
    // //     out_res << t_result[i] << " " << res.q[i] << " " << res.h[i] << " " << res.Tg[i] << '\n';
    // // }

    // // OTAP::Material mat1(121, 2800, 960, 0.3, 5000, 5000);

    // // Layer l1(mat1, 1, 10);
    // // Layer l2(mat1, 1, 10);

    // // TimeSeries initT(3);

    // // initT[0] = 300;
    // // initT[1] = 300;
    // // initT[2] = 300;

    // // LayerStack LH;
    // // LH.Push(l1);
    // // LH.Push(l2);

    // // LH.InitTemperature(initT);
    // // LH.CreateNodes();

    // // // std::cout << LH.Nodes[1].k() << "\n";

    // // ResponseSolverParams params;

    // // params.tInit = 0.0;
    // // params.tFinal;
    // // params.timestep = {{0, 0}};
    // // params.innerRadius;
    // // params.massRateChar = 0.0;
    // // params.massRatePyro = 0.0;
    // // params.propellantMass = {{0, 0}};

    // // // Frontwall radiation
    // // params.qRad = {{0, 0}};
    // // params.TRad = 300.0;
    // // // Table qGen;

    // // // Backwall radiation and convection
    // // params.hBackwall = 0.0; // FIXME: Give Options for h and characteristic length for natural convection
    // // params.TgBackwall = 300.0;
    // // params.TRadBackwall = 300.0;

    // // ResponseSolverOptions options_1;

    // // options_1.coordinateType = CoordinateType::Axisym;
    // // options_1.Ttol = 1.0;

    // // auto fluid = make_fluid(FluidType::Hansen_Air);

    // // auto hfs1 = make_HFSolver(HFSolverType::VanDriest,
    // //                           CoordinateType::Axisym,
    // //                           FlowType::Turbulent,
    // //                           fluid,
    // //                           geom.GetComponents()[1],
    // //                           7.6612,
    // //                           traj.GetPinf(t_result),
    // //                           traj.GetVinf(t_result),
    // //                           P0,
    // //                           T0,
    // //                           Me,
    // //                           Te,
    // //                           Pe,
    // //                           300.0,
    // //                           0000000.0,
    // //                           0);

    // // auto hfs2 = make_HFSolver(HFSolverType::VanDriest,
    // //                           CoordinateType::Axisym,
    // //                           FlowType::Turbulent,
    // //                           fluid,
    // //                           geom.GetComponents()[1],
    // //                           7.6612,
    // //                           traj.GetPinf(t_result),
    // //                           traj.GetVinf(t_result),
    // //                           P0,
    // //                           T0,
    // //                           Me,
    // //                           Te,
    // //                           Pe,
    // //                           300.0,
    // //                           0000000.0,
    // //                           0);

    // // DefaultResponseSolver Solver(
    // //     &LH, hfs1, hfs2,
    // //     initT,
    // //     params,
    // //     options_1);

    // // CreateNodes();

    // // table<double> t({1, 2, 3}, {4, 5, 6});

    // // MaterialManager manager("C:/Users/Arnab Mahanti/source/repos/otap_pp/docs/material_database.yaml");
    // // auto mat = manager.GetMaterialInstance("aluminium");

    // // std::cout << mat->k[500] << "\n";

    // // RadiationBC rbc1, rbc2;
    // // auto rbc1 = make_BC(BCType::Radiation);
    // // auto rbc2 = make_BC(BCType::Radiation);

    // // auto a = std::shared_ptr<RadiationBC>(new RadiationBC());
    // // BC bc;
    // // bc = a;

    // // a->Tambient = {{1, 100}, {300, 300}};
    // // a->location = BCLocation::front;

    // // rbc2->Tambient = {{1, 100}, {200, 200}};
    // // rbc2->location = BCLocation::back;

    // // auto conv = make_BC(BCType::Convection);
    // // conv->location = BCLocation::front;
    // // // conv->h = {{1, 100}, {200, 200}};
    // // // conv->tg = {{1, 100}, {200, 200}};

    // // BCArray bcarray = {rbc1, rbc2, conv};

    // auto frbridge = KnBridgedSolver(
    //                               1.0,
    //                               {},{},{},
    //                               HFSolverType::FayRiddell,
    //                               CoordinateType::Axisym,
    //                               FlowType::Default,
    //                               fluid,
    //                               geom.GetComponents()[1],
    //                               ptraj,
    //                               300.);
    // // auto test = KnBridgedSolver(1,{},{},{},)

    // auto res = frbridge.Solve(20., 300.);

    // TEST

    Geometry geom;
    geom.Push(GeometryPrimitiveType::Arc,
              0, 45, 0.7, 0);

    auto fluid = make_fluid(FluidType::Hansen_Air);

    auto manager = MaterialManager("C:/Users/Arnab Mahanti/source/repos/otap_pp/docs/material_database.yaml");
    auto al = manager.GetMaterialInstance("aluminium");
    auto cp = manager.GetMaterialInstance("carbon_phenolic");

    auto traj = make_trajectory(TrajectoryType::Velocity,
                                AmbientType::ISA,
                                FluidType::Hansen_Air,
                                "C:/Users/Arnab Mahanti/source/repos/otap_pp/docs/trajectory.dat");

    auto layers = std::make_shared<LayerStack>();
    layers->Emplace(cp, 0, 5);
    layers->Emplace(cp, 0.0, 5);
    layers->Emplace(al, 0.002, 5);

    TimeSeries Tinit{300, 300, 300, 300};

    auto hf = make_HFSolver(HFSolverType::KnBridged,
                            0.7,
                            std::initializer_list<double>{},
                            std::initializer_list<double>{},
                            std::initializer_list<double>{},
                            HFSolverType::FayRiddell,
                            CoordinateType::Axisym,
                            FlowType::Turbulent,
                            fluid,
                            geom.GetComponents()[1],
                            traj,
                            300);

    auto bcs = BCArray();

    ResponseSolverOptions options;
    options.coordinateType = CoordinateType::Cartesian;

    ResponseSolverParams params;
    params.tFinal = 150;
    params.timestep = {{0, 150}, {.1, .1}};
    params.Ttol = 0.01;

    auto response = make_ResponseSolver(ResponseSolverType::Default,
                                        layers,
                                        hf,
                                        Tinit,
                                        bcs,
                                        params,
                                        options);

    auto start = std::chrono::high_resolution_clock::now();
    auto solution = response->Solve();

    std::ofstream outfile("C:/Users/Arnab Mahanti/source/repos/otap_pp/docs/solution.csv");
    outfile << " " << ", ";
    for (size_t i = 0; i < solution.solution_T[0].size(); i++)
    {
        outfile << i << ", ";
    }
    outfile << '\n';

    for (size_t i = 0; i < solution.solution_T.size(); i++)
    {
        outfile << solution.solution_t[i] << ", ";
        for (auto &&j : solution.solution_T[i])
        {
            outfile << j << ", ";
        }
        outfile << '\n';
    }

    std::cout << "Computed in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << "ms \n";

    return 0;
}
