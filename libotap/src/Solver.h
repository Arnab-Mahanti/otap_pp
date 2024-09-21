#pragma once
#include <memory>
#include <string>

#include "OTypes.h"
#include "Fluid.h"
#include "Ambient.h"
#include "Geometry.h"
#include "MaterialLayer.h"
#include "Trajectory.h"
#include "Table.h"
#include "BC.h"
#include "Eigen/Eigen"

namespace OTAP
{

    // FIXME: Create mechanism to change the data without recreating the objects. Add methods to invalidate and validate and solve only if required (cached results)
    // FIXME: Response Solver would need an update after this change.
    // Heat Flux Solvers
    struct HFResult
    {
        TimeSeries h;
        TimeSeries Tg;
        TimeSeries q;
    };

    class HFSolverBase
    {
    protected:
        CoordinateType m_CoordinateType;
        FlowType m_FlowType;
        Fluid m_Fluid;
        bool m_hasTrajectory = false;
        Trajectory m_Trajectory = nullptr;

    public:
        HFSolverBase(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid)
            : m_CoordinateType(coordinateType), m_FlowType(flowType), m_Fluid(fluid) {}
        HFSolverBase(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, Trajectory trajectory)
            : m_CoordinateType(coordinateType), m_FlowType(flowType), m_Fluid(fluid), m_Trajectory(trajectory) { m_hasTrajectory = true; }
        virtual ~HFSolverBase() = default;
        Trajectory GetTrajectory() const { return m_Trajectory; }
        Fluid GetFluid() const { return m_Fluid; }
        virtual HFResult Solve(double time = 0.0) { return {{0},{0},{0}};}
        virtual HFResult Solve(double time, double Twall) { return Solve(); }
    };

    class FayRiddellSolver : public HFSolverBase
    {
    private:
        TimeSeries m_Pinf;
        TimeSeries m_Vinf;

        TimeSeries m_P0;
        TimeSeries m_T0;
        GeometryPrimitive m_Geom;

        double m_Tw;

        double m_currentTime = 0.0;

    public:
        FayRiddellSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, TimeSeries pinf, TimeSeries vinf, TimeSeries P0, TimeSeries T0, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_Pinf(pinf), m_Vinf(vinf), m_P0(P0), m_T0(T0), m_Tw(Tw) {}
        FayRiddellSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, Trajectory traj, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid, traj), m_Geom(geom), m_Tw(Tw)
        {
            m_hasTrajectory = true;
        }
        ~FayRiddellSolver() = default;
        virtual HFResult Solve(double time = 0.0) override;
        virtual HFResult Solve(double time, double Twall) override
        {
            m_Tw = Twall;
            m_currentTime = time;
            return Solve();
        }
    };

    class LeesSolver : public HFSolverBase
    {
    private:
        TimeSeries m_Pinf;
        TimeSeries m_Vinf;
        TimeSeries m_Tinf;
        TimeSeries m_Minf;
        TimeSeries m_P0;
        TimeSeries m_T0;
        TimeSeries m_qstag;
        GeometryPrimitive m_Geom;
        double m_Tw;
        double m_RunningLength;
        double m_currentTime;

    public:
        // FIXME: Calculate qstag from Fay & Riddelle
        LeesSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double RunningLength, TimeSeries pinf, TimeSeries Tinf, TimeSeries vinf, TimeSeries minf, TimeSeries P0, TimeSeries T0, TimeSeries qstag, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_RunningLength(RunningLength), m_Pinf(pinf), m_Tinf(Tinf), m_Vinf(vinf), m_Minf(minf), m_P0(P0), m_T0(T0), m_qstag(qstag), m_Tw(Tw) {}
        LeesSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double RunningLength, Trajectory traj, TimeSeries qstag, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid, traj), m_Geom(geom), m_RunningLength(RunningLength), m_qstag(qstag), m_Tw(Tw) { m_hasTrajectory = true; }
        ~LeesSolver() = default;
        virtual HFResult Solve(double time = 0.0) override;
        virtual HFResult Solve(double time, double Twall) override
        {
            m_Tw = Twall;
            m_currentTime = time;
            return Solve();
        }
    };

    class FreeMolecularSolver : public HFSolverBase
    {
    private:
        TimeSeries m_Pinf;
        TimeSeries m_Vinf;
        TimeSeries m_Tinf;
        GeometryPrimitive m_Geom;

        double m_Tw;
        double m_currentTime;

    public:
        FreeMolecularSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, TimeSeries pinf, TimeSeries vinf, TimeSeries Tinf, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_Pinf(pinf), m_Vinf(vinf), m_Tinf(Tinf), m_Tw(Tw) {}
        FreeMolecularSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, Trajectory traj, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid, traj), m_Geom(geom), m_Tw(Tw) { m_hasTrajectory = true; }
        ~FreeMolecularSolver() = default;
        virtual HFResult Solve(double time = 0.0) override;
        virtual HFResult Solve(double time, double Twall) override
        {
            m_Tw = Twall;
            m_currentTime = time;
            return Solve();
        }
    };

    class BeckwithGallagherSolver : public HFSolverBase
    {
    private:
        TimeSeries m_Pinf;
        TimeSeries m_Vinf;
        TimeSeries m_Tinf;
        TimeSeries m_Minf;
        TimeSeries m_P0;
        TimeSeries m_T0;
        GeometryPrimitive m_Geom;
        double m_Tw;
        double m_currentTime;

    public:
        BeckwithGallagherSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, TimeSeries pinf, TimeSeries Tinf, TimeSeries vinf, TimeSeries Minf, TimeSeries P0, TimeSeries T0, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_Pinf(pinf), m_Tinf(Tinf), m_Vinf(vinf), m_Minf(Minf), m_P0(P0), m_T0(T0), m_Tw(Tw) {}
        BeckwithGallagherSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, Trajectory traj, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid, traj), m_Geom(geom), m_Tw(Tw) { m_hasTrajectory = true; }
        ~BeckwithGallagherSolver() = default;
        virtual HFResult Solve(double time = 0.0) override;
        virtual HFResult Solve(double time, double Twall) override
        {
            m_Tw = Twall;
            m_currentTime = time;
            return Solve();
        }
    };

    class VanDriestSolver : public HFSolverBase
    {
    private:
        TimeSeries m_Pinf;
        TimeSeries m_Vinf;

        TimeSeries m_P0;
        TimeSeries m_T0;
        GeometryPrimitive m_Geom;
        double m_RunningLength;
        double m_Tw;
        TimeSeries m_Me;
        TimeSeries m_Pe;
        TimeSeries m_Te;
        double m_Rec;
        double m_IsCone;
        double m_currentTime;
        CpData m_CpData;
        double m_Loc;

    public:
        VanDriestSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double runningLength, TimeSeries pinf, TimeSeries vinf, TimeSeries P0, TimeSeries T0, TimeSeries Me, TimeSeries Te, TimeSeries Pe, double Tw = 300, double Rec = 0, double Is_Cone = 0)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_RunningLength(runningLength), m_Pinf(pinf), m_Vinf(vinf), m_P0(P0), m_T0(T0), m_Me(Me), m_Te(Te), m_Pe(Pe), m_Tw(Tw), m_Rec(Rec), m_IsCone(Is_Cone) {}
        VanDriestSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double runningLength, Trajectory traj, CpData cp, double loc, double Tw = 300, double Rec = 0, double Is_Cone = 0)
            : HFSolverBase(coordinateType, flowType, fluid, traj), m_Geom(geom), m_RunningLength(runningLength), m_CpData(cp), m_Loc(loc), m_Tw(Tw), m_Rec(Rec), m_IsCone(Is_Cone) { m_hasTrajectory = true; }
        ~VanDriestSolver() = default;
        virtual HFResult Solve(double time = 0.0) override;
        virtual HFResult Solve(double time, double Twall) override
        {
            m_Tw = Twall;
            m_currentTime = time;
            return Solve();
        }
    };

    class EckertSolver : public HFSolverBase
    {
    private:
        TimeSeries m_Pinf;
        TimeSeries m_Vinf;
        TimeSeries m_P0;
        TimeSeries m_T0;
        GeometryPrimitive m_Geom;
        double m_RunningLength;
        double m_Tw;
        TimeSeries m_Me;
        TimeSeries m_Pe;
        TimeSeries m_Te;
        double m_Rec;
        double m_IsCone;
        double m_currentTime;
        CpData m_CpData;
        double m_Loc;

    public:
        EckertSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double runningLength, TimeSeries pinf, TimeSeries vinf, TimeSeries P0, TimeSeries T0, TimeSeries Me, TimeSeries Te, TimeSeries Pe, double Tw = 300, double Rec = 0, double Is_Cone = 0)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_RunningLength(runningLength), m_Pinf(pinf), m_Vinf(vinf), m_P0(P0), m_T0(T0), m_Me(Me), m_Te(Te), m_Pe(Pe), m_Tw(Tw), m_Rec(Rec), m_IsCone(Is_Cone) {}
        EckertSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double runningLength, Trajectory traj, CpData cp, double loc, double Tw = 300, double Rec = 0, double Is_Cone = 0)
            : HFSolverBase(coordinateType, flowType, fluid, traj), m_Geom(geom), m_RunningLength(runningLength), m_CpData(cp), m_Loc(loc), m_Tw(Tw), m_Rec(Rec), m_IsCone(Is_Cone) { m_hasTrajectory = true; }
        ~EckertSolver() = default;
        virtual HFResult Solve(double time = 0.0) override;
        virtual HFResult Solve(double time, double Twall) override
        {
            m_Tw = Twall;
            m_currentTime = time;
            return Solve();
        }
    };

    class KnBridgedSolver : public HFSolverBase
    {
    private:
        double m_Tw = 300;
        double m_currentTime = 0;
        GeometryPrimitive m_geom;
        std::shared_ptr<HFSolverBase> m_continuumSolver;
        double m_characteristicLength;
        TimeSeries m_Pinffm;
        TimeSeries m_Tinffm;
        TimeSeries m_Vinffm;

    public:
        template <typename... Args>
        KnBridgedSolver(double characteristic_length, TimeSeries Pinf_fm, TimeSeries Tinf_fm, TimeSeries Vinf_fm, HFSolverType continuumSolverType, CoordinateType coordinateType, FlowType flowType, Fluid fluid, GeometryPrimitive geom, Args &&...args)
            : HFSolverBase(coordinateType, flowType, fluid), m_geom(geom), m_characteristicLength(characteristic_length), m_Pinffm(Pinf_fm), m_Tinffm(Tinf_fm), m_Vinffm(Vinf_fm)
        {
            m_continuumSolver = make_HFSolver(continuumSolverType, coordinateType, flowType, fluid, geom, std::forward<Args>(args)...);
            if (m_continuumSolver->GetTrajectory() != nullptr)
            {
                m_hasTrajectory = true;
                m_Trajectory = m_continuumSolver->GetTrajectory();
            }
        }

        ~KnBridgedSolver() = default;
        virtual HFResult Solve(double time = 0.0) override;
        virtual HFResult Solve(double time, double Twall) override
        {
            m_Tw = Twall;
            m_currentTime = time;
            return Solve();
        }
    };

    template <typename... Args>
    std::shared_ptr<HFSolverBase> make_HFSolver(HFSolverType T, Args &&...args)
    {
        switch (T)
        {
        case HFSolverType::FayRiddell:
            return safe_make_shared<FayRiddellSolver>(std::forward<Args>(args)...);
        case HFSolverType::BeckwithGallagher:
            return safe_make_shared<BeckwithGallagherSolver>(std::forward<Args>(args)...);
        case HFSolverType::Lees:
            return safe_make_shared<LeesSolver>(std::forward<Args>(args)...);
        case HFSolverType::VanDriest:
            return safe_make_shared<VanDriestSolver>(std::forward<Args>(args)...);
        case HFSolverType::Eckert:
            return safe_make_shared<EckertSolver>(std::forward<Args>(args)...);
        case HFSolverType::FreeMolecular:
            return safe_make_shared<FreeMolecularSolver>(std::forward<Args>(args)...);
        case HFSolverType::KnBridged:
            return safe_make_shared<KnBridgedSolver>(std::forward<Args>(args)...);
        default:
            assert(false);
            return nullptr;
        }
    }

    // Alias
    using HFSolver = std::shared_ptr<HFSolverBase>;

    // Response Solver
    struct ResponseSolverParams
    {
        // Initial Conditions
        double tInit = 0.0;
        double tFinal;
        // Table timestep;
        Table<double> timestep;
        double innerRadius;
        double massRateChar = 0.0;
        double massRatePyro = 0.0;
        double Sea_Level_Pressure = 101325;
        double Ttol = 1.0;
        double Air_Layer_Length = 1;
    };

    struct ResponseSolverOptions
    {
        CoordinateType coordinateType;
    };

    // FIXME: Convert to Table
    struct ResponseResult
    {
        std::vector<double> solution_t;
        TableGrid<double, TableGridOrder::ColumnFirst, TableGridOrder::RowFirst> solution_T;
        std::vector<size_t> interface_index;
        std::string message;

        void Save(const std::string &filepath, bool all = false, bool delim_whitespace = true) const
        {
            std::ofstream o(filepath);
            Save(o, all, delim_whitespace);
        }
        void Save(std::ostream &stream, bool all = false, bool delim_whitespace = true) const
        {
            TableGrid<double, TableGridOrder::RowFirst, TableGridOrder::RowFirst> temp;
            if (all)
            {
                temp = solution_T.data();
            }
            else
            {
                auto data = solution_T.data();
                std::vector<double> selected;
                for (size_t j = 0; j < solution_T.getColumnCount(); j++)
                {
                    for (auto &&i : interface_index)
                    {
                        selected.push_back(data[j][i]);
                    }
                    temp.addColumn(selected);
                    selected.clear();
                }
            }
            temp.insertRow(0, solution_t);
            temp.Save(stream, delim_whitespace);
        }
    };

    class ResponseSolverBase // FIXME: Enforce heat flux solver to have trajectory
    {
    protected:
        std::shared_ptr<LayerMesh> m_Layers;
        BCArray m_BCs;

    public:
        ResponseSolverBase(std::shared_ptr<LayerMesh> layers, TimeSeries Tinit, BCArray BCs)
            : m_Layers(layers), m_BCs(BCs)
        {
            m_Layers->InitTemperature(Tinit);
        };
        ResponseSolverBase(std::shared_ptr<LayerMesh> layers, double Tinit, BCArray BCs)
            : m_Layers(layers), m_BCs(BCs)
        {
            m_Layers->InitTemperature(Tinit);
        };

        virtual ~ResponseSolverBase() = default;
        virtual ResponseResult Solve() = 0;
    };

    class DefaultResponseSolver : public ResponseSolverBase
    {
        struct BCS
        {
            double flux_front = 0;
            double flux_back = 0;
            double propellant_mass_front = 0;
            double propellant_mass_back = 0;
            double qgen = 0;
            TimeSeries Tambient_front;
            TimeSeries Tambient_back;
            TimeSeries h_front;
            TimeSeries h_back;
            TimeSeries Tg_front;
            TimeSeries Tg_back;
            TimeSeries l_front;
            TimeSeries l_back;
            TimeSeries l_tg_front;
            TimeSeries l_tg_back;
        };

    private:
        FlowType m_FlowType;
        Fluid m_Fluid;
        HFSolver m_HFSolver;
        size_t LOPT;
        size_t MOPT;
        double Tcon;

        ResponseSolverParams m_params;
        ResponseSolverOptions m_options;

        void CYSP2(double &A1, double &B1, double &C1, double k, double Inverse_ratio);
        void CYSP1(double &A1, double &C1, double k, double Layer_thickness, double Inverse_ratio, double Outer_radius_instantaneous);
        double FreeConvection_Surface(double T, double T_AMBIENT, double Characteristic_Length, double Sea_Level_Pressure);
        void Instantaneous_MassRate_ThicknessSolver(double &Layer_thickness_1, double &Layer_thickness_2, double &Mass_rate_Char, double &Mass_rate_Pyrolysis, double Q_convective_front, double Qradiation_front, Eigen::VectorXd T, double deltat, const BCS &bcs);

        void Air_Gap(double &h_convection, double &h_radiation, double AIR_GAP, double LENGTH, double Temp_1, double Temp_2, double Sea_Level_pressure, double Emmisivity_1, double Emmisivity_2);

        void DefaultResponseMatrix(Eigen::VectorXd &T_current_step, Eigen::VectorXd Layer_thickness, Eigen::VectorXd Initial_Layer_thickness, double, double, BCS bcs, double &Mass_rate_Char, double &Mass_rate_Pyrolysis, Eigen::VectorXd k, double delt);

    public:
        // template <typename T>
        DefaultResponseSolver(std::shared_ptr<LayerMesh> layers, HFSolver hfsolver, TimeSeries Tinit, BCArray BCs, ResponseSolverParams params, ResponseSolverOptions options, size_t LOPT = 0, size_t MOPT = 1, double Tcon = 300)
            : ResponseSolverBase(layers, Tinit, BCs), m_params(params), m_options(options), m_HFSolver(hfsolver), LOPT(LOPT), MOPT(MOPT), Tcon(Tcon) { m_Fluid = m_HFSolver->GetFluid(); }
        virtual ~DefaultResponseSolver() = default;
        virtual ResponseResult Solve() override;
    };

    template <typename... Args>
    std::shared_ptr<ResponseSolverBase> make_ResponseSolver(ResponseSolverType T, Args &&...args)
    {
        switch (T)
        {
        case ResponseSolverType::Default:
            return safe_make_shared<DefaultResponseSolver>(std::forward<Args>(args)...);
        default:
            assert(false);
            return nullptr;
        }
    }

    // alias
    using ResponseSolver = std::shared_ptr<ResponseSolverBase>;

} // namespace OTAP