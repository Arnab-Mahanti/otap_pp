#pragma once

#include "OTypes.h"
#include "Fluid.h"
#include "Ambient.h"
#include "Geometry.h"
#include "MaterialLayer.h"
#include "Trajectory.h"
#include <memory>

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

    public:
        HFSolverBase(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid)
            : m_CoordinateType(coordinateType), m_FlowType(flowType), m_Fluid(fluid) {}
        virtual ~HFSolverBase() = default;
        virtual HFResult Solve(double time = 0.0) = 0;
        virtual HFResult Solve(double time, double Twall) = 0;
    };

    class FayRiddellSolver : public HFSolverBase
    {
    private:
        TimeSeries m_Pinf;
        TimeSeries m_Vinf;

        TimeSeries m_P0;
        TimeSeries m_T0;
        GeometryPrimitive m_Geom;

        Trajectory m_Trajectory;
        bool m_hasTrajectory = false;
        double m_Tw;

        double m_currentTime = 0.0;

    public:
        FayRiddellSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, TimeSeries pinf, TimeSeries vinf, TimeSeries P0, TimeSeries T0, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_Pinf(pinf), m_Vinf(vinf), m_P0(P0), m_T0(T0), m_Tw(Tw) {}
       FayRiddellSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, Trajectory traj, double Tw = 300)
           : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_Trajectory(traj), m_Tw(Tw)
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
        Trajectory m_Trajectory;
        double m_Tw;
        double m_RunningLength;
        double m_currentTime;
        bool m_hasTrajectory = false;

    public:
        // FIXME: Calculate qstag from Fay & Riddelle
        LeesSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double RunningLength, TimeSeries pinf, TimeSeries Tinf, TimeSeries vinf, TimeSeries minf, TimeSeries P0, TimeSeries T0, TimeSeries qstag, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_RunningLength(RunningLength), m_Pinf(pinf), m_Tinf(Tinf), m_Vinf(vinf), m_Minf(minf), m_P0(P0), m_T0(T0), m_qstag(qstag), m_Tw(Tw) {}
        LeesSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double RunningLength, Trajectory traj, TimeSeries qstag, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_RunningLength(RunningLength), m_Trajectory(traj), m_qstag(qstag), m_Tw(Tw) { m_hasTrajectory = true; }
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
        Trajectory m_Trajectory;
        bool m_hasTrajectory = false;

        double m_Tw;
        double m_currentTime;

    public:
        FreeMolecularSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, TimeSeries pinf, TimeSeries vinf, TimeSeries Tinf, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_Pinf(pinf), m_Vinf(vinf), m_Tinf(Tinf), m_Tw(Tw) {}
        FreeMolecularSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, Trajectory traj, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_Trajectory(traj), m_Tw(Tw) {m_hasTrajectory = true;}
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
        Trajectory m_Trajectory;
        bool m_hasTrajectory = false;


    public:
        BeckwithGallagherSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, TimeSeries pinf, TimeSeries Tinf, TimeSeries vinf, TimeSeries Minf, TimeSeries P0, TimeSeries T0, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_Pinf(pinf), m_Tinf(Tinf), m_Vinf(vinf), m_Minf(Minf), m_P0(P0), m_T0(T0), m_Tw(Tw) {}
       BeckwithGallagherSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, Trajectory traj, double Tw = 300)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_Trajectory(traj), m_Tw(Tw) {m_hasTrajectory = true;}
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
        Trajectory m_Trajectory;
        CpData m_CpData;
        bool m_hasTrajectory = false;
        double m_Loc;

    public:
        VanDriestSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double runningLength, TimeSeries pinf, TimeSeries vinf, TimeSeries P0, TimeSeries T0, TimeSeries Me, TimeSeries Te, TimeSeries Pe, double Tw = 300, double Rec = 0, double Is_Cone = 0)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_RunningLength(runningLength), m_Pinf(pinf), m_Vinf(vinf), m_P0(P0), m_T0(T0), m_Me(Me), m_Te(Te), m_Pe(Pe), m_Tw(Tw), m_Rec(Rec), m_IsCone(Is_Cone) {}
        VanDriestSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double runningLength, Trajectory traj, CpData cp, double loc, double Tw = 300, double Rec = 0, double Is_Cone = 0)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_RunningLength(runningLength), m_Trajectory(traj), m_CpData(cp),m_Loc(loc), m_Tw(Tw), m_Rec(Rec), m_IsCone(Is_Cone) {m_hasTrajectory = true;}
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
        Trajectory m_Trajectory;
        CpData m_CpData;
        bool m_hasTrajectory = false;
        double m_Loc;

    public:
        EckertSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double runningLength, TimeSeries pinf, TimeSeries vinf, TimeSeries P0, TimeSeries T0, TimeSeries Me, TimeSeries Te, TimeSeries Pe, double Tw = 300, double Rec = 0, double Is_Cone = 0)
            : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_RunningLength(runningLength), m_Pinf(pinf), m_Vinf(vinf), m_P0(P0), m_T0(T0), m_Me(Me), m_Te(Te), m_Pe(Pe), m_Tw(Tw), m_Rec(Rec), m_IsCone(Is_Cone) {}
       // EckertSolver(const CoordinateType coordinateType, const FlowType flowType, Fluid fluid, GeometryPrimitive geom, double runningLength, Trajectory traj, CpData cp,double loc, double Tw = 300, double Rec = 0, double Is_Cone = 0)
       //     : HFSolverBase(coordinateType, flowType, fluid), m_Geom(geom), m_RunningLength(runningLength), m_Trajectory(traj), m_CpData(cp),m_Loc(loc), m_Tw(Tw), m_Rec(Rec), m_IsCone(Is_Cone) {m_hasTrajectory = true;}
        ~EckertSolver() = default;
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
        if (T == HFSolverType::FayRiddell)
            return std::make_shared<FayRiddellSolver>(std::forward<Args>(args)...);
        else if (T == HFSolverType::BeckwithGallagher)
            return std::make_shared<BeckwithGallagherSolver>(std::forward<Args>(args)...);
        else if (T == HFSolverType::Lees)
            return std::make_shared<LeesSolver>(std::forward<Args>(args)...);
        else if (T == HFSolverType::VanDriest)
            return std::make_shared<VanDriestSolver>(std::forward<Args>(args)...);
        else if (T == HFSolverType::Eckert)
            return std::make_shared<EckertSolver>(std::forward<Args>(args)...);
        else if (T == HFSolverType::FreeMolecular)
            return std::make_shared<FreeMolecularSolver>(std::forward<Args>(args)...);
        else
            assert(false);
    }

    // Alias
    using HFSolver = std::shared_ptr<HFSolverBase>;

    // Response Solver
    // FIXME: Create a proper Table implementation with inbuilt interpolation
    using Table = std::vector<std::pair<double, double>>;
    struct ResponseSolverParams
    {
        // Initial Conditions
        double tInit = 0.0;
        double tFinal;
        Table timestep;
        double innerRadius;
        double massRateChar = 0.0;
        double massRatePyro = 0.0;
        Table propellantMass;

        // Frontwall radiation
        Table qRad;
        double TRad = 300.0;
        // Table qGen;

        // Backwall radiation and convection
        double hBackwall = 0.0; // FIXME: Give Options for h and characteristic length for natural convection
        double TgBackwall = 300.0;
        double TRadBackwall = 300.0;
    };

    struct ResponseSolverOptions
    {
        CoordinateType coordinateType;

        double Ttol = 1.0;
    };

    struct ResponseResult
    {
        TimeSeries T;
    };

    class ResponseSolverBase
    {
    protected:
        std::shared_ptr<LayerStack> m_Layers;

    public:
        // FIXME: Check conversion of pointer to shared_ptr
        ResponseSolverBase(LayerStack *layers, TimeSeries Tinit)
            : m_Layers(layers)
        {
            m_Layers->InitTemperature(Tinit);
        };
        ResponseSolverBase(LayerStack *layers,  double Tinit)
            : m_Layers(layers)
        {
            m_Layers->InitTemperature(Tinit);
        };

        virtual ~ResponseSolverBase() = default;
        virtual ResponseResult Solve() = 0;
    };

    class DefaultResponseSolver : public ResponseSolverBase
    {
    public:
        HFSolver m_ContinuumSolver;
        HFSolver m_FreeMolecularSolver;

    public:
      ResponseSolverParams m_params;
      ResponseSolverOptions m_options;


//FOR CYLINDRICAL/SPHERICAL COORDINATE TRANSFORMATION (RADIUS != 0)
void CYSP2(double& A1, double& B1, double& C1, double k, double Inverse_ratio, struct ResponseSolverOptions options)
{
    double Coordinate_ratio,COORDINATE_COEF;
    if(options.coordinateType==CoordinateType::Cartesian)
    {
        COORDINATE_COEF=0;
    }
    else if(options.coordinateType==CoordinateType::Cylindrical)
    {
        COORDINATE_COEF=1;
    }
    else if(options.coordinateType==CoordinateType::Spherical)
    {
        COORDINATE_COEF=2;
    }
    else
    {
        COORDINATE_COEF=0;
    }

	
	Coordinate_ratio = (COORDINATE_COEF * k) / pow(Inverse_ratio, 2);
	A1 = A1 + Coordinate_ratio;
	B1 = B1 - 2 * Coordinate_ratio;
	C1 = C1 + Coordinate_ratio;
}

//FOR CYLINDRICAL/SPHERICAL COORDINATE TRANSFORMATION (RADIUS == 0)
void CYSP1(double& A1, double& C1, double k, double Layer_thickness, double Inverse_ratio, double Outer_radius_instantaneous, struct ResponseSolverOptions options)
{
    double Coordinate_ratio,COORDINATE_COEF;
        if(options.coordinateType==CoordinateType::Cartesian)
    {
        COORDINATE_COEF=0;
    }
    else if(options.coordinateType==CoordinateType::Cylindrical)
    {
        COORDINATE_COEF=1;
    }
    else if(options.coordinateType==CoordinateType::Spherical)
    {
        COORDINATE_COEF=2;
    }
    else
    {
        COORDINATE_COEF=0;
    }
	
	Coordinate_ratio = (COORDINATE_COEF  * k * Layer_thickness) / (2 * Inverse_ratio * Outer_radius_instantaneous);
	A1 = A1 + Coordinate_ratio;
	C1 = C1 + Coordinate_ratio;
}



        template <typename T>
        DefaultResponseSolver(LayerStack *layers, HFSolver continuumSolver, HFSolver freeMolSolver, T Tinit, ResponseSolverParams params, ResponseSolverOptions options)
            : ResponseSolverBase(layers, Tinit), m_params(params), m_options(options){}
        virtual ~DefaultResponseSolver() = default;
        virtual ResponseResult Solve() override;
    };

    // alias
    using ResponseSolver = std::shared_ptr<ResponseSolverBase>;

} // namespace OTAP