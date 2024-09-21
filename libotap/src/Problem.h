#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "OTypes.h"
#include "Geometry.h"
#include "Trajectory.h"
#include "Solver.h"
#include <memory>
#include <fstream>

namespace OTAP
{

    struct ProblemSpecification
    {
        ProblemType problemType;
        CoordinateType fluxCSYS;
        CoordinateType responseCSYS;
        FlowType flowType;
        FluidType fluidType;
    };

    class ProblemBase
    {
    private:
        ProblemSpecification m_Spec;
        Geometry m_Geometry;
        Trajectory m_Trajectory;
        CpData m_CpData;

    public:
        ProblemBase(const ProblemSpecification &spec, const Geometry &geometry, Trajectory trajectory, const CpData &cp_data)
            : m_Spec(spec), m_Geometry(geometry), m_Trajectory(trajectory), m_CpData(cp_data) {}
        virtual ~ProblemBase() = default;
        virtual void Validate() const { assert(false); };
        virtual void Solve() { assert(false); };
        virtual void Solve(double location) { assert(false); }
    };

    // Returns HF for all loc for all time
    class HFProblem : public ProblemBase
    {
    private:
        std::shared_ptr<HFSolverBase> m_solver; // shared_ptr to some unknown HF (chosen at runtime) solver

    public:
        HFProblem(const ProblemSpecification &spec, const Geometry &geometry, Trajectory trajectory, const CpData &cp_data)
            : ProblemBase(spec, geometry, trajectory, cp_data) {}
        ~HFProblem() = default;
        virtual void Solve() override;
    };

    struct ResponseSolution
    {
        std::vector<double> locations;
        std::vector<ResponseResult> solution;
    };

    // TODO: Classes of Solver can be made base and inherited from

    template <typename T>
    struct GetSet
    {
        const T &get() const { return m_data; }
        void set(const T &val) { m_data = val; }
        GetSet() = default;
        GetSet(const T &val) : m_data(val) {}

    private:
        T m_data;
    };

    // FIXME: How to link to ProblemBase
    class ResponseProblem
    {
    private:
        ProblemSpecification m_spec;
        Geometry m_geom;
        std::vector<double> m_locations;
        Trajectory m_trajectory;
        std::optional<CpData> m_cpdata;
        BCArray bcs;
        // TableGrid<double,TableGridOrder::ColumnFirst,TableGridOrder::RowFirst> Tinit;

        // Params

        // public:
        // Params
        double Tinit;
        double tInit = 0.0;
        double h_scale = 1;
        double h_bias = 0;
        double tFinal;
        Table<double> timestep;
        double massRateChar = 0.0;
        double massRatePyro = 0.0;
        double Sea_Level_Pressure = 101325;
        double Ttol = 1.0;
        double Air_Layer_Length = 1;
        std::optional<double> ReC;
        std::optional<double> TransitionLength;
        HFSolverType flatplatesolverType = HFSolverType::VanDriest;
        HFSolverType stagnationSolverType = HFSolverType::FayRiddell;
        // TODO: Transitional criteria

        // helpers
        std::optional<size_t> getGeomIdFromLoc(double loc);

        explicit ResponseProblem() = default;

    public:
        ResponseProblem(const ProblemSpecification &spec, const Geometry &geom, Trajectory traj, BCArray bcs, double Tinit, size_t numlocations)
            : m_spec(spec), m_geom(geom) {}
        ResponseProblem(const ProblemSpecification &spec, const Geometry &geom, Trajectory traj, BCArray bcs, double Tinit, std::vector<double> locations)
            : m_spec(spec), m_geom(geom) {}
        virtual ~ResponseProblem() = default;
        // virtual bool Validate() const;
        virtual ResponseSolution Solve();
        virtual ResponseSolution Solve(double location, TimeSeries Tinit);
    };
} // namespace OTAP

#endif // _PROBLEM_H