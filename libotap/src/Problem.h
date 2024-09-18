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

    // FIXME: How to link to ProblemBase
    // class ResponseProblem
    // {
    // private:
    //     ProblemSpecification m_spec;
    //     Geometry m_geom;
    //     std::vector<double> m_locations;
    //     std::unordered_map<double, std::shared_ptr<LayerMesh>> m_analysis;

    // public:
    //     ResponseProblem(const ProblemSpecification &spec, const Geometry &geom, Trajectory traj, BCArray bcs, TimeSeries?? Tinit, size_t numlocations)
    //         : m_spec(spec), m_geom(geom) {}
    //     ResponseProblem(const ProblemSpecification &spec, const Geometry &geom, Trajectory traj, BCArray bcs, TimeSeries?? Tinit, std::vector<double> locations, )
    //         : m_spec(spec), m_geom(geom) {}
    //     virtual ~ResponseProblem() = default;
    //     virtual bool Validate() const;
    //     virtual ResponseSolution Solve();
    //     virtual ResponseSolution Solve(double location, TimeSeries Tinit);
    // };
} // namespace OTAP

#endif // _PROBLEM_H