#include "Problem.h"
#include "mlinterp/mlinterp.hpp"
#include <sstream>
#include <iostream>

// OTAP::Problem::Problem(const ProblemSpecification &spec, const Geometry &geometry, const TrajectoryBase &trajectory)
// {
//     return 0;
// }

namespace OTAP
{
    // helpers
    std::optional<size_t> ResponseProblem::getGeomIdFromLoc(double loc)
    {
        for (size_t i = 0; i < m_geom.Count(); i++)
        {
            if (m_geom.GetAxialLength(i) >= loc)
                return i;
        }
        return {};
    }

    void OTAP::HFProblem::Solve()
    {
        // Dispatches solvers for each geometry components based of nodal discreetization.
        // Each solver call gives history at that location.
    }

    ResponseSolution OTAP::ResponseProblem::Solve()
    {
        ResponseSolution solution;

        for (auto &&loc : m_locations)
        {
            auto gid = getGeomIdFromLoc(loc);
            assert(gid.has_value());
            auto component = m_geom.GetComponents()[gid.value()];
            HFSolver hf = safe_make_shared<HFSolverBase>();
            if (m_trajectory != nullptr)
            {
                switch (component.type)
                {
                case GeometryPrimitiveType::StagnationPoint:
                    if (component.lambda != 0 && m_spec.fluxCSYS != CoordinateType::Axisym)
                    {
                        assert(false); // beckweth
                    }
                    else
                    {
                        hf = make_HFSolver(HFSolverType::KnBridged,
                                           m_geom.GetComponents()[1].radius == 0 ? m_geom.GetComponents()[1].radius : m_geom.GetComponents()[1].radius,
                                           HFSolverType::FayRiddell,
                                           m_spec.fluxCSYS,
                                           m_spec.flowType,
                                           make_fluid(m_spec.fluidType),
                                           component,
                                           m_trajectory,
                                           Tinit);
                    }
                    break;
                case GeometryPrimitiveType::Arc:
                    if (component.length == 0)
                    {
                        // lees
                    }
                    else
                    {
                        // flat plate
                    }
                    break;
                case GeometryPrimitiveType::Line:
                    // hf = make_HFSolver(HFSolverType::KnBridged,
                    //    m_geom.GetComponents()[1].radius == 0 ? m_geom.GetComponents()[1].radius : m_geom.GetComponents()[1].radius,
                    //    flatplatesolverType,
                    //    m_spec.fluxCSYS,
                    //    m_spec.flowType,
                    //    make_fluid(m_spec.fluidType),
                    //    component,
                    //    m_trajectory,
                    //    Tinit);
                    break;
                default:
                    break;
                }
            }

            ResponseSolverParams params;
            params.tInit = tInit;
            params.tFinal = tFinal;
            params.Ttol = Ttol;
            params.timestep = timestep;
            params.Sea_Level_Pressure = Sea_Level_Pressure;
            params.massRatePyro = massRatePyro;
            params.massRateChar = massRateChar;
            params.Air_Layer_Length = Air_Layer_Length;
            // params.innerradius ===== //FIXME:

            ResponseSolverOptions options;
            options.coordinateType = m_spec.responseCSYS;

            auto mesh = std::make_shared<LayerMesh>(component.layers);

            auto response = make_ResponseSolver(ResponseSolverType::Default,
                                                mesh,
                                                hf,
                                                std::vector<double>(component.layers.GetCount() + 1, Tinit),
                                                bcs,
                                                params,
                                                options);
            auto result = response->Solve();
            solution.locations.push_back(loc);
            solution.solution.push_back(result);
        }
        return solution;
    }
}
