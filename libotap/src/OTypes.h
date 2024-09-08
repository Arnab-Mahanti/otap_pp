#ifndef _OTYPES_H
#define _OTYPES_H

#include <iostream>
#include <iterator>
#include <valarray>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <memory>
#include <functional>
#include "mlinterp/mlinterp.hpp"
#include "magic_enum/magic_enum.h"
#include "rapidcsv/rapidcsv.h"

#define DEBUG

#ifdef DEBUG

#define O_WARN(x) (std::cout << "[Warn] " << (x) << "\n")

#else

#define O_WARN(x)

#endif

namespace OTAP
{
    using TimePoints = std::vector<double>;
    using TimeSeries = std::vector<double>;

    enum class GeometryPrimitiveType
    {
        StagnationPoint,
        Line,
        Arc
    };

    enum class ProblemType
    {
        HeatFlux,
        Response,
        Optimization
    };

    // NOTE: Define SolverType corresponding to each ProblemType

    enum class HFSolverType
    {
        FayRiddell,
        VanDriest,
        BeckwithGallagher,
        Eckert,
        Lees,
        FreeMolecular,
        LeewardHeating
    };

    enum class ResponseSolverType
    {
        Default
    };

    enum class OptimizationSolverType
    {

    };

    enum class CoordinateType
    {
        Cartesian,
        Axisym,
        Spherical
    };

    enum class FlowType
    {
        Default,
        Inviscid,
        Laminar,
        Turbulent
    };

    enum class FluidType
    {
        Perfect,
        Hansen_Air,
    };

    enum class AmbientType
    {
        ISA
    };

    enum class TrajectoryType
    {
        Velocity,
        // Mach,
        // Flux,
        // Wind_Tunnel,
        // CFD_Data,
        // New_traj
    };

    enum class SourceType
    {
        RadFlux
    };

    // Helper functions

    template <typename T>
    constexpr std::string GetTypes()
    {
        auto types = magic_enum::enum_names<T>();
        std::string s = "";
        for (auto &&i : types)
        {
            s.append(i);
            s += '\0';
        }

        return s;
    }

    template <typename T>
    auto GetTypeAtIndex(int i)
    {
        return magic_enum::enum_values<T>()[i];
    }

    template <typename T>
    auto GetNameFromType(T t)
    {
        return std::string(magic_enum::enum_name<T>(t));
    }

    class CpData
    {
        std::vector<double> m_Mach;
        std::vector<double> m_Positions;
        std::vector<double> m_Cp;

    public:
        CpData(const std::vector<double> &Mach, const std::vector<double> &Pos, const std::vector<double> &Cp)
            : m_Mach(Mach), m_Positions(Pos), m_Cp(Cp) {}
        CpData(const std::string &filename, const bool delim_whitespace = true);
        CpData(const CpData &) = default; // copy constructor
        explicit CpData() = default;      // Explicit default constructor

        double GetCp(double Mach, double Pos) const;
        std::vector<double> GetCp(const std::vector<double> &Mach, const std::vector<double> &Pos) const;
    };

    template <typename T, typename... Args>
    std::shared_ptr<T> safe_make_shared(Args &&...args)
    {
        if constexpr (std::is_constructible_v<T, Args &&...>)
            return std::make_shared<T>(std::forward<Args>(args)...);
        else
            return nullptr;
    }

} // namespace OTAP

#endif // _OTYPES_H