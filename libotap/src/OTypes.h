#ifndef _OTYPES_H
#define _OTYPES_H

#include <iostream>
#include <iterator>
#include <valarray>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <functional>
#include "magic_enum/magic_enum.h"
#include "mlinterp/mlinterp.hpp"

#define DEBUG

#ifdef DEBUG

#define O_WARN(x) (std::cout << "[Warn] " << (x) << "\n")

#else

#define O_WARN(x)

#endif

namespace OTAP
{
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

    };

    enum class OptimizationSolverType
    {

    };

    enum class CoordinateType
    {
        Cartesian,
        Cylindrical,
        Spherical,
        Axisym
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

    // FIXME: Test this class for issues
    class vector : public std::vector<double>
    {
    private:
        /* data */
    public:
        vector() : std::vector<double>() {}
        vector(std::size_t s) : std::vector<double>(s) {}
        vector(std::size_t s, double val) : std::vector<double>(s, val) {}
        vector(std::initializer_list<double> v) : std::vector<double>(v) {}
        vector(const vector &_x) = default;
        vector(const std::vector<double> &v) : std::vector<double>(v) {}
        vector(vector &&) noexcept = default;
        ~vector() = default;
        operator std::valarray<double>()
        {
            return std::valarray<double>(this->data(), this->size());
        }
    };

    class CpData
    {
        std::vector<double> m_Mach;
        std::vector<double> m_Cp;
        std::vector<double> m_Positions;

    public:
        CpData(const std::vector<double> &Mach, const std::vector<double> &Pos, const std::vector<double> &Cp)
            : m_Mach(Mach), m_Positions(Pos), m_Cp(Cp) {}
        CpData(const std::string &filename);
        CpData(const CpData &) = default; // copy constructor
        explicit CpData() = default;      // Explicit default constructor

        inline double GetCp(double Mach, double Pos) const;
        inline std::vector<double> GetCp(const std::vector<double> &Mach, const std::vector<double> &Pos) const;
    };

    // FIXME: Make parser robust
    inline OTAP::CpData::CpData(const std::string &filename)
    {
        std::ifstream file(filename);
        std::string m = "";
        std::string pos = "";
        std::string line = "";
        std::string cp = "";
        std::getline(file, line);

        // Parse Mach nos.
        std::istringstream ss(line);
        while (!ss.eof())
        {
            ss >> m;
            m_Mach.push_back(std::stod(m));
            m = "";
        }

        for (auto &&i : m_Mach)
        {
            //   std::cout << i << " ";
        }
        //  std::cout << "\n";

        // Parse Cp data
        while (!file.eof())
        {
            std::getline(file, line);
            std::istringstream ss(line);
            ss >> pos;
            m_Positions.push_back(std::stod(pos));
            pos = "";
            while (!ss.eof())
            {
                ss >> m;
                m_Cp.push_back(std::stod(m));
                m = "";
            }
        }
    }

    inline double OTAP::CpData::GetCp(double Mach, double Pos) const
    {
        double cp = 0.0;
        size_t dim[] = {m_Mach.size(), m_Positions.size()};
        mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_Cp.data(), &cp, m_Mach.data(), &Mach, m_Positions.data(), &Pos);
        return cp;
    }

    std::vector<double> OTAP::CpData::GetCp(const std::vector<double> &Mach, const std::vector<double> &Pos) const
    {
        assert(Mach.size() == Pos.size());
        std::vector<double> cp(Mach.size(), 0.0);

        size_t dim[] = {m_Mach.size(), m_Positions.size()};
        mlinterp::interp<mlinterp::rnatord>(dim, cp.size(), m_Cp.data(), cp.data(), m_Mach.data(), Mach.data(), m_Positions.data(), Pos.data());

        return cp;
    }

    // std::function<double(double)> OTAP::CpData::GetCpVsM(double pos)
    // {
    //     // size_t dim[] = {m_Mach.size(), m_Positions.size()};
    //     // mlinterp::interp<mlinterp::rnatord>(dim, cp.size(), m_Cp.data(), cp.data(), m_Mach.data(), Mach.data(), m_Positions.data(), Pos.data());
    // }

    using TimePoints = std::vector<double>;
    using TimeSeries = std::vector<double>;
} // namespace OTAP

#endif // _OTYPES_H