#include "Trajectory.h"
#include <map>
#include <fstream>
#include <iomanip>
#include "mlinterp/mlinterp.hpp"
#include "Fluid.h"
namespace OTAP
{
    // FIXME: Make parser robust
    void VelocityTrajectory::ReadTrajectory(const std::string &filename)
    {
        std::ifstream file(filename);
        std::map<double, std::vector<double>> data;
        std::vector<double> line(6, 0.0);  // 6?
        std::string t;
        std::getline(file, t); // skip first line
        while (!file.eof())
        {
            // time
            file >> t;
            std::string s;
            if (t == "")
                break;

            for (size_t i = 0; i < 6; i++)
            {
                file >> s;
                line[i] = std::stod(s);
            }
            data[std::stod(t)] = std::vector(line.begin(), line.begin() + 3);
            t = "";
        }
        for (auto &&i : data)
        {
            m_timepoints.push_back(i.first);
            m_altitude.push_back(i.second[0]);
            m_velocity.push_back(i.second[1]);
        }
    }

    VelocityTrajectory::VelocityTrajectory(const AmbientType &ambientType, const FluidType &fluidType, const std::string &filename)
        : m_ambient(make_ambient(ambientType)), m_fluid(make_fluid(fluidType))
    {
        ReadTrajectory(filename);
    }

    double VelocityTrajectory::GetFreeStream() const
    {
        return 0.0;
    }

    TimeSeries VelocityTrajectory::GetPinf(TimePoints t) const
    {
        std::vector<double> altitude(t.size(), 0.0);
        size_t dim[] = {m_timepoints.size()};
        mlinterp::interp(dim, t.size(), m_altitude.data(), altitude.data(), m_timepoints.data(), t.data());
        return m_ambient->Pressure(altitude);
    }

    double VelocityTrajectory::GetPinf(double t) const
    {
        double altitude;
        size_t dim[] = {m_timepoints.size()};
        mlinterp::interp(dim, size_t(1), m_altitude.data(), &altitude, m_timepoints.data(), &t);
        return m_ambient->Pressure(altitude);
    }

    TimeSeries VelocityTrajectory::GetTinf(TimePoints t) const
    {
        std::vector<double> altitude(t.size(), 0.0);
        size_t dim[] = {m_timepoints.size()};
        mlinterp::interp(dim, t.size(), m_altitude.data(), altitude.data(), m_timepoints.data(), t.data());
        return m_ambient->Temperature(altitude);
    }

    double VelocityTrajectory::GetTinf(double t) const
    {
        double altitude;
        size_t dim[] = {m_timepoints.size()};
        mlinterp::interp(dim, size_t(1), m_altitude.data(), &altitude, m_timepoints.data(), &t);
        return m_ambient->Temperature(altitude);
    }

    TimeSeries VelocityTrajectory::GetRhoinf(TimePoints t) const
    {
        std::vector<double> altitude(t.size(), 0.0);
        size_t dim[] = {m_timepoints.size()};
        mlinterp::interp(dim, t.size(), m_altitude.data(), altitude.data(), m_timepoints.data(), t.data());
        return m_ambient->Density(altitude);
    }

    double VelocityTrajectory::GetRhoinf(double t) const
    {
        double altitude;
        size_t dim[] = {m_timepoints.size()};
        mlinterp::interp(dim, size_t(1), m_altitude.data(), &altitude, m_timepoints.data(), &t);
        return m_ambient->Density(altitude);
    }

    TimeSeries VelocityTrajectory::GetMinf(TimePoints t) const
    {
        std::vector<double> mach(t.size(), 0.0);
        auto velocity = GetVinf(t);
        auto press = GetPinf(t);
        auto temp = GetTinf(t);
        auto a = m_fluid->a(press, temp);
        std::transform(velocity.begin(), velocity.end(), a.begin(), mach.begin(), std::divides<double>());
        return mach;
    }

    double VelocityTrajectory::GetMinf(double t) const
    {
        return GetVinf(t) / m_fluid->a(GetPinf(t), GetTinf(t));
    }

    TimeSeries VelocityTrajectory::GetVinf(TimePoints t) const
    {
        std::vector<double> velocity(t.size(), 0.0);
        size_t dim[] = {m_timepoints.size()};
        mlinterp::interp(dim, t.size(), m_velocity.data(), velocity.data(), m_timepoints.data(), t.data());
        return velocity;
    }

    TimePoints VelocityTrajectory::GetTimePoints() const
    {
        return m_timepoints;
    }

    double VelocityTrajectory::GetVinf(double t) const
    {
        double velocity;
        size_t dim[] = {m_timepoints.size()};
        mlinterp::interp(dim, size_t(1), m_velocity.data(), &velocity, m_timepoints.data(), &t);
        return velocity;
    }

}