#include "Ambient.h"
#include "mlinterp/mlinterp.hpp"
#include <valarray>
namespace OTAP
{
    std::vector<double> ISA_Ambient::Pressure(std::vector<double> altitude) const
    {
        std::vector<double> pressure(altitude.size(), 0.0);
        std::for_each(altitude.begin(), altitude.end(), [](auto &&a)
                      { a /= 1000.0; });
        size_t dim[] = {m_altitudeData.size()};
        mlinterp::interp(dim, altitude.size(), m_pressureData.data(), pressure.data(), m_altitudeData.data(), altitude.data());
        return pressure;
    }
    std::vector<double> ISA_Ambient::Temperature(std::vector<double> altitude) const
    {
        std::vector<double> temperature(altitude.size(), 0.0);
        std::for_each(altitude.begin(), altitude.end(), [](auto &&a)
                      { a /= 1000.0; });
        size_t dim[] = {m_altitudeData.size()};
        mlinterp::interp(dim, altitude.size(), m_temperatureData.data(), temperature.data(), m_altitudeData.data(), altitude.data());
        return temperature;
    }
    std::vector<double> ISA_Ambient::Density(std::vector<double> altitude) const
    {
        std::vector<double> density(altitude.size(), 0.0);
        std::for_each(altitude.begin(), altitude.end(), [](auto &&a)
                      { a /= 1000.0; });
        size_t dim[] = {m_altitudeData.size()};
        mlinterp::interp(dim, altitude.size(), m_densityData.data(), density.data(), m_altitudeData.data(), altitude.data());
        return density;
    }
    // std::valarray<double> ISA_Ambient::Pressure(std::valarray<double> altitude) const
    // {
    //     vector _alt(std::begin(altitude), std::end(altitude));
    //     auto _p = Pressure(_alt);
    //     return std::valarray<double>(_p.data(), _p.size());
    // }
    // std::valarray<double> ISA_Ambient::Temperature(std::valarray<double> altitude) const
    // {
    //     vector _alt(std::begin(altitude), std::end(altitude));
    //     auto _t = Temperature(_alt);
    //     return std::valarray<double>(_t.data(), _t.size());
    // }
    // std::valarray<double> ISA_Ambient::Density(std::valarray<double> altitude) const
    // {
    //     vector _alt(std::begin(altitude), std::end(altitude));
    //     auto _d = Density(_alt);
    //     return std::valarray<double>(_d.data(), _d.size());
    // }
    double ISA_Ambient::Pressure(double altitude) const
    {
        double pressure;
        altitude /= 1000;
        size_t dim[] = {m_altitudeData.size()};
        mlinterp::interp(dim, size_t(1), m_pressureData.data(), &pressure, m_altitudeData.data(), &altitude);
        return pressure;
    }
    double ISA_Ambient::Temperature(double altitude) const
    {
        double temperature;
        altitude /= 1000;
        size_t dim[] = {m_altitudeData.size()};
        mlinterp::interp(dim, size_t(1), m_temperatureData.data(), &temperature, m_altitudeData.data(), &altitude);
        return temperature;
    }
    double ISA_Ambient::Density(double altitude) const
    {
        double density;
        altitude /= 1000;
        size_t dim[] = {m_altitudeData.size()};
        mlinterp::interp(dim, size_t(1), m_densityData.data(), &density, m_altitudeData.data(), &altitude);
        return density;
    }

}
