#include "Fluid.h"
#include <algorithm>
#include <cmath>
#include <cassert>
#include "mlinterp/mlinterp.hpp"

// NOTE: execution policy par ?
// TODO: implement in terms of valarray and boost/transform_if or alternatively c++20 range::filter
std::vector<double> OTAP::Hansen::Rho(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());

    std::vector<double> rho(P.size(), 0.0);
    auto Z_value = Z(P, T);

    std::transform(P.begin(), P.end(), T.begin(), rho.begin(), [](double p, double t)
                   { return p / (287 * t); });
    std::transform(rho.begin(), rho.end(), Z_value.begin(), rho.begin(), [](double rho, double z)
                   { return rho / z; });

    return rho;
}

// std::valarray<double> OTAP::Hansen::Rho(std::valarray<double> P, std::valarray<double> T) const
// {
//     auto Z_value = Z(P, T);
//     return P / (Z_value * 287 * T);
// }

double OTAP::Hansen::Rho(double P, double T) const
{
    auto Z_value = Z(P, T);
    return P / (Z_value * 287 * T);
}

std::vector<double> OTAP::Hansen::Z(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::transform(P.begin(), P.end(), P.begin(), [](double p)
                   { return std::log10(p / 101325); });
    std::vector<double> Z(P.size(), 0.0);
    size_t dim[] = {m_logPData.size(), m_TData.size()};
    mlinterp::interp<mlinterp::rnatord>(dim, P.size(), m_ZData.data(), Z.data(), m_logPData.data(), P.data(), m_TData.data(), T.data());

    std::transform(Z.begin(), Z.end(), T.begin(), Z.begin(),
                   [](double z, double t)
                   {
                    if(t<500)
                        return 1.0;
                    else
                        return z; });

    return Z;
}

// std::valarray<double> OTAP::Hansen::Z(std::valarray<double> P, std::valarray<double> T) const
// {
//     vector _P(std::begin(P), std::end(P));
//     vector _T(std::begin(T), std::end(T));

//     auto Z = OTAP::Hansen::Z(_P, _T);
//     return std::valarray<double>(Z.data(), Z.size());
// }

double OTAP::Hansen::Z(double P, double T) const
{
    double Z;
    if (T >= 500)
    {
        P = std::log10(P / 101325);
        size_t dim[] = {m_logPData.size(), m_TData.size()};
        mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_ZData.data(), &Z, m_logPData.data(), &P, m_TData.data(), &T);
    }
    else
    {
        Z = 1;
    }
    return Z;
}

std::vector<double> OTAP::Hansen::Pr(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::transform(P.begin(), P.end(), P.begin(), [](double p)
                   { return std::log10(p / 101325); });
    std::vector<double> Pr(P.size(), 0.0);
    size_t dim[] = {m_logPData.size(), m_TData.size()};
    mlinterp::interp<mlinterp::rnatord>(dim, P.size(), m_PrData.data(), Pr.data(), m_logPData.data(), P.data(), m_TData.data(), T.data());

    std::transform(Pr.begin(), Pr.end(), T.begin(), Pr.begin(),
                   [](double pr, double t)
                   {
                    if(t<500)
                        return 0.738;
                    else
                        return pr; });

    return Pr;
}

// std::valarray<double> OTAP::Hansen::Pr(std::valarray<double> P, std::valarray<double> T) const
// {
//     vector _P(std::begin(P), std::end(P));
//     vector _T(std::begin(T), std::end(T));

//     auto Pr = OTAP::Hansen::Pr(_P, _T);
//     return std::valarray<double>(Pr.data(), Pr.size());
// }

double OTAP::Hansen::Pr(double P, double T) const
{
    double Z;
    if (T >= 500)
    {
        P = std::log10(P / 101325);
        size_t dim[] = {m_logPData.size(), m_TData.size()};
        mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_PrData.data(), &Z, m_logPData.data(), &P, m_TData.data(), &T);
    }
    else
    {
        Z = 0.738;
    }
    return Z;
}

std::vector<double> OTAP::Hansen::H(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::transform(P.begin(), P.end(), P.begin(), [](double p)
                   { return std::log10(p / 101325); });
    std::vector<double> H(P.size(), 0.0);
    size_t dim[] = {m_logPData.size(), m_TData.size()};

    std::transform(P.begin(), P.end(), T.begin(), H.begin(),
                   [this, &P, &T, &dim](double p, double t)
                   {
                       if (t < 500)
                           return Cp(p, t) * t;
                       else
                       {
                           double h;
                           mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_HData.data(), &h, m_logPData.data(), &p, m_TData.data(), &t);
                           return h;
                       }
                   });

    return H;
}

double OTAP::Hansen::H(double P, double T) const
{

    double Z;
    if (T >= 500)
    {
        P = std::log10(P / 101325);
        size_t dim[] = {m_logPData.size(), m_TData.size()};
        mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_HData.data(), &Z, m_logPData.data(), &P, m_TData.data(), &T);
    }
    else
    {
        Z = Cp(P, T) * T;
    }
    return Z;
}

// std::valarray<double> OTAP::Hansen::H(std::valarray<double> P, std::valarray<double> T) const
// {
//     vector _P(std::begin(P), std::end(P));
//     vector _T(std::begin(T), std::end(T));

//     auto H = OTAP::Hansen::H(_P, _T);
//     return std::valarray<double>(H.data(), H.size());
// }

std::vector<double> OTAP::Hansen::Gamma(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::transform(P.begin(), P.end(), P.begin(), [](double p)
                   { return std::log10(p / 101325); });
    std::vector<double> Gamma(P.size(), 0.0);
    size_t dim[] = {m_logPData.size(), m_TData.size()};
    mlinterp::interp<mlinterp::rnatord>(dim, P.size(), m_GammaData.data(), Gamma.data(), m_logPData.data(), P.data(), m_TData.data(), T.data());

    std::transform(Gamma.begin(), Gamma.end(), T.begin(), Gamma.begin(),
                   [](double g, double t)
                   {
                    if(t<500)
                        return 1.4;
                    else
                        return g; });

    return Gamma;
}

// std::valarray<double> OTAP::Hansen::Gamma(std::valarray<double> P, std::valarray<double> T) const
// {
//     vector _P(std::begin(P), std::end(P));
//     vector _T(std::begin(T), std::end(T));

//     auto gamma = OTAP::Hansen::Gamma(_P, _T);
//     return std::valarray<double>(gamma.data(), gamma.size());
// }

double OTAP::Hansen::Gamma(double P, double T) const
{

    double Z;
    if (T >= 500)
    {
        P = std::log10(P / 101325);
        size_t dim[] = {m_logPData.size(), m_TData.size()};
        mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_GammaData.data(), &Z, m_logPData.data(), &P, m_TData.data(), &T);
    }
    else
    {
        Z = 1.4;
    }
    return Z;
}

std::vector<double> OTAP::Hansen::Cp(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::transform(P.begin(), P.end(), P.begin(), [](double p)
                   { return std::log10(p / 101325); });
    std::vector<double> Cp(P.size(), 0.0);
    size_t dim[] = {m_logPData.size(), m_TData.size()};
    mlinterp::interp<mlinterp::rnatord>(dim, P.size(), m_CpData.data(), Cp.data(), m_logPData.data(), P.data(), m_TData.data(), T.data());

    std::transform(Cp.begin(), Cp.end(), T.begin(), Cp.begin(),
                   [](double cp, double t)
                   {
                    if(t<500)
                        return 1004.64;
                    else
                        return cp; });

    return Cp;
}

// std::valarray<double> OTAP::Hansen::Cp(std::valarray<double> P, std::valarray<double> T) const
// {
//     vector _P(std::begin(P), std::end(P));
//     vector _T(std::begin(T), std::end(T));

//     auto cp = OTAP::Hansen::Cp(_P, _T);
//     return std::valarray<double>(cp.data(), cp.size());
// }

double OTAP::Hansen::Cp(double P, double T) const
{

    double Z;
    if (T >= 500)
    {
        P = std::log10(P / 101325);
        size_t dim[] = {m_logPData.size(), m_TData.size()};
        mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_CpData.data(), &Z, m_logPData.data(), &P, m_TData.data(), &T);
    }
    else
    {
        Z = 1004.64;
    }
    return Z;
}

std::vector<double> OTAP::Hansen::k(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::transform(P.begin(), P.end(), P.begin(), [](double p)
                   { return std::log10(p / 101325); });
    std::vector<double> k(P.size(), 0.0);
    size_t dim[] = {m_logPData.size(), m_TData.size()};
    mlinterp::interp<mlinterp::rnatord>(dim, P.size(), m_kData.data(), k.data(), m_logPData.data(), P.data(), m_TData.data(), T.data());

    std::transform(P.begin(), P.end(), T.begin(), k.begin(),
                   [this, &P, &T, &dim](double p, double t)
                   {
                       if (t < 500)
                           return 1.364 * mu(p, t) * 1000;
                       else
                       {
                           double k;
                           mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_kData.data(), &k, m_logPData.data(), P.data(), m_TData.data(), T.data());
                           return k;
                       }
                   });

    return k;
}

// std::valarray<double> OTAP::Hansen::k(std::valarray<double> P, std::valarray<double> T) const
// {
//     vector _P(std::begin(P), std::end(P));
//     vector _T(std::begin(T), std::end(T));

//     auto k = OTAP::Hansen::k(_P, _T);
//     return std::valarray<double>(k.data(), k.size());
// }

double OTAP::Hansen::k(double P, double T) const
{
    double Z;
    if (T >= 500)
    {
        P = std::log10(P / 101325);
        size_t dim[] = {m_logPData.size(), m_TData.size()};
        mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_kData.data(), &Z, m_logPData.data(), &P, m_TData.data(), &T);
    }
    else
    {
        Z = 1.364 * mu(P, T) * 1000;
    }
    return Z;
}

std::vector<double> OTAP::Hansen::Cv(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());

    std::vector<double> Cv(P.size(), 0.0);
    auto Gamma_value = Gamma(P, T);
    auto Cp_value = Cp(P, T);

    std::transform(Cp_value.begin(), Cp_value.end(), Gamma_value.begin(), Cv.begin(), [](double cp, double gamma)
                   { return cp / gamma; });

    return Cv;
}

// std::valarray<double> OTAP::Hansen::Cv(std::valarray<double> P, std::valarray<double> T) const
// {
//     auto Gamma_value = Gamma(P, T);
//     auto Cp_value = Cp(P, T);
//     return Cp_value / Gamma_value;
// }

double OTAP::Hansen::Cv(double P, double T) const
{
    auto Gamma_value = Gamma(P, T);
    auto Cp_value = Cp(P, T);
    return Cp_value / Gamma_value;
}

std::vector<double> OTAP::Hansen::mu(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::transform(P.begin(), P.end(), P.begin(), [](double p)
                   { return std::log10(p / 101325); });
    std::vector<double> mu(P.size(), 0.0);
    size_t dim[] = {m_logPData.size(), m_TData.size()};
    mlinterp::interp<mlinterp::rnatord>(dim, P.size(), m_muData.data(), mu.data(), m_logPData.data(), P.data(), m_TData.data(), T.data());

    std::transform(P.begin(), P.end(), T.begin(), mu.begin(),
                   [this, &P, &T, &dim](double p, double t)
                   {
                       if (t < 500)
                           return (0.1462 * 0.00001) * ((std::pow(t, 0.5)) / (1 + 112 / t));
                       else
                       {
                           double mu;
                           mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_muData.data(), &mu, m_logPData.data(), &p, m_TData.data(), &t);
                           return mu;
                       }
                   });

    return mu;
}

// std::valarray<double> OTAP::Hansen::mu(std::valarray<double> P, std::valarray<double> T) const
// {
//     vector _P(std::begin(P), std::end(P));
//     vector _T(std::begin(T), std::end(T));

//     auto mu = OTAP::Hansen::mu(_P, _T);
//     return std::valarray<double>(mu.data(), mu.size());
// }

double OTAP::Hansen::mu(double P, double T) const
{
    double Z;
    if (T >= 500)
    {
        P = std::log10(P / 101325);
        size_t dim[] = {m_logPData.size(), m_TData.size()};
        mlinterp::interp<mlinterp::rnatord>(dim, size_t(1), m_muData.data(), &Z, m_logPData.data(), &P, m_TData.data(), &T);
    }
    else
    {
        Z = (0.1462 * 0.00001) * ((std::pow(T, 0.5)) / (1 + 112 / T));
    }
    return Z;
}

std::vector<double> OTAP::Hansen::a(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());

    std::vector<double> a(P.size(), 0.0);
    auto Gamma_value = Gamma(P, T);
    auto rho_value = Rho(P, T);

    std::transform(P.begin(), P.end(), Gamma_value.begin(), a.begin(), [](double p, double gamma)
                   { return gamma * p; });
    std::transform(rho_value.begin(), rho_value.end(), a.begin(), a.begin(), [](double rho, double a)
                   { return std::sqrt(a / rho); });

    return a;
}

// std::valarray<double> OTAP::Hansen::a(std::valarray<double> P, std::valarray<double> T) const
// {
//     auto Gamma_value = Gamma(P, T);
//     auto rho_value = Rho(P, T);
//     return std::sqrt(Gamma_value * P / rho_value);
// }

double OTAP::Hansen::a(double P, double T) const
{
    auto Gamma_value = Gamma(P, T);
    auto rho_value = Rho(P, T);
    return std::sqrt(Gamma_value * P / rho_value);
}

std::vector<double> OTAP::Hansen::Lambda(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());

    std::vector<double> lambda(P.size(), 0.0);
    auto mu_value = mu(P, T);
    auto rho_value = Rho(P, T);

    std::transform(mu_value.begin(), mu_value.end(), rho_value.begin(), lambda.begin(), [](double mu, double dens)
                   { return mu / dens; });
    std::transform(T.begin(), T.end(), lambda.begin(), lambda.begin(), [](double lambda, double temp)
                   { return lambda * std::sqrt(M_PI / (2 * 287 * temp)); });

    return lambda;
}

// std::valarray<double> OTAP::Hansen::Lambda(std::valarray<double> P, std::valarray<double> T) const
// {
//     auto mu_value = mu(P, T);
//     auto rho_value = Rho(P, T);
//     return mu_value / rho_value * std::sqrt(M_PI / (2 * 287 * T));
// }

double OTAP::Hansen::Lambda(double P, double T) const
{
    auto mu_value = mu(P, T);
    auto rho_value = Rho(P, T);
    return mu_value / rho_value * std::sqrt(M_PI / (2 * 287 * T));
}

OTAP::Downstream OTAP::PerfectGas::GetShockDownStream(const OTAP::Upstream &upstream) const
{

    Downstream d;
    using namespace std;
    if (upstream.M >= 1)
    {
        d.M = std::sqrt(((((1.4 - 1) * std::pow(upstream.M, 2) + 2) / (2 * 1.4 * std::pow(upstream.M, 2) - (1.4 - 1)))));
        double T_ratio, P_ratio, rho_ratio;
        P_ratio = 1 + ((2 * 1.4) / (1.4 + 1)) * (upstream.M * upstream.M - 1);
        rho_ratio = (((1.4 + 1) * std::pow(upstream.M, 2)) / ((1.4 - 1) * std::pow(upstream.M, 2) + 2));
        T_ratio = P_ratio / rho_ratio;
        d.P = upstream.P * P_ratio;
        d.T = upstream.T * T_ratio;
        d.T0 = upstream.T * (1 + 0.5 * (1.4 - 1) * std::pow(upstream.M, 2));
        d.P0 = d.P * (std::pow((1 + 0.5 * (1.4 - 1) * std::pow(d.M, 2)), (1.4 / (1.4 - 1))));
    }
    else
    {
        d.M = upstream.M;
        d.P = upstream.P;
        d.T = upstream.T;
        d.T0 = upstream.T * (1 + 0.5 * (1.4 - 1) * std::pow(upstream.M, 2));
        d.P0 = d.P * (std::pow((1 + 0.5 * (1.4 - 1) * std::pow(d.M, 2)), (1.4 / (1.4 - 1))));
    }
    return d;
}

// FIXME: Check gamma !!!
OTAP::EdgeProp OTAP::PerfectGas::GetEdgeProperties(const Upstream &upstream, double cp) const
{
    auto v = a(upstream.P, upstream.T) * upstream.M;
    auto ps = upstream.P + cp * 0.5 * Rho(upstream.P, upstream.T)*v*v;
  
   

    auto ds = GetShockDownStream(upstream);
    auto Me = std::sqrt(2 * (std::pow(ds.P0 / ps, 0.4 / 1.4) - 1) / 0.4);

    OTAP::EdgeProp e;
    e.M = Me;
    e.P0 = ds.P0;
    e.P = ps;
    e.T0 = ds.T0;
    e.T = ds.T0 / (1 + 0.2 * Me * Me);
   
    return e;
}

OTAP::Downstream OTAP::Hansen::GetShockDownStream(const OTAP::Upstream &upstream) const
{
    PerfectGas IS;
    Downstream d;

    auto velocity_fs = upstream.M * a(upstream.P, upstream.T);
    auto rho_fs = Rho(upstream.P, upstream.T);
    // Freestream.Knudsen_number = (Freestream.Mach * sqrt((0.5 * M_PI) / (287 * Freestream.Temperature))) / (Reference_length * Freestream.Density);

    if (upstream.M < 1)
    {
        d.T = upstream.T;
        d.P = upstream.P;
        d.M = upstream.M;
        d.T0 = upstream.T * (1 + 0.5 * (Gamma(upstream.P, upstream.T) - 1) * std::pow(upstream.M, 2));
        d.P0 = d.P * (std::pow((1 + 0.5 * (Gamma(upstream.P, upstream.T) - 1) * std::pow(d.M, 2)), (Gamma(upstream.P, upstream.T) / (Gamma(upstream.P, upstream.T) - 1))));
    }

    else
    {

        double Momentum_Total, MassFlow_Rate;
        MassFlow_Rate = rho_fs * velocity_fs;
        Momentum_Total = upstream.P + rho_fs * std::pow(velocity_fs, 2);
        // Compute_Downstream_Properties_Ideal(upstream, Downstream_Ideal, Theta, Downstream_Total, Reference_length);
        auto d_guess = IS.GetShockDownStream(upstream);
        auto rho_ds = Rho(d_guess.P, d_guess.T);
        d.T = d_guess.T;

        auto H0 = H(upstream.P, upstream.T) + 0.5 * std::pow(velocity_fs, 2);
        double iter_outer, iter_inner, del_Temp, del_H, H_old, rho_old, del_rho;
        iter_outer = 0;

        del_rho = 1;

        // cout << upstream.Mach << endl;

        auto velocity_ds = 0.0;
        auto H_ds = 0.0;
        while (del_rho > 0.02 && iter_outer < 11)
        {
            // while (iter_outer < 11) {
            velocity_ds = MassFlow_Rate / rho_ds;
            d.P = Momentum_Total - rho_ds * std::pow(velocity_ds, 2);
            H_ds = H0 - 0.5 * std::pow(velocity_ds, 2);
            H_old = H_ds;
            rho_old = rho_ds;
            iter_inner = 0;
            // cout << d.Temperature <<  " " << d.Enthalpy << " " << d.Density << endl;
            // cout << d.Pressure  << " " << d.Velocity << " " << Downstream_Total.Enthalpy << endl;
            // cout << Downstream_Total.Pressure << " " << upstream.Enthalpy << " " << velocity_fs <<  endl;
            del_H = 1;
            while ((abs(del_H) > 0.01) && (iter_inner < 51))
            {
                // while (iter_inner < 51)
                H_ds = H(d.P, d.T);
                del_Temp = ((H_ds - H_old) / Cp(d.P, d.T));
                if (del_Temp / d.T > 0.3)
                {
                    del_Temp = 0.3 * d.T;
                }
                del_H = ((H_old - H_ds) / H_old);
                d.T = d.T - del_Temp;

                iter_inner = iter_inner + 1;
            }
            rho_ds = (Rho(d.P, d.T) + rho_old) / 2;
            del_rho = std::abs((Rho(d.P, d.T) - rho_old) / rho_old);
            iter_outer = iter_outer + 1;
        }

        d.M = velocity_ds / a(d.P, d.T);
        // Downstream_Real.Knudsen_number = upstream.Knudsen_number * (Downstream_Real.Viscosity / upstream.Viscosity) * (upstream.SoundSpeed / Downstream_Real.SoundSpeed) * (upstream.Density / Downstream_Real.Density);
    }
    auto gamma_ds = Gamma(d.P, d.T);
    d.T0 = d.T * (1 + 0.5 * (gamma_ds - 1) * std::pow(d.M, 2));
    d.P0 = d.P * (std::pow((1 + 0.5 * (gamma_ds - 1) * std::pow(d.M, 2)), (gamma_ds / (gamma_ds - 1))));
    // Hansen(Downstream_Total);
    // Downstream_Total.Knudsen_number = Freestream.Knudsen_number * (Downstream_Total.Viscosity / Freestream.Viscosity) * (Freestream.SoundSpeed / Downstream_Total.SoundSpeed) * (Freestream.Density / Downstream_Total.Density);
    // cout << Downstream_Total.Temperature << endl;

    return d;
}

// FIXME: Check gamma!!!
OTAP::EdgeProp OTAP::Hansen::GetEdgeProperties(const Upstream &upstream, double cp) const
{
    auto v = a(upstream.P, upstream.T) * upstream.M;
    auto ps = upstream.P + cp * 0.5 * Rho(upstream.P, upstream.T)*v*v;

    auto ds = GetShockDownStream(upstream);
    auto Me = std::sqrt(2 * (std::pow(ds.P0 / ps, 0.4 / 1.4) - 1) / 0.4);

    OTAP::EdgeProp e;
    e.M = Me;
    e.P0 = ds.P0;
    e.P = ps;
    e.T0 = ds.T0;
    e.T = ds.T0 / (1 + 0.2 * Me * Me);

    return e;
}

double OTAP::PerfectGas::Rho(double P, double T) const
{
    return P / (m_R * T);
}

double OTAP::PerfectGas::Z(double P, double T) const
{
    return 1.0;
}

double OTAP::PerfectGas::Pr(double P, double T) const
{
    return mu(P, T) * Cp(P, T) / k(P, T);
}

double OTAP::PerfectGas::H(double P, double T) const
{
    return Cp(P, T) * T;
}

double OTAP::PerfectGas::Gamma(double P, double T) const
{
    return 1.4;
}

double OTAP::PerfectGas::Cp(double P, double T) const
{
    return 1004.6;
}
// Sutherland-Eucken Conductivity for Air!
double OTAP::PerfectGas::k(double P, double T) const
{
    return 1.364 * mu(P, T) * 1000;
}

double OTAP::PerfectGas::Cv(double P, double T) const
{
    return Cp(P, T) / Gamma(P, T);
}
// Sutherland Viscosity for Air!
double OTAP::PerfectGas::mu(double P, double T) const
{
    return (0.1462 * 0.00001) * ((std::pow(T, 0.5)) / (1 + 112 / T));
}

double OTAP::PerfectGas::a(double P, double T) const
{
    return std::sqrt(Gamma(P, T) * P / Rho(P, T));
}

double OTAP::PerfectGas::Lambda(double P, double T) const
{
    return mu(P, T) / Rho(P, T) * std::sqrt(M_PI / (2 * 287 * T));
}

// std::valarray<double> OTAP::PerfectGas::Rho(std::valarray<double> P, std::valarray<double> T) const
// {
//     return P / (m_R * T);
// }

// std::valarray<double> OTAP::PerfectGas::Z(std::valarray<double> P, std::valarray<double> T) const
// {
//     return std::valarray(1.0, P.size());
// }

// std::valarray<double> OTAP::PerfectGas::Pr(std::valarray<double> P, std::valarray<double> T) const
// {
//     return mu(P, T) * Cp(P, T) / k(P, T);
// }

// std::valarray<double> OTAP::PerfectGas::H(std::valarray<double> P, std::valarray<double> T) const
// {
//     return Cp(P, T) * T;
// }

// std::valarray<double> OTAP::PerfectGas::Gamma(std::valarray<double> P, std::valarray<double> T) const
// {
//     return std::valarray(1.4, P.size());
// }

// std::valarray<double> OTAP::PerfectGas::Cp(std::valarray<double> P, std::valarray<double> T) const
// {
//     return std::valarray(1004.6, P.size());
// }
// // Sutherland-Eucken Conductivity for Air!
// std::valarray<double> OTAP::PerfectGas::k(std::valarray<double> P, std::valarray<double> T) const
// {
//     return 1.364 * mu(P, T) * 1000;
// }

// std::valarray<double> OTAP::PerfectGas::Cv(std::valarray<double> P, std::valarray<double> T) const
// {
//     return Cp(P, T) / Gamma(P, T);
// }
// // Sutherland Viscosity for Air!
// std::valarray<double> OTAP::PerfectGas::mu(std::valarray<double> P, std::valarray<double> T) const
// {
//     return (0.1462 * 0.00001) * ((std::pow(T, 0.5)) / (1 + 112 / T));
// }

// std::valarray<double> OTAP::PerfectGas::a(std::valarray<double> P, std::valarray<double> T) const
// {
//     return std::sqrt(Gamma(P, T) * P / Rho(P, T));
// }

// std::valarray<double> OTAP::PerfectGas::Lambda(std::valarray<double> P, std::valarray<double> T) const
// {
//     return mu(P, T) / Rho(P, T) * std::sqrt(M_PI / (2 * 287 * T));
// }

std::vector<double> OTAP::PerfectGas::Rho(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::vector<double> rho(P.size(), 0.0);
    std::transform(P.begin(), P.end(), T.begin(), rho.begin(),
                   [this](double p, double t)
                   {
                       return Rho(p, t);
                   });
    return rho;
}

std::vector<double> OTAP::PerfectGas::Z(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    return std::vector<double>(P.size(), 1.0);
}

std::vector<double> OTAP::PerfectGas::Pr(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::vector<double> pr(P.size(), 0.0);
    std::transform(P.begin(), P.end(), T.begin(), pr.begin(),
                   [this](double p, double t)
                   {
                       return Pr(p, t);
                   });
    return pr;
}

std::vector<double> OTAP::PerfectGas::H(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::vector<double> h(P.size(), 0.0);
    std::transform(P.begin(), P.end(), T.begin(), h.begin(),
                   [this](double p, double t)
                   {
                       return H(p, t);
                   });
    return h;
}

std::vector<double> OTAP::PerfectGas::Gamma(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    return std::vector<double>(P.size(), 1.4);
}

std::vector<double> OTAP::PerfectGas::Cp(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    return std::vector<double>(P.size(), 1004.6);
}

std::vector<double> OTAP::PerfectGas::k(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::vector<double> k_value(P.size(), 0.0);
    std::transform(P.begin(), P.end(), T.begin(), k_value.begin(),
                   [this](double p, double t)
                   {
                       return k(p, t);
                   });
    return k_value;
}

std::vector<double> OTAP::PerfectGas::Cv(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::vector<double> cv(P.size(), 0.0);
    std::transform(P.begin(), P.end(), T.begin(), cv.begin(),
                   [this](double p, double t)
                   {
                       return Cv(p, t);
                   });
    return cv;
}

std::vector<double> OTAP::PerfectGas::mu(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::vector<double> mu_value(P.size(), 0.0);
    std::transform(P.begin(), P.end(), T.begin(), mu_value.begin(),
                   [this](double p, double t)
                   {
                       return mu(p, t);
                   });
    return mu_value;
}

std::vector<double> OTAP::PerfectGas::a(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::vector<double> a_value(P.size(), 0.0);
    std::transform(P.begin(), P.end(), T.begin(), a_value.begin(),
                   [this](double p, double t)
                   {
                       return a(p, t);
                   });
    return a_value;
}

std::vector<double> OTAP::PerfectGas::Lambda(std::vector<double> P, std::vector<double> T) const
{
    assert(P.size() == T.size());
    std::vector<double> lambda(P.size(), 0.0);
    std::transform(P.begin(), P.end(), T.begin(), lambda.begin(),
                   [this](double p, double t)
                   {
                       return Lambda(p, t);
                   });
    return lambda;
}
