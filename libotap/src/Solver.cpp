#include "Solver.h"
#include "Eigen/Dense"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include <cmath>
#include <algorithm>
#include <numeric>

OTAP::HFResult OTAP::FayRiddellSolver::Solve(double time)
{
    const double P_N2 = 4.5 * std::pow(10, 7);
    const double P_O2 = 2.3 * std::pow(10, 7);
    const double T_N2 = 113200;
    const double T_O2 = 59000;
    const double beta = 0.63;
    const double Lewis_No = 1.4;
    const double E_N2 = 4183;
    const double E_O2 = 7023;

    if (m_hasTrajectory)
    {
        if (time != 0.0)
        {
            m_currentTime = time;
        }
        auto Pinf = m_Trajectory->GetPinf(m_currentTime);
        auto Vinf = m_Trajectory->GetVinf(m_currentTime);
        Upstream us;
        us.P = Pinf;
        us.T = m_Trajectory->GetTinf(m_currentTime);
        us.M = m_Trajectory->GetMinf(m_currentTime);
        auto ds = m_Fluid->GetShockDownStream(us);

        m_Pinf = std::vector({Pinf});
        m_Vinf = std::vector({Vinf});
        m_P0 = std::vector({ds.P0});
        m_T0 = std::vector({ds.T0});
    }

    TimeSeries F_N2(m_P0.size());
    std::transform(m_P0.begin(), m_P0.end(), m_T0.begin(), F_N2.begin(), [&](double p, double t)
                   { return (P_N2 * 101325 / p) * (t / T_N2) * std::exp(-T_N2 / t); });

    TimeSeries F_O2(m_P0.size());
    std::transform(m_P0.begin(), m_P0.end(), m_T0.begin(), F_O2.begin(), [&](double p, double t)
                   { return (P_O2 * 101325 / p) * (t / T_O2) * std::exp(-T_O2 / t); });

    TimeSeries Dissociation_Enthalpy(m_P0.size());
    std::transform(F_N2.begin(), F_N2.end(), F_O2.begin(), Dissociation_Enthalpy.begin(), [&](double fn2, double fo2)
                   { return 4186 * (0.79 * E_N2 * std::sqrt(fn2 / (fn2 + 1)) + 0.21 * E_O2 * sqrt(fo2 / (fo2 + 1))); });

    TimeSeries Fay_Factor(m_P0.size());
    auto pr = m_Fluid->Pr(m_P0, m_T0);
    auto mu0 = m_Fluid->mu(m_P0, m_T0);
    auto muw = m_Fluid->mu(m_P0, std::vector(m_P0.size(), m_Tw));
    auto rho0 = m_Fluid->Rho(m_P0, m_T0);
    auto rhow = m_Fluid->Rho(m_P0, std::vector(m_P0.size(), m_Tw));
    auto H0 = m_Fluid->H(m_P0, m_T0);
    auto Hw = m_Fluid->H(m_P0, std::vector(m_P0.size(), m_Tw));
    for (size_t i = 0; i < Fay_Factor.size(); i++)
    {
        Fay_Factor[i] = std::pow(pr[i], -0.6) * std::pow((rhow[i] * muw[i]), 0.1) * std::pow((rho0[i] * mu0[i]), 0.4) * (H0[i] - Hw[i]) * (1 + (std::pow(Lewis_No, beta) - 1) * (Dissociation_Enthalpy[i] / H0[i]));
    }

    // FIXME: Check validity!!!!

    TimeSeries DuDx(m_P0.size());
    TimeSeries q(m_P0.size());
    TimeSeries h(m_P0.size());
    TimeSeries Tg(m_P0.size());

    switch (m_CoordinateType)
    {
    case CoordinateType::Cartesian:
        if (m_Geom.isPlate() && m_Geom.angle == 90.0) // Flat head 2D
        {
            auto geom = m_Geom;
            std::transform(m_Vinf.begin(), m_Vinf.end(), DuDx.begin(), [&](double v)
                           { return (0.3 * v) / (2 * geom.length); });

            std::transform(DuDx.begin(), DuDx.end(), Fay_Factor.begin(), q.begin(),
                           [&](double dudx, double ff)
                           { return 0.57 * ff * std::sqrt(dudx); });
        }
        else // Cylinder 2D
        {
            auto geom = m_Geom;
            for (size_t i = 0; i < DuDx.size(); i++)
            {
                DuDx[i] = (1 / m_Geom.radius) * std::sqrt((2 * (m_P0[i] - m_Pinf[i])) / rho0[i]);
            }

            std::transform(DuDx.begin(), DuDx.end(), Fay_Factor.begin(), q.begin(),
                           [&](double dudx, double ff)
                           { return 0.57 * ff * std::sqrt(dudx); });
        }
        break;
    case CoordinateType::Axisym:
        if (m_Geom.isPlate() && m_Geom.angle == 90.0) // flathead cylinder
        {
            auto geom = m_Geom;
            std::transform(m_Vinf.begin(), m_Vinf.end(), DuDx.begin(), [&](double v)
                           { return (0.3 * v) / (2 * geom.length); });

            std::transform(DuDx.begin(), DuDx.end(), Fay_Factor.begin(), q.begin(),
                           [&](double dudx, double ff)
                           { return 0.763 * ff * std::sqrt(dudx); });
        }
        else // sphere
        {
            auto geom = m_Geom;
            for (size_t i = 0; i < DuDx.size(); i++)
            {
                DuDx[i] = (1 / m_Geom.radius) * std::sqrt((2 * (m_P0[i] - m_Pinf[i])) / rho0[i]);
            }

            std::transform(DuDx.begin(), DuDx.end(), Fay_Factor.begin(), q.begin(),
                           [&](double dudx, double ff)
                           { return 0.763 * ff * std::sqrt(dudx); });
        }

        break;
    default:
        assert("Coordinate System not implemented");
        break;
    }

    auto tw = m_Tw;
    std::transform(q.begin(), q.end(), m_T0.begin(), h.begin(), [&tw](double q, double t)
                   { return q / (t - tw); });
    Tg = m_T0;

    HFResult Result;
    Result.h = h;
    Result.Tg = Tg;
    Result.q = q;

    return Result;
}

OTAP::HFResult OTAP::LeesSolver::Solve(double time)
{
    if (m_hasTrajectory)
    {
        if (time != 0.0)
        {
            m_currentTime = time;
        }
        auto Pinf = m_Trajectory->GetPinf(m_currentTime);
        auto Vinf = m_Trajectory->GetVinf(m_currentTime);
        auto Tinf = m_Trajectory->GetTinf(m_currentTime);
        auto Minf = m_Trajectory->GetMinf(m_currentTime);

        Upstream us;
        us.P = Pinf;
        us.T = Tinf;
        us.M = Minf;
        auto ds = m_Fluid->GetShockDownStream(us);

        m_Pinf = std::vector({Pinf});
        m_Vinf = std::vector({Vinf});
        m_Tinf = std::vector({Tinf});
        m_Minf = std::vector({Minf});
        m_P0 = std::vector({ds.P0});
        m_T0 = std::vector({ds.T0});
    }

    double Angle_Axis = m_RunningLength / m_Geom.radius;
    TimeSeries D_theta(m_Pinf.size());
    TimeSeries Ratio_spherical_lees(m_Pinf.size());
    double s_dash_by_Ro;
    TimeSeries D_Thetac(m_Pinf.size());
    TimeSeries A_Thetac(m_Pinf.size());
    TimeSeries B_Thetac(m_Pinf.size());
    TimeSeries Ratio_conical_lees(m_Pinf.size());
    TimeSeries q_Lees_spherical(m_Pinf.size());
    TimeSeries q_Lees_Conical(m_Pinf.size());
    TimeSeries h(m_Pinf.size());
    TimeSeries q(m_Pinf.size());

    s_dash_by_Ro = (cos(m_Geom.angle) / sin(m_Geom.angle)) + ((m_RunningLength / m_Geom.radius) - (M_PI / 2 - m_Geom.angle));

    auto Gamma_inf = m_Fluid->Gamma(m_Pinf, m_Tinf);

    for (size_t i = 0; i < m_Pinf.size(); i++)
    {
        D_theta[i] = (1 - (1 / (Gamma_inf[i] * m_Minf[i] * m_Minf[i]))) * (pow(Angle_Axis, 2) - ((Angle_Axis * sin(4 * Angle_Axis)) / 2) + ((1 - cos(4 * Angle_Axis)) / 8)) + (4 / (Gamma_inf[i] * m_Minf[i] * m_Minf[i])) * (pow(Angle_Axis, 2) - ((Angle_Axis * sin(2 * Angle_Axis))) + ((1 - cos(2 * Angle_Axis)) / 2));
        Ratio_spherical_lees[i] = (2 * Angle_Axis * sin(Angle_Axis) * (pow(cos(Angle_Axis), 2) + (pow(sin(Angle_Axis), 2)) * (1 / (Gamma_inf[i] * m_Minf[i] * m_Minf[i])))) / sqrt(D_theta[i]);
        q_Lees_spherical[i] = m_qstag[i] * Ratio_spherical_lees[i];

        D_Thetac[i] = (1 - (1 / (Gamma_inf[i] * m_Minf[i] * m_Minf[i]))) * (pow(M_PI / 2 - m_Geom.angle, 2) - (((M_PI / 2 - m_Geom.angle) * sin(4 * (M_PI / 2 - m_Geom.angle))) / 2) + ((1 - cos(4 * (M_PI / 2 - m_Geom.angle))) / 8)) + (4 / (Gamma_inf[i] * m_Minf[i] * m_Minf[i])) * (pow((M_PI / 2 - m_Geom.angle), 2) - (((M_PI / 2 - m_Geom.angle) * sin(2 * (M_PI / 2 - m_Geom.angle)))) + ((1 - cos(2 * (M_PI / 2 - m_Geom.angle))) / 2));

        A_Thetac[i] = ((sqrt(3) / 2) * sqrt((1 - (1 / (Gamma_inf[i] * m_Minf[i] * m_Minf[i]))) * pow(sin(m_Geom.angle), 2) + (1 / (Gamma_inf[i] * m_Minf[i] * m_Minf[i])))) * sqrt(M_PI - m_Geom.angle);
        B_Thetac[i] = (3 / 16) / (pow(sin(m_Geom.angle), 4) + (pow(sin(2 * m_Geom.angle), 2) / 4) / ((1 / (Gamma_inf[i] * m_Minf[i] * m_Minf[i])))) * (D_Thetac[i] / (M_PI / 2 - m_Geom.angle)) - (pow(cos(m_Geom.angle) / sin(m_Geom.angle), 3));

        Ratio_conical_lees[i] = (A_Thetac[i] * (s_dash_by_Ro)) / (sqrt(B_Thetac[i] + pow(s_dash_by_Ro, 3)));
        q_Lees_Conical[i] = m_qstag[i] * Ratio_spherical_lees[i];

        q[i] = q_Lees_spherical[i];
    }

    auto tw = m_Tw;
    std::transform(q.begin(), q.end(), m_T0.begin(), h.begin(), [&tw](double q, double t)
                   { return q / (t - tw); });
    HFResult Result;

    Result.h = h;
    Result.Tg = m_T0;
    Result.q = q;

    return Result;
}

OTAP::HFResult OTAP::BeckwithGallagherSolver::Solve(double time)
{
    if (m_hasTrajectory)
    {
        if (time != 0.0)
        {
            m_currentTime = time;
        }
        auto Pinf = m_Trajectory->GetPinf(m_currentTime);
        auto Vinf = m_Trajectory->GetVinf(m_currentTime);
        auto Tinf = m_Trajectory->GetTinf(m_currentTime);
        auto Minf = m_Trajectory->GetMinf(m_currentTime);

        Upstream us;
        us.P = Pinf;
        us.T = Tinf;
        us.M = Minf;
        auto ds = m_Fluid->GetShockDownStream(us);

        m_Pinf = std::vector({Pinf});
        m_Vinf = std::vector({Vinf});
        m_Tinf = std::vector({Tinf});
        m_Minf = std::vector({Minf});
        m_P0 = std::vector({ds.P0});
        m_T0 = std::vector({ds.T0});
    }

    TimeSeries Rey_No(m_Pinf.size());
    TimeSeries h_laminar(m_Pinf.size());
    TimeSeries TAW(m_Pinf.size());
    TimeSeries TREC(m_Pinf.size());
    TimeSeries q_laminar(m_Pinf.size());
    TimeSeries rec_laminar(m_Pinf.size());

    TimeSeries q_turbulent(m_Pinf.size());
    TimeSeries rec_turbulent(m_Pinf.size());
    TimeSeries h_turbulent(m_Pinf.size());

    TimeSeries Down_Mach(m_Pinf.size());
    TimeSeries Down_P(m_Pinf.size());
    TimeSeries Down_T(m_Pinf.size());
    TimeSeries Down_P0(m_Pinf.size());
    TimeSeries Down_T0(m_Pinf.size());
    TimeSeries Down_mu0(m_Pinf.size());
    TimeSeries Down_Gamma0(m_Pinf.size());
    TimeSeries Down_Pr0(m_Pinf.size());
    TimeSeries Down_Rho0(m_Pinf.size());

    TimeSeries DuDx(m_Pinf.size());
    TimeSeries B_Vis(m_Pinf.size());
    TimeSeries wall_Viscosity(m_Pinf.size());

    HFResult Result;

    auto m_Rhoinf = m_Fluid->Rho(m_Pinf, m_Tinf);
    auto mu_inf = m_Fluid->mu(m_Pinf, m_Tinf);
    auto Gamma_inf = m_Fluid->Gamma(m_Pinf, m_Tinf);
    auto k_inf = m_Fluid->k(m_Pinf, m_Tinf);

    switch (m_FlowType)
    {
    case FlowType::Laminar:

        for (size_t i = 0; i < m_Pinf.size(); i++)
        {

            Rey_No[i] = (m_Rhoinf[i] * m_Vinf[i] * 2 * m_Geom.radius) / mu_inf[i];

            if (m_Minf[i] < 1)
            {
                Down_Mach[i] = m_Minf[i] * std::cos(m_Geom.lambda);
                Down_P[i] = m_Pinf[i];
                Down_T0[i] = m_Tinf[i] * (1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2));
                Down_P0[i] = Down_P[i] * (std::pow((1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2)), (Gamma_inf[i]) / (Gamma_inf[i] - 1)));
            }
            else
            {
                Down_Mach[i] = (1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2)) / (Gamma_inf[i] * std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2) - 0.5 * (Gamma_inf[i] - 1));
                Down_Mach[i] = std::sqrt(Down_Mach[i]);
                Down_P[i] = m_Pinf[i] * (1 + ((2 * Gamma_inf[i]) / (Gamma_inf[i] + 1)) * (std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2) - 1));
                Down_T0[i] = m_Tinf[i] * (1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2));
                Down_P0[i] = Down_P[i] * (std::pow((1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(Down_Mach[i], 2)), (Gamma_inf[i]) / (Gamma_inf[i] - 1)));
            }

            Down_mu0[i] = m_Fluid->mu(Down_P0[i], Down_T0[i]);
            Down_Gamma0[i] = m_Fluid->mu(Down_P0[i], Down_T0[i]);
            Down_Pr0[i] = m_Fluid->mu(Down_P0[i], Down_T0[i]);

            //	h_laminar[i] = 0.5 * std::sqrt(2 * (Rey_No[i] / m_Minf[i]) * (Down_mu0[i] / mu_inf[i])) * (std::pow((2 / Gamma_inf[i]) * (m_Tinf[i] / Down_T0[i]) * (Down_P0[i] / m_Pinf[i]) * ((Down_P0[i] / m_Pinf[i]) - 1), 0.25)) * cos(m_Geom.lambda);

            h_laminar[i] = 0.5 * std::sqrt(2 * (Rey_No[i] / m_Minf[i]) * (Down_mu0[i] / mu_inf[i])) * (std::pow((2 / Down_Gamma0[i]) * (m_Tinf[i] / Down_T0[i]) * (Down_P0[i] / m_Pinf[i]) * ((Down_P0[i] / m_Pinf[i]) - 1), 0.25));

            h_laminar[i] = h_laminar[i] * (k_inf[i] / (2 * m_Geom.radius));
            rec_laminar[i] = std::sqrt(Down_Pr0[i]);

            TAW[i] = m_Tinf[i] * (1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(m_Minf[i], 2));
            TREC[i] = rec_laminar[i] * (TAW[i] - Down_T0[i]) + Down_T0[i];
            q_laminar[i] = h_laminar[i] * (TREC[i] - m_Tw);
        }

        Result.h = h_laminar;
        Result.Tg = TREC;
        Result.q = q_laminar;
        break;

    case FlowType::Turbulent:

        for (size_t i = 0; i < m_Pinf.size(); i++)
        {

            Rey_No[i] = (m_Rhoinf[i] * m_Vinf[i] * 2 * m_Geom.radius) / mu_inf[i];

            if (m_Minf[i] < 1)
            {
                Down_Mach[i] = m_Minf[i] * std::cos(m_Geom.lambda);
                Down_P[i] = m_Pinf[i];
                Down_T0[i] = m_Tinf[i] * (1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2));
                Down_P0[i] = Down_P[i] * (std::pow((1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2)), (Gamma_inf[i]) / (Gamma_inf[i] - 1)));
            }
            else
            {
                Down_Mach[i] = (1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2)) / (Gamma_inf[i] * std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2) - 0.5 * (Gamma_inf[i] - 1));
                Down_Mach[i] = std::sqrt(Down_Mach[i]);
                Down_P[i] = m_Pinf[i] * (1 + ((2 * Gamma_inf[i]) / (Gamma_inf[i] + 1)) * (std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2) - 1));
                Down_T0[i] = m_Tinf[i] * (1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(m_Minf[i] * std::cos(m_Geom.lambda), 2));
                Down_P0[i] = Down_P[i] * (std::pow((1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(Down_Mach[i], 2)), (Gamma_inf[i]) / (Gamma_inf[i] - 1)));
            }

            Down_Rho0[i] = m_Fluid->Rho(Down_P0[i], Down_T0[i]);

            DuDx[i] = (1 / m_Geom.radius) * std::sqrt((2 * (Down_P0[i] - m_Pinf[i])) / Down_Rho0[i]);

            TAW[i] = m_Tinf[i] * (1 + 0.5 * (Gamma_inf[i] - 1) * std::pow(m_Minf[i], 2));

            B_Vis[i] = m_Fluid->mu(Down_P0[i], Down_T0[i]);

            Down_Pr0[i] = m_Fluid->mu(Down_P0[i], Down_T0[i]);

            wall_Viscosity[i] = m_Fluid->mu(m_P0[i], m_Tw);

            rec_turbulent[i] = std::pow(Down_Pr0[i], 0.333333333333333);
            h_turbulent[i] = rec_turbulent[i] * std::pow((Rey_No[i] * 0.0228 * (Down_P0[i] / m_Pinf[i]) * (wall_Viscosity[i] / B_Vis[i]) * (m_Tinf[i] / m_Tw)), 0.8) * std::pow(std::sin(m_Geom.lambda), 0.6) * std::pow((0.130319149) * (B_Vis[i] / mu_inf[i]) * std::cos(m_Geom.lambda) * DuDx[i] * ((2 * m_Geom.radius) / m_Vinf[i]), 0.2);
            h_turbulent[i] = h_turbulent[i] * (k_inf[i] / (2 * m_Geom.radius));
            //////additional line as per OTAP
            // rec_turbulent=max(rec_turbulent,rec_laminar);
            // h_turbulent=max(h_turbulent,h_laminar);
            //////////////////////////////
            TREC[i] = rec_turbulent[i] * (TAW[i] - Down_T0[i]) + Down_T0[i];
            q_turbulent[i] = h_turbulent[i] * (TREC[i] - m_Tw);
        }

        Result.h = h_turbulent;
        Result.Tg = TREC;
        Result.q = q_turbulent;

        break;

    default:
        assert("Invalid Flow_type");
        break;
    }

    return Result;
}

OTAP::HFResult OTAP::VanDriestSolver::Solve(double time)
{
    if (m_hasTrajectory)
    {
        if (time != 0.0)
        {
            m_currentTime = time;
        }
        auto Pinf = m_Trajectory->GetPinf(m_currentTime);
        auto Vinf = m_Trajectory->GetVinf(m_currentTime);
        Upstream us;
        us.P = Pinf;
        us.T = m_Trajectory->GetTinf(m_currentTime);
        us.M = m_Trajectory->GetMinf(m_currentTime);
        auto ds = m_Fluid->GetShockDownStream(us);
        auto ep = m_Fluid->GetEdgeProperties(us, m_CpData.GetCp(us.M, m_Loc));

        m_Pinf = std::vector({Pinf});
        m_Vinf = std::vector({Vinf});
        m_P0 = std::vector({ds.P0});
        m_T0 = std::vector({ds.T0});
        m_Te = std::vector({ep.T});
        m_Me = std::vector({ep.M});
        m_Pe = std::vector({ep.P});
    }

    TimeSeries q(m_Pe.size());
    TimeSeries h(m_Pe.size());
    TimeSeries Tg(m_Pe.size());
    TimeSeries Edge_Velocity(m_Pe.size());
    TimeSeries Edge_Re(m_Pe.size());
    TimeSeries Edge_Re_corrected(m_Pe.size());
    HFResult Result;
    TimeSeries Temperature_Recovery(m_Pe.size());
    TimeSeries Rec_Temperature(m_Pe.size());
    TimeSeries Rec_Pressure(m_Pe.size());
    TimeSeries DELT_RECOVERY(m_Pe.size());
    TimeSeries Hrecovery(m_Pe.size());

    int ITER_TRECOVERY;

    auto Edge_Sound_Speed = m_Fluid->a(m_Pe, m_Te);
    auto Edge_rho = m_Fluid->Rho(m_Pe, m_Te);
    auto Edge_Viscosity = m_Fluid->mu(m_Pe, m_Te);
    auto Edge_Pr = m_Fluid->Pr(m_Pe, m_Te);
    auto Edge_Enthalpy = m_Fluid->H(m_Pe, m_Te);
    auto wall_Enthalpy = m_Fluid->H(m_Pe, std::vector(m_P0.size(), m_Tw));
    auto Edge_Gamma = m_Fluid->Gamma(m_Pe, m_Te);

    std::transform(m_Me.begin(), m_Me.end(), Edge_Sound_Speed.begin(), Edge_Velocity.begin(), [&](double Me, double ae)
                   { return Me * ae; });

    for (size_t i = 0; i < Edge_Velocity.size(); i++)
    {
        Edge_Re[i] = (Edge_rho[i] * Edge_Velocity[i] * m_RunningLength) / Edge_Viscosity[i];
    }

    TimeSeries F1(m_Pe.size());
    TimeSeries F2(m_Pe.size());

    TimeSeries Ch_laminar(m_Pe.size());
    TimeSeries rec_laminar(m_Pe.size());
    TimeSeries RAF_laminar(m_Pe.size());
    TimeSeries Cf_laminar(m_Pe.size());

    TimeSeries q_laminar(m_Pe.size());
    TimeSeries Shear(m_Pe.size());

    TimeSeries Edge_Re_corrected1(m_Pe.size());
    TimeSeries A(m_Pe.size());
    TimeSeries B(m_Pe.size());
    TimeSeries ALFA(m_Pe.size());
    TimeSeries BETA(m_Pe.size());
    TimeSeries C1(m_Pe.size());
    TimeSeries C2(m_Pe.size());
    TimeSeries C4(m_Pe.size());
    TimeSeries F(m_Pe.size());
    TimeSeries FD(m_Pe.size());
    TimeSeries Cf_turbulent(m_Pe.size());
    TimeSeries Q1(m_Pe.size());
    TimeSeries Q4(m_Pe.size());
    TimeSeries Q5(m_Pe.size());
    TimeSeries Ch_turbulent(m_Pe.size());
    TimeSeries rec_turbulent(m_Pe.size());
    TimeSeries RAF_turbulent(m_Pe.size());
    TimeSeries Cf_net(m_Pe.size());
    TimeSeries Ch(m_Pe.size());
    TimeSeries q_turbulent(m_Pe.size());

    TimeSeries AA1(m_Pe.size());
    TimeSeries ARG(m_Pe.size());
    TimeSeries Intermittency_Factor(m_Pe.size());

    switch (m_FlowType)
    {

    case FlowType::Laminar:

        for (size_t i = 0; i < Edge_Velocity.size(); i++)
        {

            Edge_Re_corrected[i] = Edge_Re[i] / (2 * m_IsCone + 1);
            F1[i] = 0.416594 - 0.246733 * 0.01 * m_Me[i] - 0.817489 * 0.001 * std::pow(m_Me[i], 2) + 0.2734033 * 0.0001 * std::pow(m_Me[i], 3);
            F2[i] = -0.134671 * 0.1 + 0.2635807 * 0.001 * m_Me[i] + 0.581944 * .0001 * std::pow(m_Me[i], 2) - 0.2173257 * 0.00001 * std::pow(m_Me[i], 3);
            Ch_laminar[i] = (F1[i] + F2[i] * (m_Tw / m_Te[i])) / (std::sqrt(Edge_Re_corrected[i]));
            rec_laminar[i] = std::sqrt(Edge_Pr[i]);
            RAF_laminar[i] = std::pow(Edge_Pr[i], 2 / 3);
            Cf_laminar[i] = 2 * RAF_laminar[i] * Ch_laminar[i];
            Hrecovery[i] = Edge_Enthalpy[i] + 0.5 * Edge_Velocity[i] * Edge_Velocity[i] * rec_laminar[i];
            q_laminar[i] = Ch_laminar[i] * Edge_rho[i] * Edge_Velocity[i] * (Hrecovery[i] - wall_Enthalpy[i]);
            Shear[i] = 0.5 * Edge_rho[i] * Edge_Velocity[i] * Edge_Velocity[i] * Cf_laminar[i];
        }
        Result.q = q_laminar;

        break;

    case FlowType::Turbulent:

        double n, C3, CFS, DELTA_CFS;

        for (size_t i = 0; i < m_Pe.size(); i++)
        {

            Edge_Re_corrected[i] = Edge_Re[i] / (2 * m_IsCone + 1);
            F1[i] = 0.416594 - 0.246733 * 0.01 * m_Me[i] - 0.817489 * 0.001 * std::pow(m_Me[i], 2) + 0.2734033 * 0.0001 * std::pow(m_Me[i], 3);
            F2[i] = -0.134671 * 0.1 + 0.2635807 * 0.001 * m_Me[i] + 0.581944 * .0001 * std::pow(m_Me[i], 2) - 0.2173257 * 0.00001 * std::pow(m_Me[i], 3);
            Ch_laminar[i] = (F1[i] + F2[i] * (m_Tw / m_Te[i])) / (std::sqrt(Edge_Re_corrected[i]));
            rec_laminar[i] = std::sqrt(Edge_Pr[i]);
            RAF_laminar[i] = std::pow(Edge_Pr[i], 2 / 3);
            Cf_laminar[i] = 2 * RAF_laminar[i] * Ch_laminar[i];

            double RETS = 0;
            Edge_Re_corrected1[i] = (Edge_Re[i] - RETS) / (m_IsCone + 1);
            A[i] = std::sqrt(((Edge_Gamma[i] - 1) * std::pow(m_Me[i], 2)) / (2 * (m_Tw / m_Te[i])));
            B[i] = ((1 + 0.5 * (Edge_Gamma[i] - 1) * std::pow(m_Me[i], 2)) / (m_Tw / m_Te[i])) - 1;
            ALFA[i] = (2 * A[i] * A[i] - B[i]) / std::sqrt(B[i] * B[i] + 4 * A[i] * A[i]);
            BETA[i] = B[i] / std::sqrt(B[i] * B[i] + 4 * A[i] * A[i]);

            n = 0.707;
            C1[i] = (n + 0.5) * std::log10(m_Tw / m_Te[i]);
            C2[i] = (0.242 / (std::sqrt(A[i] * A[i] * (m_Tw / m_Te[i])))) * (std::asin(ALFA[i]) + std::asin(BETA[i]));
            C3 = 0.41;
            DELTA_CFS = 1;
            CFS = 0.001;

            while (abs(DELTA_CFS) > 0.00001)
            {
                C4[i] = std::log10(Edge_Re_corrected1[i] * CFS);
                F[i] = C2[i] / std::sqrt(CFS) - C3 - C4[i] + C1[i];
                FD[i] = -(C2[i] / 2) * std::pow(CFS, -1.5) - std::log10(exp(1)) / CFS;
                DELTA_CFS = F[i] / FD[i];
                CFS = CFS - DELTA_CFS;
            }
            Cf_turbulent[i] = CFS;
            double Pr_t, K;
            Pr_t = 0.86;
            K = 0.4;
            double Q2, Q3;
            Q1[i] = (Edge_Pr[i] / Pr_t) - 1;
            Q2 = 1 - Pr_t;
            Q3 = ((std::sqrt(2 * CFS) * Q2) / K) * (std::pow(M_PI, 2) / 6 + 1.5 * Q2);
            Q4[i] = 12.5 * CFS * (Q1[i] + 2 * std::log(1 + (5 / 6) * Q1[i]) + std::log(6) * std::log(1 + (7 / 8) * Q1[i]) - (std::log(1 + (1 / 4) * Q1[i])));
            Q5[i] = std::log(1 + (5 / 6) * Q1[i]);
            rec_turbulent[i] = Pr_t * (1 + Q3 + Q4[i]);
            RAF_turbulent[i] = Pr_t * (1 + 0.5 * Q3 + 5 * std::sqrt(0.5 * CFS) * (Q1[i] + Q5[i]));
            Ch_turbulent[i] = Cf_turbulent[i] / (2 * RAF_turbulent[i]);

            // FIXME : Re transition based on Re theta
            // if(Transition_based_on_Retheta==1)
            // {
            //     Re_transition = (0.664)/((2*Re_transition*Edge.Viscosity)/(Edge.Density*Edge.Velocity));
            // }

            if (Edge_Re[i] <= m_Rec)
            {
                Intermittency_Factor[i] = 0;
                Cf_net[i] = Cf_laminar[i] + (Intermittency_Factor[i]) * (Cf_turbulent[i] - Cf_laminar[i]);
                Ch[i] = Ch_laminar[i];
                Hrecovery[i] = Edge_Enthalpy[i] + 0.5 * Edge_Velocity[i] * Edge_Velocity[i] * rec_laminar[i];
                q_turbulent[i] = Ch[i] * Edge_rho[i] * Edge_Velocity[i] * (Hrecovery[i] - wall_Enthalpy[i]);
                Shear[i] = 0.5 * Edge_rho[i] * Edge_Velocity[i] * Edge_Velocity[i] * Cf_net[i];
            }

            else
            {
                if (m_IsCone == 0)
                {
                    if (m_Rec == 0)
                    {
                        ARG[i] = -20;
                        Intermittency_Factor[i] = 1;
                    }
                    else
                    {
                        AA1[i] = 60 + 4.86 * std::pow(m_Me[i], 1.92);
                        ARG[i] = (3 / (AA1[i] * AA1[i])) * std::pow((Edge_Re[i] - m_Rec), 2) * std::pow(m_Rec, -1.34);
                    }
                }
                else
                {
                    if (m_Rec == 0)
                    {
                        ARG[i] = -20;
                        Intermittency_Factor[i] = 1;
                    }
                    else
                    {
                        AA1[i] = 60 + 4.86 * std::pow(m_Me[i], 1.92);
                        ARG[i] = (3 / (AA1[i] * AA1[i])) * (Edge_Re[i] - m_Rec) * std::pow(m_Rec, -0.34) * (std::log(Edge_Re[i] / m_Rec));
                    }
                }

                if (ARG[i] < -10)
                {
                    Intermittency_Factor[i] = 1;
                }
                else
                {
                    Intermittency_Factor[i] = 1 - std::exp(-ARG[i]);
                }
                Cf_net[i] = Cf_laminar[i] + (Intermittency_Factor[i]) * (Cf_turbulent[i] - Cf_laminar[i]);
                Ch[i] = Cf_net[i] / (2 * RAF_turbulent[i]);
                Hrecovery[i] = Edge_Enthalpy[i] + 0.5 * Edge_Velocity[i] * Edge_Velocity[i] * rec_turbulent[i];
                q_turbulent[i] = Ch[i] * Edge_rho[i] * Edge_Velocity[i] * (Hrecovery[i] - wall_Enthalpy[i]);
                Shear[i] = 0.5 * Edge_rho[i] * Edge_Velocity[i] * Edge_Velocity[i] * Cf_net[i];
            }
        }

        //////
        Result.q = q_turbulent;
        break;

    default:
        assert("Invalid Flow_type");
        break;
    }

    TimeSeries FINDING_TRECOVERY_Enthalpy(m_Pe.size());
    TimeSeries FINDING_TRECOVERY_Cp(m_Pe.size());

    for (size_t i = 0; i < Edge_Velocity.size(); i++)
    {
        // COMPUTING RECOVERY TEMPERATURE
        Temperature_Recovery[i] = m_T0[i];
        ITER_TRECOVERY = 0;

        Rec_Pressure[i] = m_Pe[i];
        Rec_Temperature[i] = Temperature_Recovery[i];
        FINDING_TRECOVERY_Enthalpy[i] = m_Fluid->H(Rec_Pressure[i], Rec_Temperature[i]);
        FINDING_TRECOVERY_Cp[i] = m_Fluid->Cp(Rec_Pressure[i], Rec_Temperature[i]);
        DELT_RECOVERY[i] = (Hrecovery[i] - FINDING_TRECOVERY_Enthalpy[i]) / FINDING_TRECOVERY_Cp[i];

        while (ITER_TRECOVERY <= 10 && abs(DELT_RECOVERY[i] / Temperature_Recovery[i]) < 0.01)
        {
            Temperature_Recovery[i] = Temperature_Recovery[i] + DELT_RECOVERY[i];
            Rec_Temperature[i] = Temperature_Recovery[i];
            FINDING_TRECOVERY_Enthalpy[i] = m_Fluid->H(Rec_Pressure[i], Rec_Temperature[i]);
            FINDING_TRECOVERY_Cp[i] = m_Fluid->Cp(Rec_Pressure[i], Rec_Temperature[i]);
            DELT_RECOVERY[i] = (Hrecovery[i] - FINDING_TRECOVERY_Enthalpy[i]) / FINDING_TRECOVERY_Cp[i];

            ITER_TRECOVERY = ITER_TRECOVERY + 1;
        }
    }

    auto tw = m_Tw;
    std::transform(Temperature_Recovery.begin(), Temperature_Recovery.end(), Result.q.begin(), h.begin(), [&tw](double Trec, double q)
                   { return q / (Trec - tw); });
    Result.h = h;
    Result.Tg = Temperature_Recovery;

    return Result;
}

OTAP::HFResult OTAP::EckertSolver::Solve(double time)
{
    if (m_hasTrajectory)
    {
        if (time != 0.0)
        {
            m_currentTime = time;
        }
        auto Pinf = m_Trajectory->GetPinf(m_currentTime);
        auto Vinf = m_Trajectory->GetVinf(m_currentTime);
        Upstream us;
        us.P = Pinf;
        us.T = m_Trajectory->GetTinf(m_currentTime);
        us.M = m_Trajectory->GetMinf(m_currentTime);
        auto ds = m_Fluid->GetShockDownStream(us);
        auto ep = m_Fluid->GetEdgeProperties(us, m_CpData.GetCp(us.M, m_Loc));

        m_Pinf = std::vector({Pinf});
        m_Vinf = std::vector({Vinf});
        m_P0 = std::vector({ds.P0});
        m_T0 = std::vector({ds.T0});
        m_Te = std::vector({ep.T});
        m_Me = std::vector({ep.M});
        m_Pe = std::vector({ep.P});
    }

    TimeSeries q(m_Pe.size());
    TimeSeries h(m_Pe.size());
    TimeSeries Tg(m_Pe.size());
    TimeSeries Edge_Velocity(m_Pe.size());
    TimeSeries Edge_Re(m_Pe.size());
    TimeSeries Edge_Re_corrected(m_Pe.size());
    HFResult Result;

    TimeSeries Temperature_Recovery(m_Pe.size());

    auto Edge_Sound_Speed = m_Fluid->a(m_Pe, m_Te);
    auto Edge_rho = m_Fluid->Rho(m_Pe, m_Te);
    auto Edge_Viscosity = m_Fluid->mu(m_Pe, m_Te);
    auto Edge_Pr = m_Fluid->Pr(m_Pe, m_Te);
    auto Edge_Z = m_Fluid->Z(m_Pe, m_Te);
    auto Edge_Cp = m_Fluid->Cp(m_Pe, m_Te);
    auto Edge_Enthalpy = m_Fluid->H(m_Pe, m_Te);
    auto wall_Enthalpy = m_Fluid->H(m_Pe, std::vector(m_P0.size(), m_Tw));
    auto Edge_Gamma = m_Fluid->Gamma(m_Pe, m_Te);

    std::transform(m_Me.begin(), m_Me.end(), Edge_Sound_Speed.begin(), Edge_Velocity.begin(), [&](double Me, double ae)
                   { return Me * ae; });

    for (size_t i = 0; i < Edge_Velocity.size(); i++)
    {
        Edge_Re[i] = (Edge_rho[i] * Edge_Velocity[i] * m_RunningLength) / Edge_Viscosity[i];
    }

    TimeSeries rec_laminar(m_Pe.size());
    TimeSeries TREC(m_Pe.size());
    TimeSeries TREC1(m_Pe.size());
    TimeSeries HREC(m_Pe.size());
    TimeSeries RECOVERY_Enthalpy(m_Pe.size());
    TimeSeries RECOVERY_Cp(m_Pe.size());
    TimeSeries HSTAR(m_Pe.size());
    TimeSeries TSTAR(m_Pe.size());
    TimeSeries TSTAR1(m_Pe.size());
    TimeSeries STAR_Enthalpy(m_Pe.size());
    TimeSeries Ch_laminar(m_Pe.size());
    TimeSeries RAF_laminar(m_Pe.size());
    TimeSeries Cf_laminar(m_Pe.size());
    TimeSeries Hrecovery(m_Pe.size());
    TimeSeries Shear(m_Pe.size());
    TimeSeries q_laminar(m_Pe.size());
    TimeSeries STAR_Rey_No(m_Pe.size());
    TimeSeries STAR_rho(m_Pe.size());
    TimeSeries STAR_Pr(m_Pe.size());
    TimeSeries STAR_Viscosity(m_Pe.size());

    TimeSeries q_turbulent(m_Pe.size());
    TimeSeries rec_turbulent(m_Pe.size());
    TimeSeries STARt_Rey_No(m_Pe.size());
    TimeSeries STAR_rhot(m_Pe.size());
    TimeSeries STAR_Prt(m_Pe.size());
    TimeSeries STAR_Viscosityt(m_Pe.size());
    TimeSeries TSTARt(m_Pe.size());
    TimeSeries TSTARt1(m_Pe.size());
    TimeSeries HSTARt(m_Pe.size());
    TimeSeries STAR_Enthalpyt(m_Pe.size());
    TimeSeries Ch_turbulent(m_Pe.size());
    TimeSeries RAF_turbulent(m_Pe.size());
    TimeSeries Cf_turbulent(m_Pe.size());
    TimeSeries Intermittency_Factor(m_Pe.size());
    TimeSeries Cf_net(m_Pe.size());
    TimeSeries Ch(m_Pe.size());
    TimeSeries ARG(m_Pe.size());
    TimeSeries AA1(m_Pe.size());

    double ITER_1, CONVERGENCE_CRITERIA, DELTA_TEMP;
    double DELTA_TEMP_STAR, CONVERGENCE_CRITERIA2, ITER_2;

    switch (m_FlowType)
    {

    case FlowType::Laminar:

        for (size_t i = 0; i < Edge_Velocity.size(); i++)
        {

            rec_laminar[i] = std::sqrt(Edge_Pr[i]);

            TREC[i] = m_Te[i] * (1 + (Edge_Gamma[i] - 1) * 0.5 * rec_laminar[i] * m_Me[i] * m_Me[i]);
            TREC1[i] = TREC[i];
            HREC[i] = Edge_Enthalpy[i] + 0.5 * Edge_Velocity[i] * Edge_Velocity[i] * rec_laminar[i];
            ITER_1 = 0;
            CONVERGENCE_CRITERIA = 1;
            while (CONVERGENCE_CRITERIA > 0.00001 && ITER_1 < 41)
            {
                TREC[i] = 0.5 * (TREC[i] + TREC1[i]);
                RECOVERY_Enthalpy[i] = m_Fluid->H(m_Pe[i], TREC[i]);
                RECOVERY_Cp[i] = m_Fluid->Cp(m_Pe[i], TREC[i]);
                TREC1[i] = TREC[i];
                DELTA_TEMP = (HREC[i] - RECOVERY_Enthalpy[i]) / (RECOVERY_Cp[i] + (0.5 * m_Me[i] * m_Me[i] * Edge_Gamma[i] * Edge_Z[i] * 287.0));
                CONVERGENCE_CRITERIA = abs(DELTA_TEMP / TREC[i]);
                if (CONVERGENCE_CRITERIA > 0.00001)
                {
                    TREC[i] = TREC[i] + DELTA_TEMP;
                }
                else
                {
                    break;
                }
                ITER_1 = ITER_1 + 1;
            }

            ITER_2 = 0;
            CONVERGENCE_CRITERIA2 = 1;
            TSTAR[i] = 0.5 * (m_Tw + m_Te[i]) + 0.22 * (TREC[i] - m_Te[i]);
            HSTAR[i] = 0.5 * (wall_Enthalpy[i] + Edge_Enthalpy[i]) + 0.22 * (HREC[i] - Edge_Enthalpy[i]);
            while (CONVERGENCE_CRITERIA2 > 0.0001 && ITER_2 < 31)
            {

                STAR_Enthalpy[i] = m_Fluid->H(m_Pe[i], TSTAR[i]);
                DELTA_TEMP_STAR = (HSTAR[i] - STAR_Enthalpy[i]) / (Edge_Cp[i] + (0.5 * m_Me[i] * m_Me[i] * Edge_Gamma[i] * Edge_Z[i] * 287.0));
                TSTAR1[i] = TSTAR[i];
                TSTAR[i] = TSTAR[i] + DELTA_TEMP_STAR;
                TSTAR[i] = 0.5 * (TSTAR[i] + TSTAR1[i]);
                CONVERGENCE_CRITERIA2 = abs(DELTA_TEMP_STAR / TSTAR[i]);
                ITER_2 = ITER_2 + 1;
            }

            STAR_Pr[i] = m_Fluid->Pr(m_Pe[i], TSTAR[i]);
            STAR_rho[i] = m_Fluid->Rho(m_Pe[i], TSTAR[i]);
            STAR_Viscosity[i] = m_Fluid->mu(m_Pe[i], TSTAR[i]);
            RAF_laminar[i] = std::pow(STAR_Pr[i], 0.666666666666667);
            STAR_Rey_No[i] = (Edge_Velocity[i] * m_RunningLength * STAR_rho[i]) / STAR_Viscosity[i];
            STAR_Rey_No[i] = STAR_Rey_No[i] / (2 * m_IsCone + 1);
            STAR_Rey_No[i] = STAR_Rey_No[i] / 2;
            Cf_laminar[i] = 0.664 / (std::sqrt(STAR_Rey_No[i]));
            Ch_laminar[i] = Cf_laminar[i] / (2 * RAF_laminar[i]);
            Hrecovery[i] = Edge_Enthalpy[i] + 0.5 * Edge_Velocity[i] * Edge_Velocity[i] * rec_laminar[i];
            q_laminar[i] = Ch_laminar[i] * STAR_rho[i] * Edge_Velocity[i] * (Hrecovery[i] - wall_Enthalpy[i]);
            Shear[i] = 0.5 * Edge_rho[i] * Edge_Velocity[i] * Edge_Velocity[i] * Cf_laminar[i];
        }

        Result.q = q_laminar;

        break;

    case FlowType::Turbulent:
        for (size_t i = 0; i < Edge_Velocity.size(); i++)
        {

            rec_laminar[i] = std::sqrt(Edge_Pr[i]);

            TREC[i] = m_Te[i] * (1 + (Edge_Gamma[i] - 1) * 0.5 * rec_laminar[i] * m_Me[i] * m_Me[i]);
            TREC1[i] = TREC[i];
            HREC[i] = Edge_Enthalpy[i] + 0.5 * Edge_Velocity[i] * Edge_Velocity[i] * rec_laminar[i];
            ITER_1 = 0;
            CONVERGENCE_CRITERIA = 1;
            while (CONVERGENCE_CRITERIA > 0.00001 && ITER_1 < 41)
            {
                TREC[i] = 0.5 * (TREC[i] + TREC1[i]);
                RECOVERY_Enthalpy[i] = m_Fluid->H(m_Pe[i], TREC[i]);
                RECOVERY_Cp[i] = m_Fluid->Cp(m_Pe[i], TREC[i]);
                TREC1[i] = TREC[i];
                DELTA_TEMP = (HREC[i] - RECOVERY_Enthalpy[i]) / (RECOVERY_Cp[i] + (0.5 * m_Me[i] * m_Me[i] * Edge_Gamma[i] * Edge_Z[i] * 287.0));
                CONVERGENCE_CRITERIA = abs(DELTA_TEMP / TREC[i]);
                if (CONVERGENCE_CRITERIA > 0.00001)
                {
                    TREC[i] = TREC[i] + DELTA_TEMP;
                }
                else
                {
                    break;
                }
                ITER_1 = ITER_1 + 1;
            }

            ITER_2 = 0;
            CONVERGENCE_CRITERIA2 = 1;
            TSTAR[i] = 0.5 * (m_Tw + m_Te[i]) + 0.22 * (TREC[i] - m_Te[i]);
            HSTAR[i] = 0.5 * (wall_Enthalpy[i] + Edge_Enthalpy[i]) + 0.22 * (HREC[i] - Edge_Enthalpy[i]);
            while (CONVERGENCE_CRITERIA2 > 0.0001 && ITER_2 < 31)
            {

                STAR_Enthalpy[i] = m_Fluid->H(m_Pe[i], TSTAR[i]);
                DELTA_TEMP_STAR = (HSTAR[i] - STAR_Enthalpy[i]) / (Edge_Cp[i] + (0.5 * m_Me[i] * m_Me[i] * Edge_Gamma[i] * Edge_Z[i] * 287.0));
                TSTAR1[i] = TSTAR[i];
                TSTAR[i] = TSTAR[i] + DELTA_TEMP_STAR;
                TSTAR[i] = 0.5 * (TSTAR[i] + TSTAR1[i]);
                CONVERGENCE_CRITERIA2 = abs(DELTA_TEMP_STAR / TSTAR[i]);
                ITER_2 = ITER_2 + 1;
            }

            STAR_Pr[i] = m_Fluid->Pr(m_Pe[i], TSTAR[i]);
            STAR_rho[i] = m_Fluid->Rho(m_Pe[i], TSTAR[i]);
            STAR_Viscosity[i] = m_Fluid->mu(m_Pe[i], TSTAR[i]);
            RAF_laminar[i] = std::pow(STAR_Pr[i], 0.666666666666667);
            STAR_Rey_No[i] = (Edge_Velocity[i] * m_RunningLength * STAR_rho[i]) / STAR_Viscosity[i];
            STAR_Rey_No[i] = STAR_Rey_No[i] / (2 * m_IsCone + 1);
            STAR_Rey_No[i] = STAR_Rey_No[i] / 2;
            Cf_laminar[i] = 0.664 / (std::sqrt(STAR_Rey_No[i]));

            rec_turbulent[i] = std::pow(Edge_Pr[i], 0.3333333333333333);
            TREC[i] = m_Te[i] * (1 + (Edge_Gamma[i] - 1) * 0.5 * rec_turbulent[i] * m_Me[i] * m_Me[i]);
            TREC1[i] = TREC[i];
            HREC[i] = Edge_Enthalpy[i] + 0.5 * Edge_Velocity[i] * Edge_Velocity[i] * rec_turbulent[i];
            ITER_1 = 0;
            CONVERGENCE_CRITERIA = 1;
            while (CONVERGENCE_CRITERIA > 0.00001 && ITER_1 < 41)
            {
                TREC[i] = 0.5 * (TREC[i] + TREC1[i]);
                RECOVERY_Enthalpy[i] = m_Fluid->H(m_Pe[i], TREC[i]);
                RECOVERY_Cp[i] = m_Fluid->Cp(m_Pe[i], TREC[i]);
                TREC1[i] = TREC[i];
                // DELTA_TEMP=(HREC[i]-RECOVERY_Enthalpy[i])/(RECOVERY_Cp[i]+(0.5*m_Me[i]*m_Me[i]*Edge_Gamma[i]*Edge_Z[i]*287.0));
                DELTA_TEMP = (HREC[i] - RECOVERY_Enthalpy[i]) / (RECOVERY_Cp[i]);
                TREC[i] = TREC[i] + DELTA_TEMP;
                ITER_1 = ITER_1 + 1;
                CONVERGENCE_CRITERIA = abs(DELTA_TEMP / TREC[i]);
            }

            ITER_2 = 0;
            CONVERGENCE_CRITERIA2 = 1;
            TSTARt[i] = 0.5 * (m_Tw + m_Te[i]) + 0.22 * (TREC[i] - m_Te[i]);
            HSTARt[i] = 0.5 * (wall_Enthalpy[i] + Edge_Enthalpy[i]) + 0.22 * (HREC[i] - Edge_Enthalpy[i]);
            while (CONVERGENCE_CRITERIA2 > 0.0001 && ITER_2 < 31)
            {

                STAR_Enthalpyt[i] = m_Fluid->H(m_Pe[i], TSTARt[i]);
                DELTA_TEMP_STAR = (HSTARt[i] - STAR_Enthalpyt[i]) / (Edge_Cp[i] + (0.5 * m_Me[i] * m_Me[i] * Edge_Gamma[i] * Edge_Z[i] * 287.0));
                TSTARt1[i] = TSTARt[i];
                TSTARt[i] = TSTARt[i] + DELTA_TEMP_STAR;
                TSTARt[i] = 0.5 * (TSTARt[i] + TSTARt1[i]);
                CONVERGENCE_CRITERIA2 = abs(DELTA_TEMP_STAR / TSTARt[i]);
                ITER_2 = ITER_2 + 1;
            }
            STAR_Prt[i] = m_Fluid->Pr(m_Pe[i], TSTARt[i]);
            STAR_rhot[i] = m_Fluid->Rho(m_Pe[i], TSTARt[i]);
            STAR_Viscosityt[i] = m_Fluid->mu(m_Pe[i], TSTARt[i]);
            RAF_turbulent[i] = std::pow(STAR_Prt[i], 0.66666666666667);
            double RETS = 0;
            STARt_Rey_No[i] = (Edge_Velocity[i] * m_RunningLength * STAR_rhot[i]) / STAR_Viscosityt[i];
            STARt_Rey_No[i] = (STARt_Rey_No[i] - RETS) / (1 * m_IsCone + 1);
            Cf_turbulent[i] = 0.37 / (std::pow(std::log10(STARt_Rey_No[i]), 2.584));

            if (Edge_Re[i] <= m_Rec)
            {
                Intermittency_Factor[i] = 0;
                Cf_net[i] = Cf_laminar[i];
                Ch[i] = Cf_laminar[i] / (2 * RAF_laminar[i]);
                Hrecovery[i] = Edge_Enthalpy[i] + 0.5 * Edge_Velocity[i] * Edge_Velocity[i] * rec_laminar[i];
                q_turbulent[i] = Ch[i] * STAR_rho[i] * Edge_Velocity[i] * (Hrecovery[i] - wall_Enthalpy[i]);
                Shear[i] = 0.5 * Edge_rho[i] * Edge_Velocity[i] * Edge_Velocity[i] * Cf_net[i];
            }

            else
            {
                if (m_IsCone == 0)
                {
                    if (m_Rec == 0)
                    {
                        ARG[i] = -20;
                        Intermittency_Factor[i] = 1;
                    }
                    else
                    {
                        AA1[i] = 60 + 4.86 * std::pow(m_Me[i], 1.92);
                        ARG[i] = (3 / (AA1[i] * AA1[i])) * std::pow((Edge_Re[i] - m_Rec), 2) * std::pow(m_Rec, -1.34);
                    }
                }
                else
                {
                    if (m_Rec == 0)
                    {
                        ARG[i] = -20;
                        Intermittency_Factor[i] = 1;
                    }
                    else
                    {
                        AA1[i] = 60 + 4.86 * std::pow(m_Me[i], 1.92);
                        ARG[i] = (3 / (AA1[i] * AA1[i])) * (Edge_Re[i] - m_Rec) * std::pow(m_Rec, -0.34) * (std::log(Edge_Re[i] / m_Rec));
                    }
                }

                if (ARG[i] < -10)
                {
                    Intermittency_Factor[i] = 1;
                }
                else
                {
                    Intermittency_Factor[i] = 1 - std::exp(-ARG[i]);
                }
                Cf_net[i] = Cf_laminar[i] + (Intermittency_Factor[i]) * (Cf_turbulent[i] - Cf_laminar[i]);
                Ch[i] = Cf_net[i] / (2 * RAF_turbulent[i]);
                Hrecovery[i] = Edge_Enthalpy[i] + 0.5 * Edge_Velocity[i] * Edge_Velocity[i] * rec_turbulent[i];
                q_turbulent[i] = Ch[i] * STAR_rhot[i] * Edge_Velocity[i] * (Hrecovery[i] - wall_Enthalpy[i]);
                Shear[i] = 0.5 * Edge_rho[i] * Edge_Velocity[i] * Edge_Velocity[i] * Cf_net[i];
            }
        }

        Result.q = q_turbulent;

        break;

    default:
        assert("Invalid Flow_type");
        break;
    }

    TimeSeries FINDING_TRECOVERY_Enthalpy(m_Pe.size());
    TimeSeries FINDING_TRECOVERY_Cp(m_Pe.size());
    TimeSeries Rec_Pressure(m_Pe.size());
    TimeSeries Rec_Temperature(m_Pe.size());
    TimeSeries DELT_RECOVERY(m_Pe.size());

    int ITER_TRECOVERY;

    for (size_t i = 0; i < Edge_Velocity.size(); i++)
    {
        // COMPUTING RECOVERY TEMPERATURE
        Temperature_Recovery[i] = m_T0[i];
        ITER_TRECOVERY = 0;

        Rec_Pressure[i] = m_Pe[i];
        Rec_Temperature[i] = Temperature_Recovery[i];
        FINDING_TRECOVERY_Enthalpy[i] = m_Fluid->H(Rec_Pressure[i], Rec_Temperature[i]);
        FINDING_TRECOVERY_Cp[i] = m_Fluid->Cp(Rec_Pressure[i], Rec_Temperature[i]);
        DELT_RECOVERY[i] = (Hrecovery[i] - FINDING_TRECOVERY_Enthalpy[i]) / FINDING_TRECOVERY_Cp[i];

        while (ITER_TRECOVERY <= 10 && abs(DELT_RECOVERY[i] / Temperature_Recovery[i]) < 0.01)
        {
            Temperature_Recovery[i] = Temperature_Recovery[i] + DELT_RECOVERY[i];
            Rec_Temperature[i] = Temperature_Recovery[i];
            FINDING_TRECOVERY_Enthalpy[i] = m_Fluid->H(Rec_Pressure[i], Rec_Temperature[i]);
            FINDING_TRECOVERY_Cp[i] = m_Fluid->Cp(Rec_Pressure[i], Rec_Temperature[i]);
            DELT_RECOVERY[i] = (Hrecovery[i] - FINDING_TRECOVERY_Enthalpy[i]) / FINDING_TRECOVERY_Cp[i];

            ITER_TRECOVERY = ITER_TRECOVERY + 1;
        }
    }

    auto tw = m_Tw;
    std::transform(Temperature_Recovery.begin(), Temperature_Recovery.end(), Result.q.begin(), h.begin(), [&tw](double Trec, double q)
                   { return q / (Trec - tw); });
    Result.h = h;
    Result.Tg = Temperature_Recovery;

    return Result;
}

OTAP::HFResult OTAP::FreeMolecularSolver::Solve(double time)
{
    if (m_hasTrajectory)
    {
        if (time != 0.0)
        {
            m_currentTime = time;
        }
        auto Pinf = m_Trajectory->GetPinf(m_currentTime);
        auto Vinf = m_Trajectory->GetVinf(m_currentTime);
        Upstream us;
        us.P = Pinf;
        us.T = m_Trajectory->GetTinf(m_currentTime);
        us.M = m_Trajectory->GetMinf(m_currentTime);
        auto ds = m_Fluid->GetShockDownStream(us);

        m_Pinf = std::vector({Pinf});
        m_Vinf = std::vector({Vinf});
    }

    TimeSeries Molecular_Speed_ratio(m_Pinf.size());
    TimeSeries q_FreeMolecular(m_Pinf.size());
    TimeSeries Freestream_Density(m_Pinf.size());

    TimeSeries FF1(m_Pinf.size());
    TimeSeries FF2(m_Pinf.size());
    TimeSeries FF3(m_Pinf.size());
    TimeSeries FF4(m_Pinf.size());
    TimeSeries h(m_Pinf.size());

    auto Freestream_Gamma = m_Fluid->Gamma(m_Pinf, m_Tinf);

    double Flow_angle_of_incidence = M_PI / 2;

    for (size_t i = 0; i < m_Pinf.size(); i++)
    {
        Molecular_Speed_ratio[i] = m_Vinf[i] / std::sqrt(2 * 287 * m_Tinf[i]);

        FF1[i] = 0.9 * Freestream_Density[i] * 287 * m_Tinf[i] * sqrt((287 * m_Tinf[i]) / (2 * M_PI));
        FF2[i] = std::pow(Molecular_Speed_ratio[i], 2) + (Freestream_Gamma[i] / (Freestream_Gamma[i] - 1)) - ((Freestream_Gamma[i] + 1) / (2 * (Freestream_Gamma[i] - 1))) * (m_Tw / m_Tinf[i]);
        FF3[i] = std::exp(-std::pow(Molecular_Speed_ratio[i] * sin(Flow_angle_of_incidence), 2)) + std::sqrt(M_PI) * Molecular_Speed_ratio[i] * sin(Flow_angle_of_incidence) * (1 + std::erf(Molecular_Speed_ratio[i] * std::sin(Flow_angle_of_incidence)));
        FF4[i] = 0.5 * exp(-pow(Molecular_Speed_ratio[i] * sin(Flow_angle_of_incidence), 2));
        q_FreeMolecular[i] = FF1[i] * (FF2[i] * FF3[i] - FF4[i]);
    }

    HFResult Result;

    Result.q = q_FreeMolecular;

    auto tw = m_Tw;
    std::transform(m_Tinf.begin(), m_Tinf.end(), Result.q.begin(), h.begin(), [&tw](double T0, double q)
                   { return q / (T0 - tw); });

    Result.h = h;
    Result.Tg = m_Tinf;

    return Result;
}

OTAP::HFResult OTAP::KnBridgedSolver::Solve(double time)
{
    TimeSeries Pinf;
    TimeSeries Vinf;
    TimeSeries Tinf;
    TimeSeries Kninf;

    if (m_hasTrajectory)
    {
        if (time != 0.0)
        {
            m_currentTime = time;
        }
        Pinf.push_back(m_Trajectory->GetPinf(m_currentTime));
        Tinf.push_back(m_Trajectory->GetTinf(m_currentTime));
        Vinf.push_back(m_Trajectory->GetVinf(m_currentTime));
        Kninf.push_back(m_Trajectory->GetLambdainf(m_currentTime) / m_characteristicLength);
        m_Pinffm = Pinf;
        m_Tinffm = Tinf;
        m_Vinffm = Vinf;
    }
    else
    {
        auto lambda = m_Fluid->Lambda(m_Pinffm, m_Tinffm);
        std::transform(lambda.begin(), lambda.end(), std::back_inserter(Kninf), [this](double l)
                       { return l / m_characteristicLength; });
    }

    auto conres = m_continuumSolver->Solve(m_currentTime, m_Tw);
    FreeMolecularSolver fm(m_CoordinateType, m_FlowType, m_Fluid, m_geom, m_Pinffm, m_Vinffm, m_Tinffm, m_Tw);
    auto fmres = fm.Solve(m_currentTime, m_Tw);

    HFResult res;
    TimeSeries factor;
    std::transform(Kninf.begin(), Kninf.end(), std::back_inserter(factor), [](double kn)
                   { return (1 - std::exp(-kn)); });

    for (size_t i = 0; i < factor.size(); i++)
    {
        auto f = (1 - factor[i]) * conres.q[i] + factor[i] * fmres.q[i];
        res.q.push_back(f);
        res.Tg.push_back(conres.Tg[i]);
        res.h.push_back(f / (conres.Tg[i] - m_Tw));
    }

    return res;
}

/// FOR CYLINDRICAL/SPHERICAL COORDINATE TRANSFORMATION (RADIUS != 0)
void OTAP::DefaultResponseSolver::CYSP2(double &A1, double &B1, double &C1, double k, double Inverse_ratio)
{
    double Coordinate_ratio, COORDINATE_COEF;
    if (m_options.coordinateType == CoordinateType::Cartesian)
    {
        COORDINATE_COEF = 0;
    }
    else if (m_options.coordinateType == CoordinateType::Axisym)
    {
        COORDINATE_COEF = 1;
    }
    else if (m_options.coordinateType == CoordinateType::Spherical)
    {
        COORDINATE_COEF = 2;
    }
    else
    {
        COORDINATE_COEF = 0;
    }

    Coordinate_ratio = (COORDINATE_COEF * k) / pow(Inverse_ratio, 2);
    A1 = A1 + Coordinate_ratio;
    B1 = B1 - 2 * Coordinate_ratio;
    C1 = C1 + Coordinate_ratio;
}

// FOR CYLINDRICAL/SPHERICAL COORDINATE TRANSFORMATION (RADIUS == 0)
void OTAP::DefaultResponseSolver::CYSP1(double &A1, double &C1, double k, double Layer_thickness, double Inverse_ratio, double Outer_radius_instantaneous)
{
    double Coordinate_ratio, COORDINATE_COEF;
    if (m_options.coordinateType == CoordinateType::Cartesian)
    {
        COORDINATE_COEF = 0;
    }
    else if (m_options.coordinateType == CoordinateType::Axisym)
    {
        COORDINATE_COEF = 1;
    }
    else if (m_options.coordinateType == CoordinateType::Spherical)
    {
        COORDINATE_COEF = 2;
    }
    else
    {
        COORDINATE_COEF = 0;
    }

    Coordinate_ratio = (COORDINATE_COEF * k * Layer_thickness) / (2 * Inverse_ratio * Outer_radius_instantaneous);
    A1 = A1 + Coordinate_ratio;
    C1 = C1 + Coordinate_ratio;
}

double OTAP::DefaultResponseSolver::FreeConvection_Surface(double T, double T_AMBIENT, double Characteristic_Length, double Sea_Level_Pressure)
{
    double Tbar, Acceleration_due_to_Gravity;
    Tbar = 0.5 * (T + T_AMBIENT);
    Acceleration_due_to_Gravity = 9.80665;
    auto Density_flow = m_Fluid->Rho(Sea_Level_Pressure, Tbar);
    auto Viscosity_flow = m_Fluid->mu(Sea_Level_Pressure, Tbar);
    auto Conductivity_flow = m_Fluid->k(Sea_Level_Pressure, Tbar);
    auto Pr_flow = m_Fluid->Pr(Sea_Level_Pressure, Tbar);

    double delta_Temp;
    delta_Temp = Tbar - T_AMBIENT;
    double Gr_No, Gr_No_Pr, Nu;
    Gr_No = (Acceleration_due_to_Gravity * delta_Temp * pow(Characteristic_Length, 3) * std::pow(Density_flow, 2)) / (Tbar * std::pow(Viscosity_flow, 2));
    Gr_No_Pr = Gr_No * Pr_flow;

    if (Gr_No_Pr > 1000000000)
    {
        Nu = 0.13 * std::pow(Gr_No_Pr, 0.333);
    }
    else
    {
        Nu = 0.59 * std::pow(Gr_No_Pr, 0.25);
    }

    return Nu * (Conductivity_flow / Characteristic_Length);
}

// Get h_convection & h_radiation incase of air gap
void OTAP::DefaultResponseSolver::Air_Gap(double &h_convection, double &h_radiation, double AIR_GAP, double LENGTH, double Temp_1, double Temp_2, double Sea_Level_pressure, double Emmisivity_1, double Emmisivity_2)
{
    double Sigma;
    Sigma = 5.67 * pow(10, -8);
    int FACTOR;
    double Tbar, beta, Diffusivity, Dynamic_Viscosity, Acceleration_due_to_gravity;
    h_convection = 5;
    h_radiation = 5;
    auto air_fluid = make_fluid(FluidType::Hansen_Air);

    if (Temp_1 == Temp_2)
    {
        h_convection = 0;
        h_radiation = 0;
    }
    else
    {
        if (Temp_1 > Temp_2)
        {
            FACTOR = 1;
        }
        else
        {
            FACTOR = -1;
        }
        Tbar = 0.5 * (Temp_1 + Temp_2);
        beta = 1 / Tbar;

        auto Viscosity_air_gap = air_fluid->mu(Sea_Level_pressure, Tbar);
        auto Conductivity_air_gap = air_fluid->k(Sea_Level_pressure, Tbar);
        auto Cp_air_gap = air_fluid->Cp(Sea_Level_pressure, Tbar);
        auto Density_air_gap = air_fluid->Rho(Sea_Level_pressure, Tbar);
        auto Pr_air_gap = air_fluid->Pr(Sea_Level_pressure, Tbar);

        Diffusivity = Conductivity_air_gap / (Density_air_gap * Cp_air_gap);
        Dynamic_Viscosity = Viscosity_air_gap / Density_air_gap;
        if (m_options.coordinateType == CoordinateType::Cartesian)
        {
            // Parallel plate
            double AIR_GAP_HALF;
            AIR_GAP_HALF = AIR_GAP / 2;

            double Gr_No, Gr_No_s, CL, CL_BAR, A74, C3, C3RA, Nu_1, S, A1, A2, A3, A4, i, DS, AI, Nu_R, h_convection1;
            Gr_No = ((Acceleration_due_to_gravity * beta * (Temp_1 - Temp_2) * (pow(AIR_GAP, 3))) / (Dynamic_Viscosity * Diffusivity)) * FACTOR;

            Gr_No_s = Gr_No * (AIR_GAP / LENGTH);
            CL = 0.48 * pow((Pr_air_gap / (Pr_air_gap + 0.861)), 0.25);
            CL_BAR = 1.3333 * CL;
            A74 = 0.919062526849;
            C3 = pow(12 * A74 * CL_BAR, 1.3333);
            C3RA = C3 / Gr_No_s;
            Nu_1 = CL_BAR * pow(Gr_No_s, 0.25);

            S = 0;
            A1 = -1;
            i = 1;
            A4 = 0;

            if (C3RA <= 20)
            {
                while (abs(DS) > 0.0001)
                {
                    AI = i;
                    A1 = A1 * (-1);
                    A2 = pow(C3RA, (AI - 1));
                    A3 = 3 / (4 * AI - 1);
                    A4 = A4 * (AI - 1);
                    if (A4 == 0)
                    {
                        A4 = 1;
                    }
                    DS = (A1 * A2 * A3) / A4;
                    S = S + DS;
                    i = i + 1;
                }
                Nu_R = Nu_1 * S;
                h_convection1 = (Nu_R * Conductivity_air_gap) / AIR_GAP;
            }
        }
        else if (m_options.coordinateType == CoordinateType::Axisym)
        {
            // Coaxial Cylinders

            double Diameter_Outer, Diameter_Inner, Thickness, Gr_No, Gr_No_s, A1, A2, EKEF, Nu_R, h_convection1;
            Diameter_Outer = LENGTH;
            Diameter_Inner = AIR_GAP;
            Thickness = 0.5 * (Diameter_Outer - Diameter_Inner);
            Gr_No = ((Acceleration_due_to_gravity * beta * (Temp_1 - Temp_2) * (pow(Thickness, 3))) / (Dynamic_Viscosity * Diffusivity)) * FACTOR;
            A1 = log(Diameter_Outer / Diameter_Inner);
            A2 = pow(Thickness, 3) * pow(((1 / pow(Diameter_Outer, 0.6)) + (1 / pow(Diameter_Inner, 0.6))), 5);
            Gr_No_s = (Gr_No * pow(A1, 4)) / A2;
            EKEF = 0.386 * pow((Pr_air_gap / (Pr_air_gap + 0.861)), 0.25) * pow(Gr_No_s, 0.25);
            if (EKEF = 1)
            {
                EKEF = 1;
            }

            h_convection1 = (2 * Conductivity_air_gap * EKEF) / (A1 * Diameter_Inner);
            Nu_R = (h_convection1 * Thickness) / Conductivity_air_gap;
        }
    }

    h_convection = h_convection * FACTOR;

    double ELH, F, AR1, AR2, AR3, FAC, QRAD;
    ELH = LENGTH / AIR_GAP;
    F = sqrt(1 + 1 / (ELH * ELH)) - 1 / ELH;
    if (m_options.coordinateType == CoordinateType::Axisym)
    {
        F = 1 / ELH;
    }
    AR1 = (1 - Emmisivity_1) / Emmisivity_1;
    AR2 = (1 - Emmisivity_2) / Emmisivity_2;
    AR3 = 1 / F;
    FAC = 1 / (AR1 + AR2 + AR3);
    QRAD = FAC * Sigma * (pow(Temp_1, 4) - pow(Temp_2, 4));
    h_radiation = QRAD / (Temp_1 - Temp_2);
}

void OTAP::DefaultResponseSolver::Instantaneous_MassRate_ThicknessSolver(double &Layer_thickness_1, double &Layer_thickness_2, double &Mass_rate_Char, double &Mass_rate_Pyrolysis, double Q_convective_front, double Qradiation_front, Eigen::VectorXd T, double deltat, const BCS &bcs)

{
    using namespace Eigen;
    double sigma, MPD, MCD;
    sigma = 5.67 * std::pow(10, -8);
    auto Layers = *m_Layers;
    VectorXd Inverse_nodes;

    for (size_t i = 0; i < Layers.GetCount(); i++)
    {
        Inverse_nodes(i) = 1. / Layers[i].numNodes;
    }

    auto nodes = m_Layers->Nodes;
    auto numnodes = m_Layers->NumNodes;

    MPD = Mass_rate_Pyrolysis;
    MCD = Mass_rate_Char;
    if (T(0) < nodes[0].Tabl() - m_params.Ttol || Layer_thickness_1 == 0 || m_options.Sublime == 1)
    {
        MCD = 0;
    }
    else
    {
        double TT, Choice_1;
        TT = (4 * T(1) - 3 * nodes[0].Tabl() - T(2)) / (2 * Inverse_nodes(0));
        // if (Layer_thickness(0) < X1_Thin)
        // {
        // 	TT = T(1) - Temp_Ablation;
        // }
        Qradiation_front = -std::transform_reduce(bcs.Tambient_front.begin(),
                                                  bcs.Tambient_front.end(),
                                                  0.0,
                                                  std::plus<double>(),
                                                  [&](double t)
                                                  {
                                                      return nodes[0].emissivity(T(0)) * sigma * (std::pow(nodes[0].Tabl(), 4) - std::pow(t, 4));
                                                  });
        Choice_1 = Q_convective_front + Qradiation_front + (nodes[0].k(T(0)) * TT) / Layer_thickness_1;
        MCD = std::max(Choice_1, 10.0);
        MCD = MCD / nodes[0].Habl();
    }

    if (((T(numnodes[1] - 1) < nodes[numnodes[1] - 1].Tpyro() - m_params.Ttol) && MPD == 0) || Layer_thickness_1 == 0)
    {
        MPD = 0;
    }

    else
    {
        if (Layer_thickness_1 == 0)
        {
            double TT, Choice_1;
            Qradiation_front = -std::transform_reduce(bcs.Tambient_front.begin(),
                                                      bcs.Tambient_front.end(),
                                                      0.0,
                                                      std::plus<double>(),
                                                      [&](double t)
                                                      {
                                                          return nodes[numnodes[1]].emissivity(T(numnodes[1] - 1)) * sigma * (std::pow(nodes[numnodes[1] - 1].Tpyro(), 4) - std::pow(t, 4));
                                                      });
            MPD = Q_convective_front + Qradiation_front;
            TT = (-4 * T(numnodes[1]) + 3 * nodes[numnodes[1] - 1].Tpyro() + T(numnodes[1] + 1)) / (2 * Inverse_nodes(1));
            Choice_1 = nodes[numnodes[1]].k(T[numnodes[1] - 1]) * TT;
            MPD = MPD - std::min(MPD, Choice_1 / Layer_thickness_1);
            MPD = MPD / nodes[numnodes[1] - 1].Hpyro();
        }

        else
        {

            double TT, Choice_1;
            TT = (4 * T(numnodes[1] - 2) - 3 * nodes[numnodes[1] - 1].Tpyro() - T(numnodes[1] - 3)) / (2 * Inverse_nodes(0));
            MPD = (nodes[numnodes[1] - 1].k(T[numnodes[1] - 1]) * TT) / Layer_thickness_1;
            TT = (-4 * T(numnodes[1]) + 3 * nodes[numnodes[1] - 1].Tpyro() + T(numnodes[1] + 1)) / (2 * Inverse_nodes(1));
            Choice_1 = nodes[numnodes[1]].k(T[numnodes[1] - 1]) * TT;
            MPD = MPD - std::min(MPD, Choice_1 / Layer_thickness_2);
            MPD = MPD / nodes[numnodes[1] - 1].Hpyro();
        }

        if (m_options.Sublime)
        {
            double VAL;
            Layer_thickness_1 = 0;
            VAL = ((Mass_rate_Pyrolysis + MPD) * deltat) / (2 * (nodes[0].rhovirgin() - nodes[0].rhochar()));
            Layer_thickness_1 = Layer_thickness_1 - VAL;
        }
        else
        {
            double VAL, VAL1;
            VAL = ((Mass_rate_Pyrolysis + MPD) * deltat) / (2 * (nodes[0].rhovirgin() - nodes[0].rhochar()));
            Layer_thickness_1 = Layer_thickness_1 - VAL;
            VAL1 = ((Mass_rate_Char + MCD) * deltat) / (2 * (nodes[0].rhopyrogas()));
            Layer_thickness_1 = Layer_thickness_1 + VAL - VAL1;
        }
    }

    Mass_rate_Pyrolysis = MPD;
    Mass_rate_Char = MCD;
}

void OTAP::DefaultResponseSolver::DefaultResponseMatrix(Eigen::VectorXd &T, Eigen::VectorXd Layer_thickness, Eigen::VectorXd Initial_Layer_thickness, double heat_transfer_coefficient_outer, double Tg, BCS bcs, double &Mass_rate_Char, double &Mass_rate_Pyrolysis, Eigen::VectorXd k, double delt)
{
    using namespace Eigen;

    // TODO: propellant_front, l_front, l_tg_front [[unused]]
    auto &[Flux_front,
           Flux_back,
           propmass_front,
           Propellant_Mass,
           qgen,
           Tamb_front,
           Tamb_back,
           h_front,
           h_back,
           Tg_front,
           Tg_back,
           l_front,
           l_back,
           l_tg_front,
           l_tg_back] = bcs;

    auto nodes = m_Layers->Nodes;
    auto numnodes = m_Layers->NumNodes;
    auto Layers = *m_Layers;
    // Layers[0].thickness

    double Total_thickness;
    Total_thickness = 0;
    for (size_t i = 0; i < Layers.GetCount(); i++)
    {
        Total_thickness = Total_thickness + Layers[i].thickness;
    }

    MatrixXd A, B;
    VectorXd SS;
    A = MatrixXd::Zero(numnodes[Layers.GetCount()], numnodes[Layers.GetCount()]);
    B = MatrixXd::Zero(numnodes[Layers.GetCount()], 1);
    // T = VectorXd::Constant(numnodes[Layers.GetCount()]-1, 0);
    SS = VectorXd::Constant(numnodes[Layers.GetCount()], qgen);


    VectorXd Inverse_nodes(Layers.GetCount());

    for (size_t i = 0; i < Layers.GetCount(); i++)
    {
        Inverse_nodes(i) = 1. / Layers[i].numNodes;
    }

    double MPD, MCD, CPG, Outer_radius_instantaneous, Outer_Radius, Layer_thickness_0, Layer_thickness_1, Initial_Layer_thickness_1, Initial_Layer_thickness_0;

    Outer_Radius = m_params.innerRadius + Total_thickness;
    // VectorXd Layer_thickness, Initial_Layer_thickness;

    // for (size_t i = 0; i < Layers.GetCount(); i++)
    // {

    //     Layer_thickness(i) = Layers[i].thickness;
    //     Initial_Layer_thickness(i) = Layers[i].thickness;
    // }

    if (Layers.GetCount() < 2)
    {
        Layer_thickness_0 = Layer_thickness(0);
        Layer_thickness_1 = 0;
        Initial_Layer_thickness_0 = Layers[0].thickness;
        Initial_Layer_thickness_1 = 0;
    }
    else
    {
        Layer_thickness_0 = Layer_thickness(0);
        Layer_thickness_1 = Layer_thickness(1);
        Initial_Layer_thickness_0 = Layers[0].thickness;
        Initial_Layer_thickness_1 = Layers[1].thickness;
    }

    if (m_options.coordinateType != CoordinateType::Cartesian)
    {
        Outer_radius_instantaneous = Outer_Radius - Initial_Layer_thickness_0 - Initial_Layer_thickness_1 + Layer_thickness(0) + Layer_thickness(1);
    }

    double AA1_s, AA2_s, BB1_s, BB2_s, CC1_s, CC2_s, DD1_s, DD2_s;
    const auto sigma = 5.67 * std::pow(10, -8);
    MPD = Mass_rate_Pyrolysis;
    MCD = Mass_rate_Char;
    CPG = nodes[0].Cppyrogas(T(0));

    if (Layer_thickness(0) != 0)
    {
        if (MCD != 0)
        {
            A(0, 0) = 1;
            A(0, 1) = 0;
            B(0, 0) = nodes[0].Tabl();
        }

        else if (MCD == 0)
        {
            AA2_s = k(0) / std::pow(Inverse_nodes(0), 2) - (MPD * CPG * Layer_thickness(0)) / (2 * Inverse_nodes(0));
            BB2_s = -((3 * k(0) + k(1)) / (2 * (pow(Inverse_nodes(0), 2)))) - (((nodes[0].rho(T(0)) * nodes[0].Cp(T(0))) / delt) * Layer_thickness(0) * Layer_thickness(0));
            CC2_s = ((k(0) + k(1)) / (2 * (pow(Inverse_nodes(0), 2))) + (MPD * CPG * Layer_thickness(0)) / (2 * Inverse_nodes(0)));
            DD2_s = -(((nodes[0].rho(T(0)) * nodes[0].Cp(T(0))) / delt) * T(0) * Layer_thickness(0) * Layer_thickness(0)) - SS(0) * Layer_thickness(0) * Layer_thickness(0);

            if (m_options.coordinateType != CoordinateType::Cartesian)
            {

                if (Outer_radius_instantaneous <= 0.0000001)
                {
                    CYSP2(AA2_s, BB2_s, CC2_s, k(0), Inverse_nodes(0));
                }
                else if (Outer_radius_instantaneous > 0.0000001)
                {
                    CYSP1(AA2_s, CC2_s, k(0), Layer_thickness(0), Inverse_nodes(0), Outer_radius_instantaneous);
                }
            }

            AA1_s = k(0) / (2 * Inverse_nodes(0));
            BB1_s = std::reduce(h_front.begin(), h_front.end(), heat_transfer_coefficient_outer) * Layer_thickness(0);
            CC1_s = -k(0) / (2 * Inverse_nodes(0));
            DD1_s = (std::inner_product(h_front.begin(), h_front.end(), Tg_front.begin(), heat_transfer_coefficient_outer * Tg) - nodes[0].emissivity(T(0)) * sigma * std::transform_reduce(Tamb_front.begin(), Tamb_front.end(), 0., std::plus<double>(), [&T](double t)
                                                                                                                                                                                            { return std::pow(T(0), 4) - std::pow(t, 4); }) +
                     Flux_front) *
                    Layer_thickness(0);

            A(0, 0) = BB1_s * AA2_s - AA1_s * BB2_s;
            A(0, 1) = AA2_s * CC1_s - AA1_s * CC2_s;
            B(0, 0) = AA2_s * DD1_s - AA1_s * DD2_s;
        }

        //		//X1THIN OPTION TO ADD

        double Factor1, Factor2;
        for (int i = 1; i < numnodes[1] - 1; i++)
        {

            A(i, i - 1) = (k(i - 1) + k(i)) / (2 * (std::pow(Inverse_nodes(0), 2)));
            A(i, i) = -((k(i - 1) + 2 * k(i) + k(i + 1)) / (2 * (std::pow(Inverse_nodes(0), 2)))) - ((nodes[i].rho(T(i)) * nodes[i].Cp(T(i))) / delt) * Layer_thickness(0) * Layer_thickness(0);
            A(i, i + 1) = (k(i) + k(i + 1)) / (2 * (std::pow(Inverse_nodes(0), 2)));
            Factor1 = (Layer_thickness(0) / Inverse_nodes(0)) * (MPD * CPG + (MCD * nodes[i].rho(T(i)) * nodes[i].Cp(T(i))) / nodes[i].rhopyrogas() + ((i) / (numnodes[1] - 1)) * ((MPD * nodes[i].rho(T(i)) * nodes[i].Cp(T(i))) / (nodes[i].rhovirgin() - nodes[i].rhochar()) - (MCD * nodes[i].rho(T(i)) * nodes[i].Cp(T(i))) / nodes[i].rhopyrogas()));
            Factor2 = Factor1 / 2;

            if (A(i, i - 1) < Factor2)
            {
                A(i, i) = A(i, i) - Factor1;
                A(i, i + 1) = A(i, i + 1) + Factor1;
            }
            else
            {
                A(i, i - 1) = A(i, i - 1) - Factor2;
                A(i, i + 1) = A(i, i + 1) + Factor2;
            }

            B(i, 0) = -(((nodes[i].rho(T(i)) * nodes[i].Cp(T(i))) / delt) * T(i) * Layer_thickness(0) * Layer_thickness(0)) - SS(i) * Layer_thickness(0) * Layer_thickness(0);

            if (m_options.coordinateType != CoordinateType::Cartesian)
            {
                Outer_radius_instantaneous = Outer_radius_instantaneous - Layer_thickness(0) * Inverse_nodes(0);
                if (Outer_radius_instantaneous <= 0.0000001)
                {
                    CYSP2(A(i, i - 1), A(i, i), A(i, i + 1), k(i), Inverse_nodes(0));
                }
                else if (Outer_radius_instantaneous > 0.0000001)
                {
                    CYSP1(A(i, i - 1), A(i, i + 1), k(i), Layer_thickness(0), Inverse_nodes(0), Outer_radius_instantaneous);
                }
            }
        }

        if (m_options.coordinateType != CoordinateType::Cartesian)
        {
            Outer_radius_instantaneous = Outer_radius_instantaneous - Layer_thickness(0) * Inverse_nodes(0);
        }

        if (Layers.GetCount() >= 2)
        {
            if (MPD != 0)
            {
                A(numnodes[1] - 1, numnodes[1] - 2) = 0;
                A(numnodes[1] - 1, numnodes[1] - 1) = 1;
                A(numnodes[1] - 1, numnodes[1]) = 0;
                B(numnodes[1] - 1, 0) = nodes[numnodes[1] - 1].Tpyro();
            }

            else if (MPD == 0 && Layer_thickness(1) != 0)
            {

                double AA1_i, BB1_i, CC1_i, DD1_i, AA3_i, BB3_i, CC3_i, DD3_i, AA2_i, AA2D_i, CC2_i, CC2D_i, BB2, DD2;
                AA1_i = (k(numnodes[1] - 1) + k(numnodes[1] - 2)) / (2 * std::pow(Inverse_nodes(0), 2));
                BB1_i = -(k(numnodes[1] - 2) + 3 * k(numnodes[1] - 1)) / (2 * std::pow(Inverse_nodes(0), 2)) - (nodes[numnodes[1] - 1].rho(T(numnodes[1] - 1)) * nodes[numnodes[1] - 1].Cp(T(numnodes[1] - 1)) / delt) * Layer_thickness(0) * Layer_thickness(0);
                CC1_i = k(numnodes[1] - 1) / (1 * std::pow(Inverse_nodes(0), 2));
                DD1_i = -((nodes[numnodes[1] - 1].rho(T(numnodes[1] - 1)) * nodes[numnodes[1] - 1].Cp(T(numnodes[1] - 1))) / delt) * (T(numnodes[1] - 1) * Layer_thickness(0) * Layer_thickness(0)) - SS(numnodes[1] - 1) * Layer_thickness(0) * Layer_thickness(0);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {

                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(AA1_i, BB1_i, CC1_i, k(numnodes[1] - 1), Inverse_nodes(0));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(AA1_i, CC1_i, k(numnodes[1] - 1), Layer_thickness(0), Inverse_nodes(0), Outer_radius_instantaneous);
                    }
                }

                AA3_i = k(numnodes[1]) / (1 * std::pow(Inverse_nodes(1), 2));
                BB3_i = -(k(numnodes[1] + 1) + 3 * k(numnodes[1])) / (2 * std::pow(Inverse_nodes(1), 2)) - (nodes[numnodes[1]].rho(T(numnodes[1] - 1)) * nodes[numnodes[1]].Cp(T(numnodes[1] - 1)) / delt) * Layer_thickness(1) * Layer_thickness(1);
                CC3_i = (k(numnodes[1]) + k(numnodes[1] + 1)) / (2 * std::pow(Inverse_nodes(1), 2));
                DD3_i = -((nodes[numnodes[1]].rho(T(numnodes[1] - 1)) * nodes[numnodes[1]].Cp(T(numnodes[1] - 1)) / delt) * T(numnodes[1] - 1) * Layer_thickness(1) * Layer_thickness(1)) - SS(numnodes[1] - 1) * Layer_thickness(1) * Layer_thickness(1);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {

                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(AA3_i, BB3_i, CC3_i, k(numnodes[1]), Inverse_nodes(1));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(AA3_i, CC3_i, k(numnodes[1]), Layer_thickness(1), Inverse_nodes(1), Outer_radius_instantaneous);
                    }
                }

                AA2_i = -k(numnodes[1] - 1) / (2 * Inverse_nodes(0) * Layer_thickness(0));
                AA2D_i = k(numnodes[1]) / (2 * Inverse_nodes(1) * Layer_thickness(1));
                CC2_i = -AA2D_i;
                CC2D_i = -AA2_i;
                BB2 = 0;
                DD2 = MPD * nodes[numnodes[1] - 1].Hpyro();

                A(numnodes[1] - 1, numnodes[1] - 2) = AA1_i * CC2D_i * AA3_i - CC1_i * AA2_i * AA3_i;
                A(numnodes[1] - 1, numnodes[1] - 1) = BB3_i * CC1_i * AA2D_i + BB1_i * CC2D_i * AA3_i - BB2 * CC1_i * AA3_i;
                A(numnodes[1] - 1, numnodes[1]) = CC1_i * AA2D_i * CC3_i - AA3_i * CC1_i * CC2_i;
                B(numnodes[1] - 1, 0) = DD3_i * CC1_i * AA2D_i + AA3_i * CC2D_i * DD1_i - DD2 * CC1_i * AA3_i;
            }

            if ((Layers.GetCount() == 1) || (MPD == 0 && Layer_thickness(1) == 0 && Layers.GetCount() < 3))
            {
                double AA2_b, BB2_b, CC2_b, DD2_b, AA1_b, BB1_b, CC1_b, DD1_b;

                AA2_b = (k(numnodes[1] - 1) + k(numnodes[1] - 2)) / (2 * std::pow(Inverse_nodes(0), 2));
                BB2_b = -(k(numnodes[1] - 2) + 3 * k(numnodes[1] - 1)) / (2 * pow(Inverse_nodes(0), 2)) - (nodes[numnodes[1] - 1].rho(T(numnodes[1] - 1)) * nodes[numnodes[1] - 1].Cp(T(numnodes[1] - 1)) / delt) * Layer_thickness(0) * Layer_thickness(0);
                CC2_b = k(numnodes[1] - 1) / (1 * std::pow(Inverse_nodes(0), 2));
                DD2_b = -((nodes[numnodes[1] - 1].rho(T(numnodes[1] - 1)) * nodes[numnodes[1] - 1].Cp(T(numnodes[1] - 1))) / delt) * T(numnodes[1] - 1) * Layer_thickness(0) * Layer_thickness(0) - SS(numnodes[1] - 1) * Layer_thickness(0) * Layer_thickness(0);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {

                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(AA2_b, BB2_b, CC2_b, k(numnodes[1] - 1), Inverse_nodes(0));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(AA2_b, CC2_b, k(numnodes[1] - 1), Layer_thickness(0), Inverse_nodes(0), Outer_radius_instantaneous);
                    }
                }

                for (size_t i = 0; i < l_back.size(); i++)
                {
                    h_back.push_back(FreeConvection_Surface(T(numnodes[1] - 1), l_tg_back[i], l_back[i], m_params.Sea_Level_Pressure));
                    Tg_back.push_back(l_tg_back[i]);
                }

                std::vector<double> h_radiation_loss_backwall;
                for (size_t i = 0; i < Tamb_back.size(); i++)
                {
                    h_radiation_loss_backwall.push_back(sigma * nodes[numnodes[1] - 1].emissivity(numnodes[1] - 1) * (std::pow(T(numnodes[1] - 1), 3) + std::pow(T(numnodes[1] - 1), 2) * Tamb_back[i] + T(numnodes[1] - 1) * std::pow(Tamb_back[i], 2) + std::pow(Tamb_back[i], 3)));
                }

                AA1_b = k(numnodes[1] - 1) / (2 * Layer_thickness(0) * Inverse_nodes(0));
                BB1_b = std::reduce(h_back.begin(), h_back.end(), 0.) + std::reduce(h_radiation_loss_backwall.begin(), h_radiation_loss_backwall.end(), 0.) + Propellant_Mass / delt;
                CC1_b = -k(numnodes[1] - 1) / (2 * Layer_thickness(0) * Inverse_nodes(0));
                DD1_b = Flux_back + std::inner_product(h_back.begin(), h_back.end(), Tg_back.begin(), 0.) + std::inner_product(h_radiation_loss_backwall.begin(), h_radiation_loss_backwall.end(), Tamb_back.begin(), 0.) +
                        (Propellant_Mass * T(numnodes[1] - 1) / delt);

                A(numnodes[1] - 1, numnodes[1] - 2) = CC1_b * AA2_b - CC2_b * AA1_b;
                A(numnodes[1] - 1, numnodes[1] - 1) = CC1_b * BB2_b - CC2_b * BB1_b;
                B(numnodes[1] - 1, 0) = CC1_b * DD2_b - CC2_b * DD1_b;
            }
        }

        if (MPD == 0 && Layer_thickness(1) == 0 && Layers.GetCount() >= 3)
        {
            double AA1_i, BB1_i, CC1_i, DD1_i, AA3_i, BB3_i, CC3_i, DD3_i, AA2_i, AA2D_i, CC2_i, CC2D_i, BB2, DD2;

            AA1_i = (k(numnodes[1] - 1) + k(numnodes[1] - 2)) / (2 * std::pow(Inverse_nodes(0), 2));
            BB1_i = -(k(numnodes[1] - 2) + 3 * k(numnodes[1] - 1)) / (2 * std::pow(Inverse_nodes(0), 2)) - (nodes[numnodes[1] - 1].rho(T(numnodes[1] - 1)) * nodes[numnodes[1] - 1].Cp(T(numnodes[1] - 1)) / delt) * Layer_thickness(0) * Layer_thickness(0);
            CC1_i = k(numnodes[1] - 1) / (1 * std::pow(Inverse_nodes(0), 2));
            DD1_i = -((nodes[numnodes[1] - 1].rho(T(numnodes[1] - 1)) * nodes[numnodes[1] - 1].Cp(T(numnodes[1] - 1))) / delt) * (T(numnodes[1] - 1) * Layer_thickness(0) * Layer_thickness(0)) - SS(numnodes[1] - 1) * Layer_thickness(0) * Layer_thickness(0);

            if (m_options.coordinateType != CoordinateType::Cartesian)
            {

                if (Outer_radius_instantaneous <= 0.0000001)
                {
                    CYSP2(AA1_i, BB1_i, CC1_i, k(numnodes[1] - 1), Inverse_nodes(0));
                }
                else if (Outer_radius_instantaneous > 0.0000001)
                {
                    CYSP1(AA1_i, CC1_i, k(numnodes[1] - 1), Layer_thickness(0), Inverse_nodes(0), Outer_radius_instantaneous);
                }
            }

            AA3_i = k(numnodes[2] + 1) / (1 * std::pow(Inverse_nodes(2), 2));
            BB3_i = -((k(numnodes[2] + 2)) + 3 * k(numnodes[2] + 1)) / (2 * std::pow(Inverse_nodes(2), 2)) - (nodes[numnodes[2]].rho(T(numnodes[2] - 1)) * nodes[numnodes[2]].Cp(T(numnodes[2] - 1)) / delt) * Layer_thickness(2) * Layer_thickness(2);
            CC3_i = ((k(numnodes[2] + 1) + k(numnodes[2] + 2))) / (2 * std::pow(Inverse_nodes(2), 2));
            DD3_i = -((nodes[numnodes[2]].rho(T(numnodes[2] - 1)) * nodes[numnodes[2]].Cp(T(numnodes[2] - 1))) / delt) * T(numnodes[2] - 1) * Layer_thickness(2) * Layer_thickness(2) - SS(numnodes[2] - 1) * Layer_thickness(2) * Layer_thickness(2);

            if (m_options.coordinateType != CoordinateType::Cartesian)
            {

                if (Outer_radius_instantaneous <= 0.0000001)
                {
                    CYSP2(AA3_i, BB3_i, CC3_i, k(numnodes[2] + 1), Inverse_nodes(2));
                }
                else if (Outer_radius_instantaneous > 0.0000001)
                {
                    CYSP1(AA3_i, CC3_i, k(numnodes[2] + 1), Layer_thickness(2), Inverse_nodes(2), Outer_radius_instantaneous);
                }
            }

            AA2_i = -k(numnodes[1] - 1) / (2 * Inverse_nodes(0) * Layer_thickness(0));
            AA2D_i = k(numnodes[2] + 1) / (2 * Inverse_nodes(2) * Layer_thickness(2));
            CC2_i = -AA2D_i;
            CC2D_i = -AA2_i;
            BB2 = 0;
            DD2 = 0;

            A(numnodes[1] - 1, numnodes[1] - 2) = AA1_i * CC2D_i * AA3_i - CC1_i * AA2_i * AA3_i;
            A(numnodes[1] - 1, numnodes[1] - 1) = BB3_i * CC1_i * AA2D_i + BB1_i * CC2D_i * AA3_i - BB2 * CC1_i * AA3_i;
            A(numnodes[1] - 1, numnodes[1]) = CC1_i * AA2D_i * CC3_i - AA3_i * CC1_i * CC2_i;
            B(numnodes[1] - 1, 0) = DD3_i * CC1_i * AA2D_i + AA3_i * CC2D_i * DD1_i - DD2 * CC1_i * AA3_i;

            for (int io = numnodes[1]; io < numnodes[2]; io++)
            {
                A(io, io - 1) = 1;
                A(io, io) = -1;
                A(io, io + 1) = 0;
                B(io, 0) = 0;
            }
        }
    }

    if (Layers.GetCount() >= 2)
    {

        if (Layer_thickness(1) != 0)
        {
            if (MPD != 0)
            {
                A(numnodes[1] - 1, numnodes[1] - 2) = 0;
                A(numnodes[1] - 1, numnodes[1] - 1) = 1;
                A(numnodes[1] - 1, numnodes[1]) = 0;
                B(numnodes[1] - 1, 0) = nodes[numnodes[1] - 1].Tpyro();
            }
            else if (MPD == 0 && Layer_thickness(0) == 0)
            {
                AA2_s = k(numnodes[1]) / (std::pow(Inverse_nodes(1), 2));
                BB2_s = -(3 * k(numnodes[1]) + k(numnodes[1] + 1)) / (2 * std::pow(Inverse_nodes(1), 2)) - ((nodes[numnodes[1]].rho(T(numnodes[1] - 1)) * nodes[numnodes[1]].Cp(T(numnodes[1] - 1))) / delt) * Layer_thickness(1) * Layer_thickness(1);
                CC2_s = (k(numnodes[1]) + k(numnodes[1] + 1)) / (2 * std::pow(Inverse_nodes(1), 2));
                DD2_s = -((nodes[numnodes[1]].rho(T(numnodes[1] - 1)) * nodes[numnodes[1]].Cp(T(numnodes[1] - 1))) / delt) * T(numnodes[1] - 1) * Layer_thickness(1) * Layer_thickness(1) - SS(numnodes[1] - 1) * Layer_thickness(1) * Layer_thickness(1);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {

                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(AA2_s, BB2_s, CC2_s, k(numnodes[1]), Inverse_nodes(1));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(AA2_s, CC2_s, k(numnodes[1]), Layer_thickness(1), Inverse_nodes(1), Outer_radius_instantaneous);
                    }
                }

                AA1_s = k(numnodes[1]) / (2 * Inverse_nodes(1) * Layer_thickness(1));
                BB1_s = std::reduce(h_front.begin(), h_front.end(), heat_transfer_coefficient_outer);
                CC1_s = -k(numnodes[1] + 1) / (2 * Inverse_nodes(1) * Layer_thickness(1));
                DD1_s = (std::inner_product(h_front.begin(), h_front.end(), Tg_front.begin(), heat_transfer_coefficient_outer * Tg) - nodes[numnodes[1]].emissivity(T(numnodes[1] - 1)) * sigma * std::transform_reduce(Tamb_front.begin(), Tamb_front.end(), 0., std::plus<double>(), [&](double t)
                                                                                                                                                                                                                        { return std::pow(T(numnodes[1] - 1), 4) - std::pow(t, 4); }) +
                         +Flux_front);

                A(numnodes[1] - 1, numnodes[1] - 2) = 0;
                A(numnodes[1] - 1, numnodes[1] - 1) = BB1_s * AA2_s - AA1_s * BB2_s;
                A(numnodes[1] - 1, numnodes[1]) = AA2_s * CC1_s - AA1_s * CC2_s;
                B(numnodes[1] - 1, 0) = AA2_s * DD1_s - AA1_s * DD2_s;
            }

            for (int jojo = numnodes[1]; jojo < numnodes[2] - 1; jojo++)
            {

                double Factor_1, Factor_2;
                A(jojo, jojo - 1) = (k(jojo) + k(jojo + 1)) / (2 * std::pow(Inverse_nodes(1), 2));
                A(jojo, jojo) = -(k(jojo + 2) + 2 * k(jojo + 1) + k(jojo)) / (2 * std::pow(Inverse_nodes(1), 2)) - ((nodes[jojo].rho(T(jojo)) * nodes[jojo].Cp(T(jojo))) / delt) * Layer_thickness(1) * Layer_thickness(1);
                A(jojo, jojo + 1) = (k(jojo + 1) + k(jojo + 2)) / (2 * std::pow(Inverse_nodes(1), 2));
                Factor_1 = ((Layer_thickness(1)) * (MPD * (numnodes[2] - 1 - jojo) * nodes[jojo].rho(T(jojo)) * nodes[jojo].Cp(T(jojo)))) / (nodes[jojo].rhovirgin() - nodes[jojo].rhochar());
                Factor_2 = Factor_1 / 2;
                if (A(jojo, jojo - 1) < Factor_2)
                {
                    A(jojo, jojo) = A(jojo, jojo) - Factor_1;
                    A(jojo, jojo + 1) = A(jojo, jojo + 1) + Factor_1;
                }

                else
                {
                    A(jojo, jojo - 1) = A(jojo, jojo - 1) - Factor_2;
                    A(jojo, jojo + 1) = A(jojo, jojo + 1) + Factor_2;
                }

                B(jojo, 0) = -((nodes[jojo].rho(T(jojo)) * nodes[jojo].Cp(T(jojo))) / delt) * T(jojo) * Layer_thickness(1) * Layer_thickness(1) - SS(jojo) * Layer_thickness(1) * Layer_thickness(1);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {
                    Outer_radius_instantaneous = Outer_radius_instantaneous - Layer_thickness(1) * Inverse_nodes(1);
                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(A(jojo, jojo - 1), A(jojo, jojo), A(jojo, jojo + 1), k(jojo + 1), Inverse_nodes(1));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(A(jojo, jojo - 1), A(jojo, jojo + 1), k(jojo + 1), Layer_thickness(1), Inverse_nodes(1), Outer_radius_instantaneous);
                    }
                }
            }

            Outer_radius_instantaneous = Outer_radius_instantaneous - Layer_thickness(1) * Inverse_nodes(1);

            if (Layers.GetCount() < 3)
            {
                double AA2_b, BB2_b, CC2_b, DD2_b, AA1_b, BB1_b, CC1_b, DD1_b;
                AA2_b = (k(numnodes[2] - 1) + k(numnodes[2])) / (2 * std::pow(Inverse_nodes(1), 2));
                BB2_b = -(k(numnodes[2] - 1) + 3 * k(numnodes[2])) / (2 * std::pow(Inverse_nodes(1), 2)) - (nodes[numnodes[2] - 1].rho(T(numnodes[2] - 1)) * nodes[numnodes[2] - 1].Cp(T(numnodes[2] - 1)) / delt) * Layer_thickness(1) * Layer_thickness(1);
                CC2_b = k(numnodes[2]) / (1 * pow(Inverse_nodes(1), 2));
                DD2_b = -(nodes[numnodes[2] - 1].rho(T(numnodes[2] - 1)) * nodes[numnodes[2] - 1].Cp(T(numnodes[2] - 1)) / delt) * T(numnodes[2] - 1) * Layer_thickness(1) * Layer_thickness(1) - SS(numnodes[2] - 1) * Layer_thickness(1) * Layer_thickness(1);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {

                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(AA2_b, BB2_b, CC2_b, k(numnodes[2]), Inverse_nodes(1));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(AA2_b, CC2_b, k(numnodes[2]), Layer_thickness(1), Inverse_nodes(1), Outer_radius_instantaneous);
                    }
                }

                for (size_t i = 0; i < l_back.size(); i++)
                {
                    h_back.push_back(FreeConvection_Surface(T(numnodes[2] - 1), l_tg_back[i], l_back[i], m_params.Sea_Level_Pressure));
                    Tg_back.push_back(l_tg_back[i]);
                }

                std::vector<double> h_radiation_loss_backwall;
                for (size_t i = 0; i < Tamb_back.size(); i++)
                {
                    h_radiation_loss_backwall.push_back(sigma * nodes[numnodes[2] - 1].emissivity(numnodes[2] - 1) * (std::pow(T(numnodes[2] - 1), 3) + std::pow(T(numnodes[2] - 1), 2) * Tamb_back[i] + T(numnodes[2] - 1) * std::pow(Tamb_back[i], 2) + std::pow(Tamb_back[i], 3)));
                }

                AA1_b = k(numnodes[2]) / (2 * Layer_thickness(1) * Inverse_nodes(1));
                BB1_b = std::reduce(h_back.begin(), h_back.end(), 0.) + std::reduce(h_radiation_loss_backwall.begin(), h_radiation_loss_backwall.end(), 0.) + Propellant_Mass / delt;
                CC1_b = -k(numnodes[2]) / (2 * Layer_thickness(1) * Inverse_nodes(1));
                DD1_b = Flux_back + std::inner_product(h_back.begin(), h_back.end(), Tg_back.begin(), 0.) + std::inner_product(h_radiation_loss_backwall.begin(), h_radiation_loss_backwall.end(), Tamb_back.begin(), 0.) +
                        (Propellant_Mass * T(numnodes[2] - 1) / delt);

                A(numnodes[2] - 1, numnodes[2] - 2) = CC1_b * AA2_b - CC2_b * AA1_b;
                A(numnodes[2] - 1, numnodes[2] - 1) = CC1_b * BB2_b - CC2_b * BB1_b;
                B(numnodes[2] - 1, 0) = CC1_b * DD2_b - CC2_b * DD1_b;
            }

            else if (Layers.GetCount() >= 3)
            {
                double AA1_i, BB1_i, CC1_i, DD1_i, AA3_i, BB3_i, CC3_i, DD3_i, AA2_i, AA2D_i, CC2_i, CC2D_i, BB2, DD2;

                AA1_i = (k(numnodes[2] - 1) + k(numnodes[2])) / (2 * std::pow(Inverse_nodes(1), 2));
                BB1_i = -(k(numnodes[2] - 1) + 3 * k(numnodes[2])) / (2 * std::pow(Inverse_nodes(1), 2)) - (nodes[numnodes[2] - 1].rho(T(numnodes[2] - 1)) * nodes[numnodes[2] - 1].Cp(T(numnodes[2] - 1)) / delt) * Layer_thickness(1) * Layer_thickness(1);
                CC1_i = k(numnodes[2]) / (1 * pow(Inverse_nodes(1), 2));
                DD1_i = -(nodes[numnodes[2] - 1].rho(T(numnodes[2] - 1)) * nodes[numnodes[2] - 1].Cp(T(numnodes[2] - 1)) / delt) * T(numnodes[2] - 1) * Layer_thickness(1) * Layer_thickness(1) - SS(numnodes[2] - 1) * Layer_thickness(1) * Layer_thickness(1);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {

                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(AA1_i, BB1_i, CC1_i, k(numnodes[2]), Inverse_nodes(1));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(AA1_i, CC1_i, k(numnodes[2]), Layer_thickness(1), Inverse_nodes(1), Outer_radius_instantaneous);
                    }
                }

                AA3_i = (k(numnodes[2] + 1)) / (1 * std::pow(Inverse_nodes(2), 2));
                BB3_i = -(k(numnodes[2] + 2) + 3 * k(numnodes[2] + 1)) / (2 * std::pow(Inverse_nodes(2), 2)) - (nodes[numnodes[2]].rho(T(numnodes[2] - 1)) * nodes[numnodes[2]].Cp(T(numnodes[2] - 1)) / delt) * Layer_thickness(2) * Layer_thickness(2);
                CC3_i = (k(numnodes[2] + 1) + k(numnodes[2] + 2)) / (2 * std::pow(Inverse_nodes(2), 2));
                DD3_i = -((nodes[numnodes[2]].rho(T(numnodes[2] - 1)) * nodes[numnodes[2]].Cp(T(numnodes[2] - 1)) / delt) * T(numnodes[2] - 1) * Layer_thickness(2) * Layer_thickness(2)) - SS(numnodes[2] - 1) * Layer_thickness(2) * Layer_thickness(2);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {

                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(AA3_i, BB3_i, CC3_i, k(numnodes[2] + 1), Inverse_nodes(2));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(AA3_i, CC3_i, k(numnodes[2] + 1), Layer_thickness(2), Inverse_nodes(2), Outer_radius_instantaneous);
                    }
                }

                AA2_i = -(k(numnodes[2])) / (2 * Inverse_nodes(1) * Layer_thickness(1));
                AA2D_i = k(numnodes[2] + 1) / (2 * Inverse_nodes(2) * Layer_thickness(2));
                CC2_i = -AA2D_i;
                CC2D_i = -AA2_i;

                A(numnodes[2] - 1, numnodes[2] - 2) = AA1_i * CC2D_i * AA3_i - CC1_i * AA2_i * AA3_i;
                A(numnodes[2] - 1, numnodes[2] - 1) = BB3_i * CC1_i * AA2D_i + BB1_i * CC2D_i * AA3_i;
                A(numnodes[2] - 1, numnodes[2]) = CC1_i * AA2D_i * CC3_i - AA3_i * CC1_i * CC2_i;
                B(numnodes[2] - 1, 0) = DD3_i * CC1_i * AA2D_i + AA3_i * CC2D_i * DD1_i;
            }
        }
    }

    if (Layers.GetCount() >= 3)
    {
        if (MPD == 0 && Layer_thickness(0) == 0 && Layer_thickness(1) == 0 && Layers.GetCount() >= 3)
        {
            AA2_s = k(numnodes[2] + 1) / (std::pow(Inverse_nodes(2), 2));
            BB2_s = -(3 * k(numnodes[2] + 1) + k(numnodes[2] + 2)) / (2 * std::pow(Inverse_nodes(2), 2)) - (((nodes[numnodes[2]].rho(T(numnodes[2] - 1)) * nodes[numnodes[2]].Cp(T(numnodes[2] - 1))) / delt) * Layer_thickness(2) * Layer_thickness(2));
            CC2_s = (k(numnodes[2] + 1) + k(numnodes[2] + 2)) / (2 * std::pow(Inverse_nodes(2), 2));
            DD2_s = -((nodes[numnodes[2]].rho(T(numnodes[2] - 1)) * nodes[numnodes[2]].Cp(T(numnodes[2] - 1)) / delt) * T(numnodes[2] - 1) * Layer_thickness(2) * Layer_thickness(2)) - SS(numnodes[2] - 1) * Layer_thickness(2) * Layer_thickness(2);

            if (m_options.coordinateType != CoordinateType::Cartesian)
            {

                if (Outer_radius_instantaneous <= 0.0000001)
                {
                    CYSP2(AA2_s, BB2_s, CC2_s, k(numnodes[2] + 1), Inverse_nodes(2));
                }
                else if (Outer_radius_instantaneous > 0.0000001)
                {
                    CYSP1(AA2_s, CC2_s, k(numnodes[2] + 1), Layer_thickness(2), Inverse_nodes(2), Outer_radius_instantaneous);
                }
            }

            AA1_s = (k(numnodes[2] + 1)) / (2 * Inverse_nodes(2) * Layer_thickness(2));
            BB1_s = std::reduce(h_front.begin(), h_front.end(), heat_transfer_coefficient_outer);
            CC1_s = -(k(numnodes[2] + 1)) / (2 * Inverse_nodes(2) * Layer_thickness(2));
            DD1_s = (std::inner_product(h_front.begin(), h_front.end(), Tg_front.begin(), heat_transfer_coefficient_outer * Tg) - nodes[numnodes[2]].emissivity(T(numnodes[2] - 1)) * sigma * std::transform_reduce(Tamb_front.begin(), Tamb_front.end(), 0., std::plus<double>(), [&](double t)
                                                                                                                                                                                                                    { return std::pow(T(numnodes[2] - 1), 4) - std::pow(t, 4); }) +
                     +Flux_front);
            A(numnodes[2] - 1, numnodes[2] - 2) = 0;
            A(numnodes[2] - 1, numnodes[2] - 1) = BB1_s * AA2_s - AA1_s * BB2_s;
            A(numnodes[2] - 1, numnodes[2]) = AA2_s * CC1_s - AA1_s * CC2_s;
            B(numnodes[2] - 1, 0) = AA2_s * DD1_s - AA1_s * DD2_s;
        }

        ////// Before this correct numnodes usage

        // correct j+1 here
        for (int j = 3; j <= Layers.GetCount(); j++)
        {

            for (int ioi = numnodes[j - 1]; ioi < numnodes[j] - 1; ioi++)
            {
                //			A(ioi, ioi - 1) = (k(ioi + j - 3) + k(ioi + j - 2)) / (2 * pow(Inverse_nodes(j - 1) * Layer_thickness(j - 1), 2));
                //			A(ioi, ioi) = -(k(ioi + j - 3) + 2 * k(ioi + j - 2) + k(ioi + j - 1)) / (2 * pow(Inverse_nodes(j - 1) * Layer_thickness(j - 1), 2)) - ((rho(ioi + j - 2) * Cp(ioi + j - 2)) / delt);
                //			A(ioi, ioi + 1) = (k(ioi + j - 2) + k(ioi + j - 1)) / (2 * pow(Inverse_nodes(j - 1) * Layer_thickness(j - 1), 2));
                //			B(ioi, 0) = -((rho(ioi + j - 2) * Cp(ioi + j - 2)) / delt) * T(ioi) - SS(ioi);

                A(ioi, ioi - 1) = (k(ioi + j - 2) + k(ioi + j - 1)) / (2 * std::pow(Inverse_nodes(j - 1), 2));
                A(ioi, ioi) = -(k(ioi + j - 2) + 2 * k(ioi + j - 1) + k(ioi + j)) / (2 * std::pow(Inverse_nodes(j - 1), 2)) - ((nodes[ioi].rho(T(ioi)) * nodes[ioi].Cp(T(ioi))) / delt) * std::pow(Layer_thickness(j - 1), 2);
                A(ioi, ioi + 1) = (k(ioi + j - 1) + k(ioi + j)) / (2 * std::pow(Inverse_nodes(j - 1), 2));
                B(ioi, 0) = -((nodes[ioi].rho(T(ioi)) * nodes[ioi].Cp(T(ioi))) / delt) * std::pow(Layer_thickness(j - 1), 2) * T(ioi) - SS(ioi) * std::pow(Layer_thickness(j - 1), 2);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {
                    Outer_radius_instantaneous = Outer_radius_instantaneous - Layer_thickness(j - 1) * Inverse_nodes(j - 1);
                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(A(ioi, ioi - 1), A(ioi, ioi), A(ioi, ioi + 1), k(ioi + j - 1), Inverse_nodes(j - 1));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(A(ioi, ioi - 1), A(ioi, ioi + 1), k(ioi + j - 1), Layer_thickness(j - 1), Inverse_nodes(j - 1), Outer_radius_instantaneous);
                    }
                }
                //
            }

            Outer_radius_instantaneous = Outer_radius_instantaneous - Layer_thickness(j - 1) * Inverse_nodes(j - 1);

            if (j != Layers.GetCount())
            {
                double AA1_i, BB1_i, CC1_i, DD1_i, AA3_i, BB3_i, CC3_i, DD3_i, AA2_i, AA2D_i, CC2_i, CC2D_i, BB2, DD2;

                AA1_i = (k(numnodes[j] + j - 2) + k(numnodes[j] + j - 3)) / (2 * std::pow(Inverse_nodes(j - 1), 2));
                BB1_i = -(k(numnodes[j] + j - 3) + 3 * k(numnodes[j] + j - 2)) / (2 * std::pow(Inverse_nodes(j - 1), 2)) - ((nodes[numnodes[j] - 1].rho(T(numnodes[j] - 1)) * nodes[numnodes[j] - 1].Cp(T(numnodes[j] - 1))) / delt) * std::pow(Layer_thickness(j - 1), 2);
                CC1_i = k(numnodes[j] + j - 2) / (std::pow(Inverse_nodes(j - 1), 2));
                DD1_i = -((nodes[numnodes[j] - 1].rho(T(numnodes[j] - 1)) * nodes[numnodes[j] - 1].Cp(T(numnodes[j] - 1))) / delt) * std::pow(Layer_thickness(j - 1), 2) * T(numnodes[j] - 1) - SS(numnodes[j] - 1) * std::pow(Layer_thickness(j - 1), 2);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {

                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(AA1_i, BB1_i, CC1_i, k(numnodes[j] + j - 2), Inverse_nodes(j - 1));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(AA1_i, CC1_i, k(numnodes[j] + j - 2), Layer_thickness(j - 1), Inverse_nodes(j - 1), Outer_radius_instantaneous);
                    }
                }

                AA3_i = k(numnodes[j] + j - 1) / (std::pow(Inverse_nodes(j), 2));
                BB3_i = -(k(numnodes[j] + j) + 3 * k(numnodes[j] + j - 1)) / (2 * (std::pow(Inverse_nodes(j), 2))) - ((nodes[numnodes[j]].rho(T(numnodes[j] - 1)) * nodes[numnodes[j]].Cp(T(numnodes[j] - 1))) / delt) * std::pow(Layer_thickness(j), 2);
                CC3_i = (k(numnodes[j] + j - 1) + k(numnodes[j] + j)) / (2 * (std::pow(Inverse_nodes(j), 2)));
                DD3_i = -((nodes[numnodes[j]].rho(T(numnodes[j] - 1)) * nodes[numnodes[j]].Cp(T(numnodes[j] - 1))) / delt) * std::pow(Layer_thickness(j), 2) * T(numnodes[j] - 1) - SS(numnodes[j] - 1) * std::pow(Layer_thickness(j), 2);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {

                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(AA3_i, BB3_i, CC3_i, k(numnodes[j] + j - 1), Inverse_nodes(j));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(AA3_i, CC3_i, k(numnodes[j] + j - 1), Layer_thickness(j), Inverse_nodes(j), Outer_radius_instantaneous);
                    }
                }

                AA2_i = -(k(numnodes[j] + j - 2) / (2 * Inverse_nodes(j - 1) * Layer_thickness(j - 1)));
                AA2D_i = k(numnodes[j] + j - 1) / (2 * Inverse_nodes(j) * Layer_thickness(j));
                CC2_i = -AA2D_i;
                CC2D_i = -AA2_i;

                A(numnodes[j] - 1, numnodes[j] - 2) = AA1_i * CC2D_i * AA3_i - CC1_i * AA2_i * AA3_i;
                A(numnodes[j] - 1, numnodes[j] - 1) = BB3_i * CC1_i * AA2D_i + BB1_i * CC2D_i * AA3_i;
                A(numnodes[j] - 1, numnodes[j]) = CC1_i * AA2D_i * CC3_i - AA3_i * CC1_i * CC2_i;
                B(numnodes[j] - 1, 0) = DD3_i * CC1_i * AA2D_i + AA3_i * CC2D_i * DD1_i;
            }

            if (j == Layers.GetCount())
            {
                double AA2_b, BB2_b, CC2_b, DD2_b, AA1_b, BB1_b, CC1_b, DD1_b;
                //			AA2_b = (k(No_of_nodes.sum() + No_of_Layers - 1) + k(No_of_nodes.sum() + No_of_Layers - 2)) / (2 * pow(Inverse_nodes(No_of_Layers - 1) * Layer_thickness(No_of_Layers - 1), 2));
                //			BB2_b = -(3 * k(No_of_nodes.sum() + No_of_Layers - 1) + k(No_of_nodes.sum() + No_of_Layers - 2)) / (2 * pow(Inverse_nodes(No_of_Layers - 1) * Layer_thickness(No_of_Layers - 1), 2)) - (rho(No_of_nodes.sum() + No_of_Layers - 1) * Cp(No_of_nodes.sum() + No_of_Layers - 1)) / delt;
                //			CC2_b = k(No_of_nodes.sum() + No_of_Layers - 1) / (1 * pow(Inverse_nodes(No_of_Layers - 1) * Layer_thickness(No_of_Layers - 1), 2));
                //			DD2_b = -(rho(No_of_nodes.sum() + No_of_Layers - 1) * Cp(No_of_nodes.sum() + No_of_Layers - 1) / delt) * T(No_of_nodes.sum()) - SS(No_of_nodes.sum());

                AA2_b = (k(numnodes[j] + j - 2) + k(numnodes[j] + j - 3)) / (2 * std::pow(Inverse_nodes(j - 1), 2));
                BB2_b = -(k(numnodes[j] + j - 3) + 3 * k(numnodes[j] + j - 2)) / (2 * std::pow(Inverse_nodes(j - 1), 2)) - ((nodes[numnodes[j] - 1].rho(T(numnodes[j] - 1)) * nodes[numnodes[j] - 1].Cp(T(numnodes[j] - 1))) / delt) * std::pow(Layer_thickness(j - 1), 2);
                CC2_b = k(numnodes[j] + j - 2) / (std::pow(Inverse_nodes(j - 1), 2));
                DD2_b = -((nodes[numnodes[j] - 1].rho(T(numnodes[j] - 1)) * nodes[numnodes[j] - 1].Cp(T(numnodes[j] - 1))) / delt) * std::pow(Layer_thickness(j - 1), 2) * T(numnodes[j] - 1) - SS(numnodes[j] - 1) * std::pow(Layer_thickness(j - 1), 2);

                if (m_options.coordinateType != CoordinateType::Cartesian)
                {

                    if (Outer_radius_instantaneous <= 0.0000001)
                    {
                        CYSP2(AA2_b, BB2_b, CC2_b, k(numnodes[j] + j - 2), Inverse_nodes(j - 1));
                    }
                    else if (Outer_radius_instantaneous > 0.0000001)
                    {
                        CYSP1(AA2_b, CC2_b, k(numnodes[j] + j - 2), Layer_thickness(j - 1), Inverse_nodes(j - 1), Outer_radius_instantaneous);
                    }
                }

                for (size_t i = 0; i < l_back.size(); i++)
                {
                    h_back.push_back(FreeConvection_Surface(T(numnodes[j] - 1), l_tg_back[i], l_back[i], m_params.Sea_Level_Pressure));
                    Tg_back.push_back(l_tg_back[i]);
                }

                std::vector<double> h_radiation_loss_backwall;
                for (size_t i = 0; i < Tamb_back.size(); i++)
                {
                    h_radiation_loss_backwall.push_back(sigma * nodes[numnodes[j] - 1].emissivity(numnodes[j] - 1) * (std::pow(T(numnodes[j] - 1), 3) + std::pow(T(numnodes[j] - 1), 2) * Tamb_back[i] + T(numnodes[j] - 1) * std::pow(Tamb_back[i], 2) + std::pow(Tamb_back[i], 3)));
                }

                AA1_b = k(numnodes[j] + Layers.GetCount() - 2) / (2 * Layer_thickness(Layers.GetCount() - 1) * Inverse_nodes(Layers.GetCount() - 1));
                BB1_b = std::reduce(h_back.begin(), h_back.end(), 0.) + std::reduce(h_radiation_loss_backwall.begin(), h_radiation_loss_backwall.end(), 0.) + Propellant_Mass / delt;
                CC1_b = -k(numnodes[j] + Layers.GetCount() - 2) / (2 * Layer_thickness(Layers.GetCount() - 1) * Inverse_nodes(Layers.GetCount() - 1));
                DD1_b = Flux_back + std::inner_product(h_back.begin(), h_back.end(), Tg_back.begin(), 0.) + std::inner_product(h_radiation_loss_backwall.begin(), h_radiation_loss_backwall.end(), Tamb_back.begin(), 0.) +
                        (Propellant_Mass * T(numnodes[j] - 1) / delt);

                A(numnodes[j] - 1, numnodes[j] - 2) = CC1_b * AA2_b - CC2_b * AA1_b;
                A(numnodes[j] - 1, numnodes[j] - 1) = CC1_b * BB2_b - CC2_b * BB1_b;
                B(numnodes[j] - 1, 0) = CC1_b * DD2_b - CC2_b * DD1_b;
            }
        }
    }
// #define DEBUG_PRINT
#ifdef DEBUG_PRINT
    std::ofstream outfile("C:/Users/Arnab Mahanti/source/repos/otap_pp/docs/matrix.csv");
    outfile << "A:\n";
    outfile << A;
    outfile << "\n\nB:\n";
    outfile << B;
    outfile << "\n\nT:\n";
    outfile << T;
#endif

    if (Layers.GetCount() < 2)
    {
        if (Layer_thickness(0) != 0)
        {
            MatrixXd A1, B1;
            VectorXd T1;
            A1 = MatrixXd::Zero(numnodes[Layers.GetCount()], Layers.GetCount());
            B1 = MatrixXd::Zero(numnodes[Layers.GetCount()], 1);
            T1 = VectorXd::Zero(numnodes[Layers.GetCount()]);
            A1 = A;
            B1 = B;
            T1 = A1.colPivHouseholderQr().solve(B1);
            // T1=A1.lu().solve(B1);
            T = T1;
            // cout << A1 << endl;
            // cout << T1 << endl;
        }
    }
    if (Layers.GetCount() >= 2)
    {
        if (Layer_thickness(0) != 0)
        {
            MatrixXd A1, B1;
            VectorXd T1;
            A1 = MatrixXd::Zero(numnodes[Layers.GetCount()], numnodes[Layers.GetCount()]);
            B1 = MatrixXd::Zero(numnodes[Layers.GetCount()], 1);
            T1 = VectorXd::Zero(numnodes[Layers.GetCount()]);
            A1 = A;
            B1 = B;
            T1 = A1.colPivHouseholderQr().solve(B1);

            T = T1;
        }
        else if (Layer_thickness(0) == 0 && Layer_thickness(1) != 0)
        {
            MatrixXd A1, B1;
            VectorXd T1;
            A1 = MatrixXd::Zero(numnodes[Layers.GetCount()] - numnodes[1] + 1, numnodes[Layers.GetCount()] - numnodes[1] + 1);
            B1 = MatrixXd::Zero(numnodes[Layers.GetCount()] - numnodes[1] + 1, 1);

            T1 = VectorXd::Zero(numnodes[Layers.GetCount()] - numnodes[1] + 1, 1);
            A1 = A.block(numnodes[1] - 1, numnodes[1] - 1, numnodes[Layers.GetCount()] - numnodes[1] + 1, numnodes[Layers.GetCount()] - numnodes[1] + 1);
            B1 = B.block(numnodes[1] - 1, 0, numnodes[Layers.GetCount()] - numnodes[1] + 1, 1);

            T1 = A1.colPivHouseholderQr().solve(B1);
            // cout << T1 << endl;
            T.segment(numnodes[1] - 1, numnodes[Layers.GetCount()] - numnodes[1] + 1) = T1;
        }
        else
        {
            MatrixXd A1, B1;
            VectorXd T1;
            A1 = MatrixXd::Zero(numnodes[Layers.GetCount()] - numnodes[2] + 1, numnodes[Layers.GetCount()] - numnodes[2] + 1);
            B1 = MatrixXd::Zero(numnodes[Layers.GetCount()] - numnodes[2] + 1, 1);
            T1 = VectorXd::Zero(numnodes[Layers.GetCount()] - numnodes[2] + 1, 1);
            A1 = A.block(numnodes[2] - 1, numnodes[2] - 1, numnodes[Layers.GetCount()] - numnodes[2] + 1, numnodes[Layers.GetCount()] - numnodes[2] + 1);
            B1 = B.block(numnodes[2] - 1, 0, numnodes[Layers.GetCount()] - numnodes[2] + 1, 1);
            T1 = A1.colPivHouseholderQr().solve(B1);
            T.segment(numnodes[2] - 1, numnodes[Layers.GetCount()] - numnodes[2] + 1) = T1;

            // std::ofstream outfile("C:/Users/Arnab Mahanti/source/repos/otap_pp/docs/matrix.csv", std::ios::app);
            // outfile << "A1:\n";
            // outfile << A1;
            // outfile << "\n\nB1:\n";
            // outfile << B1;
            // outfile << "\n\nT1:\n";
            // outfile << T1;
        }
    }
    // m_ContinuumSolver->Solve(300);
    // return ResponseResult();
}

OTAP::ResponseResult OTAP::DefaultResponseSolver::Solve()
{
    using namespace Eigen;

    auto Layers = *m_Layers;
    int TIME_STEP_REDUCED_ITERATION = 0;
    int INNER_ITER = 0;
    int CONVERGENCE_CONDITIONS_MET = 0;
    int IXFIST = 0;
    int GONE_IN_TIME_STEP_REDUCTION = 0;
    int NODES_1;
    VectorXd Layer_thickness, Initial_Layer_thickness;
    Layer_thickness = VectorXd::Zero(Layers.GetCount());
    Initial_Layer_thickness = VectorXd::Zero(Layers.GetCount());
    ResponseResult result;

    for (size_t i = 0; i < Layers.GetCount(); i++)
    {

        Layer_thickness(i) = Layers[i].thickness;
        Initial_Layer_thickness(i) = Layers[i].thickness;
    }

    auto nodes = m_Layers->Nodes;
    auto numnodes = m_Layers->NumNodes;
    VectorXd T(numnodes[Layers.GetCount()]);
    VectorXd k(numnodes[Layers.GetCount()] + Layers.GetCount() - 1);
    // Assign initial temperature;
    for (int g = 0; g < numnodes[Layers.GetCount()]; g++)
    {
        T(g) = nodes[g].T;
    }

    double time;
    double time_store;
    double delt, delt_store;
    time = m_params.tInit;

    double Mass_rate_Pyrolysis = m_params.massRatePyro;
    double Mass_rate_Char = m_params.massRateChar;
    std::vector<double> solution_t;
    std::vector<std::vector<double>> solution_T;
    // Time Loop
    solution_t.push_back(time);
    solution_T.emplace_back(T.begin(), T.end());
    while (time < m_params.tFinal)
    {
        // delt=m_params[time];  //When Table option is incorporated
        delt = m_params.timestep[time];
        time = time + delt;
        delt_store = delt;

        TIME_STEP_REDUCED_ITERATION = 0;
        INNER_ITER = 0;
        CONVERGENCE_CONDITIONS_MET = 0;

        while (TIME_STEP_REDUCED_ITERATION <= 10 && (INNER_ITER > 10 || CONVERGENCE_CONDITIONS_MET == 0))
        {
            //				//INNER_ITER = 0;
            //
            while (INNER_ITER <= 10 && CONVERGENCE_CONDITIONS_MET == 0)
            //
            {
                INNER_ITER = INNER_ITER + 1;
                //
                if ((T(numnodes[1] - 1) > nodes[numnodes[1] - 1].Tpyro() - m_params.Ttol) && IXFIST == 0)
                {
                    if (!m_options.Sublime)
                    {
                        for (int ijk = 0; ijk < numnodes[1] - 1; ijk++)
                        {
                            NODES_1 = numnodes[1] - 1 - 1 - ijk;
                            T(NODES_1) = T(NODES_1 + 1);
                        }
                        IXFIST = 1;
                    }

                    else
                    {
                        IXFIST = 1;
                    }
                }
                // Air gap

                double h_convection_air_gap = 0;
                double h_radiation_air_gap = 0;

                for (int i = 0; i < Layers.GetCount(); i++)
                {
                    if (Layers[i].material->name == "air") // Ask
                    {
                        Air_Gap(h_convection_air_gap, h_radiation_air_gap, Layers[i].thickness, m_params.Air_Layer_Length, T(numnodes[i] - 1), T(numnodes[i + 1] - 1), m_params.Air_Layer_Length, nodes[numnodes[i] - 1].emissivity(T(numnodes[i] - 1)), nodes[numnodes[i + 1] - 1].emissivity(T(numnodes[i + 1] - 1)));
                    }
                }
                for (int nodes0_j = 0; nodes0_j <= numnodes[1] - 1; nodes0_j++)
                {
                    k(nodes0_j) = nodes[nodes0_j].k(T(nodes0_j));
                }

                if (Layers.GetCount() > 1)
                {

                    for (int mat_j = 1; mat_j < Layers.GetCount(); mat_j++)
                    {

                        for (int nodes_j = numnodes[mat_j] - 1 + mat_j; nodes_j <= numnodes[mat_j + 1] - 1 + mat_j; nodes_j++)
                        {
                            if (nodes_j == numnodes[mat_j] - 1 + mat_j)
                            {
                                k(nodes_j) = nodes[nodes_j - mat_j + 1].k(T(nodes_j - mat_j));
                            }
                            else
                            {
                                k(nodes_j) = nodes[nodes_j - mat_j].k(T(nodes_j - mat_j));
                            }

                            if (Layers[mat_j].material->name == "air")
                            {
                                k(nodes_j) = k(nodes_j) + Layer_thickness(mat_j - 1) * (abs(h_convection_air_gap + h_radiation_air_gap));
                            }
                        }
                        // add perturbation terms also here
                    }
                }

                double Twall;
                int Boundary_Node;

                Twall = T(0);
                Boundary_Node = 0;
                if (Layer_thickness(0) == 0)
                {
                    Twall = T(numnodes[1] - 1);
                    Boundary_Node = numnodes[1] - 1;
                }
                else if (Layer_thickness(0) == 0 && Layer_thickness(1) == 0)
                {
                    Twall = T(numnodes[2] - 1);
                    Boundary_Node = numnodes[2] - 1;
                }

                assert(m_HFSolver->GetTrajectory() != nullptr);
                auto res = m_HFSolver->Solve(time, Twall);
                double heat_transfer_coefficient_outer, Tg;
                heat_transfer_coefficient_outer = res.h[0];
                Tg = res.Tg[0];

                BCS BoundaryConditions;

                for (auto &&i : m_BCs)
                {
                    if (i->type == BCType::Radiation)
                    {
                        auto temp = std::dynamic_pointer_cast<RadiationBC>(i);
                        switch (temp->location)
                        {
                        case BCLocation::front:
                            BoundaryConditions.Tambient_front.push_back(temp->Tambient[time]);
                            break;
                        case BCLocation::back:
                            BoundaryConditions.Tambient_back.push_back(temp->Tambient[time]);
                            break;
                        case BCLocation::both:
                            BoundaryConditions.Tambient_front.push_back(temp->Tambient[time]);
                            BoundaryConditions.Tambient_back.push_back(temp->Tambient[time]);
                            break;
                        default:
                            break;
                        }
                    }
                    else if (i->type == BCType::Flux)
                    {
                        auto temp = std::dynamic_pointer_cast<FluxBC>(i);
                        switch (temp->location)
                        {
                        case BCLocation::front:
                            BoundaryConditions.flux_front += temp->flux[time];
                            break;
                        case BCLocation::back:
                            BoundaryConditions.flux_back += temp->flux[time];
                            break;
                        case BCLocation::both:
                            BoundaryConditions.flux_back += temp->flux[time];
                            BoundaryConditions.flux_front += temp->flux[time];
                            break;
                        default:
                            break;
                        }
                    }
                    else if (i->type == BCType::Convection)
                    {
                        auto temp = std::dynamic_pointer_cast<ConvectionBC>(i);
                        switch (temp->location)
                        {
                        case BCLocation::front:
                            BoundaryConditions.h_front.push_back(temp->h[time]);
                            BoundaryConditions.Tg_front.push_back(temp->tg[time]);
                            break;
                        case BCLocation::back:
                            BoundaryConditions.h_back.push_back(temp->h[time]);
                            BoundaryConditions.Tg_back.push_back(temp->tg[time]);
                            break;
                        case BCLocation::both:
                            BoundaryConditions.h_back.push_back(temp->h[time]);
                            BoundaryConditions.h_front.push_back(temp->h[time]);
                            BoundaryConditions.Tg_back.push_back(temp->tg[time]);
                            BoundaryConditions.Tg_front.push_back(temp->tg[time]);
                            break;
                        default:
                            break;
                        }
                    }
                    else if (i->type == BCType::NaturalConvection)
                    {
                        auto temp = std::dynamic_pointer_cast<NaturalConvectionBC>(i);
                        switch (temp->location)
                        {
                        case BCLocation::front:
                            BoundaryConditions.l_front.push_back(temp->l);
                            BoundaryConditions.l_tg_front.push_back(temp->tg[time]);
                            break;
                        case BCLocation::back:
                            BoundaryConditions.l_back.push_back(temp->l);
                            BoundaryConditions.l_tg_back.push_back(temp->tg[time]);
                            break;
                        case BCLocation::both:
                            BoundaryConditions.l_back.push_back(temp->l);
                            BoundaryConditions.l_front.push_back(temp->l);
                            BoundaryConditions.l_tg_back.push_back(temp->tg[time]);
                            BoundaryConditions.l_tg_front.push_back(temp->tg[time]);
                            break;
                        default:
                            break;
                        }
                    }
                    else if (i->type == BCType::PropellantMass)
                    {
                        auto temp = std::dynamic_pointer_cast<PropellantMassBC>(i);
                        switch (temp->location)
                        {
                        case BCLocation::front:
                            BoundaryConditions.propellant_mass_front += temp->propmass[time];
                            break;
                        case BCLocation::back:
                            BoundaryConditions.propellant_mass_back += temp->propmass[time];
                            break;
                        case BCLocation::both:
                            BoundaryConditions.propellant_mass_front += temp->propmass[time];
                            BoundaryConditions.propellant_mass_back += temp->propmass[time];
                            break;
                        default:
                            break;
                        }
                    }
                    else if (i->type == BCType::HeatGeneration)
                    {
                        auto temp = std::dynamic_pointer_cast<HeatGenBC>(i);
                        switch (temp->location)
                        {
                        case BCLocation::front:
                            BoundaryConditions.qgen += temp->qdot[time];
                            break;
                        case BCLocation::back:
                            BoundaryConditions.qgen += temp->qdot[time];
                            break;
                        case BCLocation::both:
                            BoundaryConditions.qgen += temp->qdot[time];
                            BoundaryConditions.qgen += temp->qdot[time];
                            break;
                        default:
                            break;
                        }
                    }
                }

                double Layer_thickness_old_char, Layer_thickness_old_virgin;
                Layer_thickness_old_char = Layer_thickness(0);
                Layer_thickness_old_virgin = Layer_thickness(1);

                double Q_convection_front, Q_radiation_front;
                Q_convection_front = std::transform_reduce(BoundaryConditions.h_front.begin(),
                                                           BoundaryConditions.h_front.end(),
                                                           BoundaryConditions.Tg_front.begin(),
                                                           BoundaryConditions.flux_front + res.q[0],
                                                           std::plus<double>(),
                                                           [&Twall](double h, double tg)
                                                           {
                                                               return h * (tg - Twall);
                                                           });
                Q_radiation_front = 0.;

                if (((Layers.GetCount()) >= 2 && Layer_thickness(0) != 0 && Layer_thickness(1) != 0) || (Mass_rate_Pyrolysis != 0 || T(numnodes[1] - 1) > nodes[numnodes[1] - 1].Tpyro() - 1))
                {

                    Instantaneous_MassRate_ThicknessSolver(Layer_thickness(0), Layer_thickness(1), Mass_rate_Char, Mass_rate_Pyrolysis, Q_convection_front, Q_radiation_front, T, delt, BoundaryConditions);
                }

                if (Layer_thickness(1) < 0.000001)
                {
                    Layer_thickness(1) = 0;
                }

                DefaultResponseMatrix(T, Layer_thickness, Initial_Layer_thickness, heat_transfer_coefficient_outer, Tg, BoundaryConditions, Mass_rate_Char, Mass_rate_Pyrolysis, k, delt);

                if (abs(T(Boundary_Node) - Twall) < m_params.Ttol)
                {
                    CONVERGENCE_CONDITIONS_MET = 1;
                }
                double DX1, DX2, PDX1, PDX2;
                DX1 = abs(Layer_thickness(0) - Layer_thickness_old_char);
                DX2 = abs(Layer_thickness(1) - Layer_thickness_old_virgin);
                PDX1 = 0;
                PDX2 = 0;
                if (Layer_thickness(0) != 0)
                {
                    PDX1 = DX1 / Layer_thickness(0);
                }
                if (Layer_thickness(1) != 0)
                {
                    PDX2 = DX2 / Layer_thickness(1);
                }
                if ((T(0) > nodes[0].Tabl() + 1 && Mass_rate_Char > 0.00001) || (T(numnodes[1] - 1) > nodes[numnodes[1] - 1].Tpyro() + 1 && Mass_rate_Pyrolysis > 0.00001) || (DX1 > 0.000001 || PDX1 > 0.001) || (DX2 > 0.000001 || PDX2 > 0.001))
                {
                    CONVERGENCE_CONDITIONS_MET = 0;
                }
                else
                {
                    CONVERGENCE_CONDITIONS_MET = 1;
                }
                //
                INNER_ITER = INNER_ITER + 1;
            }

            if (INNER_ITER > 10 && CONVERGENCE_CONDITIONS_MET == 0)
            {
                INNER_ITER = 0;
                delt = delt / 2;
                TIME_STEP_REDUCED_ITERATION = TIME_STEP_REDUCED_ITERATION + 1;
                time = time_store + delt;
                GONE_IN_TIME_STEP_REDUCTION = 1;
            }
        }
        if (GONE_IN_TIME_STEP_REDUCTION == 1)
        {
            GONE_IN_TIME_STEP_REDUCTION = 0;
            time_store = time;
            delt = delt_store - delt;
            delt_store = delt;
            time = time_store + delt;
        }
        solution_t.push_back(time);
        solution_T.emplace_back(T.begin(), T.end());
    }

    result.solution_T = solution_T;
    result.solution_t = solution_t;
    for (auto &&i : numnodes)
        i -= 1;
    numnodes[0] = 0;
    result.interface_index = numnodes;
    result.message = "Will see..\n";
    return result;
}
