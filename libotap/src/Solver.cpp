#include "Solver.h"
#include "Eigen/Dense"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include <cmath>

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

// FOR CYLINDRICAL/SPHERICAL COORDINATE TRANSFORMATION (RADIUS != 0)
void OTAP::DefaultResponseSolver::CYSP2(double &A1, double &B1, double &C1, double k, double Inverse_ratio)
{
    double Coordinate_ratio, COORDINATE_COEF;
    if (m_options.coordinateType == CoordinateType::Cartesian)
    {
        COORDINATE_COEF = 0;
    }
    else if (m_options.coordinateType == CoordinateType::Cylindrical)
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
    else if (m_options.coordinateType == CoordinateType::Cylindrical)
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

OTAP::ResponseResult OTAP::DefaultResponseSolver::Solve()

{
    using namespace Eigen;
    double heat_transfer_coefficient_outer, Tg, delt, Qradiation;
    heat_transfer_coefficient_outer = 0;
    Tg = 0;
    delt = 0;
    Qradiation = 0;

    auto nodes = m_Layers->Nodes;
    auto numnodes = m_Layers->NumNodes;

    // (*m_Layers)[0].thickness

    double Total_thickness;
    Total_thickness = 0;
    for (size_t i = 0; i < (*m_Layers).GetCount(); i++)
    {
        Total_thickness = Total_thickness + (*m_Layers)[i].thickness;
    }

    MatrixXd A, B;
    VectorXd T, SS;
    A = MatrixXd::Zero(numnodes[numnodes.size()], numnodes[numnodes.size()]);
    B = MatrixXd::Zero(numnodes[numnodes.size()], 1);
    T = VectorXd::Constant(numnodes[numnodes.size()], 0);
    SS = VectorXd::Constant(numnodes[numnodes.size()], 0);

    for (size_t i = 0; i < numnodes[numnodes.size() - 1]; i++)
    {
        T(i) = nodes[i].T;
    }

    VectorXd Inverse_nodes;

    for (size_t i = 0; i < (*m_Layers).GetCount(); i++)
    {
        Inverse_nodes(i) = 1 / (*m_Layers)[i].numNodes;
    }

    double MPD, MCD, CPG, Outer_radius_instantaneous, Outer_Radius, Layer_thickness_0, Layer_thickness_1, Initial_Layer_thickness_1, Initial_Layer_thickness_0;

    Outer_Radius = m_params.innerRadius + Total_thickness;
    VectorXd Layer_thickness, Initial_Layer_thickness;

    for (size_t i = 0; i < (*m_Layers).GetCount(); i++)
    {

        Layer_thickness(i) = (*m_Layers)[i].thickness;
        Initial_Layer_thickness(i) = (*m_Layers)[i].thickness;
    }

    if (m_options.coordinateType != CoordinateType::Cartesian)
    {
        Outer_radius_instantaneous = Outer_Radius - Initial_Layer_thickness_0 - Initial_Layer_thickness_1 + Layer_thickness(0) + Layer_thickness(1);
    }

    double AA1_s, AA2_s, BB1_s, BB2_s, CC1_s, CC2_s, DD1_s, DD2_s;
    double sigma;
    sigma = 5.67 * pow(10, -8);
    MPD = m_params.massRatePyro;
    MCD = m_params.massRateChar;
    CPG = nodes[0].Cppyrogas(nodes[0].T);

    if (Layer_thickness_0 != 0)
    {
        if (MCD != 0)
        {
            A(0, 0) = 1;
            A(0, 1) = 0;
            B(0, 0) = nodes[0].Tabl();
        }

        else if (MCD == 0)
        {
            AA2_s = nodes[0].k() / std::pow(Inverse_nodes(0), 2) - (MPD * CPG * Layer_thickness(0)) / (2 * Inverse_nodes(0));
            BB2_s = -((3 * nodes[0].k() + nodes[1].k()) / (2 * (pow(Inverse_nodes(0), 2)))) - (((nodes[0].rho() * nodes[0].Cp()) / delt) * Layer_thickness(0) * Layer_thickness(0));
            CC2_s = ((nodes[0].k() + nodes[1].k()) / (2 * (pow(Inverse_nodes(0), 2))) + (MPD * CPG * Layer_thickness(0)) / (2 * Inverse_nodes(0)));
            DD2_s = -(((nodes[0].rho() * nodes[0].Cp()) / delt) * T(0) * Layer_thickness(0) * Layer_thickness(0)) - SS(0) * Layer_thickness(0) * Layer_thickness(0);

            if (m_options.coordinateType != CoordinateType::Cartesian)
            {

                if (Outer_radius_instantaneous <= 0.0000001)
                {
                    CYSP2(AA2_s, BB2_s, CC2_s, nodes[0].k(), Inverse_nodes(0));
                }
                else if (Outer_radius_instantaneous > 0.0000001)
                {
                    CYSP1(AA2_s, CC2_s, nodes[0].k(), Layer_thickness(0), Inverse_nodes(0), Outer_radius_instantaneous);
                }
            }

            AA1_s = nodes[0].k() / (2 * Inverse_nodes(0));
            BB1_s = heat_transfer_coefficient_outer * Layer_thickness(0);
            CC1_s = -nodes[0].k() / (2 * Inverse_nodes(0));
            DD1_s = (heat_transfer_coefficient_outer * Tg - nodes[0].emissivity() * sigma * (pow(T(0), 4) - pow(m_params.TRad, 4)) + Qradiation) * Layer_thickness(0);

            A(0, 0) = BB1_s * AA2_s - AA1_s * BB2_s;
            A(0, 1) = AA2_s * CC1_s - AA1_s * CC2_s;
            B(0, 0) = AA2_s * DD1_s - AA1_s * DD2_s;
        }
    }
    m_ContinuumSolver->Solve(300);
    return ResponseResult();
}
