#pragma once

#include <ryml.hpp>
#include <ryml_std.hpp>

namespace OTAP
{
    // FIXME: parse yaml
    class Material
    {
    private:
        std::string m_Name;
        double m_k,m_kchar,m_rho,m_rhochar,m_rhopyrogas,m_Cp,m_Cpchar,m_emissivity,m_Tabl,m_Tpyro,m_Habl,m_Hpyro;
    public:
        Material(double k, double rho, double Cp, double emissivity, double Tabl, double Tpyro)
        :m_k(k),m_kchar(k),m_rho(rho),m_rhochar(rho),m_Cp(Cp), m_Cpchar(Cp), m_emissivity(emissivity), m_Tabl(Tabl), m_Tpyro(Tpyro)
        {};
        double k(double T){return m_k;};
        double kchar(double T){return m_kchar;};
        double rho(double T){return m_rho;};
        double rhochar(double T){return m_rhochar;};
        double rhopyrogas(double T){return 0.0;};
        double Cp(double T){return m_Cp;};
        double Cpchar(double T){return m_Cp;};
        double Cppyrogas(double T){return 0.0;};
        double emissivity(double T){return m_emissivity;};
        double Tabl(){return m_Tabl;};
        double Tpyro(){return m_Tpyro;};
        double Habl(double T){return 0.0;};
        double Hpyro(double T){return 0.0;};
    };

} // namespace OTAP
