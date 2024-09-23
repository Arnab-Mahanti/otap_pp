#pragma once

#include <memory>
#include <string>
#include <ryml.hpp>
#include <ryml_std.hpp>
#include "Table.h"

namespace OTAP
{
    struct Material
    {
        std::string name;
        Table<double> rho;
        Table<double> rhochar;
        Table<double> rhopyrogas;
        Table<double> k;
        Table<double> kchar;
        Table<double> Cp;
        Table<double> Cpchar;
        Table<double> Cppyrogas;
        Table<double> emissivity;
        double Tpyro;
        double Tabl;
        double Hpyro;
        double Habl;
        bool sublime = false;
    };

    // FIXME: parse yaml
    class MaterialManager
    {
    private:
        std::unordered_map<std::string, std::shared_ptr<Material>> m_cache;
        ryml::Tree m_database;

    public:
        MaterialManager(const std::string &databasePath = std::string());
        ~MaterialManager() = default;

        std::shared_ptr<Material> GetMaterialInstance(const std::string &name, bool sublime = false);
        std::vector<std::string> GetMaterialList();
        std::shared_ptr<Material> GetMaterialInstance(const size_t &index, bool sublime = false);
    };

} // namespace OTAP
