#include "OTypes.h"
#include "Material.h"
#include <ryml.hpp>
#include <ryml_std.hpp>
#include "FileUtils.h"

OTAP::MaterialManager::MaterialManager(const std::string &databasePath)
{
    auto contents = OTAP::file_get_contents(databasePath.c_str());
    m_database = ryml::parse_in_arena(ryml::to_csubstr(contents));
}

std::shared_ptr<OTAP::Material> OTAP::MaterialManager::GetMaterialInstance(const std::string &name)
{
    if (m_cache.find(name) != m_cache.end())
        return m_cache[name];
    else
    {
        auto dataid = m_database["data"].id();
        if (m_database.find_child(dataid, ryml::to_csubstr(name)) != ryml::NONE)
        {
            //TODO: Read directly from schema
            constexpr char *tabularkeys[] =
                {
                    "density",
                    "char_density",
                    "pyrogas_density",
                    "conductivity",
                    "char_conductivity",
                    "specific_heat",
                    "char_specific_heat",
                    "pyrogas_specific_heat",
                    "emissivity",
                };

            constexpr char *scalarkeys[] =
                {
                    "pyrolysis_temperature",
                    "ablation_temperature",
                    "pyrolysis_enthalpy",
                    "ablation_enthalpy",
                };

            auto data = m_database["data"][name.c_str()];
            std::vector<double> col1, col2;
            double scalar = 0.0;
            std::vector<Table<double>> tabulardata;
            std::vector<double> scalardata;

            for (auto &&i : tabularkeys)
            {
                data[i]["T"] >> col1;
                data[i]["val"] >> col2;
                assert(!col1.empty() && !col2.empty());
                assert(col1.size() == col2.size());
                tabulardata.push_back(Table<double>({col1, col2}));
            }

            for (auto &&i : scalarkeys)
            {
                data[i] >> scalar;
                scalardata.push_back(scalar);
            }

            return std::shared_ptr<Material>(new Material{
                name,
                tabulardata[0],
                tabulardata[1],
                tabulardata[2],
                tabulardata[3],
                tabulardata[4],
                tabulardata[5],
                tabulardata[6],
                tabulardata[7],
                tabulardata[8],
                scalardata[0],
                scalardata[1],
                scalardata[2],
                scalardata[3],
            });
        }
    }
    assert(false);
    return nullptr;
}

std::shared_ptr<OTAP::Material> OTAP::MaterialManager::GetMaterialInstance(const size_t &index)
{

    return nullptr;
}
