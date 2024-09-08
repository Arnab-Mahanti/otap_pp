#include "OTypes.h"

OTAP::CpData::CpData(const std::string &filename, const bool delim_whitespace)
{
    rapidcsv::Document cpdata(filename,
                              rapidcsv::LabelParams(0, 0),
                              delim_whitespace ? rapidcsv::SeparatorParams() : rapidcsv::SeparatorParams(false, ',', true),
                              rapidcsv::ConverterParams(),
                              rapidcsv::LineReaderParams(true, '#', true));
    auto machs = cpdata.GetColumnNames();
    auto pos = cpdata.GetRowNames();
    std::transform(machs.begin(), machs.end(), std::back_inserter(m_Mach), [](auto &&a)
                   { return std::stod(a); });
    std::transform(pos.begin(), pos.end(), std::back_inserter(m_Positions), [](auto &&a)
                   { return std::stod(a); });

    for (size_t i = 0; i < cpdata.GetRowCount(); i++)
        for (size_t j = 0; j < cpdata.GetColumnCount(); j++)
            m_Cp.push_back(cpdata.GetCell<double>(j, i));
}

double OTAP::CpData::GetCp(double Mach, double Pos) const
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
