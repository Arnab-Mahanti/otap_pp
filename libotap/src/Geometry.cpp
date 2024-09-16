#include <cassert>
#include "Geometry.h"

// TODO: Implement sweep functionality for geometry

// NOTE: Disabled template system

// template <OTAP::GeometryPrimitiveType T>
// void OTAP::Geometry::Push(double length, double angle, double radius, double sweep)
// {
//     static_assert(false);
// }

// template <>
// void OTAP::Geometry::Push<OTAP::GeometryPrimitiveType::StagnationPoint>(double length, double angle, double radius, double sweep)
// {
//     m_geometryComponents.push_back({OTAP::GeometryPrimitiveType::StagnationPoint, 0.0, 0.0, 0.0, 0.0});
// }

// template <>
// void OTAP::Geometry::Push<OTAP::GeometryPrimitiveType::Line>(double length, double angle, double radius, double sweep)
// {
//     m_geometryComponents.push_back({OTAP::GeometryPrimitiveType::Line, length, angle, 0.0, 0.0});
// }

// template <>
// void OTAP::Geometry::Push<OTAP::GeometryPrimitiveType::Arc>(double length, double angle, double radius, double sweep)
// {
//     m_geometryComponents.push_back({OTAP::GeometryPrimitiveType::Arc, length, angle, radius, 0.0});
// }

void OTAP::Geometry::Push(OTAP::GeometryPrimitiveType T, double length, double angle, double radius, double sweep)
{
    m_geometryComponents.push_back({T, length, angle, radius, sweep, LayerStack()});
}

void OTAP::Geometry::Push(OTAP::GeometryPrimitiveType T, double length, double angle, double radius, double sweep, OTAP::LayerStack ls)
{
    m_geometryComponents.push_back({T, length, angle, radius, sweep, ls});
}

double OTAP::Geometry::GetRunningLength(size_t index) const
{
    double runningLength = 0.0;
    index = std::min(m_geometryComponents.size() - 1, index);
    if (index == 0)
    {
        return 0;
    }
    for (size_t i = 0; i <= index - 1; i++)
    {
        runningLength += m_geometryComponents[i].RunningLength();
    }
    return runningLength;
}