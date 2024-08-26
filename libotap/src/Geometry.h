#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include <vector>
#include "GeometryPrimitive.h"

namespace OTAP
{

    class Geometry
    {
    private:
        std::vector<GeometryPrimitive> m_geometryComponents;
        
    public:
        double m_characteristicLength=1.0;

        Geometry()
        {
            m_geometryComponents.push_back({OTAP::GeometryPrimitiveType::StagnationPoint, 0.0, 0.0, 0.0});
        }
        
        void Push(OTAP::GeometryPrimitiveType T, double length, double angle, double radius, double sweep);

        //NOTE: Template disabled

        // template <GeometryPrimitiveType T>
        // void Push(double length = 0.0, double angle = 0.0, double radius = 0.0, double sweep = 0.0);

        inline void Delete(size_t index){m_geometryComponents.erase(m_geometryComponents.begin()+index);}
        inline const std::vector<GeometryPrimitive> &GetComponents() const { return m_geometryComponents; }
        inline size_t Count() const { return m_geometryComponents.size(); }

        double GetRunningLength(size_t index) const; 
    };

} // namespace OTAP

#endif // _GEOMETRY_H