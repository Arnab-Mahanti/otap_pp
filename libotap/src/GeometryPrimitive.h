#ifndef _GEOMETRY_PRIMITIVE_H
#define _GEOMETRY_PRIMITIVE_H

#include "OTypes.h"
#include "MaterialLayer.h"
#include <cmath>

namespace OTAP
{

    struct GeometryPrimitive
    {
        GeometryPrimitiveType type;
        double length;
        double angle;
        double radius;
        double lambda;
        LayerStack layers = LayerStack();

        double Total_Length();

        bool isPlate() const
        {
            return type == GeometryPrimitiveType::Line;
        }

        bool isArc() const
        {
            return type == GeometryPrimitiveType::Arc;
        }

        double RunningLength() const
        {
            switch (type)
            {
            case GeometryPrimitiveType::StagnationPoint:
                return 0.0;
                break;
            case GeometryPrimitiveType::Line:
                return length;
                break;
            case GeometryPrimitiveType::Arc:
                return radius * angle * M_PI / 180.0 + length;
                break;
            default:
                break;
            }
        }
    };

} // namespace OTAP

#endif // _GEOMETRY_PRIMITIVE_H