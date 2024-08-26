#include "GeometryViewport.h"
#include <iostream>
#include <cmath>

void GeometryViewport::OnUIRender()
{

    ImGui::Begin("Geometry Viewer");
    auto p = ImGui::GetCursorScreenPos();
    auto sz = ImGui::GetContentRegionAvail();
    auto drawList = ImGui::GetWindowDrawList();

    double x0 = p.x;
    double y0 = p.y;

    x0 += (sz.x / 2);
    y0 += (sz.y / 10);

    // TODO: Shift drawing functionality to OTAP geometry
    // draw

    for (auto &&geom : m_state.geometry.GetComponents())
    {
        switch (geom.type)
        {
        case OTAP::GeometryPrimitiveType::StagnationPoint:
            drawList->AddCircleFilled({x0, y0}, 3.0f, IM_COL32(179, 51, 77, 255));
            break;
        case OTAP::GeometryPrimitiveType::Line:
            drawList->AddLine({x0, y0},
                              {x0 + geom.length * std::sin(geom.angle), y0 + geom.length * std::cos(geom.angle)},
                              IM_COL32(200, 200, 200, 255), 2.0);
            x0 += geom.length * std::sin(geom.angle);
            y0 += geom.length * std::cos(geom.angle);
            break;
        case OTAP::GeometryPrimitiveType::Arc:
            if (geom.length == 0)
            {
                drawList->PathArcTo({x0, y0 + geom.radius}, geom.radius, -M_PI_2, -M_PI_2 + geom.angle);
                drawList->PathStroke(IM_COL32(200, 200, 200, 255), 0, 2);
                x0 += geom.radius * std::sin(geom.angle);
                y0 += geom.radius * (1-std::cos(geom.angle));
            }
            else
            {
                drawList->PathArcTo({x0, y0 + geom.radius}, geom.radius, -M_PI_2, -M_PI_2 + geom.angle);
                drawList->PathStroke(IM_COL32(200, 200, 200, 255), 0, 2);
                x0 += geom.radius * std::sin(geom.angle);
                y0 += geom.radius * (1-std::cos(geom.angle));
                drawList->AddLine({x0, y0},
                                  {x0 + geom.length * std::cos(geom.angle), y0 + geom.length * std::sin(geom.angle)},
                                  IM_COL32(200, 200, 200, 255), 2.0);
                x0 += geom.length * std::cos(geom.angle);
                y0 += geom.length * std::sin(geom.angle);
            }
            break;
        default:
            break;
        }
    }

    // ImGui::Image()

    ImGui::End();
}
