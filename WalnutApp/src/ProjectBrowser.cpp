#include "imrad.h"
#include "ProjectBrowser.h"
// #include "GeometryAddPopup.h"
#include <iostream>
#include "nfd.hpp"

// ProjectBrowser projectBrowser;

static ImGuiTableFlags flags = ImGuiTableFlags_BordersV | ImGuiTableFlags_BordersOuterH | ImGuiTableFlags_Resizable | ImGuiTableFlags_RowBg | ImGuiTableFlags_NoBordersInBody | ImGuiTableFlags_SizingStretchProp;

void ProjectBrowser::OnAttach()
{
    isOpen = true;
}

void ProjectBrowser::OnDetach()
{
    isOpen = false;
}

// void ProjectBrowser::OnUIRender()
// {
//     /// @style test
//     /// @unit px
//     /// @begin TopWindow
//     auto *ioUserData = (ImRad::IOUserData *)ImGui::GetIO().UserData;
//     ImGui::SetNextWindowSize({300, 800}, ImGuiCond_FirstUseEver);
//     // TODO: Remove close buttons
//     if (isOpen && ImGui::Begin("Project Browser###ProjectBrowser", nullptr, ImGuiWindowFlags_NoCollapse))
//     {
//         /// @separator

//         /// @begin CollapsingHeader
//         ImGui::SetNextItemOpen(false, ImGuiCond_Appearing);
//         if (ImGui::CollapsingHeader("Geometry"))
//         {
//             /// @separator

//             /// @begin Table
//             if (ImGui::BeginTable("table1", 2,
//                                   ImGuiTableFlags_BordersInnerH |
//                                       ImGuiTableFlags_BordersInnerV |
//                                       ImGuiTableFlags_BordersOuterH |
//                                       ImGuiTableFlags_BordersOuterV |
//                                       ImGuiTableFlags_Hideable |
//                                       ImGuiTableFlags_RowBg,
//                                   {0, 0}))
//             {
//                 ImGui::TableSetupColumn("A", ImGuiTableColumnFlags_None, 0);
//                 ImGui::TableSetupColumn("B", ImGuiTableColumnFlags_None, 0);

//                 /// @separator

//                 auto &component = m_state.geometry.GetComponents();
//                 /// @begin Text
//                 ImGui::TableNextRow(0, 0);
//                 ImGui::TableSetColumnIndex(0);
//                 ImGui::TextUnformatted(
//                     OTAP::GetNameFromType(component[0].type).c_str());
//                 for (size_t i = 1; i < m_state.geometry.Count(); i++)
//                 {
//                     ImGui::TableNextRow(0, 0);
//                     ImGui::TableSetColumnIndex(0);
//                     ImGui::TextUnformatted(
//                         OTAP::GetNameFromType(component[i].type).c_str());
//                     ImGui::TableSetColumnIndex(1);
//                     ImGui::PushItemWidth(-std::numeric_limits<float>::epsilon());
//                     ImGui::PushID(i);
//                     if (ImGui::Button("Delete"))
//                     {
//                         OnGeometryDelete(i);
//                     }
//                     ImGui::PopID();
//                 }
//                 /// @end Text

//                 /// @separator
//                 ImGui::EndTable();
//             }
//             /// @end Table

//             /// @begin Button
//             ImGui::Indent(1 * ImGui::GetStyle().IndentSpacing / 2);
//             ImGui::Button("Add", {0, 0});
//             if (ImGui::IsItemHovered(ImGuiHoveredFlags_None))
//                 ImGui::SetTooltip("Adds geometry to the stack");
//             if (ImGui::IsItemClicked())
//                 geometryAddPopup.OpenPopup();
//             geometryAddPopup.Draw();
//             /// @end Button

//             /// @separator
//         }
//         /// @end CollapsingHeader

//         /// @begin CollapsingHeader
//         ImGui::SetNextItemOpen(false, ImGuiCond_Appearing);
//         if (ImGui::CollapsingHeader("Trajectory"))
//         {
//             /// @separator

//             /// @begin Table
//             if (ImGui::BeginTable("table2", 3, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_SizingFixedSame | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoPadOuterX | ImGuiTableFlags_NoPadInnerX, {0, 0}))
//             {
//                 ImGui::TableSetupColumn("item1", ImGuiTableColumnFlags_WidthStretch, 0);
//                 ImGui::TableSetupColumn("spacing1", ImGuiTableColumnFlags_WidthFixed, 8);
//                 ImGui::TableSetupColumn("item2", ImGuiTableColumnFlags_WidthStretch, 0);
//                 ImGui::TableNextRow(0, 0);
//                 ImGui::TableSetColumnIndex(0);
//                 /// @separator

//                 /// @begin Text
//                 ImGui::TextUnformatted("Type");
//                 /// @end Text

//                 /// @begin Combo
//                 ImRad::TableNextColumn(2);
//                 ImGui::SetNextItemWidth(200);
//                 ImGui::Combo("##SelectedTrajectory", &m_SelectedTrajectoryType, OTAP::GetTypes<OTAP::TrajectoryType>().c_str());
//                 /// @end Combo

//                 /// @separator

//                 ImGui::TableNextRow(0, 0);
//                 ImGui::TableSetColumnIndex(0);
//                 ImGui::TextUnformatted("File Name");
//                 ImRad::TableNextColumn(2);
//                 ImGui::SetNextItemWidth(200);
//                 ImGui::BeginDisabled();
//                 ImGui::TextUnformatted(m_TrajectoryFile.c_str());
//                 ImGui::EndDisabled();
//                 ImGui::SameLine();
//                 if (ImGui::Button("..."))
//                 {
//                     NFD::Guard nfdGuard;
//                     NFD::UniquePath outPath;
//                     nfdfilteritem_t filterItem[1] = {{"Text", "txt"}};

//                     nfdresult_t result = NFD::OpenDialog(outPath, filterItem, 1);
//                     if (result == NFD_OKAY)
//                     {
//                         m_TrajectoryFile = outPath.get();
//                     }
//                 }
//                 if (ImGui::Button("OK"))
//                 {
//                     OnTrajectoryAdd();
//                 }

//                 ImGui::EndTable();
//             }
//             /// @end Table

//             /// @separator
//         }
//         /// @end CollapsingHeader

//         /// @separator
//         ImGui::End();
//     }
//     /// @end TopWindow
// }

void ProjectBrowser::OnUIRender()
{
    auto *ioUserData = (ImRad::IOUserData *)ImGui::GetIO().UserData;
    ImGui::SetNextWindowSize({375, 800}, ImGuiCond_FirstUseEver);

    const float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;
    const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();

    if (isOpen && ImGui::Begin("Project Browser###ProjectBrowser", nullptr, ImGuiWindowFlags_NoCollapse))
    {
        if (ImGui::BeginTable("3way", 3, flags))
        {
            ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_NoHide);
            ImGui::TableSetupColumn("Status", ImGuiTableColumnFlags_WidthFixed, TEXT_BASE_WIDTH * 12.0f);
            ImGui::TableSetupColumn("Description", ImGuiTableColumnFlags_WidthFixed, TEXT_BASE_WIDTH * 18.0f);
            ImGui::TableHeadersRow();

            // Show various sections
            ShowGeometry();
            ShowTrajectory();

            ImGui::EndTable();
        }
        ImGui::End();
    }
}

void ProjectBrowser::ShowGeometry()
{
    // Geometry
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    if (ImGui::TreeNodeEx("Geometry", ImGuiTreeNodeFlags_SpanFullWidth))
    {
        // List all geometries
        for (auto &&geometry : m_state.geometries)
        {
            if (ImGui::TreeNodeEx(geometry.first.c_str(), ImGuiTreeNodeFlags_SpanFullWidth))
            {
                // List all geometry componenets
                auto components = geometry.second->GetComponents();
                for (size_t i = 0; i < components.size(); i++)
                {
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TreeNodeEx(
                        OTAP::GetNameFromType(components[i].type).c_str(),
                        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_Bullet | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_SpanFullWidth);
                    if (ImGui::IsItemHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left))
                    {
                        OnGeometryPrimitiveAdd(*geometry.second, true, i);
                    }
                }
                ImGui::TreePop();
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::Button(" + "))
                    ImGui::OpenPopup("Add Geometry Primitive");
                OnGeometryPrimitiveAdd(*(geometry.second));
            }
        }
        ImGui::TreePop();
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::Button(" + "))
            ImGui::OpenPopup("Add Geometry");
        OnGeometryAdd();
    }
}

void ProjectBrowser::OnGeometryPrimitiveAdd(OTAP::Geometry &geometry, bool edit, size_t editindex)
{
    if (edit)
    {
        m_SelectedGeometryType = (int)geometry.GetComponents()[editindex].type;
        m_GeometryLength = geometry.GetComponents()[editindex].length;
        m_GeometryAngle = geometry.GetComponents()[editindex].angle;
        m_GeometryRadius = geometry.GetComponents()[editindex].radius;
        m_GeometrySweep = geometry.GetComponents()[editindex].lambda;
    }

    const float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;
    const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
    ImGui::SetNextWindowSize(ImVec2(300, 225));
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (ImGui::BeginPopup("Add Geometry Primitive", ImGuiWindowFlags_NoResize ^ ImGuiWindowFlags_NoTitleBar))
    {
        if (ImGui::BeginTable("add_primitive", 2, flags ^ ImGuiTableFlags_Resizable))
        {
            ImGui::TableSetupColumn("data", ImGuiTableColumnFlags_NoHide);
            ImGui::TableSetupColumn("value", ImGuiTableColumnFlags_NoHide);

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Type");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::Combo("##m_SelectedGeometryType",
                         &m_SelectedGeometryType,
                         OTAP::GetTypes<OTAP::GeometryPrimitiveType>().c_str());

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Length");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::InputDouble("##m_GeometryLength", &m_GeometryLength, 1, 0.0, "%.3f");

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Radius");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::InputDouble("##m_GeometryRadius", &m_GeometryRadius, 1, 0.0, "%.3f");

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Angle");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderAngle("##m_GeometryAngle", &m_GeometryAngle, -180, 180, "%.3f");

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Sweep");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderAngle("##m_GeometrySweep", &m_GeometrySweep, -180, 180, "%.3f");

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::Button(" OK "))
            {
                if (edit)
                {
                    geometry.GetComponents()[editindex].type = (OTAP::GeometryPrimitiveType)m_SelectedGeometryType;
                    geometry.GetComponents()[editindex].angle = m_GeometryAngle;
                    geometry.GetComponents()[editindex].length = m_GeometryLength;
                    geometry.GetComponents()[editindex].radius = m_GeometryRadius;
                    geometry.GetComponents()[editindex].lambda = m_GeometrySweep;
                }
                else
                {
                    geometry.Push(OTAP::GetTypeAtIndex<OTAP::GeometryPrimitiveType>(m_SelectedGeometryType),
                                  m_GeometryLength,
                                  m_GeometryAngle,
                                  m_GeometryRadius,
                                  m_GeometrySweep);
                }
                ImGui::CloseCurrentPopup();
                m_SelectedGeometryType = -1;
                m_GeometryLength = 0;
                m_GeometryAngle = 0;
                m_GeometryRadius = 0;
                m_GeometrySweep = 0;
            }
            ImGui::EndTable();
        }
        ImGui::EndPopup();
    }
}

void ProjectBrowser::OnGeometryAdd()
{
    const float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;
    const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
    if (ImGui::BeginPopup("Add Geometry"))
    {
        if (ImGui::BeginTable("2way", 2, flags))
        {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Name");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(15 * TEXT_BASE_WIDTH);
            ImGui::InputText("##geom_name", m_geometryName, IM_ARRAYSIZE(m_geometryName), ImGuiInputTextFlags_AutoSelectAll);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            if (ImGui::Button(" OK "))
            {
                m_state.geometries[m_geometryName] = std::make_shared<OTAP::Geometry>(); 
                const auto count = m_state.geometries.size();
                const auto name = "Geom {" + std::to_string(count) + "}";
                std::strncpy(m_geometryName, name.c_str(), std::min(name.size(), size_t(IM_ARRAYSIZE(m_geometryName))));
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndTable();
        }
        ImGui::EndPopup();
    }
}

void ProjectBrowser::OnGeometryDelete(size_t index)
{
    // m_state.geometry.Delete(index);
}

void ProjectBrowser::ShowTrajectory()
{
}

void ProjectBrowser::OnTrajectoryAdd()
{
    // m_state.trajectory.reset(OTAP::make_trajectory(
    //     OTAP::GetTypeAtIndex<OTAP::TrajectoryType>(m_SelectedTrajectoryType),
    //     m_TrajectoryFile));
}

void ProjectBrowser::OnUpdate(float ts)
{
}
