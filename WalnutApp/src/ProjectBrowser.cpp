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
            bool geometry_node_open = ImGui::TreeNodeEx(geometry.first.c_str(), ImGuiTreeNodeFlags_SpanFullWidth);
            ImGui::PushID(geometry.first.c_str());
            if (ImGui::BeginPopupContextItem("##geom_delete"))
            {
                if (ImGui::Button("Delete"))
                {
                    m_state.geometries.erase(geometry.first);
                    ImGui::CloseCurrentPopup();
                }
                ImGui::EndPopup();
            }
            ImGui::PopID();
            if (geometry_node_open)
            {
                // List all geometry componenets
                auto &components = geometry.second->GetComponents();
                for (size_t i = 0; i < components.size(); i++)
                {
                    ImGui::PushID(i);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TreeNodeEx(
                        OTAP::GetNameFromType(components[i].type).c_str(),
                        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_Bullet | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_SpanFullWidth);
                    if (i > 0)
                    {
                        if (ImGui::IsItemHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left))
                        {
                            ImGui::OpenPopup("Add Geometry Primitive");
                        }
                        OnGeometryPrimitiveAdd(*geometry.second, i);
                        ImGui::TableNextColumn();
                        if (ImGui::SmallButton(" - "))
                        {
                            components.erase(components.begin() + i);
                        }
                    }
                    ImGui::PopID();
                }
                ImGui::TreePop(); // individual geometry
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::PushID(geometry.first.c_str());
                if (ImGui::Button(" + "))
                    ImGui::OpenPopup("Add Geometry Primitive");
                OnGeometryPrimitiveAdd(*(geometry.second));
                ImGui::PopID();
            }
        }
        ImGui::TreePop(); // all geometry
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::Button(" + "))
            ImGui::OpenPopup("Add Geometry");
        OnGeometryAdd();
    }
}

void ProjectBrowser::OnGeometryPrimitiveAdd(OTAP::Geometry &geometry, size_t editindex)
{
    static bool editing = false;
    const float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;
    const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
    ImGui::SetNextWindowSize(ImVec2(300, 225));
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (ImGui::BeginPopupModal("Add Geometry Primitive", &isOpen, ImGuiWindowFlags_NoResize ^ ImGuiWindowFlags_NoTitleBar))
    {
        if (editindex && !editing)
        {
            m_SelectedGeometryType = (int)geometry.GetComponents()[editindex].type;
            m_GeometryLength = geometry.GetComponents()[editindex].length;
            m_GeometryAngle = geometry.GetComponents()[editindex].angle;
            m_GeometryRadius = geometry.GetComponents()[editindex].radius;
            m_GeometrySweep = geometry.GetComponents()[editindex].lambda;

            editing = true;
        }
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
                if (m_SelectedGeometryType == (int)OTAP::GeometryPrimitiveType::Arc && editindex == 1)
                {
                    geometry.GetComponents()[editindex].type = (OTAP::GeometryPrimitiveType)m_SelectedGeometryType;
                    geometry.GetComponents()[editindex].angle = m_GeometryAngle;
                    geometry.GetComponents()[editindex].length = m_GeometryLength;
                    geometry.GetComponents()[editindex].radius = m_GeometryRadius;
                    geometry.GetComponents()[editindex].lambda = m_GeometrySweep;
                }
                else if (editindex && m_SelectedGeometryType != (int)OTAP::GeometryPrimitiveType::Arc)
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
                editing = false;
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
                static size_t count = 0;
                m_state.geometries[m_geometryName] = std::make_shared<OTAP::Geometry>();
                count++;
                const auto name = "Geom {" + std::to_string(count) + "}";
                std::strncpy(m_geometryName, name.c_str(), std::min(name.size(), size_t(IM_ARRAYSIZE(m_geometryName))));
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndTable();
        }
        ImGui::EndPopup();
    }
}

void ProjectBrowser::ShowTrajectory()
{
    // Trajectory
    ImGui::TableNextRow();
    ImGui::TableNextColumn();

    if (ImGui::TreeNodeEx("Trajectory", ImGuiTreeNodeFlags_SpanFullWidth))
    {
        // List all trajectories
        for (auto &&trajectory : m_state.trajectories)
        {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            bool trajectory_node_open = ImGui::TreeNodeEx(trajectory.first.c_str(), ImGuiTreeNodeFlags_SpanFullWidth);

            ImGui::PushID(trajectory.first.c_str());
            if (ImGui::BeginPopupContextItem("##traj_delete"))
            {
                if (ImGui::Button("Delete"))
                {
                    m_state.trajectories.erase(trajectory.first);
                    ImGui::CloseCurrentPopup();
                }

                if (ImGui::Button("Edit"))
                {
                    auto traj = trajectory.second.second;

                    m_SelectedTrajectoryType = traj ? (int)traj->trajectoryType : -1;
                    m_SelectedFluidType = traj ? (int)traj->fluidType : -1;
                    m_SelectedAmbientType = traj ? (int)traj->ambientType : -1;

                    ImGui::OpenPopup("Edit Trajectory");
                }
                OnTrajectoryEdit(trajectory);
                ImGui::EndPopup();
            }
            ImGui::PopID();

            if (trajectory_node_open)
            {
                // List all trajectory settings

                {
                    ImGui::PushID(trajectory.first.c_str());

                    auto traj = trajectory.second.second;

                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted("Type");
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    ImGui::TextUnformatted(
                        traj ? OTAP::GetNameFromType<OTAP::TrajectoryType>(traj->trajectoryType).c_str() : "");

                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted("Fluid");
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    ImGui::TextUnformatted(
                        traj ? OTAP::GetNameFromType<OTAP::FluidType>(traj->fluidType).c_str() : "");

                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted("Ambient");
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    ImGui::TextUnformatted(
                        traj ? OTAP::GetNameFromType<OTAP::AmbientType>(traj->ambientType).c_str() : "");

                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted("File");
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    ImGui::TextUnformatted(
                        traj ? trajectory.second.first.c_str() : "");

                    ImGui::PopID();
                }

                ImGui::TreePop(); // individual geometry
            }
        }
        ImGui::TreePop(); // all geometry
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::Button(" + "))
            ImGui::OpenPopup("Add Trajectory");
        OnTrajectoryAdd();
    }
}

void ProjectBrowser::OnTrajectoryAdd()
{
    const float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;
    const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
    if (ImGui::BeginPopup("Add Trajectory"))
    {
        if (ImGui::BeginTable("2way", 2, flags))
        {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Name");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(15 * TEXT_BASE_WIDTH);
            ImGui::InputText("##traj_name", m_trajectoryName, IM_ARRAYSIZE(m_trajectoryName), ImGuiInputTextFlags_AutoSelectAll);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            if (ImGui::Button(" OK "))
            {
                static size_t count = 0;
                m_state.trajectories[m_trajectoryName] = std::make_pair("", OTAP::Trajectory());
                count++;
                const auto name = "Traj {" + std::to_string(count) + "}";
                std::strncpy(m_trajectoryName, name.c_str(), std::min(name.size(), size_t(IM_ARRAYSIZE(m_trajectoryName))));
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndTable();
        }
        ImGui::EndPopup();
    }
}

void ProjectBrowser::OnTrajectoryEdit(std::pair<const std::string, std::pair<std::string, OTAP::Trajectory>> &trajectory)
{
    const float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;
    const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
    // ImGui::PushID(trajectory.first.c_str());
    if (ImGui::BeginPopupModal("Edit Trajectory"))
    {
        if (ImGui::BeginTable("3way", 3, flags ^ ImGuiTableFlags_Resizable))
        {
            {
                // m_TrajectoryFile = traj ? (int)trajectory.second.first.c_str() : "Enter Trajectory File";

                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TextUnformatted("Type");
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::Combo("##m_SelectedTrajectoryType",
                             &m_SelectedTrajectoryType,
                             OTAP::GetTypes<OTAP::TrajectoryType>().c_str());

                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TextUnformatted("Fluid");
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::Combo("##m_SelectedFluidType",
                             &m_SelectedFluidType,
                             OTAP::GetTypes<OTAP::FluidType>().c_str());

                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TextUnformatted("Ambient");
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::Combo("##m_SelectedAmbientType",
                             &m_SelectedAmbientType,
                             OTAP::GetTypes<OTAP::AmbientType>().c_str());

                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TextUnformatted("File");
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::InputText("##traj_name", m_TrajectoryFile, IM_ARRAYSIZE(m_trajectoryName), ImGuiInputTextFlags_AutoSelectAll);
                ImGui::TableNextColumn();
                if (ImGui::SmallButton("..."))
                {
                    NFD::Guard nfdGuard;
                    NFD::UniquePath outPath;

                    nfdfilteritem_t filterItem[1] = {{"Datafile", "dat,txt,csv"}};

                    nfdresult_t result = NFD::OpenDialog(outPath, filterItem, 1);
                    if (result == NFD_OKAY)
                    {
                        std::strncpy(m_TrajectoryFile, outPath.get(), std::min(std::strlen(outPath.get()), size_t(IM_ARRAYSIZE(m_TrajectoryFile))));
                    }
                }

                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                if (ImGui::Button("OK"))
                {
                    trajectory.second.first = m_TrajectoryFile;
                    trajectory.second.second = OTAP::make_trajectory(
                        OTAP::GetTypeAtIndex<OTAP::TrajectoryType>(m_SelectedTrajectoryType),
                        OTAP::GetTypeAtIndex<OTAP::AmbientType>(m_SelectedAmbientType),
                        OTAP::GetTypeAtIndex<OTAP::FluidType>(m_SelectedFluidType),
                        trajectory.second.first);
                    ImGui::CloseCurrentPopup(); // this modal popup
                }
            }

            ImGui::EndTable();
        }
        ImGui::EndPopup();
    }
    // ImGui::PopID();
}

void ProjectBrowser::OnUpdate(float ts)
{
}
