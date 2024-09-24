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
            ShowLayerStack();
            ShowBC();

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
        if (ImGui::Button(" + ##geom"))
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
    if (ImGui::BeginPopupModal("Add Geometry Primitive", &isOpen, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoTitleBar))
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
        if (ImGui::BeginTable("2way##geom", 2, flags))
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

                ImGui::TreePop(); // individual trajectory
            }
        }
        ImGui::TreePop(); // all trajectory
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::Button(" + ##traj"))
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
        if (ImGui::BeginTable("2way##traj", 2, flags))
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

void ProjectBrowser::ShowLayerStack()
{
    // LayerStack
    ImGui::TableNextRow();
    ImGui::TableNextColumn();

    if (ImGui::TreeNodeEx("LayerStacks", ImGuiTreeNodeFlags_SpanFullWidth))
    {
        // List all LayerStacks
        for (auto &&layerstack : m_state.layerstacks)
        {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            bool layerstack_node_open = ImGui::TreeNodeEx(layerstack.first.c_str(), ImGuiTreeNodeFlags_SpanFullWidth);

            ImGui::PushID(layerstack.first.c_str());
            if (ImGui::BeginPopupContextItem("##ls_delete"))
            {
                if (ImGui::Button("Delete"))
                {
                    m_state.layerstacks.erase(layerstack.first);
                    ImGui::CloseCurrentPopup();
                }
                ImGui::EndPopup();
            }

            if (layerstack_node_open)
            {
                // List all layers
                for (size_t i = 0; i < layerstack.second->GetCount(); i++)
                {
                    ImGui::PushID(i);

                    {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        ImGui::TreeNodeEx(
                            ("Layer " + std::to_string(i)).c_str(),
                            ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_Bullet | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_SpanFullWidth);

                        if (ImGui::IsItemHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left))
                        {
                            ImGui::OpenPopup("Add Layer");
                        }
                        OnLayerAdd(*layerstack.second, true, i);

                        ImGui::TableNextColumn();
                        if (ImGui::SmallButton(" - "))
                        {
                            layerstack.second->Delete(i);
                        }
                    }

                    ImGui::PopID();
                }

                ImGui::TreePop(); // individual layerstack

                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::Button(" + ##layer"))
                    ImGui::OpenPopup("Add Layer");
                OnLayerAdd(*layerstack.second);
            }
            ImGui::PopID();
        }
        ImGui::TreePop(); // all layerstacks
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::Button(" + ##layerstack"))
            ImGui::OpenPopup("Add LayerStack");
        OnLayerStackAdd();
    }
}

void ProjectBrowser::OnLayerStackAdd()
{
    const float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;
    const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
    if (ImGui::BeginPopup("Add LayerStack"))
    {
        if (ImGui::BeginTable("2way##ls", 2, flags))
        {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Name");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(15 * TEXT_BASE_WIDTH);
            ImGui::InputText("##ls_name", m_layerstackName, IM_ARRAYSIZE(m_layerstackName), ImGuiInputTextFlags_AutoSelectAll);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            if (ImGui::Button(" OK "))
            {
                static size_t count = 0;
                m_state.layerstacks[m_layerstackName] = std::make_shared<OTAP::LayerStack>();
                count++;
                const auto name = "LayerStack {" + std::to_string(count) + "}";
                std::strncpy(m_layerstackName, name.c_str(), std::min(name.size(), size_t(IM_ARRAYSIZE(m_layerstackName))));
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndTable();
        }
        ImGui::EndPopup();
    }
}

void ProjectBrowser::OnLayerAdd(OTAP::LayerStack &layerstack, bool editing, size_t index)
{
    const float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;
    const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
    ImGui::SetNextWindowSize(ImVec2(300, 175));
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    static auto matnames = m_state.materialManager.GetMaterialList();

    if (ImGui::BeginPopupModal("Add Layer", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoTitleBar))
    {
        static bool toedit = false;
        static bool sublime = false;
        static int item_current_idx = 0;                                      // Here we store our selection data as an index.
        const char *combo_preview_value = matnames[item_current_idx].c_str(); // Pass in the preview value visible before opening the combo (it could be anything)

        if (editing && !toedit)
        {
            m_layerThickness = layerstack[index].thickness;
            m_layerNumnodes = layerstack[index].numNodes;
            item_current_idx = std::distance(matnames.begin(), std::find(matnames.begin(), matnames.end(), layerstack[index].material->name));
            toedit = true;
            sublime = layerstack[index].material->sublime;
        }

        if (ImGui::BeginTable("add_layer", 2, flags ^ ImGuiTableFlags_Resizable))
        {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Material");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);

            if (ImGui::BeginCombo("##Select Material", combo_preview_value))
            {
                for (int n = 0; n < matnames.size(); n++)
                {
                    const bool is_selected = (item_current_idx == n);
                    if (ImGui::Selectable(matnames[n].c_str(), is_selected))
                        item_current_idx = n;

                    // Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
                    if (is_selected)
                        ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::Checkbox("Sublime", &sublime);

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Thickness");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::InputFloat("##layer_thickness", &m_layerThickness, 0.0001, 0.001, "%.6f");

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Nodes");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::InputInt("##layer_numnodes", &m_layerNumnodes, 1, 10);
            ImGui::EndTable();
        }

        ImGui::Separator();

        if (ImGui::Button("OK ##add_layer"))
        {
            if (toedit)
            {
                layerstack[index].thickness = m_layerThickness;
                layerstack[index].numNodes = m_layerNumnodes;
                layerstack[index].material = m_state.materialManager.GetMaterialInstance(matnames[item_current_idx], sublime);

                toedit = false;
            }
            else
            {
                layerstack.Emplace(m_state.materialManager.GetMaterialInstance(matnames[item_current_idx], sublime),
                                   std::max(0.0, (double)m_layerThickness),
                                   (size_t)m_layerNumnodes);
                m_layerThickness = 0;
                m_layerNumnodes = 0;
            }

            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}

void ProjectBrowser::ShowBC()
{
    // BC
    ImGui::TableNextRow();
    ImGui::TableNextColumn();

    if (ImGui::TreeNodeEx("Boundary Conditions", ImGuiTreeNodeFlags_SpanFullWidth))
    {
        size_t i = 0;
        for (auto &&bc : m_state.bcs)
        {
            ImGui::PushID(i);
            {
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TreeNodeEx(
                    bc.first.c_str(),
                    ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_Bullet | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_SpanFullWidth);

                static bool edit = false;
                if (ImGui::IsItemHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left))
                {
                    edit = true;
                }
                // EditBC(bc, edit);

                ImGui::TableNextColumn();
                if (ImGui::SmallButton(" - "))
                {
                    m_state.bcs.erase(bc.first);
                }
            }
            ImGui::PopID();
            i++;
        }

        ImGui::TreePop();

        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        if (ImGui::Button(" + ##bc"))
        {
            ImGui::OpenPopup("Choose BC Type");
        }
        ShowBCchoices();
    }
}

void ProjectBrowser::ShowBCchoices()
{
    if (ImGui::BeginPopup("Choose BC Type"))
    {
        auto bcChoices = OTAP::GetTypesAsStrings<OTAP::BCType>();
        static int selected_bc = -1;
        for (int i = 0; i < bcChoices.size(); i++)
            if (ImGui::Selectable(bcChoices[i].c_str(), &selected_bc, ImGuiSelectableFlags_DontClosePopups))
                selected_bc = i;

        if (selected_bc != -1)
        {
            auto selectedbctype = OTAP::GetTypeAtIndex<OTAP::BCType>(selected_bc);
            bool editting = false;
            switch (selectedbctype)
            {
            case OTAP::BCType::Convection:
                ImGui::OpenPopup("ConvectionBC");
                if (AddConvectionBC(selected_bc, editting))
                    ImGui::CloseCurrentPopup();
                break;
            case OTAP::BCType::Flux:
                ImGui::OpenPopup("FluxBC");
                if (AddFluxBC(selected_bc, editting))
                    ImGui::CloseCurrentPopup();
                break;
            case OTAP::BCType::HeatGeneration:
                ImGui::OpenPopup("HeatGenBC");
                if (AddHeatGenBC(selected_bc, editting))
                    ImGui::CloseCurrentPopup();
                break;
            case OTAP::BCType::NaturalConvection:
                ImGui::OpenPopup("NatConvBC");
                if (AddNaturalConvectionBC(selected_bc, editting))
                    ImGui::CloseCurrentPopup();
                break;
            case OTAP::BCType::PropellantMass:
                ImGui::OpenPopup("PropMassBC");
                if (AddPropellantMassBC(selected_bc, editting))
                    ImGui::CloseCurrentPopup();
                break;
            case OTAP::BCType::Radiation:
                ImGui::OpenPopup("RadBC");
                if (AddRadiationBC(selected_bc, editting))
                    ImGui::CloseCurrentPopup();
                break;
            default:
                break;
            }
        }
        ImGui::EndPopup();
    }
}

void ProjectBrowser::EditBC(const std::pair<const std::string, OTAP::BC> &bc, bool &peditting)
{
    int done = -1;
    static bool editting = false;
    if (peditting && !editting)
    {
        editting = true;
    }
    
    if (editting)
    {

        switch (bc.second->type)
        {
        case OTAP::BCType::Convection:
            ImGui::OpenPopup("ConvectionBC");
            AddConvectionBC(done, editting, bc);
            break;
        case OTAP::BCType::Flux:
            ImGui::OpenPopup("FluxBC");
            AddFluxBC(done, editting, bc);
            break;
        case OTAP::BCType::HeatGeneration:
            ImGui::OpenPopup("HeatGenBC");
            AddHeatGenBC(done, editting, bc);
            break;
        case OTAP::BCType::NaturalConvection:
            ImGui::OpenPopup("NatConvBC");
            AddNaturalConvectionBC(done, editting, bc);
            break;
        case OTAP::BCType::PropellantMass:
            ImGui::OpenPopup("PropMassBC");
            AddPropellantMassBC(done, editting, bc);
            break;
        case OTAP::BCType::Radiation:
            ImGui::OpenPopup("RadBC");
            AddRadiationBC(done, editting, bc);
            break;
        default:
            break;
        }
    }
}

bool ProjectBrowser::AddConvectionBC(int &done, bool &editting, const std::pair<const std::string, OTAP::BC> &pbc)
{
    const float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;
    const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
    static bool toedit = false;
    ImGui::SetNextWindowSize(ImVec2(475, 225));
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    static bool h_delim = false;
    static bool tg_delim = false;
    static bool front = false;
    static bool back = false;
    auto bc = const_cast<std::pair<const std::string, OTAP::BC> &>(pbc);
    auto cbc = std::dynamic_pointer_cast<OTAP::ConvectionBC>(bc.second);

    if (editting && !toedit)
    {
        std::strncpy(m_bcName, bc.first.c_str(), std::min(bc.first.size(), size_t(IM_ARRAYSIZE(m_bcName))));
        front = (cbc->location == OTAP::BCLocation::both || cbc->location == OTAP::BCLocation::front) ? true : false;
        back = (cbc->location == OTAP::BCLocation::both || cbc->location == OTAP::BCLocation::back) ? true : false;

        toedit = true;
        editting = false;
    }

    if (ImGui::BeginPopupModal("ConvectionBC", nullptr, ImGuiWindowFlags_NoResize))
    {
        if (ImGui::BeginTable("4way##convbc", 4, flags | ImGuiTableColumnFlags_NoResize | ImGuiWindowFlags_NoTitleBar))
        {
            ImGui::TableSetupColumn("Parameter", ImGuiTableColumnFlags_NoHide);
            ImGui::TableSetupColumn("File", ImGuiTableColumnFlags_NoHide);
            ImGui::TableSetupColumn("", ImGuiTableColumnFlags_NoHide);
            ImGui::TableSetupColumn("csv", ImGuiTableColumnFlags_NoHide);

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Name");
            ImGui::TableNextColumn();
            if (toedit)
            {
                ImGui::TextDisabled(bc.first.c_str());
            }
            else
            {
                ImGui::InputText("##conv_bcname", m_bcName, IM_ARRAYSIZE(m_bcName), ImGuiInputTextFlags_AutoSelectAll);
            }
            ImGui::TableHeadersRow();

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Convection Coefficient");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(15 * TEXT_BASE_WIDTH);
            ImGui::InputText("##h_filename", m_hfilename, IM_ARRAYSIZE(m_hfilename), ImGuiInputTextFlags_AutoSelectAll);
            ImGui::TableNextColumn();
            if (ImGui::SmallButton("...##h"))
            {
                NFD::Guard nfdGuard;
                NFD::UniquePath outPath;

                nfdfilteritem_t filterItem[1] = {{"Datafile", "dat,txt,csv"}};

                nfdresult_t result = NFD::OpenDialog(outPath, filterItem, 1);
                if (result == NFD_OKAY)
                {
                    std::strncpy(m_hfilename, outPath.get(), std::min(std::strlen(outPath.get()), size_t(IM_ARRAYSIZE(m_hfilename))));
                }
            }

            ImGui::TableNextColumn();
            ImGui::Checkbox("", &h_delim);

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Ambient Temperature");
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(15 * TEXT_BASE_WIDTH);
            ImGui::InputText("##tg_filename", m_tgfilename, IM_ARRAYSIZE(m_tgfilename), ImGuiInputTextFlags_AutoSelectAll);
            ImGui::TableNextColumn();
            if (ImGui::SmallButton("...##tg"))
            {
                NFD::Guard nfdGuard;
                NFD::UniquePath outPath;

                nfdfilteritem_t filterItem[1] = {{"Datafile", "dat,txt,csv"}};

                nfdresult_t result = NFD::OpenDialog(outPath, filterItem, 1);
                if (result == NFD_OKAY)
                {
                    std::strncpy(m_tgfilename, outPath.get(), std::min(std::strlen(outPath.get()), size_t(IM_ARRAYSIZE(m_tgfilename))));
                }
            }

            ImGui::TableNextColumn();
            ImGui::Checkbox("", &tg_delim);

            ImGui::TableNextRow();
            ImGui::TableNextColumn();

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Apply to");
            ImGui::TableNextColumn();
            ImGui::Checkbox("Front", &front);
            ImGui::TableNextColumn();
            ImGui::Checkbox("Back", &back);
            ImGui::EndTable();
        }

        ImGui::Separator();
        if (ImGui::Button(" OK "))
        {
            if (toedit)
            {
                // bc.first = m_bcName;
                if (std::strcmp(m_hfilename, "") != 0)
                    cbc->h.Load(m_hfilename, h_delim);
                if (std::strcmp(m_tgfilename, "") != 0)
                    cbc->tg.Load(m_tgfilename, tg_delim);
                bc.second->location = (front && back) ? OTAP::BCLocation::both : (front ? OTAP::BCLocation::front : OTAP::BCLocation::back);

                toedit = false;
                editting = false;
            }
            else
            {

                auto bc = std::make_shared<OTAP::ConvectionBC>();
                bc->h.Load(m_hfilename, !h_delim);
                bc->tg.Load(m_tgfilename, !tg_delim);
                bc->location = (front && back) ? OTAP::BCLocation::both : (front ? OTAP::BCLocation::front : OTAP::BCLocation::back);
                m_state.bcs[m_bcName] = bc;

                bcnamecount++;
                const auto name = "BC {" + std::to_string(bcnamecount) + "}";
                std::strncpy(m_bcName, name.c_str(), std::min(name.size(), size_t(IM_ARRAYSIZE(m_bcName))));
                done = -1;
            }

            std::strcpy(m_hfilename, "");
            std::strcpy(m_tgfilename, "");
            front = false;
            back = false;
            h_delim = false;
            tg_delim = false;
            ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
            return true;
        }
    }

    ImGui::EndPopup();
    return false;
}

bool ProjectBrowser::AddFluxBC(int &done, bool &editting, const std::pair<const std::string, OTAP::BC> &bc)
{
}

bool ProjectBrowser::AddHeatGenBC(int &done, bool &editting, const std::pair<const std::string, OTAP::BC> &bc)
{
}

bool ProjectBrowser::AddNaturalConvectionBC(int &done, bool &editting, const std::pair<const std::string, OTAP::BC> &bc)
{
}

bool ProjectBrowser::AddRadiationBC(int &done, bool &editting, const std::pair<const std::string, OTAP::BC> &bc)
{
}

bool ProjectBrowser::AddPropellantMassBC(int &done, bool &editting, const std::pair<const std::string, OTAP::BC> &bc)
{
}

void ProjectBrowser::OnUpdate(float ts)
{
}
