#include "imrad.h"
#include "ProjectBrowser.h"
#include "GeometryAddPopup.h"
#include <iostream>
#include "nfd.hpp"

// ProjectBrowser projectBrowser;

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
    /// @style test
    /// @unit px
    /// @begin TopWindow
    auto *ioUserData = (ImRad::IOUserData *)ImGui::GetIO().UserData;
    ImGui::SetNextWindowSize({300, 800}, ImGuiCond_FirstUseEver);
    // TODO: Remove close buttons
    if (isOpen && ImGui::Begin("Project Browser###ProjectBrowser", nullptr, ImGuiWindowFlags_NoCollapse))
    {
        /// @separator

        /// @begin CollapsingHeader
        ImGui::SetNextItemOpen(false, ImGuiCond_Appearing);
        if (ImGui::CollapsingHeader("Geometry"))
        {
            /// @separator

            /// @begin Table
            if (ImGui::BeginTable("table1", 2,
                                  ImGuiTableFlags_BordersInnerH |
                                      ImGuiTableFlags_BordersInnerV |
                                      ImGuiTableFlags_BordersOuterH |
                                      ImGuiTableFlags_BordersOuterV |
                                      ImGuiTableFlags_Hideable |
                                      ImGuiTableFlags_RowBg,
                                  {0, 0}))
            {
                ImGui::TableSetupColumn("A", ImGuiTableColumnFlags_None, 0);
                ImGui::TableSetupColumn("B", ImGuiTableColumnFlags_None, 0);

                /// @separator

                auto &component = m_state.geometry.GetComponents();
                /// @begin Text
                ImGui::TableNextRow(0, 0);
                ImGui::TableSetColumnIndex(0);
                ImGui::TextUnformatted(
                    OTAP::GetNameFromType(component[0].type).c_str());
                for (size_t i = 1; i < m_state.geometry.Count(); i++)
                {
                    ImGui::TableNextRow(0, 0);
                    ImGui::TableSetColumnIndex(0);
                    ImGui::TextUnformatted(
                        OTAP::GetNameFromType(component[i].type).c_str());
                    ImGui::TableSetColumnIndex(1);
                    ImGui::PushItemWidth(-std::numeric_limits<float>::epsilon());
                    ImGui::PushID(i);
                    if (ImGui::Button("Delete"))
                    {
                        OnGeometryDelete(i);
                    }
                    ImGui::PopID();
                }
                /// @end Text

                /// @separator
                ImGui::EndTable();
            }
            /// @end Table

            /// @begin Button
            ImGui::Indent(1 * ImGui::GetStyle().IndentSpacing / 2);
            ImGui::Button("Add", {0, 0});
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_None))
                ImGui::SetTooltip("Adds geometry to the stack");
            if (ImGui::IsItemClicked())
                geometryAddPopup.OpenPopup();
            geometryAddPopup.Draw();
            /// @end Button

            /// @separator
        }
        /// @end CollapsingHeader

        /// @begin CollapsingHeader
        ImGui::SetNextItemOpen(false, ImGuiCond_Appearing);
        if (ImGui::CollapsingHeader("Trajectory"))
        {
            /// @separator

            /// @begin Table
            if (ImGui::BeginTable("table2", 3, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_SizingFixedSame | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoPadOuterX | ImGuiTableFlags_NoPadInnerX, {0, 0}))
            {
                ImGui::TableSetupColumn("item1", ImGuiTableColumnFlags_WidthStretch, 0);
                ImGui::TableSetupColumn("spacing1", ImGuiTableColumnFlags_WidthFixed, 8);
                ImGui::TableSetupColumn("item2", ImGuiTableColumnFlags_WidthStretch, 0);
                ImGui::TableNextRow(0, 0);
                ImGui::TableSetColumnIndex(0);
                /// @separator

                /// @begin Text
                ImGui::TextUnformatted("Type");
                /// @end Text

                /// @begin Combo
                ImRad::TableNextColumn(2);
                ImGui::SetNextItemWidth(200);
                ImGui::Combo("##SelectedTrajectory", &m_SelectedTrajectoryType, OTAP::GetTypes<OTAP::TrajectoryType>().c_str());
                /// @end Combo

                /// @separator

                ImGui::TableNextRow(0, 0);
                ImGui::TableSetColumnIndex(0);
                ImGui::TextUnformatted("File Name");
                ImRad::TableNextColumn(2);
                ImGui::SetNextItemWidth(200);
                ImGui::BeginDisabled();
                ImGui::TextUnformatted(m_TrajectoryFile.c_str());
                ImGui::EndDisabled();
                ImGui::SameLine();
                if (ImGui::Button("..."))
                {
                    NFD::Guard nfdGuard;
                    NFD::UniquePath outPath;
                    nfdfilteritem_t filterItem[1] = {{"Text", "txt"}};

                    nfdresult_t result = NFD::OpenDialog(outPath, filterItem, 1);
                    if (result == NFD_OKAY)
                    {
                        m_TrajectoryFile = outPath.get();
                    }
                }
                if (ImGui::Button("OK"))
                {
                    OnTrajectoryAdd();
                }

                ImGui::EndTable();
            }
            /// @end Table

            /// @separator
        }
        /// @end CollapsingHeader

        /// @separator
        ImGui::End();
    }
    /// @end TopWindow
}

void ProjectBrowser::OnGeometryDelete(size_t index)
{
    m_state.geometry.Delete(index);
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
