#pragma once

#include "imrad.h"
#include "Walnut/Layer.h"
#include "otap.h"

class ProjectBrowser : public Walnut::Layer
{
public:
    /// @begin interface
    virtual void OnAttach() override;
    virtual void OnDetach() override;

    virtual void OnUpdate(float ts) override;
    virtual void OnUIRender() override;

    /// @end interface

private:
    /// @begin impl

    bool isOpen = true;
    bool p_open = false;
    int m_SelectedTrajectoryType = -1;
    std::string m_TrajectoryFile = "Enter Trajectory File";
    OTAP::State &m_state = OTAP::State::GetInstance();

    void OnGeometryDelete(size_t index);
    void OnTrajectoryAdd();
    /// @end impl
};

// extern ProjectBrowser projectBrowser;
