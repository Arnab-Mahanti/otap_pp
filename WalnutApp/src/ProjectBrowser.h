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

    ProjectBrowser()
    {
    }

    /// @end interface

private:
    // General

    bool isOpen = true;
    bool p_open = false;
    OTAP::State &m_state = OTAP::State::GetInstance();

    // Geometry
    bool m_AddGeom = false;
    int m_SelectedGeometryType = -1;
    double m_GeometryLength = 0.0;
    double m_GeometryRadius = 0.0;
    float m_GeometryAngle = 0.0f;
    float m_GeometrySweep = 0.0f;
    char m_geometryName[50] = "Geom {0}";

    void OnGeometryAdd();
    void OnGeometryDelete(size_t index);
    void OnGeometryPrimitiveAdd(OTAP::Geometry &geometry, size_t editindex = 0);
    void ShowGeometry();

    // Trajectory
    char m_trajectoryName[50] = "Traj {0}";
    char m_TrajectoryFile[500] = "Enter Trajectory File";
    int m_SelectedTrajectoryType = -1;
    int m_SelectedFluidType = -1;
    int m_SelectedAmbientType = -1;

    void OnTrajectoryAdd();
    void OnTrajectoryEdit(std::pair<const std::string, std::pair<std::string, OTAP::Trajectory>> &trajectory);
    void OnTrajectoryDelete(size_t index);
    void ShowTrajectory();
};

// extern ProjectBrowser projectBrowser;
