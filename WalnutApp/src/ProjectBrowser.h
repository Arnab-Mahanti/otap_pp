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
    int m_SelectedTrajectoryType = -1;
    int m_SelectedGeometryType = -1;
    double m_GeometryLength = 0.0;
    double m_GeometryRadius = 0.0;
    float m_GeometryAngle = 0.0f;
    float m_GeometrySweep = 0.0f;
    char m_geometryName[50] = "Geom {0}";

    void OnGeometryAdd();
    void OnGeometryDelete(size_t index);
    void OnGeometryPrimitiveAdd(OTAP::Geometry &geometry, bool edit = false, size_t editindex = 0);
    void ShowGeometry();

    // Trajectory
    std::unordered_map<std::string, OTAP::Trajectory> m_Trajectories;
    std::string m_TrajectoryFile = "Enter Trajectory File";

    void OnTrajectoryAdd();
    void OnTrajectoryDelete(size_t index);
    void ShowTrajectory();
};

// extern ProjectBrowser projectBrowser;
