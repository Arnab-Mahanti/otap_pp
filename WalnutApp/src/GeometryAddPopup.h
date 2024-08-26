

#pragma once
#include "imrad.h"

class GeometryAddPopup
{
public:
    /// @begin interface
    void OpenPopup(std::function<void(ImRad::ModalResult)> clb = [](ImRad::ModalResult) {});
    void ClosePopup(ImRad::ModalResult mr);
    void Draw();

    int m_SelectedGeometryType = -1;
    double m_GeometryLength;
    float m_GeometryAngle;
    double m_GeometryRadius;
    float m_GeometrySweep;
    /// @end interface

private:
    /// @begin impl
    void Init();
    std::string GetGeometryTypes() const;

    ImGuiID ID = 0;
    ImRad::ModalResult modalResult;
    ImRad::Animator animator;
    std::function<void(ImRad::ModalResult)> callback;

    bool isOpen = false;

    void PushGeometry();
    /// @end impl
};

extern GeometryAddPopup geometryAddPopup;
