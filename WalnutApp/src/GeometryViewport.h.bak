#ifndef _VIEWPORT_H
#define _VIEWPORT_H

#include "imgui.h"
#include "Walnut/Layer.h"
#include "Renderer.h"
#include "otap.h"

class GeometryViewport : public Walnut::Layer
{
private:
    // Renderer m_Renderer;
    OTAP::Engine &m_state = OTAP::State::GetInstance();

public:
    virtual void OnUIRender() override;
};

#endif //_VIEWPORT_H