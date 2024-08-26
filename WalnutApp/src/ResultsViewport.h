#ifndef _RESULTS_VIEWPORT_H
#define _RESULTS_VIEWPORT_H

#include "imgui.h"
#include "Walnut/Layer.h"
#include "Renderer.h"

class ResultsViewport : public Walnut::Layer
{
private:
    Renderer m_Renderer;
public:
    virtual void OnUIRender() override;
};


#endif //_RESULTS_VIEWPORT_H