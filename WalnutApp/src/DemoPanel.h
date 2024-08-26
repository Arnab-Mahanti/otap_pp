#include "Walnut/Application.h"
#include "Walnut/EntryPoint.h"

#include "Walnut/Image.h"
#include "Renderer.h"
#include "Problem.h"
#include "Walnut/Timer.h"

using namespace Walnut;

class DemoPanel : public Walnut::Layer
{
public:
	virtual void OnUIRender() override
	{
		ImGui::ShowDemoWindow();
	}

private:
};


