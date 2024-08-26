#include "Walnut/Application.h"
#include "Walnut/EntryPoint.h"

#include "Walnut/Image.h"
#include "Renderer.h"
#include "Problem.h"
#include "Walnut/Timer.h"

#include "ProjectBrowser.h"
#include "DemoPanel.h"
#include "GeometryViewport.h"
#include "ResultsViewport.h"

using namespace Walnut;

Walnut::Application *Walnut::CreateApplication(int argc, char **argv)
{
	Walnut::ApplicationSpecification spec;
	spec.Name = "OTAP++";

	Walnut::Application *app = new Walnut::Application(spec);
	app->PushLayer<ProjectBrowser>();
	app->PushLayer<DemoPanel>();
	app->PushLayer<GeometryViewport>();
	app->PushLayer<ResultsViewport>();


	app->SetMenubarCallback([app]()
							{
		if (ImGui::BeginMenu("File"))
		{
			if(ImGui::MenuItem("Switch Theme"))
				app->SwitchTheme();
			
			if (ImGui::MenuItem("Exit"))
				app->Close();

				ImGui::EndMenu();
		} });
	return app;
}