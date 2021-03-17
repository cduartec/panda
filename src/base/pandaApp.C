#include "pandaApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
pandaApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

pandaApp::pandaApp(InputParameters parameters) : MooseApp(parameters)
{
  pandaApp::registerAll(_factory, _action_factory, _syntax);
}

pandaApp::~pandaApp() {}

void
pandaApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"pandaApp"});
  Registry::registerActionsTo(af, {"pandaApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
pandaApp::registerApps()
{
  registerApp(pandaApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
pandaApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  pandaApp::registerAll(f, af, s);
}
extern "C" void
pandaApp__registerApps()
{
  pandaApp::registerApps();
}
