#include "pandaApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"





template <>
InputParameters
validParams<pandaApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

pandaApp::pandaApp(InputParameters parameters) : MooseApp(parameters)
{
  pandaApp::registerAll(_factory, _action_factory, _syntax);
}

pandaApp::~pandaApp() {}

void
pandaApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
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
