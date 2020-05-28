/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "RandomProp.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "MooseRandom.h"
registerMooseObject("pandaApp",RandomProp);

template <>
InputParameters validParams<RandomProp>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("This material generates random values"
                             "from normal distribution");
  params.addParam<Real>("prop_value" , 2.0e-3 , "property mean value");
  params.addParam<Real>("prop_variation", 0.0,"stantdard deviation");
  return params;
}

RandomProp::RandomProp(const InputParameters & parameters)
  : Material(parameters),
    _prop_value(getParam<Real>("prop_value")),
    _prop_variation(getParam<Real>("prop_variation")),
    _normal_prop(declareProperty<Real>("normal_prop"))
{

}

void RandomProp::computeQpProperties()
{  
  _normal_prop[_qp] = 0.0;
  _normal_prop[_qp] =  _prop_value +  getRandomReal() * _prop_value;
}
