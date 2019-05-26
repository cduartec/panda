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
  params.addRequiredCoupledVar("gc","Fracture toughness");
  params.addParam<Real>("prop_value" , 2.0e-3 , "Property value");
  return params;
}

RandomProp::RandomProp(const InputParameters & parameters)
  : Material(parameters),
    _gc(coupledValue("gc")),
    _prop_value(getParam<Real>("prop_value")),
    _gc_prop(declareProperty<Real>("gc_prop"))
{
}

void RandomProp::computeQpProperties()
{  
  Real qp_gc = _gc[_qp]; 
  _gc_prop[_qp] = qp_gc;
}
