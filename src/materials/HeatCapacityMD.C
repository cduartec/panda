/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "HeatCapacityMD.h"

// libmesh includes
#include "libmesh/quadrature.h"


registerMooseObject("pandaApp", HeatCapacityMD);

template <>
InputParameters validParams<HeatCapacityMD>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("This material calculates the heat "
                             "capacity C following an einstein relation ");
  params.addRequiredCoupledVar("temp" , "Temperature");
  params.addParam<Real>("A" , 4186.79e-6 , "Units");
  params.addParam<Real>("b_3" , 1.93e-10 , "Constant");
  params.addParam<Real>("b_2" , -7.0e-7 , "Constant");
  params.addParam<Real>("b_1" , 0.001 , "Constant");
  params.addParam<Real>("b_0" , 0.0177 , "Constant");
  return params;
}

HeatCapacityMD::HeatCapacityMD(const InputParameters & parameters)
  : Material(parameters),
    _temp(coupledValue("temp")),
    _A(getParam<Real>("A")),
    _b_3(getParam<Real>("b_3")),
    _b_2(getParam<Real>("b_2")),
    _b_1(getParam<Real>("b_1")),
    _b_0(getParam<Real>("b_0")),
    _specific_heat(declareProperty<Real>("specific_heat"))
{
}

void HeatCapacityMD::computeQpProperties()
{
    Real qp_temp = _temp[_qp];
    _specific_heat[_qp] = _A *(_b_3*pow(qp_temp,3.0) + _b_2*pow(qp_temp,2.0) + _b_1*pow(qp_temp,1.0) + _b_0*pow(qp_temp,0.0));
}
