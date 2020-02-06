//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrystalPlasticitySlipResistanceGSSThermalSoft.h"

registerMooseObject("pandaApp", CrystalPlasticitySlipResistanceGSSThermalSoft);

template <>
InputParameters
validParams<CrystalPlasticitySlipResistanceGSSThermalSoft>()
{
  InputParameters params = validParams<CrystalPlasticitySlipResistance>();
  params.addParam<std::string>("uo_state_var_name",
                               "Name of state variable property: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addClassDescription("Phenomenological constitutive models' slip resistance base class.  "
                             "Override the virtual functions in your class");
  params.addRequiredCoupledVar("temp","Temperature");
  params.addRequiredParam<Real>("temp_ref", "Reference temperature");
  params.addRequiredParam<Real>("temp_melt", "Melting temperature");
  params.addRequiredParam<Real>("q", "Exponent for thermal softening");
  return params;
}

CrystalPlasticitySlipResistanceGSSThermalSoft::CrystalPlasticitySlipResistanceGSSThermalSoft(
    const InputParameters & parameters)
  : CrystalPlasticitySlipResistance(parameters),
    _mat_prop_state_var(
        getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_state_var_name"))),
    _temp(coupledValue("temp")),
    _temp_ref(getParam<Real>("temp_ref")),
    _temp_melt(getParam<Real>("temp_melt")),
    _q(getParam<Real>("q"))
{
}

bool
CrystalPlasticitySlipResistanceGSSThermalSoft::calcSlipResistance(unsigned int qp,
                                                       std::vector<Real> & val) const
{
  for (unsigned int i = 0; i < _variable_size; ++i)
    val[i] = _mat_prop_state_var[qp][i] * ( 1.0 - std::pow( (_temp[qp] - _temp_ref) / (_temp_melt - _temp_ref), _q));

  return true;
}
