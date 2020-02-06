//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CRYSTALPLASTICITYSLIPRESISTANCEGSSTHERMALSOFT_H
#define CRYSTALPLASTICITYSLIPRESISTANCEGSSTHERMALSOFT_H

#include "CrystalPlasticitySlipResistance.h"

class CrystalPlasticitySlipResistanceGSSThermalSoft;

template <>
InputParameters validParams<CrystalPlasticitySlipResistanceGSSThermalSoft>();

/**
 * Phenomenological constitutive model slip resistance userobject class.
 */
class CrystalPlasticitySlipResistanceGSSThermalSoft : public CrystalPlasticitySlipResistance
{
public:
  CrystalPlasticitySlipResistanceGSSThermalSoft(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

protected:
  const MaterialProperty<std::vector<Real>> & _mat_prop_state_var;
  const VariableValue & _temp;
  const Real _temp_ref;
  const Real _temp_melt;
  const Real _q;
};

#endif // CRYSTALPLASTICITYSLIPRESISTANCEGSSTHERMALSOFT_H
