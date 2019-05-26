/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef HEATCAPACITYMD_H
#define HEATCAPACITYMD_H

#include "Material.h"
#include "RankTwoTensor.h"

// Forward Declarations
class HeatCapacityMD;
class Function;

template <>
InputParameters validParams<HeatCapacityMD>();

/**
 * Heat capacity as a function of temperature einstein relation
 */
class HeatCapacityMD : public Material
{
public:
  HeatCapacityMD(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  const VariableValue & _temp;

  const Real _A;

  const Real _b_3;

  const Real _b_2;

  const Real _b_1;

  const Real _b_0;

  MaterialProperty<Real> & _specific_heat;

};

#endif // HEATCAPACITYMD_H
