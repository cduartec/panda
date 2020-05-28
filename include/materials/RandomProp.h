/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef RANDOMPROP_H
#define RANDOMPROP_H

#include "Material.h"
#include "RankTwoTensor.h"

// Forward Declarations
class RandomProp;
class Function;

template <>
InputParameters validParams<RandomProp>();

/**
 * Random value from a normal distribution
 */
class RandomProp : public Material
{
public:
  RandomProp(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  
  const Real _prop_value;
  const Real _prop_variation;

  MaterialProperty<Real> & _normal_prop;

};

#endif // RANDOMPROP_H
