/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTECRACKFRICTIONHEATENERGYDIENES3D_H
#define COMPUTECRACKFRICTIONHEATENERGYDIENES3D_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "DerivativeMaterialInterface.h"
/**
 * ComputeCrackFrictionHeatEnergyDienes computes the energy from crack friction
 * and, if currentlyComputingJacobian, then the derivative of this quantity wrt total strain
 */
class ComputeCrackFrictionHeatEnergyDienes3D : public DerivativeMaterialInterface<Material>
{
public:
  ComputeCrackFrictionHeatEnergyDienes3D(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// optional parameter that allows multiple mechanics materials to be defined
  std::string _base_name;

  const VariableValue & _c;
  const VariableValue & _dcdx;
  const VariableValue & _dcdy;
  const VariableValue & _dcdz;

  MaterialProperty<std::vector<Real>> & _crack_normal;
  MaterialProperty<Real> & _crack_normal_norm;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _strain_rate;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;

  const Real _friction_coefficient;
  const Real _l;

  MaterialProperty<std::vector<Real>> & _friction_force;
  MaterialProperty<Real> & _friction_normal_force;
  MaterialProperty<std::vector<Real>> & _slide_velocity;
  MaterialProperty<Real> & _slide_velocity_parallel;

  MaterialProperty<Real> & _crack_surface_density;
  MaterialProperty<Real> & _heat_source_rate;

  MaterialProperty<RankTwoTensor> & _velocity_gradient;

};

#endif // COMPUTECRACKFRICTIONHEATENERGYDIENES3D_H
