//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEPRESSUREDEPENDENTELASTICITYTENSOR_H
#define COMPUTEPRESSUREDEPENDENTELASTICITYTENSOR_H

#include "ComputeElasticityTensor.h"
#include "ElementPropertyReadFile.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"


class ComputePressureDependentElasticityTensor;

template <>
InputParameters validParams<ComputePressureDependentElasticityTensor>();

/**
 * ComputeElasticityTensor defines an elasticity tensor material object as a function of
 * pressure field.
 */
class ComputePressureDependentElasticityTensor : public ComputeElasticityTensor
{
public:
  ComputePressureDependentElasticityTensor(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor();

  virtual void assignEulerAngles();

  const ElementPropertyReadFile * _read_prop_user_object;

  MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;

  /// Crystal Rotation Matrix
  MaterialProperty<RankTwoTensor> & _crysrot;
  /// Rotation matrix
  RotationTensor _R;

  /// Elasticity tensor coefficient 0.
  RankFourTensor _Aijkl0;
  /// Elasticity tensor coefficient 1.
  RankFourTensor _Aijkl1;
  /// Elasticity tensor coefficient 2.
  RankFourTensor _Aijkl2;
  /// Pressure  variable.
  const VariableValue & _p;
  VariableName _p_name;

  /// Derivative of elasticity tensor with respect to pressure.
  MaterialProperty<RankFourTensor> & _delasticity_tensor_dp;

};

#endif // COMPUTEPRESSUREDEPENDENTELASTICITYTENSOR_H
