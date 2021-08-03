//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEFINITESTRAINELASTICSTRESSBIRCHMURNAGHANNEW_H
#define COMPUTEFINITESTRAINELASTICSTRESSBIRCHMURNAGHANNEW_H

#include "ComputeStressBase.h"
#include "GuaranteeConsumer.h"

class ComputeFiniteStrainElasticStressBirchMurnaghanNew;

template <>
InputParameters validParams<ComputeFiniteStrainElasticStressBirchMurnaghanNew>();

/**
 * ComputeFiniteStrainElasticStressBirchMurnaghanNew computes the stress following elasticity
 * theory for finite strains
 * third-order Birch Murnaghan equation of state
 * finite strain formalism as in Luscher2017
 * no free energy calculation for efficiency
 */
class ComputeFiniteStrainElasticStressBirchMurnaghanNew : public ComputeStressBase, public GuaranteeConsumer
{
public:
  ComputeFiniteStrainElasticStressBirchMurnaghanNew(const InputParameters & parameters);

  void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

  /**
   * InitialStress Deprecation: remove this method
   *
   * Rotates initial_stress via rotation_increment.
   * In large-strain scenarios this must be used before addQpInitialStress
   */
  virtual void rotateQpInitialStress();

  // temperature
  const VariableValue & _temp;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;


  // exponent in Murnaghan EOS
  const Real _n_Murnaghan;

  // reference bulk modulus
  const Real _Bulk_Modulus_Ref;

  // KT0 prime correction to the bulk modulus (Menikoff 2001)
  const Real _Bulk_Modulus_Cor;

  // Gruneisen G (or Gamma) parameter
  const Real _G_Gruneisen;

  // Us-Up slope in Mie-Gruneisen EOS
  const Real _s_UsUp;

  // Von Neumann coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

  // Mesh spacing
  const Real _L0;

  // speed of sound
  const Real _C_s;

  // volumetric thermal exapnsion coeffcient, as in Austin Barton 2015
  const Real _thermal_expansion;

  // reference temperature, as in Luscher2017
  const Real _reference_temperature;
  const Real _reference_density;
  const MaterialProperty<Real> & _logarithmic_strain;

  /// Material property declarin minimum compresibility ratio over the history of deformation 
  MaterialProperty<Real> & _Jt_minimum;

  const MaterialProperty<Real> & _specific_heat;

  const MaterialProperty<Real> & _Q;
  // second piola-kirchoff stress tensor
  MaterialProperty<RankTwoTensor> & _pk2;

  // Lagrangian strain
  MaterialProperty<RankTwoTensor> & _lag_e;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
};

#endif // COMPUTEFINITESTRAINELASTICSTRESSBIRCHMURNAGHANNEW_H
