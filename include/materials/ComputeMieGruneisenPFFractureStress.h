/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State

#pragma once

#include "ComputeStressBase.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "PiecewiseBilinear.h"
#include "RankFourTensor.h"
#include "ExternalPetscSolverApp.h"
#include "petscblaslapack.h"

class ComputeMieGruneisenPFFractureStress : public ComputeStressBase
{
public:
  static InputParameters validParams();

  ComputeMieGruneisenPFFractureStress(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  /// Coupled order parameter defining the crack
  const VariableValue & _c;

  /// Material property defining crack width, declared elsewhere
  const MaterialProperty<Real> & _l;

  /// Material property defining pressure, declared elsewhere
  const MaterialProperty<Real> & _pressure;

  /// Material property defining gc parameter, declared elsewhere
  const MaterialProperty<Real> & _gc;

  /// Use current value of history variable
  bool _use_current_hist;

  /// Use PETSc's VI (Reduced space active set solvers for variational inequalities based on Newton's method) solver
  bool _use_snes_vi_solver;

  /// History variable that prevents crack healing, declared in this material
  MaterialProperty<Real> & _H;

  /// Old value of history variable
  const MaterialProperty<Real> & _H_old;

  /// material property for fracture energy barrier
  const MaterialProperty<Real> & _barrier;

  /// Material property for elastic energy
  MaterialProperty<Real> & _E;

  /// Derivative of elastic energy w.r.t damage variable
  MaterialProperty<Real> & _dEdc;

  /// Second-order derivative of elastic energy w.r.t damage variable
  MaterialProperty<Real> & _d2Ed2c;

  /// Derivative of stress w.r.t damage variable
  MaterialProperty<RankTwoTensor> & _dstress_dc;

  /// Second-order derivative of elastic energy w.r.t damage variable and strain
  MaterialProperty<RankTwoTensor> & _d2Fdcdstrain;

  /// Material property for energetic degradation function
  const MaterialProperty<Real> & _D;

  /// Derivative of degradation function w.r.t damage variable
  const MaterialProperty<Real> & _dDdc;

  /// Second-order derivative of degradation w.r.t damage variable
  const MaterialProperty<Real> & _d2Dd2c;

  /// Material property for damage indicator function
  const MaterialProperty<Real> & _I;

  /// Derivative of damage indicator function w.r.t damage variable
  const MaterialProperty<Real> & _dIdc;

  /// Second-order derivative of damage indicator function w.r.t damage variable
  const MaterialProperty<Real> & _d2Id2c;

  const Real _Gamma;

  const MaterialProperty<Real> & _density;

  const MaterialProperty<Real> & _specific_heat;

  const VariableValue & _temperature;
  const Real _ref_temperature;

  const Real _s;

  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  /// Damping coefficients
  const Real _C0;
  const Real _C1;

  /// Maximum element size
  const Real _Le;

  /// Material property defining speed of sound in material
  const Real _sound_speed;

  MaterialProperty<Real> & _bulk_modulus;
  MaterialProperty<Real> & _pressure_eos;

  MaterialProperty<RankTwoTensor> & _stress_eos_elastic;
  MaterialProperty<RankTwoTensor> & _stress_cpl_elastic;

  virtual void computeQpStress() override;
};
