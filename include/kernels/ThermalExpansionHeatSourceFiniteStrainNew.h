#ifndef THERMALEXPANSIONHEATSOURCEFINITESTRAINNEW_H
#define THERMALEXPANSIONHEATSOURCEFINITESTRAINNEW_H

#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class ThermalExpansionHeatSourceFiniteStrainNew;

template <>
InputParameters validParams<ThermalExpansionHeatSourceFiniteStrainNew>();

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Luscher2017 model and Amor damage model
 */
class ThermalExpansionHeatSourceFiniteStrainNew : public HeatSource
{
public:
  ThermalExpansionHeatSourceFiniteStrainNew(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;


  // Reference bulk modulus in the equation of state
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

  // reference temperature with zero thermal expansion
  const Real _reference_temperature;
  const Real _reference_density;
  const MaterialProperty<Real> & _logarithmic_strain;
  const MaterialProperty<Real> & _specific_heat;
  // const MaterialProperty<Real> & _p;

  /// Pressure  variable.
  const VariableValue & _p;
  VariableName _p_name;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient; // deformation gradient
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old; // deformation gradient, previous timestep

  // const MaterialProperty<RankTwoTensor> & _fp; // plastic deformation gradient
  // const MaterialProperty<RankTwoTensor> & _fp_old; // plastic deformation gradient, previous timestep

  const MaterialProperty<RankFourTensor> & _elasticity_tensor; //elasticity tensor


  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // THERMALEXPANSIONHEATSOURCEFINITESTRAINNEW_H
