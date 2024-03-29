/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINCPPFFRACTUREMIEGRUNEISENCONSTH_H
#define FINITESTRAINCPPFFRACTUREMIEGRUNEISENCONSTH_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCPPFFractureMieGruneisenConstH uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 * Principal values of the Lagrangian strain used to calculate damage with Amor 2009 model
 * Calculation of the broken plastic energy for temperature calculation
 * Mie Gruneisen equation of state
 * finite strain formalism as in Luscher2017
 * Computes the stress and free energy derivatives for the phase field
 * Allen-Cahn formalism
 */
class FiniteStrainCPPFFractureMieGruneisenConstH;

template<>
InputParameters validParams<FiniteStrainCPPFFractureMieGruneisenConstH>();

class FiniteStrainCPPFFractureMieGruneisenConstH : public FiniteStrainCrystalPlasticity
{
public:
  FiniteStrainCPPFFractureMieGruneisenConstH(const InputParameters & parameters);

protected:
  /// Function required to initialize statefull material properties
  virtual void initQpStatefulProperties();
  /**
   * This function set variables for internal variable solve.
   */
  virtual void preSolveStatevar();

  /**
   * This function solves internal variables.
   */
  virtual void solveStatevar();

  /**
   * This function update internal variable after solve.
   */
  virtual void postSolveStatevar();


  // update slip system resistances and output slip increment
  virtual void updateGss();

  /**
   * Update elastic and plastic work
   */
  virtual void update_energies();

  /**
   * This function calculate stress residual.
   */
  virtual void calcResidual( RankTwoTensor & );

  virtual void getSlipIncrements();

  // damage phase field
  const VariableValue & _c;

  // temperature
  const VariableValue & _temp;

  /// Small number to avoid non-positive definiteness at or near complete damage
  const Real _kdamage;

  /// Use current value of history variable
  bool _use_current_hist;

  // Gruneisen G (or Gamma) parameter
  const Real _G_Gruneisen;

  // Us-Up slope in Mie-Gruneisen EOS
  const Real _s_UsUp;

  // reference bulk modulus
  const Real _Bulk_Modulus_Ref;

  // Von Neumann coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

  // Average element size
  const Real _h_e;

  // Speed of sound
  const Real _c_l;

  // reference temperature, as in Luscher2017
  const Real _reference_temperature;

  // prefactor of the plastic contribution to damage
  const Real _plastic_factor;
  
  //Specific heat
  const MaterialProperty<Real> & _specific_heat;
  
  //Density
  const MaterialProperty<Real> & _density;

  /// Material property defining crack width, declared elsewhere
  const MaterialProperty<Real> & _l;

  /// Material property defining gc parameter, declared elsewhere
  const MaterialProperty<Real> & _gc;

  /// Total energy and derivatives, declared in this material
  MaterialProperty<Real> & _F;
  MaterialProperty<Real> & _dFdc;
  MaterialProperty<Real> & _d2Fdc2;

  /// Elastic and plastic energies
  MaterialProperty<Real> & _W0e_pos;
  MaterialProperty<Real> & _W0e_neg;
  MaterialProperty<Real> & _W0p;
  const MaterialProperty<Real> & _W0p_old;
  MaterialProperty<Real> & _W0p_broken;
  const MaterialProperty<Real> & _W0p_broken_old;

  /// Total energy derivatives, declared in this material
  MaterialProperty<RankTwoTensor> & _dstress_dc;
  MaterialProperty<RankTwoTensor> & _d2Fdcdstrain;

  /// History variable that prevents crack healing, declared in this material
  MaterialProperty<Real> & _hist;

  /// Old value of history variable
  const MaterialProperty<Real> & _hist_old;

  MaterialProperty<RankTwoTensor> & _pk2_undamaged;
  MaterialProperty<RankTwoTensor> & _fe_out; // Elastic deformation gradient for output
  MaterialProperty<std::vector<Real>> & _slip_incr_out; // slip increment output
  MaterialProperty<std::vector<Real>> & _tau_out; // slip increment output

  /// Pressure  variable.
  const VariableValue & _p;
  VariableName _p_name;

  /// Heat rates
  MaterialProperty<Real> &  _heat_rate_vis; //Heat rate due to shock dissipation
  MaterialProperty<Real> &  _heat_rate_therm1; //Heat due to thermo-elastic coupling volumetric
  MaterialProperty<Real> &  _heat_rate_therm2; //Heat due to thermo-elastic coupling coupled
  MaterialProperty<Real> &  _heat_rate_p; //Heat due to plasticity
 
  /// Lagrangian strain rate
  MaterialProperty<Real> &  _Tr_E_dot; //Tr(C^{-1}_e\dot{E}_e)

  Real _W0p_tmp;
  Real _W0p_tmp_old;
  Real _W0p_broken_tmp;
  Real _W0p_broken_tmp_old;

};

#endif //FINITESTRAINCRYSTALPLASTICITYPFFRACTURESTRESSMIEGRUNEISENCONST_H
