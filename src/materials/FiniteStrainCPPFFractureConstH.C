#include "FiniteStrainCPPFFractureConstH.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"
#include "MathUtils.h"

#include <fstream>
#include <cmath>

registerMooseObject("pandaApp", FiniteStrainCPPFFractureConstH);

template<>
InputParameters validParams<FiniteStrainCPPFFractureConstH>()
{
  InputParameters params = validParams<FiniteStrainCrystalPlasticity>();
  params.addClassDescription("Crystal Plasticity class. Damage split volumetric and coupled volumetric/deviatoric."
                             "Calculation of the broken plastic energy for temperature calculation,"
                             "Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)"
                             "Finite strain formalism as in (Luscher 2017)."
                             "Computes the stress and free energy derivatives for the phase field "
                             "fracture model.");
  params.addRequiredCoupledVar("c","Damage variable");
  params.addRequiredCoupledVar("temp","Temperature");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<MaterialPropertyName>(
      "F_name", "E_elpl", "Name of material property storing the elastic and plastic energy driving damage");
  params.addRequiredParam<Real>("G_Gruneisen", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("s_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "Reference bulk modulus");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient, artificial viscosity");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient, artificial viscosity");
  params.addRequiredParam<Real>("h_e", "Average element size mesh");
  params.addRequiredParam<Real>("c_l", "Speed of sound");
  params.addRequiredParam<Real>("reference_temperature", "Reference temperature for thermal expansion");
  params.addRequiredParam<Real>("plastic_factor","Prefactor of the plastic contribution to damage");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  params.addRequiredCoupledVar("p", "Pressure");
  params.addParam<MaterialPropertyName>("delasticity_tensor_dp", "Derivative of C_ijkl  with pressure");  
  return params;
}

FiniteStrainCPPFFractureConstH::FiniteStrainCPPFFractureConstH(const InputParameters & parameters) :
    FiniteStrainCrystalPlasticity(parameters),
    _c(coupledValue("c")),
    _temp(coupledValue("temp")),
    _kdamage(getParam<Real>("kdamage")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _G_Gruneisen(getParam<Real>("G_Gruneisen")),
    _s_UsUp(getParam<Real>("s_UsUp")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _h_e(getParam<Real>("h_e")),
    _c_l(getParam<Real>("c_l")),
    _reference_temperature(getParam<Real>("reference_temperature")), 
    _plastic_factor(getParam<Real>("plastic_factor")), 
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density_name")),
    _l(getMaterialProperty<Real>("l")),
    _gc(getMaterialProperty<Real>("gc_prop")),
    _F(declareProperty<Real>(getParam<MaterialPropertyName>("F_name"))),
    _dFdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("F_name"),
                                          getVar("c", 0)->name())),
    _d2Fdc2(declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("F_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _W0e_pos(declareProperty<Real>("W0e_pos")), // positive (= damaging) elastic energy
    _W0e_neg(declareProperty<Real>("W0e_neg")), // negative (= non-damaging) elastic energy
    _W0p(declareProperty<Real>("W0p")), // plastic energy (unbroken)
    _W0p_old(getMaterialPropertyOld<Real>("W0p")), // plastic energy of previous increment (unbroken)
    _W0p_broken(declareProperty<Real>("W0p_broken")), // plastic energy (broken)
    _W0p_broken_old(getMaterialPropertyOld<Real>("W0p_broken")), // plastic energy of previous increment (broken)
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    _hist(declareProperty<Real>("hist")), // history variable = (never decreasing) positive elastic energy
    _hist_old(getMaterialPropertyOld<Real>("hist")), // history variable = (never decreasing) positive elastic energy
    _pk2_undamaged(declareProperty<RankTwoTensor>("pk2_undamaged")), // undamaged 2nd Piola Kirchoff Stress
    _fe_out(declareProperty<RankTwoTensor>("fe_out")), // Elastic deformation gradient for output
    _slip_incr_out(declareProperty<std::vector<Real>>("slip_incr_out")), // slip system strain increment for output
    _tau_out(declareProperty<std::vector<Real>>("tau_out")), // slip system strain increment for output
    _p(coupledValue("p")),
    _p_name(getVar("p", 0)->name()),
    _heat_rate_vis(declareProperty<Real>("heat_rate_vis")), //Heat rate due to plastic dissipation
    _heat_rate_therm1(declareProperty<Real>("heat_rate_therm1")), // Heat rate due to thermo-elastic coupling
    _heat_rate_therm2(declareProperty<Real>("heat_rate_therm2")), // Heat rate due to thermo-elastic coupling
    _heat_rate_p(declareProperty<Real>("heat_rate_p")), // Heat rate due to plasticity
    _Tr_E_dot(declareProperty<Real>("Tr_E_dot")) // Trace Lagrangian strain rate
{
    _err_tol = false;
}

void
FiniteStrainCPPFFractureConstH::initQpStatefulProperties()
{
  _stress[_qp].zero();

  _fp[_qp].zero();
  _fp[_qp].addIa(1.0);

  _pk2[_qp].zero();
  _acc_slip[_qp] = 0.0;
  _lag_e[_qp].zero();

  _update_rot[_qp].zero();
  _update_rot[_qp].addIa(1.0);

  _hist[_qp] = 0.0; // history variable = (never decreasing) positive elastic energy
  
  //Heat sources
  _heat_rate_vis[_qp] = 0.0;
  _heat_rate_therm1[_qp] = 0.0;
  _heat_rate_therm2[_qp] = 0.0;
  _heat_rate_p[_qp] = 0.0;

  _Tr_E_dot[_qp] = 0.0;
  
  initSlipSysProps(); // Initializes slip system related properties
  initAdditionalProps();
}

void
FiniteStrainCPPFFractureConstH::preSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _gss_tmp = _gss_old[_qp];
    _W0p_tmp = _W0p_old[_qp];
    _W0p_broken_tmp = _W0p_broken_old[_qp];
    _accslip_tmp_old = _acc_slip_old[_qp];
  }
  else
  {
    if (_first_step_iter)
    {
      _gss_tmp = _gss_tmp_old = _gss_old[_qp];
      _W0p_tmp = _W0p_tmp_old = _W0p_old[_qp];
      _W0p_broken_tmp = _W0p_broken_tmp_old = _W0p_broken_old[_qp];
      _accslip_tmp_old = _acc_slip_old[_qp];
    }
    else
    {
      _gss_tmp = _gss_tmp_old;
      _W0p_tmp = _W0p_tmp_old;
      _W0p_broken_tmp = _W0p_broken_tmp_old;
    }
  }
}

void
FiniteStrainCPPFFractureConstH::solveStatevar()
{
  Real gmax, gdiff;
  unsigned int iterg;
  std::vector<Real> gss_prev(_nss);

  gmax = 1.1 * _gtol;
  iterg = 0;

  while (gmax > _gtol && iterg < _maxiterg) // Check for slip system resistance update tolerance
  {
    preSolveStress();
    solveStress();
    if (_err_tol)
      return;

    update_energies(); // update elastic, plastic and total energies

    postSolveStress(); // Update _fp_old_inv = _fp_old

    gss_prev = _gss_tmp;

    update_slip_system_resistance(); // Update slip system resistance

    gmax = 0.0;
    for (unsigned i = 0; i < _nss; ++i)
    {
      gdiff = std::abs(gss_prev[i] - _gss_tmp[i]); // Calculate increment size

      if (gdiff > gmax)
        gmax = gdiff;
    }
    iterg++;
  }

  if (iterg == _maxiterg)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainCPPFFractureConstH: Hardness Integration error gmax", gmax, "\n");
#endif
    _err_tol = true;
  }
}

void
FiniteStrainCPPFFractureConstH::postSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _gss[_qp] = _gss_tmp;
    _W0p[_qp] = _W0p_tmp;
    _W0p_broken[_qp] = _W0p_broken_tmp;
    _acc_slip[_qp] = _accslip_tmp;
  }
  else
  {
    if (_last_step_iter)
    {
      _gss[_qp] = _gss_tmp;
      _W0p[_qp] = _W0p_tmp;
      _W0p_broken[_qp] = _W0p_broken_tmp;
      _acc_slip[_qp] = _accslip_tmp;
    }
    else
    {
      _gss_tmp_old = _gss_tmp;
      _W0p_tmp_old = _W0p_tmp;
      _W0p_broken_tmp_old = _W0p_broken_tmp;
      _accslip_tmp_old = _accslip_tmp;
    }
  }
}

/**
 * Old function to update slip system resistances.
 * Kept to avoid code break at computeQpstress
 * output slip increment
 */
void
FiniteStrainCPPFFractureConstH::updateGss()
{
  DenseVector<Real> hb(_nss);
  Real qab;

  Real a = _hprops[4]; // Kalidindi

  _slip_incr_out[_qp].resize(_nss);

  _accslip_tmp = _accslip_tmp_old;
  for (unsigned int i = 0; i < _nss; ++i)
    _accslip_tmp += std::abs(_slip_incr(i));

  for (unsigned int i = 0; i < _nss; ++i)
    _slip_incr_out[_qp][i] = _slip_incr(i);

  for (unsigned int i = 0; i < _nss; ++i)
    hb(i) = _h0 * std::pow(std::abs(1.0 - _gss_tmp[i] / _tau_sat), a) *
            copysign(1.0, 1.0 - _gss_tmp[i] / _tau_sat);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    if (_max_substep_iter == 1) // No substepping
      _gss_tmp[i] = _gss_old[_qp][i];
    else
      _gss_tmp[i] = _gss_tmp_old[i];

    for (unsigned int j = 0; j < _nss; ++j)
    {
      unsigned int iplane, jplane;
      iplane = i / 3;
      jplane = j / 3;

      if (iplane == jplane) // Kalidindi
        qab = 1.0;
      else
        qab = _r;

      _gss_tmp[i] += qab * hb(j) * std::abs(_slip_incr(j));
      _dgss_dsliprate(i, j) = qab * hb(j) * copysign(1.0, _slip_incr(j)) * _dt;
    }
  }
}

// Update slip system resistance, elastic, plastic and total work
void
FiniteStrainCPPFFractureConstH::update_energies()
{
  RankTwoTensor cauchy_stress_undamaged, cauchy_stress, WpToTrace, WpBrokenToTrace, invFe;
  Real detFe;
  Real c = _c[_qp];


  Real hist_variable = _hist_old[_qp]; // history variable = (never decreasing) positive elastic energy
  if (_use_current_hist)
    hist_variable = _hist[_qp];

  if (_max_substep_iter == 1) //No substepping
  {
    _W0p_tmp = _W0p_old[_qp];
    _W0p_broken_tmp = _W0p_broken_old[_qp];
  }
  else
  {
    _W0p_tmp = _W0p_tmp_old;
    _W0p_broken_tmp = _W0p_broken_tmp_old;
  }

  // Update plastic work
  detFe = _fe.det();
  invFe = _fe.inverse();

  // _pk2[_qp] is the updated piola-kirchoff
  cauchy_stress = _fe * _pk2[_qp] * _fe.transpose()/detFe;
  cauchy_stress_undamaged = _fe * _pk2_undamaged[_qp] * _fe.transpose()/detFe;

  WpBrokenToTrace = cauchy_stress * _fe * ( _fp_inv.inverse() - _fp_old_inv.inverse() ) * _fp_inv * invFe * detFe;
  WpToTrace = cauchy_stress_undamaged * _fe * ( _fp_inv.inverse() - _fp_old_inv.inverse() ) * _fp_inv * invFe * detFe;

  _W0p_broken_tmp += WpBrokenToTrace.trace();
  _W0p_tmp += WpToTrace.trace();

  // Total free energy density 
  _F[_qp] = (hist_variable + _plastic_factor * _W0p_tmp) * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) +
            _W0e_neg[_qp] + _gc[_qp] / (2 * _l[_qp]) * c * c;

  // derivative of total free energy density wrt c
  _dFdc[_qp] = -(hist_variable + _plastic_factor * _W0p_tmp) * 2.0 * (1.0 - c) * (1 - _kdamage) +
               _gc[_qp] / _l[_qp] * c;

  // 2nd derivative of total free energy density wrt c
  _d2Fdc2[_qp] = (hist_variable + _plastic_factor * _W0p_tmp) * 2.0 * (1 - _kdamage) + _gc[_qp] / _l[_qp];

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history variable
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = -cauchy_stress_undamaged * (1.0 - c) * (1 - _kdamage);

}

void
FiniteStrainCPPFFractureConstH::calcResidual( RankTwoTensor &resid )
{
  RankTwoTensor iden, ce, invce, ee, ce_pk2, eqv_slip_incr, pk2_new, temporal;
  RankTwoTensor ce_old, fe_old, ee_old, ee_rate, invce_ee_rate, inv_fp_old;
  RankTwoTensor thermal_eigenstrain, pk2_new_vol;
  Real trD, Kb, Je, detFe;
  Real c = _c[_qp];
  Real temp = _temp[_qp];
  Real xfac = ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage);
  Real thermal_expansion_coeff; 
  detFe = _fe.det();
  
  //Bulk modulus
  // Kb = K in Luscher2017
  // Kb = (1/9) I : C : I
  //Real Kb = 0.0;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      Kb +=  _elasticity_tensor[_qp](i, i, j, j);
  Kb = (1.0 / 9.0) * Kb;

  // Thermal expansion coefficient depends on Gruneisen parameter, bulk modulus and sound speed
  thermal_expansion_coeff = _G_Gruneisen * _density[_qp] * _specific_heat[_qp] / Kb;

  iden.zero();
  temporal.zero();
  iden.addIa(1.0);

  //Elastic part deformation gradient
  _fe = _dfgrd_tmp * _fp_prev_inv; // _fp_inv  ==> _fp_prev_inv

  //Right cauchy strain tensor 
  ce = _fe.transpose() * _fe;
  ce_pk2 = ce * _pk2_tmp;
  ce_pk2 = ce_pk2 / _fe.det();

  //Calculate Schmid tensor and resolved shear stresses

  _tau_out[_qp].resize(_nss);//Resolved shear stress for output

  for (unsigned int i = 0; i < _nss; ++i){
    _tau(i) = ce_pk2.doubleContraction(_s0[i]);
    _tau_out[_qp][i] = _tau(i);}

  getSlipIncrements(); // Calculate dslip,dslipdtau

  if (_err_tol)
    return;

  eqv_slip_incr.zero();

  for (unsigned int i = 0; i < _nss; ++i)
    eqv_slip_incr += _s0[i] * _slip_incr(i);

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;
  _fe = _dfgrd_tmp * _fp_inv;
  _fe_out[_qp] = _fe; // Elastic deformation gradient for output
  
  ce = _fe.transpose() * _fe;
  invce = ce.inverse();
  ee = 0.5 * ( ce - iden );
  Je = _fe.det(); // Jacobian = relative volume
  
  // Thermal eigenstrain (equation (18) in Luscher2017-------------------------------------------------
  // Lagrangian strain E_thermal = 1/2 (F_thermal^T F_thermal - I)
  // finite strain formula (Lubarda2002): F_thermal = exp((alpha/3)*(T-T_ref)) I
  thermal_eigenstrain = (1.0 / 2.0)
                      * (std::exp((2.0/3.0) * thermal_expansion_coeff * (temp - _reference_temperature)) - 1.0)
                      * iden;

  //Second Piola-Kirchoff stress----------------------------------------------------------------------
  //Volumetric stress + damage-------------------------------------------------------------------------
  //Formula (Luscher2017):  K*J^(2/3)*(delta - alpha_v)*C^{-1}_e 
  Real delta;
  delta = 1.5 * (std::pow(Je , 2.0/3.0) - 1.0);
  pk2_new_vol = Kb * std::pow(Je , 2.0/3.0)
                   * (delta  -  thermal_eigenstrain.trace())
                   * invce;
  if (Je > 1.0){
  pk2_new = xfac * pk2_new_vol;
  } else {
  pk2_new = pk2_new_vol;
  }
  _pk2_undamaged[_qp] = pk2_new_vol;
  
 // Deviatoric stress + damage (equation (18) in Luscher2017): 
 // Formula (Eq. (18) in Luscher2017) : C_ijkl * (E_e - alpha)-K*J^(2/3)*(delta - alpha_v)*C^{-1}_e
  pk2_new += xfac * _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);
  _pk2_undamaged[_qp] += _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);
  pk2_new +=  -1.0 * xfac * Kb * std::pow(Je , 2.0/3.0)
                   * (delta  -  thermal_eigenstrain.trace())
                   * invce;
  _pk2_undamaged[_qp] += -1.0 * Kb * std::pow(Je , 2.0/3.0)
                          * (delta  - thermal_eigenstrain.trace())
                          * invce;
  
  // Calculate bulk viscosity damping-----------------------------------------------------------
  // C0 * dot(J) / J * |dot(J) / J| + C1 * dot(J) / J
  // C0 should be chosen of the order of rho * Le^2, rho = density, Le = element size
  // C1 should be chosen of the order of rho * Le * cs, cs = sound speed
  // Maheo et al. Mechanics Research Communications 38 (2011) 81 88
  trD = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  trD /= _deformation_gradient_old[_qp].det();
  Real J_dot;
  J_dot = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  if (Je <  1.0){
  pk2_new.addIa( _C0 * trD * _h_e * _h_e * std::abs(trD) * _density[_qp] );
  pk2_new.addIa( _C1 * trD * _h_e * _density[_qp] * _c_l);
  }
  //_pk2_undamaged.addIa( _C0 * trD * _h_e * _h_e * std::abs(trD) * _density[_qp] );
 // _pk2_undamaged.addIa( _C1 * trD * _h_e * _density[_qp] * _c_l);
   

  // Free energy calculation ---------------------------------------------------------------------
  // Energy that contributes to damage \Psi^{+}
  _W0e_pos[_qp] = 0.0;

  // Energy that DOES NOT contribute to damaga \Psi^{-} 
  _W0e_neg[_qp] = 0.0;

  // Volumetric free energy = Psi_EOS in Luscher2017------------------------------------------------

  if (Je >= 1.0) {// In expansion 
     _W0e_pos[_qp] = 1.0/2.0 * Kb * delta * delta - Kb * delta * thermal_eigenstrain.trace();
  } else {// In compression 
    _W0e_neg[_qp] = 1.0/2.0 * Kb * delta * delta - Kb * delta * thermal_eigenstrain.trace();
  }

  // Volumetric-Deviatoric (coupled) free energy or Psi_cpl (Luscher2017)
  RankTwoTensor elastic_energy_tensor, thermal_coupling_tensor;
  elastic_energy_tensor = _elasticity_tensor[_qp] * ee;
  elastic_energy_tensor = 0.5 * ee * elastic_energy_tensor;
  thermal_coupling_tensor = _elasticity_tensor[_qp] * thermal_eigenstrain; 
  thermal_coupling_tensor = ee * thermal_coupling_tensor;

  _W0e_pos[_qp] += elastic_energy_tensor.trace(); // 1/2 * Ee : C : Ee
  _W0e_pos[_qp] -= thermal_coupling_tensor.trace(); // - Ee : C : alpha in equation 15 of Luscher2017
  _W0e_pos[_qp] -= 0.5 * Kb * delta * delta; // - 1/2 K delta^2 in equation 15 of Luscher2017
  _W0e_pos[_qp] += Kb * delta * thermal_eigenstrain.trace(); // + K delta alpha_v in equation 15 of Luscher2017

  // Assign history variable and derivative
  if (_W0e_pos[_qp] > _hist_old[_qp])
    _hist[_qp] = _W0e_pos[_qp];
  else
    _hist[_qp] = _hist_old[_qp];

  // Used in StressDivergencePFFracTensors Jacobian
  // same approximation as above
  if (c < 1.0)
    _dstress_dc[_qp] = -_pk2_undamaged[_qp] * 2.0 * (1.0 - c) * (1 - _kdamage);
  else
    _dstress_dc[_qp].zero();


  // Calculate heat rates-----------------------------------------------------------
  // Heat rate due to shock dissipation
  inv_fp_old = _fp_old[_qp].inverse();
  fe_old = _deformation_gradient_old[_qp] * inv_fp_old;
  ce_old = fe_old.transpose() * fe_old; // Ce old = Cauchy tensor at previous time step
  ee_old = 0.5 * (ce_old - iden);
  ee_rate = (ee - ee_old) / _dt;
  invce_ee_rate = invce * ee_rate;

  _heat_rate_vis[_qp]  = _C0 * trD * std::abs(trD) * _density[_qp] * _h_e * _h_e ;
  _heat_rate_vis[_qp] += _C1 * trD * _density[_qp] * _h_e;
  _heat_rate_vis[_qp] *= ee_rate.trace();

  //Heat rate due to thermo-elastic coupling
  //heat source = - G_Gruneisen * rho_0 * C_v * T * Tr(Ce^-1 dot(Ee))
   _heat_rate_therm1[_qp] = - _G_Gruneisen * _density[_qp] * _specific_heat[_qp] * temp
                           * invce_ee_rate.trace() * (1.0 - _c[_qp]) * (1.0 - _c[_qp]);
  //From Psi_cpl
   _heat_rate_therm2[_qp] = thermal_expansion_coeff * temp
              * std::exp((2.0/3.0) * thermal_expansion_coeff * (temp - _reference_temperature))
              * (Kb * std::pow(Je , 2.0/3.0) * invce_ee_rate.trace() - (1.0/3.0) * thermal_coupling_tensor.trace())
              * (1.0 - _c[_qp]) * (1.0 - _c[_qp]);
  
  // Change to positive sign
  _heat_rate_therm1[_qp] *= -1.0;
  _heat_rate_therm2[_qp] *= -1.0;
  // Heat rate due to plasticity
  _heat_rate_p[_qp] = 0.5 * (1.0 - _c[_qp]) * (1.0 - _c[_qp])
                    * ( _W0p_tmp - _W0p_tmp_old ) / _dt; 
  // Lagrangian strain rate trace
  _Tr_E_dot[_qp] = ee_rate.trace();

  resid = _pk2_tmp - pk2_new;
}

// Calculate slip increment,dslipdtau. Override to modify.
void
FiniteStrainCPPFFractureConstH::getSlipIncrements()
{
 // Real Je;
 // Je = _deformation_gradient[_qp].det();
  for (unsigned int i = 0; i < _nss; ++i)
  { 
      //_slip_incr_tol  = 0.025 * (10.0 - std::pow(Je,3.0)/0.1125);
      _slip_incr(i) = _a0(i) * std::pow(std::abs(_tau(i) / _gss_tmp[i]), 1.0 / _xm(i)) *
                      copysign(1.0, _tau(i)) * _dt;
      _dslipdtau(i) = _a0(i) / _xm(i) * copysign(1.0, _tau(i)) *
                      std::pow(std::abs(_tau(i) / _gss_tmp[i]), 1.0 / _xm(i) - 1.0) / _gss_tmp[i] *
                      _dt;

      if (std::abs(_slip_incr(i)) > _slip_incr_tol * _dt )
      {
        _slip_incr(i) = _slip_incr_tol * copysign(1.0, _tau(i))  * _dt ;
        _dslipdtau(i) = 0.0;
       //_err_tol = true;
#ifdef DEBUG
       mooseWarning("Maximum allowable slip increment exceeded ", std::abs(_slip_incr(i)));
#endif 
       //if ( i == (_nss - 1.0)) {
       }
  }
}
