/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"
#include "MathUtils.h"

registerMooseObject("pandaApp", FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT);

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT>()
{
  InputParameters params = validParams<FiniteStrainCrystalPlasticity>();
  params.addClassDescription("Crystal Plasticity class. Damage. Amor 2009."
                             "Calculation of the broken plastic energy for temperature calculation,"
                             "Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)"
                             "finite strain formalism as in (Luscher 2017)."
                             "Computes the stress and free energy derivatives for the phase field "
                             "fracture model.");
  params.addRequiredCoupledVar("c","Order parameter for damage");
  params.addRequiredCoupledVar("temp","Temperature");
  params.addRequiredCoupledVar("h_max","Maximum element size");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<MaterialPropertyName>(
      "F_name", "E_elpl", "Name of material property storing the elastic and plastic energy driving damage");
  params.addRequiredParam<Real>("G_Gruneisen", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("s_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus = rho_0 c_0^2 (equation 12 in Menon 2014)");
  //params.addRequiredParam<Real>("bulk_modulus_cor", "KT0 prime correction (Menikoff 2001)");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  //params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("plastic_factor","Prefactor of the plastic contribution to damage");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  params.addRequiredCoupledVar("p", "Pressure");
  params.addParam<MaterialPropertyName>("delasticity_tensor_dp", "Derivative of C_ijkl  with pressure");  
  params.addRequiredParam<Real>("temp_melt", "Melting temperature");
  params.addRequiredParam<Real>("temp_init", "Initial temperature softening");
  params.addRequiredParam<Real>("q", "Exponent for thermal softening");
  return params;
}

FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT(const InputParameters & parameters) :
    FiniteStrainCrystalPlasticity(parameters),
    _c(coupledValue("c")),
    _temp(coupledValue("temp")),
    _h_max(coupledValue("h_max")),
    _kdamage(getParam<Real>("kdamage")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _G_Gruneisen(getParam<Real>("G_Gruneisen")),
    _s_UsUp(getParam<Real>("s_UsUp")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    //_Bulk_Modulus_Cor(getParam<Real>("bulk_modulus_cor")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    //_thermal_expansion(getParam<Real>("thermal_expansion")), // volumetric thermal exapnsion coeffcient, as in Austin Barton 2015
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _temp_melt(getParam<Real>("temp_melt")), // Melting temperature
    _temp_init(getParam<Real>("temp_init")), // Initial temperature softening
    _q(getParam<Real>("q")), // Thermal softening exponent
    _plastic_factor(getParam<Real>("plastic_factor")), // prefactor of the plastic contribution to damage
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
    _sigma_eos(declareProperty<RankTwoTensor>("sigma_eos")), // Cauchy stress eos
    _sigma_dev(declareProperty<RankTwoTensor>("sigma_dev")), // Cauchy stress dev
    _slip_incr_out(declareProperty<std::vector<Real>>("slip_incr_out")), // slip system strain increment for output
    _degradation(declareProperty<Real>("degradation")),
    _tau_out(declareProperty<std::vector<Real>>("tau_out")), // slip system strain increment for output
    _gss_thermal(declareProperty<std::vector<Real>>("gss_thermal")),
    _p(coupledValue("p")),
    _p_name(getVar("p", 0)->name()),
    _delasticity_tensor_dp(getMaterialPropertyDerivativeByName<RankFourTensor>("elasticity_tensor",  _p_name)),
    _p_ev(declareProperty<Real>("p_ev")),
    _rot(declareProperty<RankTwoTensor>("rot")) // Rotation tensor considering material rotation and crystal orientation
{
    _err_tol = false;
}

void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::initQpStatefulProperties()
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
  _degradation[_qp] = 1.0; // No thermal degradation
  initSlipSysProps(); // Initializes slip system related properties
  initAdditionalProps();
}


void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::getInitSlipSysRes()
{

  if (_gprops.size() <= 0)
    mooseError("FiniteStrainCrystalPLasticity: Error in reading slip system resistance properties: "
               "Specify input in .i file or in slip_sys_res_prop_file or in slip_sys_file");

  _gss[_qp].resize(_nss, 0.0);

  unsigned int num_data_grp = 3; // Number of data per group e.g. start_slip_sys, end_slip_sys,
                                 // value

  for (unsigned int i = 0; i < _gprops.size() / num_data_grp; ++i)
  {
    Real vs, ve;
    unsigned int is, ie;

    vs = _gprops[i * num_data_grp];
    ve = _gprops[i * num_data_grp + 1];

    if (vs <= 0 || ve <= 0)
      mooseError("FiniteStrainCrystalPLasticity: Indices in gss property read must be positive "
                 "integers: is = ",
                 vs,
                 " ie = ",
                 ve);

    if (vs != floor(vs) || ve != floor(ve))
      mooseError("FiniteStrainCrystalPLasticity: Error in reading slip system resistances: Values "
                 "specifying start and end number of slip system groups should be integer");

    is = static_cast<unsigned int>(vs);
    ie = static_cast<unsigned int>(ve);

    if (is > ie)
      mooseError("FiniteStrainCrystalPLasticity: Start index is = ",
                 is,
                 " should be greater than end index ie = ",
                 ie,
                 " in slip system resistance property read");

    for (unsigned int j = is; j <= ie; ++j)
      _gss[_qp][j - 1] = _gprops[i * num_data_grp + 2];
  }

  for (unsigned int i = 0; i < _nss; ++i)
    if (_gss[_qp][i] <= 0.0)
      mooseError("FiniteStrainCrystalPLasticity: Value of resistance for slip system ",
                 i + 1,
                 " non positive");
}



// Read flow rate parameters from .i file

void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::getFlowRateParams()
{
  if (_flowprops.size() <= 0)
    mooseError("FiniteStrainCrystalPLasticity: Error in reading flow rate  properties: Specify "
               "input in .i file or a slip_sys_flow_prop_file_name");

  _a0.resize(_nss);
  _xm.resize(_nss);

  unsigned int num_data_grp = 2 + _num_slip_sys_flowrate_props; //Number of data per group e.g
                                                                //start_slip_sys, end_slip_sys,
                                                                //value1, value2, ..

  for (unsigned int i = 0; i < _flowprops.size() / num_data_grp; ++i)
  {
    Real vs, ve;
    unsigned int is, ie;

    vs = _flowprops[i * num_data_grp];
    ve = _flowprops[i * num_data_grp + 1];

    if (vs <= 0 || ve <= 0)
      mooseError("FiniteStrainCrystalPLasticity: Indices in flow rate parameter read must be "
                 "positive integers: is = ",
                 vs,
                 " ie = ",
                 ve);

    if (vs != floor(vs) || ve != floor(ve))
      mooseError("FiniteStrainCrystalPLasticity: Error in reading flow props: Values specifying "
                 "start and end number of slip system groups should be integer");

    if (vs != floor(vs) || ve != floor(ve))
      mooseError("FiniteStrainCrystalPLasticity: Error in reading flow props: Values specifying "
                 "start and end number of slip system groups should be integer");

    is = static_cast<unsigned int>(vs);
    ie = static_cast<unsigned int>(ve);

    if (is > ie)
      mooseError("FiniteStrainCrystalPLasticity: Start index is = ",
                 is,
                 " should be greater than end index ie = ",
                 ie,
                 " in flow rate parameter read");

    for (unsigned int j = is; j <= ie; ++j)
    {
      _a0(j - 1) = _flowprops[i * num_data_grp + 2];
      _xm(j - 1) = _flowprops[i * num_data_grp + 3];
    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
  {
    if (!(_a0(i) > 0.0 && _xm(i) > 0.0))
    {
      mooseWarning(
          "FiniteStrainCrystalPlasticity: Non-positive flow rate parameters ", _a0(i), ",", _xm(i));
      break;
    }
  }
}

// Read hardness parameters from .i file
void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::getHardnessParams()
{
  if (_hprops.size() <= 0)
    mooseError("FiniteStrainCrystalPLasticity: Error in reading hardness properties: Specify input "
               "in .i file or a slip_sys_hard_prop_file_name");

  _r = _hprops[0];
  _h0 = _hprops[1];
  _tau_init = _hprops[2];
  _tau_sat = _hprops[3];
}

void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::preSolveQp()
{
  //Initialize variable
  if (_first_substep)
  {
    _Jacobian_mult[_qp].zero(); // Initializes jacobian for preconditioner
    calc_schmid_tensor();
  }

  if (_max_substep_iter == 1)
    _dfgrd_tmp = _deformation_gradient[_qp]; // Without substepping
  else
    _dfgrd_tmp = _dfgrd_scale_factor * _delta_dfgrd + _dfgrd_tmp_old;
  
  _err_tol = false;

  //Computing themal softening coefficient
  Real current_temp = _temp[_qp];
  _degradation[_qp] = 1.0;
  if ( current_temp <= _temp_melt & current_temp >= _temp_init)
  {
    _degradation[_qp] = 1.0 - std::pow( (current_temp  - _temp_init) / (_temp_melt - _temp_init), _q);
  }
  else if ( current_temp  > _temp_melt)
  {
    _degradation[_qp] = 0.2;
  }

}


void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::solveQp()
{

  preSolveStatevar();
  solveStatevar();
  if (_err_tol)
     return;
  postSolveStatevar();
}


void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::preSolveStatevar()
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
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::solveStatevar()
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
    mooseWarning("FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT: Hardness Integration error gmax", gmax, "\n");
#endif
    _err_tol = true;
  }
}

void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::postSolveStatevar()
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

 // Update slip system resistance. Overide to incorporate new slip system resistance laws
void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::update_slip_system_resistance()
{
  updateGss();
}


/**
 * Old function to update slip system resistances.
 * Kept to avoid code break at computeQpstress
 * output slip increment
 */
void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::updateGss()
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


  // Real val = std::cosh(_h0 * _accslip_tmp / (_tau_sat - _tau_init)); // Karthik
  // val = _h0 * std::pow(1.0/val,2.0); // Kalidindi

  for (unsigned int i = 0; i < _nss; ++i)
    // hb(i)=val;
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
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::update_energies()
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

  // Total free energy density (- sign in front of W0e_neg? W0e_neg is positive in compression)
  _F[_qp] = (hist_variable + _plastic_factor * _W0p_tmp) * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) -
            _W0e_neg[_qp] + _gc[_qp] / (2 * _l[_qp]) * c * c;

  // derivative of total free energy density wrt c
  _dFdc[_qp] = -(hist_variable + _plastic_factor * _W0p_tmp) * 2.0 * (1.0 - c) * (1 - _kdamage) +
               _gc[_qp] / _l[_qp] * c;

  // 2nd derivative of total free energy density wrt c
  _d2Fdc2[_qp] = (hist_variable + _plastic_factor * _W0p_tmp) * 2.0 * (1 - _kdamage) + _gc[_qp] / _l[_qp];

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history variable
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = -cauchy_stress_undamaged * (1.0 - c) * (1 - _kdamage);

  //_dW0p_broken_dstrain[_qp].zero();
  //_dW0p_dstrain[_qp].zero();
  //_dstress_dc[_qp] = -cauchy_stress_undamaged * (2.0 * (1.0 - c)); defined in calcResidual
}

void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::calcResidual( RankTwoTensor &resid )
{
  RankTwoTensor iden, ce, invce, ee, ce_pk2, eqv_slip_incr, pk2_new, temporal;
  Real trD, Kb, Je, J_dot, detFe;
  Real c = _c[_qp];
  Real temp = _temp[_qp];
  Real h_max = _h_max[_qp];
  Real xfac = ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage);
  Real thermal_expansion_coeff; // thermal expansion coefficient depends on Gruneisen parameter, bulk modulus and sound speed
  Real Eev, dp_dEev, pk2_ev;

  detFe = _fe.det();
  // Anisotropic elasticity (Luscher2017)
  // Kb = K in Luscher2017
  // Kb = (1/9) I : C : I
  //Real Kb = 0.0;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      Kb +=  _elasticity_tensor[_qp](i, i, j, j);

  Kb = (1.0 / 9.0) * Kb;

  // or reference bulk modulus can be an input
  //Kb = _Bulk_Modulus_Ref;

  // thermal expansion coefficient depends on Gruneisen parameter, bulk modulus and sound speed
  // condition: for small strain and temperature variation the second Piola-Kirchoff stress has to become:
  // (Luscher 2017): C : (Ee - alpha)
  // therefore the relationship holds:
  // G_Gruneisen * rho_0 * C_v / K_0 = alpha_thermal
  thermal_expansion_coeff = _G_Gruneisen * _density[_qp] * _specific_heat[_qp] / Kb;

  iden.zero();
  temporal.zero();
  iden.addIa(1.0);

  _fe = _dfgrd_tmp * _fp_prev_inv; // _fp_inv  ==> _fp_prev_inv

  ce = _fe.transpose() * _fe;
  ce_pk2 = ce * _pk2_tmp;
  ce_pk2 = ce_pk2 / _fe.det();

  // Calculate Schmid tensor and resolved shear stresses

  _tau_out[_qp].resize(_nss);

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
  
  //Obtain rotation tensor  
  _rot[_qp] = get_current_rotation(_deformation_gradient[_qp]);

  ce = _fe.transpose() * _fe;
  invce = ce.inverse();
  ee = 0.5 * ( ce - iden );
  Je = _fe.det(); // Jacobian = relative volume
  J_dot = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  Real Peos, Peos_pos, Peos_neg, V0V, eta;

  // Cauchy pressure, Mie Gruneisen (Menon, 2014), (Zhang, 2011)
  V0V = 1.0 / Je; // relative volume
  eta = 1.0 - Je; // eta = 1 - (v / v0) (Menon, 2014)
  Real delta;
  // Pcor = correcting pressure = linearized form of the EOS
  // equation (18) in Luscher2017
  delta = 1.5 * (std::pow(Je , 2.0/3.0) - 1.0);
  RankTwoTensor thermal_eigenstrain, thermal_eigenstrain_dev;

  // thermal eigenstrain (equation (18) in Luscher2017)
  // Lagrangian strain E_thermal = 1/2 (F_thermal^T F_thermal - I)
  // finite strain formula (Lubarda2002): F_thermal = exp((alpha/3)*(T-T_ref)) I
  thermal_eigenstrain = (1.0 / 2.0)
                      * (std::exp((2.0/3.0) * thermal_expansion_coeff * (temp - _reference_temperature)) - 1.0)
                      * iden;
  thermal_eigenstrain_dev = thermal_eigenstrain- thermal_eigenstrain.trace()*iden;

//  if (Je < 1.0){


  Peos = _G_Gruneisen * _density[_qp] * _specific_heat[_qp] * (temp - _reference_temperature) * V0V;
  Peos += Kb * eta * (1.0 - (_G_Gruneisen / 2.0) * (V0V - 1.0)) / std::pow((1.0 - _s_UsUp * eta) , 2.0);
  Peos = - Peos; // negative stress in compression

  Peos_pos = (std::abs(Peos) + Peos) / 2.0; // volumetric expansion
  Peos_neg = (std::abs(Peos) - Peos) / 2.0; // volumetric compression

  // positive and negative volumetric stresses (equation 17 in Luscher2017) + damage
  pk2_new = Je * (Peos_pos * xfac - Peos_neg) * invce;

  // Undamaged second piola-kirchoff stress to calculate undamaged plastic work
  _pk2_undamaged[_qp] = Je * Peos * invce;

  _sigma_eos[_qp] =  Peos * iden;

  
  // deviatoric stress + damage (equation (18) in Luscher2017): C : (Ee - alpha)
  pk2_new += xfac * _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);

  //Undamaged second piola-kirchoff stress to calculate undamaged plastic work
  _pk2_undamaged[_qp] += _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);

  // deviatoric stress + damage (equation (18) in Luscher2017): C : (Ee - alpha)
  pk2_new += xfac * _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);

  // Undamaged second piola-kirchoff stress to calculate undamaged plastic work
  _pk2_undamaged[_qp] += _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);

  //temporal  = _elasticity_tensor[_qp] * (ee - thermal_eigenstrain) ;

  pk2_new += xfac * (-1.0/3.0 * temporal.trace() * iden + _elasticity_tensor[_qp] * (ee - thermal_eigenstrain));
         
  pk2_new +=  -1.0 * xfac * Kb * std::pow(Je , 2.0/3.0)
              * (delta  -  thermal_eigenstrain.trace())
              * invce;

  // Undamaged second piola-kirchoff stress to calculate undamaged plastic work
  //_pk2_undamaged[_qp] += -1.0/3.0 * temporal.trace() * iden + _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);

  _pk2_undamaged[_qp] += -1.0 * Kb * std::pow(Je , 2.0/3.0)
                         * (delta  - thermal_eigenstrain.trace())
                         * invce;

  _sigma_dev[_qp] = -1.0/3.0 * temporal.trace() * iden + _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);
  
  _sigma_dev[_qp] = _fe * _sigma_dev[_qp] * _fe.transpose()/detFe;
  //} 
  //else {

  // For positive deformation
  
  //pk2_new = xfac * _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);
  //_pk2_undamaged[_qp] = _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);
  //}
  // volumetric free energy = Psi_EOS in Luscher2017
  // still to implement for Mie Gruneisen
  _W0e_pos[_qp] = 0.0;
  _W0e_neg[_qp] = 0.0;
  if (Je > 1.0) { // only in expansion
    _W0e_pos[_qp] = 1.0*(_G_Gruneisen * _density[_qp] * _specific_heat[_qp] * (_reference_temperature-temp) * std::log(1/V0V)
                  + (Kb*V0V*(2*_s_UsUp-2-_G_Gruneisen))/(2*(_s_UsUp-1)*std::pow(_s_UsUp, 2.0)*(V0V+_s_UsUp*(1-V0V)))
                  - ((2*_s_UsUp-2-_G_Gruneisen)*Kb)/(2*(_s_UsUp-1)*std::pow(_s_UsUp, 2.0))
                  + (_G_Gruneisen*(1-2*_s_UsUp)+2*std::pow((_s_UsUp-1), 2.0))/(2*std::pow(_s_UsUp, 2.0)*std::pow((_s_UsUp-1), 2.0))
                  * std::log((V0V+_s_UsUp*(1-V0V))/(V0V))*Kb
                  + (Kb*_G_Gruneisen)/(2*std::pow((_s_UsUp-1), 2.0))*std::log(1/V0V));
  } else { 
     _W0e_pos[_qp] = 0.0;
  }

  // volumetric coupling free energy = Psi_cpl in Luscher2017
  RankTwoTensor elastic_energy_tensor, thermal_coupling_tensor;

  elastic_energy_tensor = _elasticity_tensor[_qp] * ee;
  elastic_energy_tensor = 0.5 * ee * elastic_energy_tensor;
  _W0e_pos[_qp] += elastic_energy_tensor.trace(); // 1/2 * Ee : C : Ee

  thermal_coupling_tensor = _elasticity_tensor[_qp] * thermal_eigenstrain; // C : alpha in equation 15 of Luscher2017
  thermal_coupling_tensor = ee * thermal_coupling_tensor;
  _W0e_pos[_qp] -= thermal_coupling_tensor.trace(); // - Ee : C : alpha in equation 15 of Luscher2017

  _W0e_pos[_qp] -= 0.5 * Kb * delta * delta; // - 1/2 K delta^2 in equation 15 of Luscher2017

  _W0e_pos[_qp] += Kb * delta * thermal_eigenstrain.trace(); // + K delta alpha_v in equation 15 of Luscher2017

  // Assign history variable and derivative
  if (_W0e_pos[_qp] > _hist_old[_qp])
    _hist[_qp] = _W0e_pos[_qp];
  else
    _hist[_qp] = _hist_old[_qp];

  // Used in PFFracBulkRate Jacobian
  // approximation: it should be only the positive part, undamaged pk2 is all
  //_dW0e_dstrain[_qp] = _pk2_undamaged[_qp];

  // Used in StressDivergencePFFracTensors Jacobian
  // same approximation as above
  if (c < 1.0)
    _dstress_dc[_qp] = -_pk2_undamaged[_qp] * 2.0 * (1.0 - c) * (1 - _kdamage);
  else
    _dstress_dc[_qp].zero();

  // Calculate bulk viscosity damping
  // C0 * dot(J) / J * |dot(J) / J| + C1 * dot(J) / J
  // C0 should be chosen of the order of rho * Le^2, rho = density, Le = element size
  // C1 should be chosen of the order of rho * Le * cs, cs = sound speed
  // Maheo et al. Mechanics Research Communications 38 (2011) 81 88
  trD = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  trD /= _deformation_gradient_old[_qp].det();

//  if (J_dot < 0.0){
  pk2_new.addIa( _C0 * trD * h_max * h_max * std::abs(trD) * _density[_qp] );
  pk2_new.addIa( _C1 * trD * h_max * _density[_qp]);
//  }

  resid = _pk2_tmp - pk2_new;
}

// Calculate slip increment,dslipdtau. Override to modify.
void
FiniteStrainCrystalPlasticityPFFractureStressMieGruneisenNewT::getSlipIncrements()
{
 // Real Je;
 // Je = _deformation_gradient[_qp].det();
  _gss_thermal[_qp].resize(_nss);

  for (unsigned i = 0; i < _nss; ++i)
  {
     _gss_thermal[_qp][i] =  _gss_tmp[i] * _degradation[_qp];
  }

  for (unsigned int i = 0; i < _nss; ++i)
  { 
      //_slip_incr_tol  = 0.025 * (10.0 - std::pow(Je,3.0)/0.1125);
      _slip_incr(i) = _a0(i) * std::pow(std::abs(_tau(i) / _gss_thermal[_qp][i]), 1.0 / _xm(i)) *
                      copysign(1.0, _tau(i)) * _dt;
      _dslipdtau(i) = _a0(i) / _xm(i) * copysign(1.0, _tau(i)) *
                      std::pow(std::abs(_tau(i) / _gss_thermal[_qp][i]), 1.0 / _xm(i) - 1.0) / _gss_thermal[_qp][i] *
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
