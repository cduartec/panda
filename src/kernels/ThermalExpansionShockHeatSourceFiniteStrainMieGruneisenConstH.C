#include "ThermalExpansionShockHeatSourceFiniteStrainMieGruneisenConstH.h"

registerMooseObject("pandaApp",ThermalExpansionShockHeatSourceFiniteStrainMieGruneisenConstH);

template <>
InputParameters
validParams<ThermalExpansionShockHeatSourceFiniteStrainMieGruneisenConstH>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Thermal expansion heat source kernel"
                             "generic kernel for finite strain"
                             "Luscher2017 model and Amor damage model"
                             "Mie Gruneisen equation of state"
                             "(Menon, 2014) (Zhang, 2011)");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addRequiredParam<Real>("G_Gruneisen", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("s_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  params.addCoupledVar("c","Phase field damage variable: Used to indicate calculation of Off Diagonal Jacobian term");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("h_e", "Element size");
  params.addRequiredParam<Real>("c_l", "Speed of sound");
  params.addRequiredParam<Real>("beta_v", "Amount of shock dissipation converted into heat");
  return params;
}

ThermalExpansionShockHeatSourceFiniteStrainMieGruneisenConstH::ThermalExpansionShockHeatSourceFiniteStrainMieGruneisenConstH(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _kdamage(getParam<Real>("kdamage")),
    _G_Gruneisen(getParam<Real>("G_Gruneisen")), // Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS
    _s_UsUp(getParam<Real>("s_UsUp")), // Us-Up slope in Mie-Gruneisen EOS
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density_name")),
    _deformation_gradient(getMaterialPropertyByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient
    _deformation_gradient_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient, previous timestep
    _fp(getMaterialPropertyByName<RankTwoTensor>(_base_name + "fp")), // Plastic deformation gradient
    _fp_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "fp")), // Plastic deformation gradient, previous timestep
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")), //elasticity tensor
    _c(coupledValue("c")),
    _c_coupled(isCoupled("c")),
    _c_var(_c_coupled ? coupled("c") : 0),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _h_e(getParam<Real>("h_e")),
    _c_l(getParam<Real>("c_l")),
    _beta_v(getParam<Real>("beta_v"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
ThermalExpansionShockHeatSourceFiniteStrainMieGruneisenConstH::computeQpResidual()
{
  Real heat_source, Je;
  Real Kb = _Bulk_Modulus_Ref;
  Real thermal_expansion_coeff; // thermal expansion coefficient depends on Gruneisen parameter, bulk modulus and sound speed
  RankTwoTensor ce, ce_old, invce, iden, ee, ee_old, ee_rate, invce_ee_rate;
  RankTwoTensor fe, fe_old, inv_fp, inv_fp_old;

  // calculate elastic deformation gradient
  inv_fp = _fp[_qp].inverse();
  inv_fp_old = _fp_old[_qp].inverse();
  fe = _deformation_gradient[_qp] * inv_fp;
  fe_old = _deformation_gradient_old[_qp] * inv_fp_old;

  // preliminary calculation of deformation-related tensors
  Je = fe.det(); // Jacobian = relative volume
  ce = fe.transpose() * fe; // Ce = Cauchy tensor
  ce_old = fe_old.transpose() * fe_old; // Ce old = Cauchy tensor at previous time step
  invce = ce.inverse();
  iden.zero();
  iden.addIa(1.0);
  ee = 0.5 * (ce - iden);
  ee_old = 0.5 * (ce_old - iden);
  ee_rate = (ee - ee_old) / _dt;
  invce_ee_rate = invce * ee_rate;

  // Equation of state contribution (Psi_EOS in Luscher2017 and P_eos in Menon 2014, equation 16)
  // heat source = - G_Gruneisen * rho_0 * C_v * T * Tr(Ce^-1 dot(Ee))
  heat_source = - _G_Gruneisen * _density[_qp] * _specific_heat[_qp] * _u[_qp]
              * invce_ee_rate.trace();

  if (Je > 1.0) {
    heat_source = heat_source * (1.0 - _c[_qp]) * (1.0 - _c[_qp]);
  }

  // Adding viscous energy contribution (q*vdot)
  // Calculate bulk viscosity damping
  // C0 * dot(J) / J * |dot(J) / J| + C1 * dot(J) / J
  // C0 should be chosen of the order of rho * Le^2, rho = density, Le = element size
  // C1 should be chosen of the order of rho * Le * cs, cs = sound speed
  // Maheo et al. Mechanics Research Communications 38 (2011) 81 88
  Real J_dot;
  J_dot = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  if (Je < 1.0){
  Real viscous_energy, trD;
  RankTwoTensor pk2_vis, pk2_vis_E_dot;
  // Real trD;
  trD = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  trD /= _deformation_gradient_old[_qp].det();
  pk2_vis.addIa( _C0 * trD * _h_e * _h_e * std::abs(trD) * _density[_qp] );
  pk2_vis.addIa( _C1 * trD * _h_e * _density[_qp] * _c_l);
  pk2_vis_E_dot = ee_rate * pk2_vis;
  viscous_energy  = _beta_v * pk2_vis_E_dot.trace() ;
  heat_source += abs(viscous_energy); 
  }
  // Coupling contribution (Psi_cpl in Luscher2017, equation 15)
  // no difference for this heat source in the case of Birch-Murnaghan or Mie-Gruneisen
  // just the thermal expansion coefficient has to be redefined using the relationship:
  // G_Gruneisen * rho_0 * C_v / K_0 = alpha_thermal
  //
  // heat source = - alpha T / 3 exp(2/3 alpha (T-T0)) (dot(Ee) : C : I)
  //               + K alpha T Je^2/3 Tr(Ce^-1 dot(Ee)) exp(2/3 alpha (T-T0))
  thermal_expansion_coeff = _G_Gruneisen * _density[_qp] * _specific_heat[_qp] / Kb;

  RankTwoTensor thermal_coupling_tensor;
  thermal_coupling_tensor = _elasticity_tensor[_qp] * iden;
  thermal_coupling_tensor = ee_rate * thermal_coupling_tensor;
  heat_source += thermal_expansion_coeff * _u[_qp]
              * std::exp((2.0/3.0) * thermal_expansion_coeff * (_u[_qp] - _reference_temperature))
              * (Kb * std::pow(Je , 2.0/3.0) * invce_ee_rate.trace() - (1.0/3.0) * thermal_coupling_tensor.trace())
              * (1.0 - _c[_qp]) * (1.0 - _c[_qp]);

  return - heat_source * _test[_i][_qp];
}

Real
ThermalExpansionShockHeatSourceFiniteStrainMieGruneisenConstH::computeQpJacobian()
{
  Real Kb = _Bulk_Modulus_Ref;
  RankTwoTensor ce, ce_old, iden, ee, ee_old, ee_rate;
  RankTwoTensor fe, fe_old, inv_fp, inv_fp_old;

  // calculate elastic deformation gradient
  inv_fp = _fp[_qp].inverse();
  inv_fp_old = _fp_old[_qp].inverse();
  fe = _deformation_gradient[_qp] * inv_fp;
  fe_old = _deformation_gradient_old[_qp] * inv_fp_old;

  ce = fe.transpose() * fe; // Ce = Cauchy tensor
  ce_old = fe_old.transpose() * fe_old; // Ce old = Cauchy tensor at previous time step
  iden.zero();
  iden.addIa(1.0);
  ee = 0.5 * (ce - iden);
  ee_old = 0.5 * (ce_old - iden);
  ee_rate = (ee - ee_old) / _dt;

  // approximate (small strain) Jacobian, first order in the strain rate and temperature
  // so that only the EOS term contribution has to be considered
  return _G_Gruneisen * _density[_qp] * _specific_heat[_qp] * ee_rate.trace() * _phi[_j][_qp] * _test[_i][_qp];
}

Real
ThermalExpansionShockHeatSourceFiniteStrainMieGruneisenConstH::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  Real Kb = _Bulk_Modulus_Ref;
  RankTwoTensor ce, ce_old, iden, ee, ee_old, ee_rate;
  RankTwoTensor fe, fe_old, inv_fp, inv_fp_old;

  // calculate elastic deformation gradient
  inv_fp = _fp[_qp].inverse();
  inv_fp_old = _fp_old[_qp].inverse();
  fe = _deformation_gradient[_qp] * inv_fp;
  fe_old = _deformation_gradient_old[_qp] * inv_fp_old;

  ce = fe.transpose() * fe; // Ce = Cauchy tensor
  ce_old = fe_old.transpose() * fe_old; // Ce old = Cauchy tensor at previous time step
  iden.zero();
  iden.addIa(1.0);
  ee = 0.5 * (ce - iden);
  ee_old = 0.5 * (ce_old - iden);
  ee_rate = (ee - ee_old) / _dt;

  val = 0.0;
  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k]) {
      // approximate (small strain) off-diagonal Jacobian, first order in the strain rate and temperature
      val = _G_Gruneisen * _density[_qp] * _specific_heat[_qp] * _u[_qp] * _grad_phi[_j][_qp](k) * _test[_i][_qp] / _dt;
    }
  }

  if (jvar == _c_var) {
    // approximate (small strain) off-diagonal Jacobian, first order in the strain rate and temperature
    // assumption that damage is always present
    val = - 2.0 * (1.0 - _c[_qp]) * _G_Gruneisen * _density[_qp] * _specific_heat[_qp]
        * _u[_qp] * ee_rate.trace() * _phi[_j][_qp] * _test[_i][_qp];
  }

  return val;
}
