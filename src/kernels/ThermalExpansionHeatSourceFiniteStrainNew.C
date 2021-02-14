#include "ThermalExpansionHeatSourceFiniteStrainNew.h"
registerMooseObject("pandaApp", ThermalExpansionHeatSourceFiniteStrainNew);
template <>
InputParameters
validParams<ThermalExpansionHeatSourceFiniteStrainNew>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Thermal expansion heat source kernel"
                             "generic kernel for finite strain"
                             "Luscher2017 model and Amor damage model");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("bulk_modulus_cor", "KT0 prime correction (Menikoff 2001, Yoo 1998)");
  params.addRequiredParam<Real>("G_Gruneisen", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("s_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("L0", "Mesh spacing");
  params.addRequiredParam<Real>("C_s", "speed of sound");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("reference_density", "reference density for MG");
  params.addParam<MaterialPropertyName>("logarithmic_strain", "reference density for MG");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  /*params.addParam<MaterialPropertyName>(
      "p", "pressure", "Property name of the pressue property");*/
  params.addRequiredCoupledVar("p", "Pressure");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  return params;
}

ThermalExpansionHeatSourceFiniteStrainNew::ThermalExpansionHeatSourceFiniteStrainNew(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _Bulk_Modulus_Cor(getParam<Real>("bulk_modulus_cor")),
    _G_Gruneisen(getParam<Real>("G_Gruneisen")),
    _s_UsUp(getParam<Real>("s_UsUp")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _L0(getParam<Real>("L0")),
    _C_s(getParam<Real>("C_s")),
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _reference_density(getParam<Real>("reference_density")),
    _logarithmic_strain(getMaterialProperty<Real>("logarithmic_strain")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    // _p(getMaterialProperty<Real>("p")),
    _p(coupledValue("p")),
    _p_name(getVar("p", 0)->name()),
    _deformation_gradient(getMaterialPropertyByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient
    _deformation_gradient_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient, previous timestep
    // _fp(getMaterialPropertyByName<RankTwoTensor>(_base_name + "fp")), // Plastic deformation gradient
    // _fp_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "fp")), // Plastic deformation gradient, previous timestep
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")), //elasticity tensor
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
ThermalExpansionHeatSourceFiniteStrainNew::computeQpResidual()
{
  Real heat_source, Je;
  Real Kb = _Bulk_Modulus_Ref;
  RankTwoTensor ce, ce_old, invce, iden, ee, ee_old, ee_rate, invce_ee_rate;
  RankTwoTensor fe, fe_old, inv_fp, inv_fp_old;

  // calculate elastic deformation gradient
  // inv_fp = _fp[_qp].inverse();
  // inv_fp_old = _fp_old[_qp].inverse();
  fe = _deformation_gradient[_qp]; // * inv_fp;
  fe_old = _deformation_gradient_old[_qp]; // * inv_fp_old;

  // preliminary calculation of deformation-related tensors
  Je = _logarithmic_strain[_qp]; //fe.det(); // Jacobian = relative volume
  ce = fe.transpose() * fe; // Ce = Cauchy tensor
  ce_old = fe_old.transpose() * fe_old; // Ce old = Cauchy tensor at previous time step
  invce = ce.inverse();
  iden.zero();
  iden.addIa(1.0);
  ee = 0.5 * (ce - iden);
  ee_old = 0.5 * (ce_old - iden);
  ee_rate = (ee - ee_old) / _dt;
  invce_ee_rate = invce * ee_rate;

  // Equation of state contribution (Psi_EOS in Luscher2017, equation 16)
  // heat source = - K alpha T exp(-n alpha (T - T0)) Je Tr(Ce^-1 dot(Ee))
  /*heat_source = - Kb * _thermal_expansion * _u[_qp]
              * std::exp(-_n_Murnaghan * _thermal_expansion * (_u[_qp] - _reference_temperature))
              * Je * invce_ee_rate.trace();*/
   heat_source = 0.0; //-_G_Gruneisen * _u[_qp] * ee_rate.trace();

   if (_deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() > 0.0) {
   // if (_deformation_gradient[_qp].det() > 0.0) {
   /*if (_deformation_gradient_old[_qp].det() > 0.0)
   heat_source -= _p[_qp] * ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) /_dt;
   else if (_deformation_gradient_old[_qp].det() < 0.0)
   heat_source -= _p[_qp] * ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) /_dt;*/
   Real V0V = _deformation_gradient[_qp].det();
   Real V0V_old = _deformation_gradient_old[_qp].det();
   Real KT0prime = _Bulk_Modulus_Cor;
   heat_source += (9.0 / 16.0 / _reference_density * Kb * (std::pow ((std::pow(V0V , 2.0/3.0) - 1.0) , 3.0) * KT0prime
                 +  (std::pow ((std::pow(V0V , 2.0/3.0) - 1.0) , 2.0) * (6.0 - 4.0 * std::pow(V0V , 2.0/3.0)))) /_dt);
   heat_source -= (9.0 / 16.0 / _reference_density * Kb * (std::pow ((std::pow(V0V_old , 2.0/3.0) - 1.0),3.0) * KT0prime
                 +  (std::pow ((std::pow(V0V_old , 2.0/3.0) - 1.0) , 2.0) * (6.0 - 4.0 * std::pow(V0V_old , 2.0/3.0)))) /_dt);
     }

   if (_deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() < 0.0) {//   {if (_deformation_gradient[_qp].det() < 0.0)
  /*heat_source += - _G_Gruneisen * _density[_qp] * _specific_heat[_qp] * _u[_qp]
              * invce_ee_rate.trace();*/
  Real viscous_energy, delta_v, P_Hugon_new, P_Hugon_old, delta_e_new, delta_e_old, chi_old, chi_new; // trD;
  /*trD = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  trD /= _deformation_gradient_old[_qp].det();
  viscous_energy  = _C0 * trD * std::abs(trD) * _density[_qp] * h_max * h_max ;
  viscous_energy += _C1 * trD * _density[_qp] *  h_max;
  viscous_energy *= invce_ee_rate.trace();*/ 
  delta_v = abs( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _reference_density;
  chi_old = 1.0 - _deformation_gradient_old[_qp].det();
  P_Hugon_old = _Bulk_Modulus_Ref * chi_old / (1.0 - _s_UsUp * chi_old) / (1.0 - _s_UsUp * chi_old); 
  chi_new = 1.0 - _deformation_gradient[_qp].det();
  P_Hugon_new = _Bulk_Modulus_Ref * chi_new / (1.0 - _s_UsUp * chi_new) / (1.0 - _s_UsUp * chi_new); 
  delta_e_new = abs (P_Hugon_new * chi_new / _reference_density) / 2.0; // - P_Hugon_old;
  delta_e_old = abs (P_Hugon_old * chi_old / _reference_density) / 2.0; 
  // if (_deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() < 0.0)
  heat_source += (_reference_density * abs(delta_e_new - delta_e_old) )/ _dt ; // - (_Bulk_Modulus_Ref * thermal_expansion_coeff * (_u[_qp] - _reference_temperature))
  /*else if (_deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() > 0.0)
  heat_source -= (_reference_density * abs(delta_e_new - delta_e_old) )/ _dt ;*/
  }

  /*if (Je > 1.0) {
    heat_source = heat_source * (1.0 - _c[_qp]) * (1.0 - _c[_qp]);
  }*/

  // Coupling contribution (Psi_cpl in Luscher2017, equation 15)
  // heat source = - alpha T / 3 exp(2/3 alpha (T-T0)) (dot(Ee) : C : I)
  //               + K alpha T Je^2/3 Tr(Ce^-1 dot(Ee)) exp(2/3 alpha (T-T0))
  /*RankTwoTensor thermal_coupling_tensor;
  thermal_coupling_tensor = _elasticity_tensor[_qp] * iden;
  thermal_coupling_tensor = ee_rate * thermal_coupling_tensor;
  heat_source += _thermal_expansion * _u[_qp]
              * std::exp((2.0/3.0) * _thermal_expansion * (_u[_qp] - _reference_temperature))
              * (Kb * std::pow(Je , 2.0/3.0) * invce_ee_rate.trace() - (1.0/3.0) * thermal_coupling_tensor.trace());*/

  return - heat_source * _test[_i][_qp];
}

Real
ThermalExpansionHeatSourceFiniteStrainNew::computeQpJacobian()
{
  Real Kb = _Bulk_Modulus_Ref;
  Real thermal_expansion = _G_Gruneisen * _reference_density * _specific_heat[_qp] / _Bulk_Modulus_Ref;
  RankTwoTensor ce, ce_old, iden, ee, ee_old, ee_rate;
  RankTwoTensor fe, fe_old, inv_fp, inv_fp_old;

  // calculate elastic deformation gradient
  // inv_fp = _fp[_qp].inverse();
  // inv_fp_old = _fp_old[_qp].inverse();
  fe = _deformation_gradient[_qp]; // * inv_fp;
  fe_old = _deformation_gradient_old[_qp]; // * inv_fp_old;

  ce = fe.transpose() * fe; // Ce = Cauchy tensor
  ce_old = fe_old.transpose() * fe_old; // Ce old = Cauchy tensor at previous time step
  iden.zero();
  iden.addIa(1.0);
  ee = 0.5 * (ce - iden);
  ee_old = 0.5 * (ce_old - iden);
  ee_rate = (ee - ee_old) / _dt;

  // approximate (small strain) Jacobian, first order in the strain rate and temperature
  return (1.0/3.0) * Kb * thermal_expansion * ee_rate.trace() * _phi[_j][_qp] * _test[_i][_qp];
}

Real
ThermalExpansionHeatSourceFiniteStrainNew::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  Real Kb = _Bulk_Modulus_Ref;
  Real thermal_expansion = _G_Gruneisen * _reference_density * _specific_heat[_qp] / _Bulk_Modulus_Ref;
  RankTwoTensor ce, ce_old, iden, ee, ee_old, ee_rate;
  RankTwoTensor fe, fe_old, inv_fp, inv_fp_old;

  // calculate elastic deformation gradient
  // inv_fp = _fp[_qp].inverse();
  // inv_fp_old = _fp_old[_qp].inverse();
  fe = _deformation_gradient[_qp]; // * inv_fp;
  fe_old = _deformation_gradient_old[_qp]; // * inv_fp_old;

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
      val = (1.0/3.0) * Kb * thermal_expansion * _u[_qp] * _grad_phi[_j][_qp](k) * _test[_i][_qp] / _dt;
    }
  }


  return val;
}
