/// Calculates heat generated due to thermal expansion

#include "ComputeThermalExpansionHeatSourceSmallStrainMieGruneisen.h"

registerMooseObject("pandaApp",ComputeThermalExpansionHeatSourceSmallStrainMieGruneisen);

InputParameters
ComputeThermalExpansionHeatSourceSmallStrainMieGruneisen::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Thermal expansion heat source kernel generic kernel for small strain"
                             "Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)");
  params.addRequiredParam<Real>("Gamma", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>("density", "density", "Property name of the density material property");
  params.addCoupledVar("c","Phase field damage variable: Used to indicate calculation of Off Diagonal Jacobian term");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addRequiredParam<Real>("Le","Maximum element size");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("sound_speed","Speed of sound in the material");
  params.addRequiredParam<Real>("beta_v","proportion of viscous energy contributing to heat generation");
  return params;
}

ComputeThermalExpansionHeatSourceSmallStrainMieGruneisen::ComputeThermalExpansionHeatSourceSmallStrainMieGruneisen(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")), //elasticity tensor
    _Gamma(getParam<Real>("Gamma")), // Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS
    _ref_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density")),
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _c(coupledValue("c")),
    _temp(coupledValue("temperature")),
    _Le(getParam<Real>("Le")),
    _c_coupled(isCoupled("c")),
    _c_var(_c_coupled ? coupled("c") : 0),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _sound_speed(getParam<Real>("sound_speed")),
    _beta_v(getParam<Real>("beta_v")),
    _shock_heat(declareProperty<Real>(_base_name + "shock_heat"))
{
}

void
ComputeThermalExpansionHeatSourceSmallStrainMieGruneisen::computeQpProperties()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  Real K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);

  // Equation of state contribution (Psi_EOS in Luscher2017 and P_eos in Menon 2014, equation 16)
  // heat source = - G_Gruneisen * rho_0 * C_v * T * Tr(Ce^-1 dot(Ee))
  Real jacob, shock_heat;
  RankTwoTensor ce, invce, ee_rate, invce_ee_rate;
  jacob = 1.0 + _mechanical_strain[_qp].trace();
  ce = 2.0 * _mechanical_strain[_qp] + I2; // Ce = Cauchy tensor
  invce = ce.inverse();
  ee_rate = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  invce_ee_rate = invce * ee_rate;
  shock_heat = - _Gamma * _density[_qp] * _specific_heat[_qp] * _temp[_qp] * invce_ee_rate.trace();
  if (jacob > 1.0) {
    shock_heat = shock_heat * (1.0 - _c[_qp]) * (1.0 - _c[_qp]);
  }
  

  // Artificial viscosity contribution (Maheo et al. Mechanics Research Communications 38 (2011))

  Real trD, viscous_energy;
  trD = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;
  viscous_energy = _C0 * _density[_qp] * trD * std::abs(trD) * std::pow(_Le,2.0) / std::pow(jacob,2.0) + _C1 * _density[_qp] * _sound_speed * trD * _Le / jacob;
  viscous_energy = _beta_v * viscous_energy * invce_ee_rate.trace();
  shock_heat += std::abs(viscous_energy);
  _shock_heat[_qp] = shock_heat;
}

