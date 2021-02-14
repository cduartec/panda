//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeFiniteStrainElasticStressBirchMurnaghanNew.h"

registerMooseObject("pandaApp", ComputeFiniteStrainElasticStressBirchMurnaghanNew);

template <>
InputParameters
validParams<ComputeFiniteStrainElasticStressBirchMurnaghanNew>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Compute stress using elasticity for finite strains"
                             "third-order Birch Murnaghan equation of state"
                             "finite strain formalism as in Luscher2017"
                             "no free energy calculation for efficiency");
  params.addRequiredCoupledVar("temp","Temperature");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("bulk_modulus_cor", "KT0 prime correction (Menikoff 2001, Yoo 1998)");
  params.addRequiredParam<Real>("G_Gruneisen", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("s_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("L0", "Mesh spacing");
  params.addRequiredParam<Real>("C_s", "speed of sound");
  params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("reference_density", "reference density for MG");
  params.addParam<MaterialPropertyName>("logarithmic_strain", "reference density for MG");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>(
      "heat_source", "Q", "Property name of the chemical heat source property");
  return params;
}

ComputeFiniteStrainElasticStressBirchMurnaghanNew::ComputeFiniteStrainElasticStressBirchMurnaghanNew(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
    GuaranteeConsumer(this),
    _temp(coupledValue("temp")),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _Bulk_Modulus_Cor(getParam<Real>("bulk_modulus_cor")),
    _G_Gruneisen(getParam<Real>("G_Gruneisen")),
    _s_UsUp(getParam<Real>("s_UsUp")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _L0(getParam<Real>("L0")),
    _C_s(getParam<Real>("C_s")),
    _thermal_expansion(getParam<Real>("thermal_expansion")), // volumetric thermal exapnsion coeffcient, as in Austin Barton 2015
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _reference_density(getParam<Real>("reference_density")),
    _logarithmic_strain(getMaterialProperty<Real>("logarithmic_strain")),
    _Jt_minimum(declareProperty<Real>("Jt_minimum")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _Q(getMaterialProperty<Real>("heat_source")),
    _pk2(declareProperty<RankTwoTensor>("pk2")), // 2nd Piola Kirchoff Stress
    _lag_e(declareProperty<RankTwoTensor>("lage")), // Lagrangian strain
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient"))
{
}

void
ComputeFiniteStrainElasticStressBirchMurnaghanNew::initialSetup()
{
}

void
ComputeFiniteStrainElasticStressBirchMurnaghanNew::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  _Jt_minimum[_qp] = 1.0; // history variable for the minimum compression ration ever
}

void
ComputeFiniteStrainElasticStressBirchMurnaghanNew::computeQpStress()
{
  Real trD, Kb, KT0prime, Je, Peos, V0V, delta;
  Real temp = _temp[_qp];
  RankTwoTensor iden, ce, invce, thermal_eigenstrain;

  // reference bulk modulus is an input
  Kb = _Bulk_Modulus_Ref;
  KT0prime = _Bulk_Modulus_Cor;

  // tensor calculation
  iden.zero();
  iden.addIa(1.0);

  ce = _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp]; // Cauchy-Green strain tensor

  _lag_e[_qp] = ce - iden; // Lagrangian strain tensor
  _lag_e[_qp] = _lag_e[_qp] * 0.5;

  invce = ce.inverse();

  Je =  _deformation_gradient[_qp].det(); // Jacobian = relative volume // 1.0 - _logarithmic_strain[_qp]; //

  if ( _deformation_gradient[_qp].det() < _Jt_minimum[_qp])
  _Jt_minimum[_qp] = _deformation_gradient[_qp].det();

  // Cauchy pressure, third-order Birch Murnaghan (Menikoff2001, Yoo1998)
  V0V = 1.0 / Je; // relative volume
  Peos = 1.5 * Kb * (std::pow(V0V , 7.0/3.0) - std::pow(V0V , 5.0/3.0));
  Peos = Peos * (1.0 + 0.75 * (KT0prime - 4.0) * (std::pow(V0V , 2.0/3.0) - 1.0));
  Peos = - Peos; // negative stress in compression
  // Peos += (Kb / _n_Murnaghan) * (std::exp(-_n_Murnaghan * _thermal_expansion * (temp - _reference_temperature)) - 1.0);

  Real eta_minimum = 1.0 - _Jt_minimum[_qp];
  Real delta_e_shock = 0.0;
  // if (eta_minimum < 1.0) delta_e_shock = std::pow(eta_minimum * _C_s,2.0) / std::pow((1.0 - _s_UsUp * eta_minimum) , 2.0) / 2.0;
  Peos -= _G_Gruneisen * _reference_density * (delta_e_shock + _Q[_qp]);

  // volumetric stress (equation 17 in Luscher2017)
  _pk2[_qp] = Je * Peos * invce;

  // thermal eigenstrain (equation (18) in Luscher2017)
  // Lagrangian strain E_thermal = 1/2 (F_thermal^T F_thermal - I)
  // finite strain formula (Lubarda2002): F_thermal = exp((alpha/3)*(T-T_ref)) I
  /*thermal_eigenstrain = (1.0 / 2.0)
                      * (std::exp((2.0/3.0) * _thermal_expansion * (temp - _reference_temperature)) - 1.0)
                      * iden;*/

  // deviatoric/isochoric stress (equation (18) in Luscher2017): C : (Ee - alpha)
  _pk2[_qp] += _elasticity_tensor[_qp] * (_lag_e[_qp]); // - thermal_eigenstrain);

  // Pcor = correcting pressure = linearized form of the EOS
  // equation (18) in Luscher2017
  delta = 1.5 * (std::pow(Je , 2.0/3.0) - 1.0);
  /*_pk2[_qp] -= Kb * std::pow(Je , 2.0/3.0)
           * (delta * iden - 3.0 * thermal_eigenstrain)
           * invce;*/

  _pk2[_qp] -= Kb * std::pow(Je , 2.0/3.0)
           * (delta * iden) * invce;
           

  // Calculate bulk viscosity damping
  // C0 * dot(J) / J * |dot(J) / J| + C1 * dot(J) / J
  // C0 should be chosen of the order of rho * Le^2, rho = density, Le = element size
  // C1 should be chosen of the order of rho * Le * cs, cs = sound speed
  // Maheo et al. Mechanics Research Communications 38 (2011) 81 88
  trD = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  trD /= _deformation_gradient_old[_qp].det();

  // _pk2[_qp].addIa( _C0 * trD * std::abs(trD) );
  // _pk2[_qp].addIa( _C1 * trD );

  if ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() < 0.0) {
  _pk2[_qp].addIa( _C0 * trD * _L0 * _L0 * std::abs(trD) * _reference_density ); // ( _C0 * trD * h_max * h_max * std::abs(trD) * _reference_density )
  _pk2[_qp].addIa( _C1 * trD * _L0 * _C_s * _reference_density); // ( _C1 * trD * h_max * _reference_density)
  }
  // Cauchy stress
  _stress[_qp] = _deformation_gradient[_qp] * _pk2[_qp] * _deformation_gradient[_qp].transpose() / Je;

  // Compute dstress_dstrain
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp]; // This is NOT the exact jacobian
}

void
ComputeFiniteStrainElasticStressBirchMurnaghanNew::rotateQpInitialStress()
{
}
