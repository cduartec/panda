/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State

#include "ComputeMieGruneisenPFFractureStress.h"

registerMooseObject("pandaApp", ComputeMieGruneisenPFFractureStress);

InputParameters
ComputeMieGruneisenPFFractureStress::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addRequiredCoupledVar("c", "Name of damage variable");
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field fracture model, with small strain"
                             "Considers Mie Gruneisen EOS and artificial viscosity damping");
  params.addParam<bool>("use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<bool>("use_snes_vi_solver",false,"Use PETSc's SNES variational inequalities solver to enforce damage "
                        "irreversibility condition and restrict damage value <= 1.");
  params.addParam<MaterialPropertyName>("barrier_energy", "Name of material property for fracture energy barrier.");
  params.addParam<MaterialPropertyName>("E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>("D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<MaterialPropertyName>("I_name", "indicator", "Name of material property for damage indicator function.");
  params.addParam<MaterialPropertyName>("F_name", "local_fracture_energy", "Name of material property for local fracture energy function.");
  params.addRequiredParam<Real>("Gamma", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addParam<MaterialPropertyName>("density", "density", "Name of Material Property that provides the density");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Name of Material Property that provides the density");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("slope_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("Le","Maximum element size");
  params.addRequiredParam<Real>("sound_speed","Speed of sound in the material");
  return params;
}

ComputeMieGruneisenPFFractureStress::ComputeMieGruneisenPFFractureStress(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _c(coupledValue("c")),
    _l(getMaterialProperty<Real>("l")),
    _pressure(getDefaultMaterialProperty<Real>("fracture_pressure")),
    _gc(getMaterialProperty<Real>("gc_prop")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _use_snes_vi_solver(getParam<bool>("use_snes_vi_solver")),
    _H(declareProperty<Real>("hist")),
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _barrier(getDefaultMaterialProperty<Real>("barrier_energy")),
    _E(declareProperty<Real>(getParam<MaterialPropertyName>("E_name"))),
    _dEdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name())),
    _d2Ed2c(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    _D(getMaterialProperty<Real>("D_name")),
    _dDdc(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name())),
    _d2Dd2c(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _I(getDefaultMaterialProperty<Real>("I_name")),
    _dIdc(getMaterialPropertyDerivative<Real>("I_name", getVar("c", 0)->name())),
    _d2Id2c(getMaterialPropertyDerivative<Real>("I_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _Gamma(getParam<Real>("Gamma")),
    _density(getMaterialProperty<Real>("density")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _temperature(coupledValue("temperature")),
    _ref_temperature(getParam<Real>("reference_temperature")),
    _s(getParam<Real>("slope_UsUp")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("Le")),
    _sound_speed(getParam<Real>("sound_speed")),
    _bulk_modulus(declareProperty<Real>("bulk_modulus")),
    _pressure_eos(declareProperty<Real>("pressure_eos")),
    _stress_eos_elastic(declareProperty<RankTwoTensor>(_base_name + "stress_eos_elastic")),
    _stress_cpl_elastic(declareProperty<RankTwoTensor>(_base_name + "stress_cpl_elastic"))
{
}

void
ComputeMieGruneisenPFFractureStress::initQpStatefulProperties()
{
  _H[_qp] = 0.0;
  const double B=1E10;
  const double PI=3.141592653589793238463;
//  if ( std::abs(_q_point[_qp](1)-100E-3) < 40E-3 ) {
//    if ( std::abs(_q_point[_qp](0)-125E-3) < (0.1*_l[_qp]) ) {
//      _H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](0)-125E-3)/(0.5*_l[_qp]) ); } }
//  if ( std::sqrt( std::pow(_q_point[_qp](0)-0.125,2.0) + std::pow(_q_point[_qp](1)-0.080,2.0)) < 0.040 ) {
//    if ( std::sqrt( std::pow(_q_point[_qp](0)-0.125,2.0) + std::pow(_q_point[_qp](1)-0.080,2.0)) > 0.035 ) {
//      if ( _q_point[_qp](1)-0.100 > 0.0 ) {
//        _H[_qp] = B*_gc[_qp]/4/(0.5*_l[_qp]); } } }
  if ( std::abs(_q_point[_qp](1)+_q_point[_qp](0)-200.0) < 50.0 ) {
    if ( std::abs(_q_point[_qp](1)-_q_point[_qp](0)) < (0.25*_l[_qp]) ) {
      _H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](1)-_q_point[_qp](0))/(0.5*_l[_qp]) ); } }
}

void
ComputeMieGruneisenPFFractureStress::computeQpStress()
{
  // Calculate pressure from Mie Gruneisen (Menon, 2014), (Zhang, 2011)
  // https://en.wikipedia.org/wiki/Mie%E2%80%93Gruneisen_equation_of_state
  Real K0, delta, eta, temperature, peos;
  RankTwoTensor stress_eos, stress, stress_cpl;
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  _bulk_modulus[_qp] = K0;
  delta = _mechanical_strain[_qp].trace();
  eta = - delta;
  temperature = _temperature[_qp];
  peos = - K0 * eta * (1.0 - (_Gamma * eta / 2.0)) / std::pow((1.0 - _s * eta), 2.0)
         - _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature);
  _pressure_eos[_qp] = peos;
  stress_eos = peos * I2;
  _stress_eos_elastic[_qp] = stress_eos;
  stress_cpl = _elasticity_tensor[_qp] * _mechanical_strain[_qp] - K0 * delta * I2;
  _stress_cpl_elastic[_qp] = stress_cpl;
  stress = stress_eos + stress_cpl;

  // Create the positive and negative projection tensors
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
  std::vector<Real> eigval;
  RankTwoTensor eigvec;
  RankFourTensor Ppos = stress_cpl.positiveProjectionEigenDecomposition(eigval, eigvec);

  // Project the positive and negative stresses
  Real peos_pos;
  RankTwoTensor stress_eos_pos, stress_eos_neg, stress_cpl_pos, stress_cpl_neg, stress0pos, stress0neg;
  peos_pos = (std::abs(peos) + peos) / 2.0;
  stress_eos_pos = peos_pos * I2;
  stress_eos_neg = stress_eos - stress_eos_pos;
  stress_cpl_pos = Ppos * stress_cpl;
  stress_cpl_neg = stress_cpl - stress_cpl_pos;
  stress0pos = stress_eos_pos + stress_cpl_pos;
  stress0neg = stress - stress0pos;

  // Compute the positive and negative elastic energies
  Real F_pos, F_neg;
  Real A, B, C;
  A = K0 * (_Gamma * (1.0/(2.0 * std::pow(_s,2.0)) + 0.5 - 1.0 / _s) - 1.0);
  B = K0 * (_Gamma * (1.0 / _s - 0.5) + 1.0);
  C = - K0 * _Gamma / (2.0 * std::pow(_s,2.0));
  if (delta>=0.0) {
    F_pos = A / (_s - std::pow(_s,2.0) * eta) - B * std::log(1.0 - _s * eta) / _s + C * eta - A / _s +
            _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature) * delta;
    F_neg = 0.0;
  }
  else {
    F_pos = 0.0;
    F_neg = A / (_s - std::pow(_s,2.0) * eta) - B * std::log(1.0 - _s * eta) / _s + C * eta - A / _s +
            _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature) * delta;
  }

  F_pos += (stress_cpl_pos).doubleContraction(_mechanical_strain[_qp]) / 2.0;
  F_neg += (stress_cpl_neg).doubleContraction(_mechanical_strain[_qp]) / 2.0;

  _stress[_qp] = stress0pos * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg;

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];

  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];

  _Jacobian_mult[_qp] = (I4sym - (1 - _D[_qp]) * Ppos) * _elasticity_tensor[_qp];

  // Calculate bulk-viscosity stress term
  Real trD, jacob, q_bv;
  // trD = _strain_increment[_qp].trace() / _dt;
  trD = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;
  // trD /= (1.0 + _mechanical_strain_old[_qp].trace());
  jacob = 1.0 + _mechanical_strain[_qp].trace();
  q_bv = 0.0;
  if (jacob < 1.0) {
    q_bv = ( _C0 * _density[_qp] * trD * std::abs(trD) * std::pow(_Le,2.0) / std::pow(jacob,2.0) ) + ( _C1 * _density[_qp] * _sound_speed * trD * _Le / jacob );
  }
  _stress[_qp] += q_bv * I2;


  // // Assign history variable
  Real hist_variable = _H_old[_qp];
  if (_use_snes_vi_solver)
  {
    _H[_qp] = F_pos;

    if (_use_current_hist)
      hist_variable = _H[_qp];
  }
  else
  {
    if (F_pos > _H_old[_qp])
      _H[_qp] = F_pos;
    else
      _H[_qp] = _H_old[_qp];

    if (_use_current_hist)
      hist_variable = _H[_qp];

    if (hist_variable < _barrier[_qp])
      hist_variable = _barrier[_qp];
  }

  // Elastic free energy density
  _E[_qp] =
      hist_variable * _D[_qp] + F_neg - _pressure[_qp] * _mechanical_strain[_qp].trace() * _I[_qp];
  _dEdc[_qp] =
      hist_variable * _dDdc[_qp] - _pressure[_qp] * _mechanical_strain[_qp].trace() * _dIdc[_qp];
  _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp] -
                 _pressure[_qp] * _mechanical_strain[_qp].trace() * _d2Id2c[_qp];
}
