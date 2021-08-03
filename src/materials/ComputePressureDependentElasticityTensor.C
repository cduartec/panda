//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputePressureDependentElasticityTensor.h"
#include "RotationTensor.h"

registerMooseObject("pandaApp", ComputePressureDependentElasticityTensor);

template <>
InputParameters
validParams<ComputePressureDependentElasticityTensor>()
{
  InputParameters params = validParams<ComputeElasticityTensor>();
  params.addClassDescription("Compute pressure dependent elasticity tensor.");
  params.addParam<UserObjectName>("read_prop_user_object",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element "
                                  "specific property values from file");
  params.addRequiredParam<std::vector<Real>>("A0_ijkl",
                                             "Coefficient A0");
  params.addRequiredParam<std::vector<Real>>("A1_ijkl",
                                             "Coefficient A1");
  params.addRequiredParam<std::vector<Real>>("A2_ijkl",
                                             "Coefficient A2");
  params.addParam<MooseEnum>(
      "fill_method0", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  params.addParam<MooseEnum>(
      "fill_method1", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  params.addParam<MooseEnum>(
      "fill_method2", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  params.addRequiredCoupledVar("p", "Pressure");
  return params;
}

ComputePressureDependentElasticityTensor::ComputePressureDependentElasticityTensor(
    const InputParameters & parameters)
  : ComputeElasticityTensor(parameters),
    _read_prop_user_object(isParamValid("read_prop_user_object")
                               ? &getUserObject<ElementPropertyReadFile>("read_prop_user_object")
                               : nullptr),
    _Euler_angles_mat_prop(declareProperty<RealVectorValue>("Euler_angles")),
    _crysrot(declareProperty<RankTwoTensor>("crysrot")),
    _R(_Euler_angles),
    _Aijkl0(getParam<std::vector<Real>>("A0_ijkl"),
            (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method0")),
    _Aijkl1(getParam<std::vector<Real>>("A1_ijkl"),
            (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method1")),
    _Aijkl2(getParam<std::vector<Real>>("A2_ijkl"),
            (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method2")),
    _p(coupledValue("p")),
    _p_name(getVar("p", 0)->name()),
    _delasticity_tensor_dp(
        declarePropertyDerivative<RankFourTensor>(_elasticity_tensor_name, _p_name))
{ 
  // the base class guarantees constant in time, but in this derived class the
  // tensor will rotate over time once plastic deformation sets in
  revokeGuarantee(_elasticity_tensor_name, Guarantee::CONSTANT_IN_TIME);
  
  // the base class performs a passive rotation, but the crystal plasticity
  // materials use active rotation: recover unrotated _Cijkl here
  _Aijkl0.rotate(_R.transpose());
  _Aijkl1.rotate(_R.transpose());
  _Aijkl2.rotate(_R.transpose());
  _Cijkl.rotate(_R.transpose());
}

void
ComputePressureDependentElasticityTensor::assignEulerAngles()
{
  if (_read_prop_user_object)
  {
    _Euler_angles_mat_prop[_qp](0) = _read_prop_user_object->getData(_current_elem, 0);
    _Euler_angles_mat_prop[_qp](1) = _read_prop_user_object->getData(_current_elem, 1);
    _Euler_angles_mat_prop[_qp](2) = _read_prop_user_object->getData(_current_elem, 2);
  }
  else
    _Euler_angles_mat_prop[_qp] = _Euler_angles;
}

void
ComputePressureDependentElasticityTensor::computeQpElasticityTensor()
{
  Real pb  =  -1.0 * _p[_qp];
  // Properties assigned at the beginning of every call to material calculation
  assignEulerAngles();
  _R.update(_Euler_angles_mat_prop[_qp]);
  _crysrot[_qp] = _R.transpose();
  // Assign elasticity tensor at a given quad point
  if ( pb >  0.0)//No substepping
  {
    _elasticity_tensor[_qp] = _Cijkl + _Aijkl1 * pb + _Aijkl2 * pb * pb;
    _delasticity_tensor_dp[_qp] = _Aijkl1 + 2.0 *_Aijkl2 * pb;
  }
  else 
  { 
    _elasticity_tensor[_qp] = _Cijkl;
    _delasticity_tensor_dp[_qp] = 0.0;
  }
  _elasticity_tensor[_qp].rotate(_crysrot[_qp]);
  // Define derivative of elasticity tensor with respect to pressure.
}
