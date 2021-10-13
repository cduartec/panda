#Shock in PBX using sylgard as binder
#Elastic
# Mie Gruneisen EOS
#Units
#-lenght um
#-time ns
#-Pressure GPa

[Mesh]
  type = FileMesh
  file = pbx_h10_o1.msh
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Modules]
  [./PhaseField]
    [./Nonconserved]
      [./c]
        free_energy = E_elpl
        kappa = kappa_op
        mobility = L
      [../]
    [../]
  [../]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = FINITE
        additional_generate_output = 'stress_xx stress_yy stress_xy stress_zz stress_zx stress_zy'
        save_in = 'resid_x resid_y'
      [../]
    [../]
  [../]
[]

[Variables]
  [./temp]
  [../]
[]

[AuxVariables]
  [./resid_x]
  [../]
  [./resid_y]
  [../]
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vel_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./accel_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./vel_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./accel_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./von_mises]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rotation_norm]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Wp]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./We]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fpxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fpyy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_incr1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_incr2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_incr3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_incr4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_incr5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_incr6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_incr7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_incr8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_incr9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_incr10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./acc_slip]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./h_max]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pressure]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dcdx]
    order = FIRST
    family = MONOMIAL
  [../]
  [./dcdy]
    order = FIRST
    family = MONOMIAL
  [../]
  [./friction_normal_force]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slide_velocity_parallel]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./crack_rho]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gc]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Q_dot_friction]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Q_dot_therm1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Q_dot_therm2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Q_dot_shock]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Q_dot_p]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./c11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./TrCEdot]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[ICs]
  [./tempIC]
    type = FunctionIC
    variable = temp
    function = tempfunc
  [../]
[]

[Functions]
  [./tfunc]
    type = ParsedFunction
    value = '0.1'
    #value = 'if(t<400,1.0,0)'
  [../]
  [./tempfunc]
    type = ParsedFunction
    value = '300'
  [../]
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = 0.3025
    gamma = 0.6
  [../]
  [./inertia_y]
    type = InertialForce
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = 0.3025
    gamma = 0.6
  [../]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = c
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = c
  [../]
  [./heat_conduction]
    type = HeatConduction
    variable = temp
  [../]
  [./hct]
    type = HeatConductionTimeDerivative
    variable = temp
    specific_heat = specific_heat
    density_name = density
  [../]
  [./thermoelastic_heat_source]
    type = ThermalExpansionShockHeatSourceFiniteStrainMieGruneisenConstH
    variable = temp
    G_Gruneisen = 0.7
    s_UsUp = 2.29 #New MD results
    bulk_modulus_ref = 15.31 # New MD results
    reference_temperature =300.0
    c = c
    specific_heat = specific_heat
    density_name = density
    beta_v = 0.4
    C0 = 0.1
    C1 = 1.0
    c_l = 2.88
    h_e = 10.0
    block = 'hmx'
  [../]
  [./thermoelastic_heat_source_binder]
    type = ThermalExpansionShockHeatSourceFiniteStrainMieGruneisenConstHNoCpl
    variable = temp
    G_Gruneisen = 1.5
    s_UsUp = 1.66 #New MD results
    bulk_modulus_ref = 0.24 # New MD results
    reference_temperature =300.0
    c = c
    specific_heat = specific_heat
    density_name = density
    beta_v = 0.4
    C0 = 0.1
    C1 = 1.0
    c_l = 1.13 #bulk speed of sound
    h_e = 10.0
    block = 'GB'
  [../]
  [./plasticheat]
    type = PlasticHeatingSourceDamage
    variable = temp
    plastic_factor = 0.5
    W0p = 'W0p'
    W0p_broken = 'W0p_broken'
    #dW0p_dstrain = '0'
    #dW0p_broken_dstrain = '0'
    c = c
    displacements = 'disp_x disp_y'
    block = 'hmx'
  [../]
  [./friction]
    type = CrackFrictionHeatSource
    variable = temp
    friction_coefficient = 0.5
    dcdx = dcdx
    dcdy = dcdy
    displacements = 'disp_x disp_y'
  [../]
[]

[UserObjects]
  [./prop_read]
    type = ElementPropertyReadFile
    prop_file_name = 'euler_angles_010.txt'
    nprop = 3
    read_type = element
  [../]
  [./prop_read2]
    type = ElementPropertyReadFile
    prop_file_name = 'euler_angles_pbx_h10_n47.txt'
    nprop = 3
    read_type = element
  [../]
[]

[AuxKernels]
  [./strain_xx]
    type = RankTwoAux
    variable = strain_xx
    rank_two_tensor = total_strain
    index_j = 0
    index_i = 0
  [../]
  [./strain_yy]
    type = RankTwoAux
    variable = strain_yy
    rank_two_tensor = total_strain
    index_j = 1
    index_i = 1
  [../]
  [./strain_xy]
    type = RankTwoAux
    variable = strain_xy
    rank_two_tensor = total_strain
    index_j = 0
    index_i = 1
  [../]
  [./accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = 0.3025
  [../]
  [./vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = 0.6
  [../]
  [./accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = 0.3025
  [../]
  [./vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = 0.6
  [../]
  [./von_mises_kernel]
    #Calculates the von mises stress and assigns it to von_mises
    type = RankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./rotation_norm]
    #Calculates the von mises stress and assigns it to von_mises
    type = RankTwoScalarAux
    variable = rotation_norm
    rank_two_tensor = update_rot
    execute_on = timestep_end
    scalar_type = L2Norm
  [../]
  [./Wp]
    type = MaterialRealAux
    variable = Wp
    property = W0p
  [../]
  [./We]
    type = MaterialRealAux
    variable = We
    property = W0e_pos
  [../]
  [./fpxx]
    type = RankTwoAux
    variable = fpxx
    rank_two_tensor = fp
    index_j = 0
    index_i = 0
  [../]
  [./fpyy]
    type = RankTwoAux
    variable = fpyy
    rank_two_tensor = fp
    index_j = 1
    index_i = 1
  [../]
  [./fpxy]
    type = RankTwoAux
    variable = fpyy
    rank_two_tensor = fp
    index_j = 0
    index_i = 1
  [../]
  [./gss1]
    type = MaterialStdVectorAux
    variable = gss1
    property = gss
    index = 0
    execute_on = timestep_end
  [../]
  [./gss2]
    type = MaterialStdVectorAux
    variable = gss2
    property = gss
    index = 1
    execute_on = timestep_end
  [../]
 [./gss3]
    type = MaterialStdVectorAux
    variable = gss3
    property = gss
    index = 2
    execute_on = timestep_end
  [../]
  [./gss4]
    type = MaterialStdVectorAux
    variable = gss4
    property = gss
    index = 3
    execute_on = timestep_end
  [../]
  [./gss5]
    type = MaterialStdVectorAux
    variable = gss5
    property = gss
    index = 4
    execute_on = timestep_end
  [../]
  [./gss6]
    type = MaterialStdVectorAux
    variable = gss6
    property = gss
    index = 5
    execute_on = timestep_end
  [../]
  [./gss7]
    type = MaterialStdVectorAux
    variable = gss7
    property = gss
    index = 6
    execute_on = timestep_end
  [../]
  [./gss8]
    type = MaterialStdVectorAux
    variable = gss8
    property = gss
    index = 7
    execute_on = timestep_end
  [../]
[./gss9]
    type = MaterialStdVectorAux
    variable = gss9
    property = gss
    index = 8
    execute_on = timestep_end
  [../]
  [./gss10]
    type = MaterialStdVectorAux
    variable = gss10
    property = gss
    index = 9
    execute_on = timestep_end
  [../]
  [./slip_incr1]
    type = MaterialStdVectorAux
    variable = slip_incr1
    property = slip_incr_out
    index = 0
    execute_on = timestep_end
  [../]
  [./slip_incr2]
    type = MaterialStdVectorAux
    variable = slip_incr2
    property = slip_incr_out
    index = 1
  [../]
  [./slip_incr3]
    type = MaterialStdVectorAux
    variable = slip_incr3
    property = slip_incr_out
    index = 2
  [../]
  [./slip_incr4]
    type = MaterialStdVectorAux
    variable = slip_incr4
    property = slip_incr_out
    index = 3
  [../]
  [./slip_incr5]
    type = MaterialStdVectorAux
    variable = slip_incr5
    property = slip_incr_out
    index = 4
  [../]
[./slip_incr6]
    type = MaterialStdVectorAux
    variable = slip_incr6
    property = slip_incr_out
    index = 5
  [../]
  [./slip_incr7]
    type = MaterialStdVectorAux
    variable = slip_incr7
    property = slip_incr_out
    index = 6
  [../]
  [./slip_incr8]
    type = MaterialStdVectorAux
    variable = slip_incr8
    property = slip_incr_out
    index = 7
  [../]
  [./slip_incr9]
    type = MaterialStdVectorAux
    variable = slip_incr9
    property = slip_incr_out
    index = 8
  [../]
  [./slip_incr10]
    type = MaterialStdVectorAux
    variable = slip_incr10
    property = slip_incr_out
    index = 9
  [../]
  [./acc_slip]
    type = MaterialRealAux
    variable = acc_slip
    property = acc_slip
  [../]
  [./tau1]
    type = MaterialStdVectorAux
    variable = tau1
    property = tau_out
    index = 0
    execute_on = timestep_end
  [../]
  [./tau2]
    type = MaterialStdVectorAux
    variable = tau2
    property = tau_out
    index = 1
    execute_on = timestep_end
  [../]
  [./tau3]
    type = MaterialStdVectorAux
    variable = tau3
    property = tau_out
    index = 2
    execute_on = timestep_end
  [../]
  [./tau4]
    type = MaterialStdVectorAux
    variable = tau4
    property = tau_out
    index = 3
    execute_on = timestep_end
  [../]
  [./tau5]
    type = MaterialStdVectorAux
    variable = tau5
    property = tau_out
    index = 4
    execute_on = timestep_end
  [../]
  [./tau6]
    type = MaterialStdVectorAux
    variable = tau6
    property = tau_out
    index = 5
    execute_on = timestep_end
  [../]
  [./tau7]
    type = MaterialStdVectorAux
    variable = tau7
    property = tau_out
    index = 6
    execute_on = timestep_end
  [../]
  [./tau8]
    type = MaterialStdVectorAux
    variable = tau8
    property = tau_out
    index = 7
    execute_on = timestep_end
  [../]
  [./tau9]
    type = MaterialStdVectorAux
    variable = tau9
    property = tau_out
    index = 8
    execute_on = timestep_end
  [../]
  [./tau10]
    type = MaterialStdVectorAux
    variable = tau10
    property = tau_out
    index = 9
    execute_on = timestep_end
  [../]
  [./h_max]
    type = ElementLengthAux
    variable = h_max
    method = min
    execute_on = initial
  [../]
  [./pressure]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = pressure
    scalar_type = Hydrostatic
  [../]
  [./dcdx]
    type = VariableGradientComponent
    variable = dcdx
    gradient_variable = c
    component = 'x'
    use_displaced_mesh = true
  [../]
  [./dcdy]
    type = VariableGradientComponent
    variable = dcdy
    gradient_variable = c
    component = 'y'
    use_displaced_mesh = true
  [../]
  [./crack_rho]
    type = MaterialRealAux
    variable = crack_rho
    property = crack_surface_density
  [../]
  [./Q_dot_friction]
    type = MaterialRealAux
    variable = Q_dot_friction
    property = heat_source_rate
  [../]
  [./f_x]
    type = MaterialStdVectorAux
    variable = f_x
    property = friction_force
    index = 0
  [../]
  [./f_y]
    type = MaterialStdVectorAux
    variable = f_y
    property = friction_force
    index = 1
  [../]
  [./s_x]
    type = MaterialStdVectorAux
    variable = s_x
    property = slide_velocity
    index = 0
  [../]
  [./s_y]
    type = MaterialStdVectorAux
    variable = s_y
    property = slide_velocity
    index = 1
  [../]
  [./friction_normal_force]
    type = MaterialRealAux
    variable = friction_normal_force
    property = friction_normal_force
  [../]
  [./slide_velocity_parallel]
    type = MaterialRealAux
    variable = slide_velocity_parallel
    property = slide_velocity_parallel
  [../]
  [./gc]
    type = MaterialRealAux
    variable = gc
    property = gc_prop
  [../]
  [./Q_dot_therm1]
    type = MaterialRealAux
    variable = Q_dot_therm1
    property = heat_rate_therm1
  [../]
  [./Q_dot_therm2]
    type = MaterialRealAux
    variable = Q_dot_therm2
    property = heat_rate_therm2
  [../]
  [./Q_dot_shock]
    type = MaterialRealAux
    variable = Q_dot_shock
    property = heat_rate_vis
  [../]
  [./Q_dot_p]
    type = MaterialRealAux
    variable = Q_dot_p
    property = heat_rate_p
  [../]
  [./c11]
    type = RankFourAux
    variable = c11
    rank_four_tensor = elasticity_tensor
    index_i = 0
    index_j = 0
    index_k = 0
    index_l = 0
    execute_on = timestep_end
  [../]
  [./TrCEdot]
    type = MaterialRealAux
    variable = TrCEdot
    property = Tr_E_dot
  [../]
[]

[BCs]
  [./BC_velocity]
    type = PresetVelocity
    variable = disp_x
    boundary = left
    function = tfunc
    velocity = 1.0
  [../] 
  [./fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'right'
    value = 0
  [../]
  [./fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'top bottom'
    value = 0
  [../]
[]

[Materials]
  [./crysp_hmx]
    type = FiniteStrainCPPFFractureMieGruneisenConstH
    gtol = 1e-2
    abs_tol = 1e-6
    slip_incr_tol = 0.0025
    maximum_substep_iteration = 1
    gen_random_stress_flag = true
    slip_sys_file_name = input_slip_sys_HMX_austin_P21n.txt
    nss = 10
    num_slip_sys_flowrate_props = 2 #Number of properties in a slip system
    #flowprops = '1 10 0.0  0.1'
    flowprops = '1 1 1.0e-3 0.1 2 3 1.46e-3 0.1 4 4 2.0e-3 0.1 5 5 5.6e-6 0.1 6 6 17.7 0.1 7 8 2.04e-3 0.1 9 10 34.9e-3 0.1'
    hprops = '1.0 9.34e-3 0.10303 0.15573 2.5'
    gprops = '1 10 0.10303'
    C0 = 0.1
    C1 = 1.0
    c_l = 2.88
    G_Gruneisen = 0.7
    s_UsUp = 2.29 #New MD results
    bulk_modulus_ref = 14.24 #New MD results
    c = c
    temp = temp
    reference_temperature = 300.0
    plastic_factor = 0.0
    kdamage = 1e-06
    F_name = E_elpl
    specific_heat = specific_heat
    density_name = density
    h_e = 10.0
    p = pressure
    min_line_search_step_size = 0.01
    block = 'hmx' 
  [../]
  [./crysp_binder]
    type = FiniteStrainCPPFFractureMieGruneisenConstH
    gtol = 1e-2
    abs_tol = 1e-6
    slip_incr_tol = 0.025
    maximum_substep_iteration = 1
    gen_random_stress_flag = true
    slip_sys_file_name = input_slip_sys_HMX_austin_P21n.txt
    nss = 10
    num_slip_sys_flowrate_props = 2 #Number of properties in a slip system
    flowprops = '1 10 0.0  0.1'
    hprops = '1.0 9.34e-3 0.10303 0.15573 2.5'
    gprops = '1 10 0.10303'
    C0 = 0.1
    C1 = 1.0
    c_l = 1.13 #Bulk speed of sound
    G_Gruneisen = 1.5
    s_UsUp = 1.66 #New MD results
    bulk_modulus_ref = 1.325 #New MD results
    c = c
    temp = temp
    reference_temperature = 300.0
    plastic_factor = 0.0
    kdamage = 1e-06
    F_name = E_elpl
    specific_heat = specific_heat
    density_name = density
    h_e = 10.0
    p = pressure
    min_line_search_step_size = 0.01
    block = 'GB'
  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '0.002 80 10'
    block = 'hmx'
  [../]
  [./pf_interface]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '0.400 80 10'
    block = 'GB'
  [../]
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = 0#'1.0/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l'
  [../]
  [./elasticity_tensor_hmx]
    type = ComputeElasticityTensorCP
    #C_ijkl = '1111 1122 1133 1123 1113 1112 2222 2233 2223 2213 2212 3333 3323 3313 3312 2323 2313 2312 1313 1312 1212'
    C_ijkl = '25.1 9.7 12.8 0.0 -1.3 0.0 22.3 11.8 0.0 4.6 0.0 21.8 0.0 1.4 0.0 9.7 0.0 3.18 11.036 0.0 8.66'
    fill_method = symmetric21
    read_prop_user_object = prop_read2
    block = 'hmx' 
  [../]
  [./elasticity_tensor_binder]
    type = ComputeElasticityTensorCP
    #C_ijkl = '1111 1122 1133 1123 1113 1112 2222 2233 2223 2213 2212 3333 3323 3313 3312 2323 2313 2312 1313 1312 1212'
    C_ijkl = '1.397 1.289 1.289 0.000 0.000 0.000 1.397 1.289 0.000 0.000 0.000 1.397 0.000 0.000 0.000 0.054 0.000 0.000 0.054 0.000 0.054'
    fill_method = symmetric21
    read_prop_user_object = prop_read
    block = 'GB'
  [../]
  [./density_hmx]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '1.85'
    block = 'hmx' 
  [../]
  [./density_binder]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '1.03'
    block = 'GB' 
  [../]
  [./thermal_conductivity_hmx]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '0.31e-6'
    block = 'hmx' 
  [../]
  [./thermal_conductivity_binder]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '0.27e-6'
    block = 'GB' 
  [../]
  [./specific_heat_hmx]
    type = HeatCapacityMD
    temp = temp
    block = 'hmx' 
  [../]
  [./specific_heat_binder]
    type = GenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '1624.0e-6'
    block = 'GB' 
  [../]
  [./crackfrictionheatenergy]
    type = ComputeCrackFrictionHeatEnergyDienes
    friction_coefficient = 0.5
    dcdx = dcdx
    dcdy = dcdy
    c = c
    l = 40
  [../]
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./Volume]
    type = VolumePostprocessor
    execute_on = 'initial'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201 hypre boomeramg 20'
  line_search = 'none'

  #l_max_its = 100
  #nl_max_its = 100

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-7

  start_time = 0.0
  end_time = 20000.0
  dtmax = 0.01 #tmax = hmin/(C0 = 3.0 um/ns)
  dtmin = 1.0e-6
[]

[Outputs]
  exodus = true
  interval = 100
  checkpoint = True
[]
