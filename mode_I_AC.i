#This input uses PhaseField-Nonconserved Action to add phase field fracture bulk rate kernels
# Mie Gruneisen EOS
#Units
#-lenght um
#-time ns
#-Pressure GPa

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmin = 0.0
  ymin = 0.0
  xmax = 1.0
  ymax = 1.0
  elem_type = QUAD4
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
    family = LAGRANGE
    order = FIRST
  [../]
  [./disp_y]
    family = LAGRANGE
    order = FIRST
  [../]
  [./c]
    family = LAGRANGE
    order = FIRST
  [../]
  [./temp]
  [../]
[]

[AuxVariables]
  [./resid_x]
  [../]
  [./resid_y]
  [../]
  [./resid_z]
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
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
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
  [./We]
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
#  [./cIC]
#    type = FunctionIC
#    variable = c
#    function = xyfunc_c
#  [../]
[]

[Functions]
  [./tfunc]
    type = ParsedFunction
    value = '0.01*t'
  [../]
  [./tfunc2]
    type = ParsedFunction
    value = '--0.01*t'
  [../]
  [./tempfunc]
    type = ParsedFunction
    value = '300'
  [../]
  [./xyfunc_c]
    type = ParsedFunction
    value = 'min(1.0*exp(-abs(y-0.5)/l)*if((x)/0.5,0,1.0),1.0)'
    vals = 4.0e-2
    vars = l
  [../]
[]

[Kernels]
  [./DynamicTensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
  [./ACbulk]
    type = AllenCahn
    variable = c
    f_name = E_elpl
  [../]
  [./ACinterface]
    type = ACInterface
    variable = c
  [../]
  [./dcdt]
    type = TimeDerivative
    variable = c
  [../]
  [./inertia_x]
    type = InertialForce
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = 0.3025
    gamma = 0.6
    use_displaced_mesh = False
  [../]
  [./inertia_y]
    type = InertialForce
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = 0.3025
    gamma = 0.6
    use_displaced_mesh = False
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
  [./ACOffDiag]
    type = AllenCahnElasticEnergyOffDiag
    variable = c
    displacements = 'disp_x disp_y'
  [../]
  [./heat_conduction]
    type = HeatConduction
    variable = temp
    use_displaced_mesh = False
  [../]
  [./hct]
    type = HeatConductionTimeDerivative
    variable = temp
    specific_heat = specific_heat
    density_name = density
    use_displaced_mesh = False
  [../]
  [./thermoelastic_heat_source]
    type = ThermalExpansionHeatSourceFiniteStrainMieGruneisen
    variable = temp
    G_Gruneisen = 0.7
    s_UsUp = 1.79 
    bulk_modulus_ref = 17.822 # Menon 2014: K_0 = rho_0 * c_0 * c_0
    reference_temperature =300.0
    c = c
    specific_heat = specific_heat
    density_name = density
    C0 = 0.955e-5
    C1 = 2.865e-2
  [../]
[]

[UserObjects]
  [./prop_read]
    type = ElementPropertyReadFile
    prop_file_name = 'euler_ang_file_010.txt'
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
  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
  [../]
  [./stress_xy]
    type = RankTwoAux
    variable = stress_xy
    rank_two_tensor = stress
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
  [./We]
    type = MaterialRealAux
    variable = We
    property = W0e_pos
  [../]
[]

[BCs]
  [./ydisp1]
    type = FunctionPresetBC
    variable = disp_y
    boundary = top
    function = tfunc
  [../]  
 # [./ydisp2]
 #   type = FunctionPresetBC
 #   variable = disp_y
 #   boundary = bottom
 #   function = tfunc2
 # [../]  
  [./fixy]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
[]

[Materials]
  [./crysp_hmx]
    type = FiniteStrainCrystalPlasticityPFFractureStressMieGruneisen
    gtol = 1e-2
    abs_tol = 1e-5
    slip_incr_tol = 0.0025
    maximum_substep_iteration = 1
    gen_random_stress_flag = true
    slip_sys_file_name = input_slip_sys_HMX_austin.txt
    nss = 10
    num_slip_sys_flowrate_props = 2 #Number of properties in a slip system
    flowprops = '1 10 0.0 0.1'
    #flowprops = '1 1 1.0e-3 0.1 2 3 1.46e-3 0.1 4 4 2.0e-3 0.1 5 5 5.6e-6 0.1 6 6 17.7 0.1 7 8 2.04e-3 0.1 9 10 34.9e-3 0.1'
    hprops = '1.0 9.34e-3 0.10303 0.15573 2.5'
    gprops = '1 10 0.10303'
    C0 = 1.91e-5
    C1 = 5.73e-2
    G_Gruneisen = 0.7
    s_UsUp = 1.79
    bulk_modulus_ref = 17.822 # Menon 2014: K_0 = rho_0 * c_0 * c_0
    c = c
    temp = temp 
    reference_temperature = 300.0
    plastic_factor = 0.0
    kdamage = 1e-06
    F_name = E_elpl
    specific_heat = specific_heat
    density_name = density
  [../]
  [./strain_hmx]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y'
  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '0.002 4.0e-2 100.0'
  [../]
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1.0/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l'
  [../]
  [./elasticity_tensor_hmx]
    type = ComputeElasticityTensorCP
    C_ijkl = '21.15 10.18 9.77 0.058 4.05 -0.18 20.34 13.35 0.23 6.96 0.14 21.27 -0.004 5.01 0.19 8.79 0.32 4.16 6.20 0.22 12.00' 
    fill_method = symmetric21
    read_prop_user_object = prop_read
  [../]
  [./HMX]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity density'
    prop_values = '0.31E-6 1.891' 
  [../]
  [./specific_heat_hmx]
    type = HeatCapacityMD
    temp = temp
  [../]
#  [./specific_heat_hmx]
#    type = GenericConstantMaterial
#    prop_names = 'specific_heat'
#    prop_values = '2357.3e-6'
#  [../]
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201 hypre boomeramg 20'
  line_search = 'none'

  #l_max_its = 10
  #nl_max_its = 10

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-5

  start_time = 0.0
  end_time = 10000.0
  dtmax = 1.0e-4
  dtmin = 1.0e-9
  [./TimeStepper]
    type = ConstantDT
    dt = 1.0e-4
  [../]
[]

[Outputs]
  exodus = true
  interval = 100
  checkpoint = true
  num_files = 4
[]
