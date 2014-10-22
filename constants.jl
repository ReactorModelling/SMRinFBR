################################################################################
#                                   Constants                                  #
################################################################################
const global Nz      = 200   # Number of collocation points in z direction
const global Ncomp   = 6    # Number of chemical components
const Z,A,B,Q        = colloc(Nz-2,1,1) # Collocation points and matrices
global Z
global A
global B
global Q
const global I       = eye(Nz) # Identity matrix
const global R       = 8.3145 # Gas constant [J/K mol]

################################################################################
#                              Reactor parameters                              #
################################################################################
const global void     = 0.528   # Void fraction
const global dInner   = 0.102   # Inner tube diameter [m]
const global Ta       = 1100    # Ambient temperature [K]
const global U        = 56.783  # Heat transfer coefficient [J/K m2 s]
const global rhoCat   = 2355.2; # Catalyst density          [kgcat/m^3]
const global dParticle = 0.0173;# Particle diameter         [m]
const global lambdaSt  = 52     # Heat coef. for tube metal[W/mK]
const global rInner = dInner/2  # Inner tube radius [m]
const global rOuter = 0.066     # Outer tube radius [m]



################################################################################
#                              Molar masses                                    #
################################################################################
const global molarMass = [
                            16.04 # CH4
                            28.01 # CO
                            44.01 # CO2
                            2.02 # H2
                            18.02 # H2O
                            28.01 # N2
                          ]*1e-3 # [kg/mol]

################################################################################
#                              Reaction rate data                              #
################################################################################

# Preexponential factors for the rate constants      [mol/kgcat s]
const global aj = [
                    4.255e15 # Reaction 1
                    1.955e6 # Reaction 2
                    1.020E15 # Reaction 3
                  ]/3.6

# Preexponential factors for the adsorbtion constants
const global ax = [
                    6.65e-9 # Factor for CH4 [Pa^-1]
                    8.23e-10 # Factor for CO            [Pa^-1]
                    6.12e-14 # Factor for H2            [Pa^-1]
                    1.77e5 # Factor for H2O           [-]
                  ]


# Activation energies for the reactions              [J/mol]
const global actEn = [
                        240.1e3 # Activation energy for rx. 1
                        67.13e3 # Activation energy for rx. 2
                        243.9e3 # Activation energy for rx. 3
                     ]
const global N = [-1     1   0   3   -1   0    # rx 1
                   0    -1   1   1   -1   0    # rx 2
                  -1     0   1   4   -2   0];  # rx 3

################################################################################
#                              Enthalpy data                                   #
################################################################################
# Reaction enthalpies at 298K                        [J/mol]
const global ent298 = [
                        206.1e3 # Enthalpy for rx. 1 
                        -41.15e3 # Enthalpy for rx. 2
                        164.9e3 # Enthalpy for rx. 3
                      ]

# Reaction enthalpies at 948K                        [J/mol]
const global ent948 = [
                        224.0e3 # Enthalpy for rx. 1 
                        -37.30e3 # Enthalpy for rx. 2
                        187.5e3 # Enthalpy for rx. 3
                      ]

# Adsorption enthalpies                              [J/mol]
const global adEnt = [
                       -38.28e3 # Enthalpy of adsorption for CH4
                       -70.65e3 # Enthalpy of adsorption for CO
                       -82.90e3 # Enthalpy of adsorption for H2
                       88.68e3 # Enthalpy of adsorption for H2O
                     ]
################################################################################
#                              Heat capacity data                              #
################################################################################

const global cpCoeff = [
                        1.925e4  5.213e1  1.197e-2 -1.132e-5 # CH4
                        3.087e4 -1.285e1  2.789e-2 -1.272e-5 # CO
                        1.980e4  7.344e1 -5.602e-2  1.715e-5 # CO2
                        2.714e4  9.274e0 -1.381e-2  7.645e-6 # H2
                        3.224e4  1.924e0  1.055e-2  3.596e-6 # H2O
                        3.115e4 -1.357e1  2.268e-2 -1.168e-5 # N2
                        ]'*1e-3

################################################################################
#                              Viscosity data                                  #
################################################################################

const global b = [
                  1.00e-6 # Coefficient for CH4      [kg/msK^0.5]
                  1.50e-6 # Coefficient for CO       [kg/msK^0.5]
                  1.50e-6 # Coefficient for CO2      [kg/msK^0.5]
                  0.65e-6 # Coefficient for H2       [kg/msK^0.5]
                  1.74e-6 # Coefficient for H2O      [kg/msK^0.5]
                  1.40e-6 # Coefficient for N2       [kg/msK^0.5]
                  ]
const global s = [
                  165 # Coefficient for CH4      [K]
                  220 # Coefficient for CO       [K]
                  220 # Coefficient for CO2      [K]
                  67  # Coefficient for H2       [K]
                  626 # Coefficient for H2O      [K]
                  108 # Coefficient for N2       [K]
                  ]
################################################################################
#                               Conductivity data                              #
################################################################################
const global lambda = [
                      -1.8690e-03   8.7270e-05   1.1790e-07  -3.6140e-11  # CH4 [W/mK^n]
                       5.0670e-04   9.1025e-05  -3.5240e-08   8.1990e-12  # CO  [W/mK^n]
                      -7.2150e-03   8.0150e-05   5.4770e-09  -1.0530e-11  # CO2 [W/mK^n]
                       8.0990e-03   6.6890e-04  -4.1580e-07   1.5620e-10  # H2  [W/mK^n]
                       7.3410e-03  -1.0130e-05   1.8010e-07  -9.1000e-11  # H2O [W/mK^n]
                       3.9190e-04   9.9660e-05  -5.0670e-08   1.5040e-11  # N2  [W/mK^n]
                      ]

################################################################################
#                               Diffusion Volumes                              #
################################################################################
const global sumV = [
                       2.5140e+01   # CH4 [-]
                       1.8010e+01   # CO  [-]
                       2.6900e+01   # CO2 [-]
                       6.1200e+00   # H2  [-]
                       1.3100e+01   # H2O [-]
                       1.8500e+01   # N2  [-]
                      ]

