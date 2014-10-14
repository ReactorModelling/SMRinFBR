################################################################################
#                                   Constants                                  #
################################################################################
const global Nz      = 20   # Number of collocation points in z direction
const global Ncomp   = 6    # Number of chemical components
const z,a,b,q        = colloc(Nz-2,1,1) # Collocation points and matrices
global z
global a
global b
global q
const global I       = eye(Nz) # Identity matrix
const global R       = 8.3145 # Gas constant [J/K mol]

################################################################################
#                              Reactor parameters                              #
################################################################################
const global void     = 0.528 # Void fraction
const global dInner   = 0.102 # Inner tube diameter [m]
const global Ta       = 1100 # Ambient temperature [K]
const global U        = 56.783 # Heat transfer coefficient [J/K m2 s]
const global rhoCat   = 2355.2;  # Catalyst density          [kgcat/m^3]
const global dParticle     = 0.0173; # Particle diameter         [m]
const global lambdaSt  = 52 # Heat coef. for tube metal[W/mK]
const global rInner = dInner/2 # Inner tube radius [m]
const global rOuter = 0.066 # Outer tube radius [m]



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

################################################################################
#                              Enthalpy data                                   #
################################################################################
#=
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

cpCoeff
CP(1,1) = 1.925E4;        % 1. coefficient for CH4   [J/kmoleK]
CP(1,2) = 5.213E1;        % 2. coefficient for CH4   [J/kmoleK^2]
CP(1,3) = 1.197E-2;       % 3. coefficient for CH4   [J/kmoleK^3]
CP(1,4) =-1.132E-5;       % 4. coefficient for CH4   [J/kmoleK^4]

CP(2,1) = 3.087E4;        % 1. coefficient for CO    [J/kmoleK]
CP(2,2) =-1.285E1;        % 2. coefficient for CO    [J/kmoleK^2]
CP(2,3) = 2.789E-2;       % 3. coefficient for CO    [J/kmoleK^3]
CP(2,4) =-1.272E-5;       % 4. coefficient for CO    [J/kmoleK^4]

CP(3,1) = 1.980E4;        % 1. coefficient for CO2   [J/kmoleK]
CP(3,2) = 7.344E1;        % 2. coefficient for CO2   [J/kmoleK^2]
CP(3,3) =-5.602E-2;       % 3. coefficient for CO2   [J/kmoleK^3]
CP(3,4) = 1.715E-5;       % 4. coefficient for CO2   [J/kmoleK^4]

CP(4,1) = 2.714E4;        % 1. coefficient for H2    [J/kmoleK]
CP(4,2) = 0.9274E1;       % 2. coefficient for H2    [J/kmoleK^2]
CP(4,3) =-1.381E-2;       % 3. coefficient for H2    [J/kmoleK^3]
CP(4,4) = 0.7645E-5;      % 4. coefficient for H2    [J/kmoleK^4]

CP(5,1) = 3.224E4;        % 1. coefficient for H2O   [J/kmoleK]
CP(5,2) = 0.1924E1;       % 2. coefficient for H2O   [J/kmoleK^2]
CP(5,3) = 1.055E-2;       % 3. coefficient for H2O   [J/kmoleK^3]
CP(5,4) = 0.3596E-5;      % 4. coefficient for H2O   [J/kmoleK^4]

CP(6,1) = 3.115E4;        % 1. coefficient for N2    [J/kmoleK]
CP(6,2) =-1.357E1;        % 2. coefficient for N2    [J/kmoleK^2]
CP(6,3) = 2.680E-2;       % 3. coefficient for N2    [J/kmoleK^3]
CP(6,4) =-1.168E-5;       % 4. coefficient for N2    [J/kmoleK^4]