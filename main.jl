workspace()
include("colloc.jl")
include("getFrictionFactor.jl")
include("getReynolds.jl")
include("getMolarFractions.jl")
include("getAvgMolarMass.jl")
include("getViscosity.jl")
include("getHeatCapacity.jl")
include("getReaction.jl")
include("continuityEquation.jl")
include("crossSectionalAverage.jl")
include("ergunEquation.jl")
#=
Main script for simulating the steam methane reforming in a fixed bed reactor
using the method of orthogonal collocation.

Inlet values:
    p       = 2.9e6 Pa
    T       = 793 K    
    uz      = 1.89 m/s
    wCH4    = 0.1911
    wCO     = 0.0001
    wCO2    = 0.0200
    wH2     = 0.0029
    wH2O    = 0.7218
    wN2     = 0.0641
=# 
const global SMALL = 10eps(Float64)
include("constants.jl")

################################################################################
#                              Inlet conditions                                #
################################################################################
const global pIn     = 2.9e6 # Inlet pressure[Pa]
const global Tin     = 793 # Inlet temperature [K]
const global uzIn    = 1.89 # Inlet velocity [m/s]
const global wCH4in  = 0.1911 # Inlet mass fraction of CH4
const global wCOin   = 0.0001 # Inlet mass fraction of CO
const global wCO2in  = 0.0200 # Inlet mass fraction of CO2
const global wH2in   = 0.0029 # Inlet mass fraction of H2
const global wH2Oin  = 0.7218 # Inlet mass fraction of H2O
const global wN2in   = 0.0641 # Inlet mass fraction of N2


# Initial guess
uz      = uzIn*ones(Nglob)
p       = pIn*ones(Nglob)
T       = [Tin*ones(Nr); 500.0*ones(Nglob-Nr)] #Tin*ones(Nglob)#
wCH4    = wCH4in*ones(Nglob)
wCO     = wCOin*ones(Nglob)
wCO2    = wCO2in*ones(Nglob)
wH2     = wH2in*ones(Nglob)
wH2O    = wH2Oin*ones(Nglob)
wN2     = wN2in*ones(Nglob)

w   = [wCH4 wCO wCO2 wH2 wH2O wN2]          # Matrix with all the mass fractions
x   = getMolarFractions(w)                 # Matrix with all the molar fractions
M   = getAvgMolarMass(x)                      # Average molar mass [kg mol^{-1}]
rho = M.*p./(R*T)
mu  = getViscosity(T,x)                                       # Viscosity [Pa s]
Re  = getReynolds(rho, uz, mu)                                 # Reynolds number
f   = getFrictionFactor(Re)                                    # Friction factor
cp  = getHeatCapacity(T,x)                    # Heat capacity [J K^{-1} kg^{-1}]
dH, reaction = getReaction(T,x,p)      # Enthalpy of reaction and reaction rates
                                       # [J kg^{-1} s^{-1}] [mol kg^{-1} s^{-1}]

A_uz = zeros(Nglob,Nglob)
b_uz = zeros(Nglob)
continuityEquation(uz, rho, A_uz, b_uz)
A_p = zeros(Nz,Nz)
b_p = zeros(Nz)
ergunEquation(p, rho, uz, f, A_p, b_p)

for i = 1:100
    p = kron(A_p\b_p, ones(Nr))
    rho = M.*p./(R*T)
    uz = A_uz\b_uz

    continuityEquation(uz, rho, A_uz)
    ergunEquation(p, rho, uz, f, b_p)



end