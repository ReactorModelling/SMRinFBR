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
include("getHeatCoefficients.jl")
include("energyEquation.jl")
include("speciesMassBalance.jl")
include("getDiffusivity.jl")
include("getDerivative.jl")
include("couple.jl")
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
const global wIn     = [wCH4in, wCOin, wCO2in, wH2in, wH2Oin, wN2in]

# Initial guess: inlet values in the whol reactor
uz      = uzIn*ones(Nglob)
p       = pIn*ones(Nglob)
T       = Tin*ones(Nglob)
wN2     = wN2in*ones(Nglob)
wCH4    = wCH4in*ones(Nglob)
wCO     = wCOin*ones(Nglob) 
wCO2    = wCO2in*ones(Nglob) 
wH2     = wH2in*ones(Nglob) 
wH2O    = wH2Oin*ones(Nglob) 

w    = [wCH4 wCO wCO2 wH2 wH2O wN2]         # Matrix with all the mass fractions
wVec = vec(w[:,[CompIndex,2]])            # A vector with all the mass fractions

x   = getMolarFractions(w)                 # Matrix with all the molar fractions
M   = getAvgMolarMass(x)                      # Average molar mass [kg mol^{-1}]
rho = M.*p./(R*T)                                          # Density [kg m^{-3}]
mu  = getViscosity(T,x)                                       # Viscosity [Pa s]
Re  = getReynolds(rho, uz, mu)                                 # Reynolds number
f   = getFrictionFactor(Re)                                    # Friction factor
cp  = getHeatCapacity(T,x)                    # Heat capacity [J K^{-1} kg^{-1}]
D   = getDiffusivity(uz)                              # Diffusivity [m^2 s^{-1}]
dH, reaction = getReaction(T,x,p)      # Enthalpy of reaction and reaction rates
                                       # [J kg^{-1} s^{-1}] [mol kg^{-1} s^{-1}]
lambdaEff, U = getHeatCoefficients(Re, T, x, mu, cp, M) # Effective conductivity
                                                        # [W m^{-1} K^{-1}]
                                                        # Heat coefficient
                                                        # [W m^{-2} K^{-1}]
# Initialize A matrix and b vector for uz
A_uz = zeros(Nglob,Nglob)
b_uz = zeros(Nglob)
# Fill in values (done by reference)
continuityEquation!(uz, rho, A_uz, b_uz)

# Initialize A matrix and b vector for p
# Note that dim(A_p) = Nz x Nz (no radial variation)
A_p = zeros(Nz,Nz)
b_p = zeros(Nz)
# Fill in values (done by reference)
ergunEquation!(p, rho, uz, f, A_p, b_p)

# Initialize A matrix and b vector for T
A_T = zeros(Nglob, Nglob)
b_T = zeros(Nglob)
# Fill in values (done by reference)
energyEquation!(T, rho, uz, cp, dH, U, lambdaEff, A_T, b_T)

# Initialize an array with all the matrices for the mass fractions to be solved
A_w = {zeros(Nglob,Nglob) for i in CompIndex}
# Initialize an array with all the vectors for the mass fractions to be solved
b_w = {zeros(Nglob) for i in CompIndex}
# Fill in values in all the vectors and matrices in the arrays (by reference)
speciesMassBalance!(w, rho, uz, reaction, D, A_w, b_w)

# Initialize total A (sparse)
A = zeros(10Nglob,10Nglob)
# Initialize total B
B = zeros(10Nglob,1)
# Fill A and B (by reference)
couple!(w, T, uz, rho, p, 
        f, cp, dH, U, lambdaEff, reaction, D, 
        A_w, b_w, A_T, b_T, A_uz, b_uz, A_p, b_p, 
        A, B)
# Initialize total variable vector
variables = [wVec, T, uz, rho, p]

# Under-relaxation factor
const gamma_w     = 5e-2                                       #  Mass fractions
const gamma_T     = 5e-2                                           # Temperature

const maxIter   = 1000                             # Max iterations in the loops
totIter         = 1                                           # Total iterations
totRes          = 1.0                           #  Initialize the total residual

while totRes > 1e-1
    # Solve for all the variables
    variables = A\B
    # Extract the variables
    wVec = gamma_w*variables[1:6Nglob] + (1-gamma_w)*wVec
    w    = reshape(wVec,Nglob,Ncomp)[:,revertIndex]
    T    = gamma_T*variables[6Nglob+1:7Nglob] + (1-gamma_T)*T
    uz   = variables[7Nglob+1:8Nglob]
    rho  = variables[8Nglob+1:9Nglob]
    p    = variables[9Nglob+1:end]

    # Update dependent variables
    x   = getMolarFractions(w)             # Matrix with all the molar fractions
    M   = getAvgMolarMass(x)                  # Average molar mass [kg mol^{-1}]
    rho = M.*p./(R*T)                                      # Density [kg m^{-3}]
    mu  = getViscosity(T,x)                                   # Viscosity [Pa s]
    Re  = getReynolds(rho, uz, mu)                             # Reynolds number
    f   = getFrictionFactor(Re)                                # Friction factor
    cp  = getHeatCapacity(T,x)                # Heat capacity [J K^{-1} kg^{-1}]
    D   = getDiffusivity(uz)                          # Diffusivity [m^2 s^{-1}]
    dH, reaction = getReaction(T,x,p)  # Enthalpy of reaction and reaction rates
                                       # [J kg^{-1} s^{-1}] [mol kg^{-1} s^{-1}]
    lambdaEff, U = getHeatCoefficients(Re, T, x, mu, cp, M)       # Conductivity
                                                             # [W m^{-1} K^{-1}]
                                                              # Heat coefficient
                                                             # [W m^{-2} K^{-1}]

    # Update matrix and vector
    couple!(w, T, uz, rho, p,
            f, cp, dH, U, lambdaEff, reaction, D, 
            A_w, b_w, A_T, b_T, A_uz, b_uz, A_p, b_p, 
            A, B)
    # Calculate the total residual
    totRes  = norm(A*variables - B)
    # Display information
    println("Total residual: $totRes")
    println("Outer loop iterations: $totIter")
    # Update iteration number
    totIter += 1
end

################################################################################
#        Convert solution vectors to matrices and save the result              #
################################################################################

T   = reshape(T, Nr, Nz)'
rho = reshape(rho, Nr, Nz)'
uz  = reshape(uz, Nr, Nz)'
p   = reshape(p, Nr, Nz)'
w   = {reshape(w[:,i],Nr,Nz)' for i in 1:Ncomp}
x   = {reshape(x[:,i],Nr,Nz)' for i in 1:Ncomp}

using HDF5, MAT

c = matopen("coupled.mat", "w") do file
    write(file, "r", r)
    write(file, "z", Z)
    write(file, "T", T)
    write(file, "rho", rho)
    write(file, "x", x)
    write(file, "w", w)
    write(file, "uz", uz)
    write(file, "p", p)
end
