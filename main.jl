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
include("coupleMassTemperature.jl")
include("coupleVelocityDensityPressure.jl")
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

# Initial guess: inlet values in the whole reactor
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
wVec = vec(w[:,[CompIndex,2]])

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
# Initialize A matrix and b vector for A*uz = b
A_uz = zeros(Nglob,Nglob)
b_uz = zeros(Nglob)
# Fill in values (done by reference)
continuityEquation!(uz, rho, A_uz, b_uz)

# Initialize A matrix and b vector for A*p = b
# Note that dim(A_p) = Nz x Nz (no radial variation)
A_p = zeros(Nz,Nz)
b_p = zeros(Nz)
# Fill in values (done by reference)
ergunEquation!(p, rho, uz, f, A_p, b_p)

# Initialize A matrix and b vector for A*T = b
A_T = zeros(Nglob, Nglob)
b_T = zeros(Nglob)
# Fill in values (done by reference)
energyEquation!(T, rho, uz, cp, dH, U, lambdaEff, A_T, b_T)

# Initialize an array with all the matrices for the mass fractions to be solved
A_w = {zeros(Nglob,Nglob) for i in CompIndex}
# Initialize an array with all the vectors for the mass fractions to be solved
b_w = {zeros(Nglob) for i in CompIndex}
# Fill in values in all the vectors and matrices in the arrays
speciesMassBalance!(w, rho, uz, reaction, D, A_w, b_w)


A_uzRhoP = zeros(3Nglob,3Nglob)
b_uzRhoP = zeros(3Nglob)
coupleVelocityDensityPressure!(uz, rho, p,
                               M, T, f, 
                               A_uz, b_uz, A_p, b_p, 
                               A_uzRhoP, b_uzRhoP)
uzRhoP = [uz, rho, p]

A_wT = zeros(7Nglob,7Nglob)
b_wT = zeros(7Nglob)
coupleMassTemperature!(w, T, 
                       rho, uz, cp, dH, U, lambdaEff, reaction, D,
                       A_w, b_w, A_T, b_T,
                       A_wT, b_wT)
wT = [wVec, T]

# Under-relaxation factor
const gamma_wT     = 5e-2                       # Mass fractions and temperature

const maxIter   = 1000                             # Max iterations in the loops
totIter         = 1                                           # Total iterations
totRes          = 1.0                            # Initialize the total residual

while totRes > 1e-2
    ############################################################################
    #                           T-w iteration loop                             #
    ############################################################################

    # Calculate the residuals
    res_wT = norm(A_wT*wT - b_wT)
    iter = 0                                      # Initialize iteration numbers
    while res_wT > 1e-2 && iter < maxIter
        # Solve for w and T and apply under-relaxation
        wT = gamma_wT*(A_wT\b_wT) + (1 - gamma_wT)*wT

        wVec = wT[1:6Nglob]
        T    = wT[6Nglob+1:end]
        w    = reshape(wVec,Nglob,Ncomp)[:,revertIndex]

        # Update dependent variables
        x   = getMolarFractions(w)         # Matrix with all the molar fractions
        cp  = getHeatCapacity(T,x)            # Heat capacity [J K^{-1} kg^{-1}]
        dH, reaction = getReaction(T,x,p) # Reaction enthalpy and reaction rates
                                       # [J kg^{-1} s^{-1}] [mol kg^{-1} s^{-1}]

        # Update the matrices and vectors
        coupleMassTemperature!(w, T,
                               rho, uz, cp, dH, U, lambdaEff, reaction, D,
                               A_w, b_w, A_T, b_T,
                               A_wT, b_wT)

        # Update the residual
        res_wT = norm(A_wT*wT - b_wT)
        # Update iteration number
        iter += 1
    end
    # Display information
    println("T-w iterations: $iter")
    println("T-w residual : $res_wT")
    println("Min T: $(minimum(T))")
    println("Min w: $(minimum(w))")

    M   = getAvgMolarMass(x)                  # Average molar mass [kg mol^{-1}]
    mu  = getViscosity(T,x)                                   # Viscosity [Pa s]
    Re  = getReynolds(rho, uz, mu)                             # Reynolds number
    f   = getFrictionFactor(Re)                                # Friction factor


    ############################################################################
    #                           uz-rho-p iteration loop                        #
    ############################################################################
    
    coupleVelocityDensityPressure!(uz, rho, p, 
                                   M, T, f, 
                                   A_uz, b_uz, A_p, b_p, 
                                   A_uzRhoP, b_uzRhoP)

    res_uzRhoP = norm(A_uzRhoP*uzRhoP - b_uzRhoP)
    iter = 0
    while res_uzRhoP > 1e-7 && iter < maxIter
        # Solve for velocity, density and pressure  
        uzRhoP = A_uzRhoP\b_uzRhoP

        # Extract variables
        uz = uzRhoP[1:Nglob]
        rho = uzRhoP[Nglob+1:2Nglob]
        p = uzRhoP[2Nglob+1:end]

        # Update dependent variables
        Re  = getReynolds(rho, uz, mu)                         # Reynolds number
        f   = getFrictionFactor(Re)                            # Friction factor

        # Update vectors and matrices
        coupleVelocityDensityPressure!(uz, rho, p, 
                                       M, T, f, 
                                       A_uz, b_uz, A_p, b_p, 
                                       A_uzRhoP, b_uzRhoP)
        # Update the residuals
        res_uzRhoP = norm(A_uzRhoP*uzRhoP - b_uzRhoP)
        # Update iterations
        iter += 1
    end
    println("uz-rho-p residual : $res_uzRhoP")
    println("uz-rho-p iterations: $iter")

    dH, reaction = getReaction(T,x,p) # Reaction enthalpy and reaction rates
                                       # [J kg^{-1} s^{-1}] [mol kg^{-1} s^{-1}]
    lambdaEff, U = getHeatCoefficients(Re, T, x, mu, cp, M)
    D            = getDiffusivity(uz)

    # Update the matrices and vectors
    coupleMassTemperature!(w, T,
                           rho, uz, cp, dH, U, lambdaEff, reaction, D,
                           A_w, b_w, A_T, b_T,
                           A_wT, b_wT)
    coupleVelocityDensityPressure!(uz, rho, p, 
                                   M, T, f, 
                                   A_uz, b_uz, A_p, b_p, 
                                   A_uzRhoP, b_uzRhoP)
    # Update the residuals
    res_uzRhoP  = norm(A_uzRhoP*uzRhoP - b_uzRhoP)
    res_wT      = norm(A_wT*wT - b_wT)
    # Calculate the total residual
    totRes  = res_uzRhoP + res_wT
    # Display information
    println("uz-rho-p residual : $res_uzRhoP")
    println("w-T residual : $res_wT")
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

c = matopen("combined.mat", "w") do file
    write(file, "r", r)
    write(file, "z", Z)
    write(file, "T", T)
    write(file, "rho", rho)
    write(file, "x", x)
    write(file, "w", w)
    write(file, "uz", uz)
    write(file, "p", p)
end
