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

w   = [wCH4 wCO wCO2 wH2 wH2O wN2]          # Matrix with all the mass fractions

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

# Under-relaxation factors
const gamma_w     = 5e-2                                        # Mass fractions
const gamma_T     = 5e-2                                           # Temperature
const gamma_uz    = 1e-0                                              # Velocity
const gamma_p     = 1e-0                                              # Pressure

const maxIter   = 500                            # Max iterations in the loops
totIter         = 1                                           # Total iterations
totRes          = 1.0                           # Initialize the total residual

while totRes > 1e-1
    ############################################################################
    #                           T-w iteration loop                             #
    ############################################################################

    # Calculate the residuals
    res_T = norm(A_T*T - b_T)                    # 2-norm of the residuals in T
    res_w = [
                norm(A_w[i]*w[:,CompIndex[i]]-b_w[i])
                for i in 1:length(CompIndex)
            ]                  # A vector with the 2-norms of the residuals in w
    iter = 0                                      # Initialize iteration numbers
    while (res_T > 1e-2 || maximum(res_w) > 1e-6) && iter < maxIter
        # Solve for T and apply under-relaxation
        T = gamma_T*(A_T\b_T) + (1 - gamma_T)*T

        for i = 1:length(CompIndex)            # Loop through the mass fractions
            c = CompIndex[i]                       # Extract the component index
            w[:,c] = gamma_w*(A_w[i]\b_w[i]) + (1-gamma_w)*w[:,c]        # Solve
        end
        # Solve for CO using 1 - sum of the other mass fractions
        w[:,2] = 1 - sum(w[:,CompIndex],2)

        # Update dependent variables
        x   = getMolarFractions(w)         # Matrix with all the molar fractions
        cp  = getHeatCapacity(T,x)            # Heat capacity [J K^{-1} kg^{-1}]
        dH, reaction = getReaction(T,x,p) # Reaction enthalpy and reaction rates
                                       # [J kg^{-1} s^{-1}] [mol kg^{-1} s^{-1}]

        # Update the matrices and vectors
        energyEquation!(T, rho, uz, cp, dH, U, lambdaEff, A_T, b_T)
        speciesMassBalance!(w, rho, uz, reaction, D, A_w, b_w)

        # Update the residuals
        res_T   = norm(A_T*T - b_T)
        res_w   = [
                    norm(A_w[i]*w[:,CompIndex[i]] - b_w[i])
                    for i in 1:length(CompIndex)
                  ]
        # Update iteration number
        iter += 1
    end
    # Display information
    println("T-w iterations: $iter")
    println("T residual : $res_T")
    println("w residual : $(maximum(res_w))")
    println("Min T: $(minimum(T))")
    println("Min w: $(minimum(w))")

    # Update dependent variables
    M   = getAvgMolarMass(x)                  # Average molar mass [kg mol^{-1}]
    rho = M.*p./(R*T)                                      # Density [kg m^{-3}]
    mu  = getViscosity(T,x)                                   # Viscosity [Pa s]
    Re  = getReynolds(rho, uz, mu)                             # Reynolds number
    f   = getFrictionFactor(Re)                                # Friction factor

    ############################################################################
    #                           uz-p iteration loop                            #
    ############################################################################

    # Update the matrices and vectors for uz and p
    ergunEquation!(p, rho, uz, f, b_p)
    continuityEquation!(uz, rho, A_uz)

    # Calculate the residuals
    res_uz = norm(A_uz*uz - b_uz)                    # 2-norm of residuals in uz
    res_p  = norm(A_p*p[1:Nr:end] - b_p)              # 2-norm of residuals in p
    iter   = 0                                     # Initialize iteration number
    while (res_uz > 1e-12 || res_p > 1e-7) && iter < maxIter
        # Solve for p and apply under-relaxation
        p  = gamma_p*kron(A_p\b_p, ones(Nr)) + (1-gamma_p)*p
        # Solve for uz and apply under-relaxation
        uz = gamma_uz*(A_uz\b_uz) + (1-gamma_uz)*uz

        # Update dependent variables
        rho = M.*p./(R*T)                                  # Density [kg m^{-3}]
        Re  = getReynolds(rho, uz, mu)                         # Reynolds number
        f   = getFrictionFactor(Re)                            # Friction factor

        # Update dependent variables
        ergunEquation!(p, rho, uz, f, b_p)
        continuityEquation!(uz, rho, A_uz)

        # Update the residuals
        res_uz = norm(A_uz*uz - b_uz)
        res_p  = norm(A_p*p[1:Nr:end] - b_p)
        # Update iteration number
        iter += 1
    end
    # Display information
    println("p-uz iterations: $iter")
    println("p residual : $res_p")
    println("uz residual : $res_uz")

    # Update dependent variables
    lambdaEff, U = getHeatCoefficients(Re, T, x, mu, cp, M)
    D            = getDiffusivity(uz)

    # Update the matrices and vectors
    ergunEquation!(p, rho, uz, f, b_p)
    continuityEquation!(uz, rho, A_uz)
    energyEquation!(T, rho, uz, cp, dH, U, lambdaEff, A_T, b_T)
    speciesMassBalance!(w, rho, uz, reaction, D, A_w, b_w)

    # Update the residuals
    res_p   = norm(A_p*p[1:Nr:end] - b_p)
    res_uz  = norm(A_uz*uz - b_uz)
    res_T   = norm(A_T*T - b_T)
    res_w   = [
                norm(A_w[i]*w[:,CompIndex[i]]-b_w[i])
                for i in 1:length(CompIndex)
              ]
    # Calculate the total residual
    totRes  = res_p + res_uz + res_T + sum(res_w)
    # Display information
    println("p residual : $res_p")
    println("uz residual : $res_uz")
    println("T residual : $res_T")
    println("w residual : $(maximum(res_w))")
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

c = matopen("data/pairwiseSegregated.mat", "w") do file
    write(file, "r", r)
    write(file, "z", Z)
    write(file, "T", T)
    write(file, "rho", rho)
    write(file, "x", x)
    write(file, "w", w)
    write(file, "uz", uz)
    write(file, "p", p)
end
