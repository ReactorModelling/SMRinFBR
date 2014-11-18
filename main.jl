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

const global SMALL = 1000eps(Float64)
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

# Initial guess
uz      = uzIn*ones(Nglob)
p       = pIn*ones(Nglob)
T       = Tin*ones(Nglob)
#=
wCH4    = [wCH4in*ones(Nr); 0.048ones(Nglob - Nr)]
wCO     = [wCOin*ones(Nr); 0.1479ones(Nglob - Nr)]
wCO2    = [wCO2in*ones(Nr); 0.1804ones(Nglob - Nr)]
wH2     = [wH2in*ones(Nr); 0.0643ones(Nglob - Nr)]
wH2O    = [wH2Oin*ones(Nr); 0.4954ones(Nglob - Nr)]
=#

wN2     = wN2in*ones(Nglob)
wCH4    = wCH4in*ones(Nglob)
wCO     = wCOin*ones(Nglob) 
wCO2    = wCO2in*ones(Nglob) 
wH2     = wH2in*ones(Nglob) 
wH2O    = wH2Oin*ones(Nglob) 

w   = [wCH4 wCO wCO2 wH2 wH2O wN2]          # Matrix with all the mass fractions
w ./= sum(w,2)
x   = getMolarFractions(w)                 # Matrix with all the molar fractions
M   = getAvgMolarMass(x)                      # Average molar mass [kg mol^{-1}]
rho = M.*p./(R*T)
mu  = getViscosity(T,x)                                       # Viscosity [Pa s]
Re  = getReynolds(rho, uz, mu)                                 # Reynolds number
f   = getFrictionFactor(Re)                                    # Friction factor
cp  = getHeatCapacity(T,x)                    # Heat capacity [J K^{-1} kg^{-1}]
dH, reaction = getReaction(T,x,p)      # Enthalpy of reaction and reaction rates
                                       # [J kg^{-1} s^{-1}] [mol kg^{-1} s^{-1}]
lambdaEff, U = getHeatCoefficients(Re, T, x, mu, cp, M)

D = getDiffusivity(uz)

A_uz = zeros(Nglob,Nglob)
b_uz = zeros(Nglob)
continuityEquation(uz, rho, A_uz, b_uz)

A_p = zeros(Nz,Nz)
b_p = zeros(Nz)
ergunEquation(p, rho, uz, f, A_p, b_p)

A_T = zeros(Nglob, Nglob)
b_T = zeros(Nglob)
energyEquation(T, rho, uz, cp, dH, U, lambdaEff, A_T, b_T)

A_w = {zeros(Nglob,Nglob) for i in CompIndex}
b_w = {zeros(Nglob) for i in CompIndex}
speciesMassBalance(w, rho, uz, reaction, D, A_w, b_w)

gamma_w = 5e-3
gamma_T = 5e-3
Gamma_T = 5e-1
Gamma_w = 5e-1
totIter = 1
maxIter = 100000
totRes  = 1.0
while totRes > 1e-2

    res_p = norm(A_p*p[1:Nr:end] - b_p)/mean(p)
    iter = 1
    while res_p > 1e-6 && iter < maxIter
        p = kron(A_p\b_p, ones(Nr))
        ergunEquation(p, rho, uz, f, b_p)
        res_p = norm(A_p*p[1:Nr:end]-b_p)/abs(mean(p))
        iter += 1
    end
    println("p iterations: $iter")
    println("p residual : $res_p")

    res_uz = norm(A_uz*uz - b_uz)/mean(uz)
    iter = 1
    while res_uz > 1e-6 && iter < maxIter
        uz = A_uz\b_uz
        continuityEquation(uz, rho, A_uz)
        res_uz = norm(A_uz*uz - b_uz)/abs(mean(uz))
        iter += 1
    end
    println("uz iterations: $iter")
    println("uz residual : $res_uz")
    
    Told = copy(T)
    res_T = norm(A_T*T - b_T)/mean(T)
    
    iter = 1
    while res_T > 1e-6 && iter < maxIter
        T = gamma_T*(A_T\b_T) + (1 - gamma_T)*T
        energyEquation(T, rho, uz, cp, dH, U, lambdaEff, A_T, b_T)
        res_T = norm(A_T*T - b_T)/abs(mean(T))
        iter += 1
    end

    T = Gamma_T*T + (1 - Gamma_T)*Told
    T = max(1.0, T)
    println("T iterations: $iter")
    println("T residual : $res_T")
    println("Min T: $(minimum(T))")

    

    for i = 1:length(CompIndex)
        c = CompIndex[i]
        wOld = copy(w[:,c])
        res_w = norm(A_w[i]*w[:,c] - b_w[i])
        iter = 1
        while res_w > 1e-6 && iter < maxIter
            w[:,c] = gamma_w*(A_w[i]\b_w[i]) + (1 - gamma_w)*w[:,c]
            speciesMassBalance(w, rho, uz, reaction, D, A_w, b_w)
            res_w = norm(A_w[i]*w[:,c] - b_w[i])/abs(mean(w[:,c]))
            iter += 1
        end
        println("$(Comp[CompIndex[i]]) iterations: $iter")
        println("$(Comp[CompIndex[i]]) residual : $res_w")
        w[:,c] = Gamma_w*w[:,c] + (1-Gamma_w)*wOld
        w[:,c] = max(SMALL, min(1.0, w[:,c]))
    end
    w ./= sum(w,2)
    w[:,2] = 1 - sum(w[:,CompIndex],2)
    w[:,2] = max(SMALL, min(1.0, w[:,2]))
    w ./= sum(w,2)

    x   = getMolarFractions(w)                 # Matrix with all the molar fractions
    M   = getAvgMolarMass(x)                      # Average molar mass [kg mol^{-1}]
    rho = M.*p./(R*T)
    mu  = getViscosity(T,x)                                       # Viscosity [Pa s]
    Re  = getReynolds(rho, uz, mu)                                 # Reynolds number
    f   = getFrictionFactor(Re)                                    # Friction factor
    cp  = getHeatCapacity(T,x)                    # Heat capacity [J K^{-1} kg^{-1}]
    dH, reaction = getReaction(T,x,p)      # Enthalpy of reaction and reaction rates
                                           # [J kg^{-1} s^{-1}] [mol kg^{-1} s^{-1}]
    lambdaEff, U = getHeatCoefficients(Re, T, x, mu, cp, M)
    D = getDiffusivity(uz)

    ergunEquation(p, rho, uz, f, b_p)
    continuityEquation(uz, rho, A_uz)
    energyEquation(T, rho, uz, cp, dH, U, lambdaEff, A_T, b_T)
    speciesMassBalance(w, rho, uz, reaction, D, A_w, b_w)

    res_p   = norm(A_p*p[1:Nr:end] - b_p)/mean(p)
    res_uz  = norm(A_uz*uz - b_uz)/mean(uz)
    res_T   = norm(A_T*T - b_T)/mean(T)
    res_w   = sum([norm(A_w[i]*w[:,CompIndex[i]] - b_w[i]) for i in 1:length(CompIndex)])
    totRes  = res_p + res_uz + res_T + res_w
    println("p residual: $res_p")
    println("uz residual $res_uz")
    println("T residual $res_T")
    println("w residual $res_w")
    println("Total resisdual: $totRes")
    println("Outer loop iterations: $totIter")
    totIter += 1
end
