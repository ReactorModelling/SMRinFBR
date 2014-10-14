include("colloc.jl")
include("getFrictionFactor.jl")
include("getReynolds.jl")
include("getMolarFractions.jl")
include("getAvgMolarMass.jl")
include("getViscosity.jl")
include("getHeatCapacity.jl")
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
const pIn     = 2.9e6 # Inlet pressure[Pa]
const Tin     = 793 # Inlet temperature [K]
const uzIn    = 1.89 # Inlet velocity [m/s]
const wCH4in  = 0.1911 # Inlet mass fraction of CH4
const wCOin   = 0.0001 # Inlet mass fraction of CO
const wCO2in  = 0.0200 # Inlet mass fraction of CO2
const wH2in   = 0.0029 # Inlet mass fraction of H2
const wH2Oin  = 0.7218 # Inlet mass fraction of H2O
const wN2in   = 0.0641 # Inlet mass fraction of N2

const rhoIn   = 0.018*pIn./(R*Tin) #Inlet density [kg/m3]



# Initial guess
uz  = uzIn*ones(Nz)
rho = rhoIn*ones(Nz)
p   = pIn*ones(Nz)
T   = Tin*ones(Nz)
wCH4  = 0.1911*ones(Nz)
wCO   = 0.0001*ones(Nz)
wCO2  = 0.0200*ones(Nz)
wH2   = 0.0029*ones(Nz)
wH2O  = 0.7218*ones(Nz)
wN2   = 0.0641*ones(Nz)
w     = [wCH4 wCO wCO2 wH2 wH2O wN2]
x         = getMolarFractions(w)
mu        = getViscosity(T,x)
Re        = getReynolds(rho, uz, mu)
f         = getFrictionFactor(Re)
cp        = getHeatCapacity(T,x)

# Define A matrices and b vectors
A_p      = [1 zeros(1,Nz-1); A[2:end,:]]
b_p      = [pIn; -1/dInner*(f.*rho.*uz.^2)[2:end,:]]
A_uz     = [1 zeros(1,Nz-1); (rho.*A + (A*rho).*I)[2:end,:]]
b_uz     = [uzIn; zeros(Nz-1)]
A_T      = [1 zeros(1,Nz-1); (rho.*cp.*uz.*A + 4*U/dInner)[2:end,:]]
b_T      = [Tin; 4*U/dInner*Ta*ones(Nz-1,1)]

gamma_p = 0.5
gamma_uz = 0.5
gamma_T = 0.5

converged = false
iter = 1

while (!converged && (iter <= 100000))
    println("Iteration number: $iter")
    # Solve ergun's equation for pressure
    p        = (1-gamma_p)*p + gamma_p*(A_p\b_p)

    # Solve the ideal gas law for density
    rho = 0.018*p./(R*T)

    # solve the continuity for velocity
    uz = (1-gamma_uz)*uz + gamma_uz*(A_uz\b_uz)

    # Solve for temperature
    T = (1-gamma_T)*T + gamma_T*(A_T\b_T)

    # Update parameters
    mu        = getViscosity(T,x)
    Re        = getReynolds(rho, uz, mu)
    f         = getFrictionFactor(Re)
    cp        = getHeatCapacity(T,x)

    # Update matrices and vectors
    b_p      = [pIn; -1/dInner*(f.*rho.*uz.^2)[2:end,:]]
    A_uz     = [1 zeros(1,Nz-1); (rho.*A + (A*rho).*I)[2:end,:]]
    A_T      = [1 zeros(1,Nz-1); (rho.*cp.*uz.*A + 4*U/dInner.*I)[2:end,:]]

    # Calculate residuals
    residual_p   = sqrt((A_p*p - b_p)'*(A_p*p - b_p))[1]/pIn
    residual_uz  = sqrt((A_uz*uz - b_uz)'*(A_uz*uz - b_uz))[1]/uzIn
    residual_T   = sqrt((A_T*T-b_T)'*(A_T*T-b_T))[1]/Tin

    println("Residual of p: $residual_p")
    println("Residual of uz: $residual_uz")
    println("Residual of T: $residual_T")


    converged = (residual_p < 1e-10) && (residual_uz < 1e-10) && (residual_T < 1e-10)

    iter += 1

end