include("colloc.jl")
include("getFrictionFactor.jl")
include("getReynolds.jl")
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

converged = false
mu        = 3e-5
Re        = getReynolds(rho, uz, mu)
f         = getFrictionFactor(Re)
cp        = getHeatCapacity(Tin, )

# Define A matrices and b vectors
A_p      = [1 zeros(1,Nz-1); a[2:end,:]]
b_p      = [pIn; -1/dInner*(f.*rho.*uz.^2)[2:end,:]]
A_uz     = [1 zeros(1,Nz-1); (rho.*a + (a*rho).*I)[2:end,:]]
b_uz     = [uzIn; zeros(Nz-1)]
A_T      = [1 zeros(1,Nz-1); (rho.*cp.*a + 4*U/dInner)[2:end,:]]
b_T      = [Tin; 4*U/dInner*Ta*ones(Nz-1,1)]

gamma_p = 1
gamma_uz = 1
gamma_

iter = 1

while (!converged && (iter <= 100))
    println("Iteration number: $iter")
    # Solve ergun's equation for pressure
    p        = (1-gamma_p)*p + gamma_p*(A_p\b_p)

    # Solve the ideal gas law for density
    rho = 0.018*p./(R*T)

    # solve the continuity for velocity
    uz = (1-gamma_uz)*uz + gamma_uz*(A_uz\b_uz)

    # Solve for temperature
    T = (1-gamma_T)*T + gamma_T*(A_T\b_T)

    # Update matrices and vectors
    Re       = getReynolds(rho, uz, mu)
    f        = getFrictionFactor(Re)
    b_p      = [pIn; -1/dInner*(f.*rho.*uz.^2)[2:end,:]]
    A_uz     = [1 zeros(1,Nz-1); (rho.*a + (a*rho).*I)[2:end,:]]
    A_T      = [1 zeros(1,Nz-1); (rho.*cp.*a + 4*U/dInner)[2:end,:]]

    # Calculate residuals
    residual_p   = sqrt((A_p*p - b_p)'*(A_p*p - b_p))[1]/pIn
    residual_uz  = sqrt((A_uz*uz - b_uz)'*(A_uz*uz - b_uz))[1]/uzIn
    residual_A

    println(residual_p)
    println(residual_uz)

    converged = (residual_p < 1e-10) && (residual_uz < 1e-10)

    iter += 1

end
