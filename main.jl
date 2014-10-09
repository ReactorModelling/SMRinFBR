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

const global Nz      = 20
const global Ncomp   = 6
const z,a,b,q        = colloc(Nz-2,1,1)
const global I       = eye(Nz)
const global R       = 8.3145
const global dInner  = 1e-4
const global void     = 0.528
const global dInner = 0.102
const pIn     = 2.9e6
const Tin     = 793
const uzIn    = 1.89
const rhoIn   = 0.018*pIn./(R*Tin)
const wCH4in  = 0.1911
const wCOin   = 0.0001
const wCO2in  = 0.0200
const wH2in   = 0.0029
const wH2Oin  = 0.7218
const wN2in   = 0.0641

# Initial guess
uz  = uzIn*ones(Nz)
rho = rhoIn*ones(Nz)
p   = pIn*ones(Nz)
T   = Tin*ones(Nz)
wCH4  = 0.1911
wCO   = 0.0001
wCO2  = 0.0200
wH2     = 0.0029
wH2O    = 0.7218
wN2     = 0.0641

converged = false
mu       = 3e-5;
Re       = getReynolds(rho, uz, mu)
f        = getFrictionFactor(Re)

# Define A matrices and b vectors
A_p      = [1 zeros(1,Nz-1); a[2:end,:]]
b_p      = [pIn; -1/dInner*(f.*rho.*uz.^2)[2:end,:]]
A_uz     = [1 zeros(1,Nz-1); (rho.*a + (a*rho).*I)[2:end,:]]
b_uz     = [uzIn; zeros(Nz-1)]
A_T      = 

gamma_p = 1
gamma_uz = 1

iter = 1

while (!converged && (iter <= 100))
    println(iter)
    # Solve ergun's equation for pressure
    p        = (1-gamma_p)*p + gamma_p*(A_p\b_p)

    # Solve the ideal gas law for density
    rho = 0.018*p./(R*T)

    # solve the continuity for velocity
    uz = (1-gamma_uz)*uz + gamma_uz*(A_uz\b_uz)

    # Update matrices and vectors
    Re       = getReynolds(rho, uz, mu)
    f        = getFrictionFactor(Re)
    b_p      = [pIn; -1/dInner*(f.*rho.*uz.^2)[2:end,:]]
    A_uz     = [1 zeros(1,Nz-1); (rho.*a + (a*rho).*I)[2:end,:]]

    # Calculate residuals
    residual_p   = sqrt((A_p*p - b_p)'*(A_p*p - b_p))[1]
    residual_uz  = sqrt((A_uz*uz - b_uz)'*(A_uz*uz - b_uz))[1]

    println(residual_p)
    println(residual_uz)

    converged = (residual_p < 1e-6) && (residual_uz < 1e-6)

    iter += 1

end
