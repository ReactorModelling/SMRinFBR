include("colloc.jl")
#=
Main script for simulating the steam methane reforming in a fixed bed reactor
using the method of orthogonal collocation.

Inlet values:
    p   = 2.9e6 Pa
    T   = 793 K    
    uz  = 1.89 m/s
=# 

const global Nz      = 20
const z,a,b,q        = colloc(Nz-2,1,1)
const global I       = eye(Nz)
const global R       = 8.3145
const global dInner  = 1e-4
const pIn     = 2.9e6
const Tin     = 793
const uzIn    = 1.89
const rhoIn   = 0.018*pIn./(R*Tin)

# Initial guess
uz  = uzIn*ones(Nz)
rho = rhoIn*ones(Nz)
p   = pIn*ones(Nz)
T   = Tin*ones(Nz)

converged = false

# Define A matrices and b vectors
A_p      = [1 zeros(1,Nz-1); a[2:end,:]]
b_p      = [pIn; -6*(rho.*uz.^2)[2:end,:]]
A_rho    = [1 zeros(1,Nz-1); a[2:end,:]]
b_rho    = [rhoIn; 0.018/R * (1./T.*(a*p) - p./T.^2.*(a*T))[2:end,:]]
A_uz     = [1 zeros(1,Nz-1); (rho.*a + (a*rho).*I)[2:end,:]]
b_uz     = [uzIn; zeros(Nz-1)]

iter = 1

while (!converged && (iter < 100))
    println(iter)
    # Solve ergun's equation for pressure
    p        = A_p\b_p

    # Solve the ideal gas law for density
    rho = A_rho\b_rho

    # solve the continuity for velocity
    u_z = A_uz\b_uz

    # Update matrices and vectors
    b_p = [pIn; -6*(rho.*uz.^2)[2:end,:]]
    b_rho    = [rhoIn; 0.018/R * (1./T.*(a*p) - p./T.^2.*(a*T))[2:end,:]]
    A_uz     = [1 zeros(1,Nz-1); (rho.*a + (a*rho).*I)[2:end,:]]

    # Calculate residuals
    residual_p   = sqrt((A_p*p - b_p)'*(A_p*p - b_p))[1]
    residual_rho = sqrt((A_rho*rho - b_rho)'*(A_rho*rho - b_rho))[1]
    residual_uz  = sqrt((A_uz*uz - b_uz)'*(A_uz*uz - b_uz))[1]

    println(residual_p)
    println(residual_uz)
    println(residual_rho)

    converged = (residual_p < 1e-5) && (residual_rho < 1e-5) && (residual_uz < 1e-5)

    iter += 1

end
