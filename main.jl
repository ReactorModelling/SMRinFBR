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
const pIn     = 2.9e6
const Tin     = 793
const uzIn    = 1.89
const rhoIn   = 0.018*pIn./(R*Tin)

# Initial guess
uz  = uzIn*ones(Nz)
rho = rhoIn*ones(Nz)
p   = pIn*ones(Nz)
T   = Tin*ones(Nz)

# Solve ergun's equation for pressure
A_p      = a
A_p[1,:] = zeros(1,Nz)
A_p[1,1] = 1
b_p      = -100 * rho.*uz.^2
b_p[1]   = pIn
p        = A_p\b_p

# Solve the ideal gas law for density
A_rho = a
A_rho[1,:] = zeros(1,Nz)
A_rho[1,1] = 1
b_rho   = 0.018/R * (1./T.*(a*p) - p./T.^2.*(a*T))
b_rho[1] = rhoIn
rho = A_rho\b_rho

# solve the continuity for velocity

A_uz         = (rho.*a + (a*rho).*I)
A_uz[1,:]    = zeros(1,Nz)
A_uz[1,1]    = 1
b_uz         = zeros(Nz)
b_uz[1]      = uzIn

u_z          = A_uz\b_uz