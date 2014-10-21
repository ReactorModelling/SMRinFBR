include("colloc.jl")
include("getFrictionFactor.jl")
include("getReynolds.jl")
include("getMolarFractions.jl")
include("getAvgMolarMass.jl")
include("getViscosity.jl")
include("getHeatCapacity.jl")
include("getReaction.jl")
using PyPlot
hold(true)
yscale("log")
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
wCH4  = wCH4in*ones(Nz)
wCO   = wCOin*ones(Nz)
wCO2  = wCO2in*ones(Nz)
wH2   = wH2in*ones(Nz)
wH2O  = wH2Oin*ones(Nz)
wN2   = wN2in*ones(Nz)

    # Calculate parameters
w   = [wCH4 wCO wCO2 wH2 wH2O wN2]          # Matrix with all the mass fractions
x   = getMolarFractions(w)                 # Matrix with all the molar fractions
mu  = getViscosity(T,x)                                       # Viscosity [Pa s]
Re  = getReynolds(rho, uz, mu)                                 # Reynolds number
f   = getFrictionFactor(Re)                                    # Friction factor
cp  = getHeatCapacity(T,x)                    # Heat capacity [J K^{-1} kg^{-1}]
M   = getAvgMolarMass(x)                      # Average molar mass [kg mol^{-1}]
dH, reaction = getReaction(T,x,p)      # Enthalpy of reaction and reaction rates
                                       # [J kg^{-1} s^{-1}] [mol kg^{-1} s^{-1}]

# Define A matrices and b vectors
A_p      = [1 zeros(1,Nz-1); A[2:end,:]]
b_p      = [pIn; -1/dInner*(f.*rho.*uz.^2)[2:end,:]]
A_uz     = [1 zeros(1,Nz-1); (rho.*A + (A*rho).*I)[2:end,:]]
b_uz     = [uzIn; zeros(Nz-1)]
A_T      = [1 zeros(1,Nz-1); (rho.*cp.*uz.*A + 4*U/dInner)[2:end,:]]
b_T      = [Tin; 4*U/dInner*Ta*ones(Nz-1,1)]
#A_wCH4   = [1 zeros(1,Nz-1); (rho.*uz.*A + uz.*(A*rho).*I + rho.*(A*uz).*I)[2:end,:]]

#b_wCH4   = [wCH4in; ((1-void)*rhoCat*molarMass[1]*(reaction*N[:,1]))[2:end,:]]

# Under-relaxation factors
gamma_p = 1
gamma_uz = 1
gamma_T = 0.50

converged = false
iter = 1

while (!converged && (iter <= 100000))
    println("Iteration number: $iter")

    # Solve ergun's equation for pressure
    p        = (1-gamma_p)*p + gamma_p*(A_p\b_p)

    # Solve the ideal gas law for density
    rho = M.*p./(R*T)

    # solve the continuity for velocity
    uz = (1-gamma_uz)*uz + gamma_uz*(A_uz\b_uz)

    # Solve for temperature
    T = (1-gamma_T)*T + gamma_T*(A_T\b_T)
    #println(size(A_wCH4))
    #wCH4 = (1-gamma_T)*wCH4 + gamma_T*(A_wCH4\b_wCH4)
    # Update parameters
    x         = getMolarFractions(w)
    mu        = getViscosity(T,x)
    Re        = getReynolds(rho, uz, mu)
    f         = getFrictionFactor(Re)
    cp        = getHeatCapacity(T,x)
    M         = getAvgMolarMass(x)
    dH, reaction = getReaction(T,x,p)


    # Update matrices and vectors
    b_p      = [pIn; -1/dInner*(f.*rho.*uz.^2)[2:end,:]]
    A_uz     = [1 zeros(1,Nz-1); (rho.*A + (A*rho).*I)[2:end,:]]
    A_T      = [1 zeros(1,Nz-1); (rho.*cp.*uz.*A + 4*U/dInner.*I)[2:end,:]]
    #A_wCH4   = [1 zeros(1,Nz-1); (rho.*uz.*A + uz.*(A*rho).*I + rho.*(A*uz).*I)[2:end,:]]
    #b_wCH4   = [wCH4in; ((1-void)*rhoCat*molarMass[1]*(reaction*N[:,1]))[2:end,:]]

    # Calculate residuals
    p_err       = vec(A_p*p - b_p)
    uz_err      = vec(A_uz*uz - b_uz)
    T_err       = vec(A_T*T - b_T)
    residual_p  = sqrt(dot(p_err, p_err))/mean(p)
    residual_uz = sqrt(dot(uz_err, uz_err))/mean(uz)
    residual_T  = sqrt(dot(T_err, T_err))/mean(T)
    #residual_wCH4= sqrt((A_wCH4*wCH4-b_wCH4)'*(A_wCH4*wCH4-b_wCH4))[1]/mean(wCH4)
    
    println("Residual of p:  $residual_p")
    println("Residual of uz: $residual_uz")
    println("Residual of T:  $residual_T")
    #println("Residual of CH4:  $residual_wCH4")

    if iter%2 == 0
        plot(iter,residual_T,"xk")
        plot(iter,residual_p,"xr")
        plot(iter,residual_uz,"xb")
        show()
    end

    converged = (residual_p + residual_uz + residual_T) < 1e-8

    iter += 1

end

figure(2)
plot(Z,T)
figure(3)
plot(Z,p)
figure(4)
plot(Z,uz)