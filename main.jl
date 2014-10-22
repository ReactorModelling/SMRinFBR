include("colloc.jl")
include("getFrictionFactor.jl")
include("getReynolds.jl")
include("getMolarFractions.jl")
include("getAvgMolarMass.jl")
include("getViscosity.jl")
include("getHeatCapacity.jl")
include("getReaction.jl")

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
const global SMALL = 10eps(Float64)
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




# Initial guess
uz  = uzIn*ones(Nz)
p   = pIn*ones(Nz)
T   = [Tin; 500*ones(Nz-1)]
wCH4 = [wCH4in; wCH4in*ones(Nz-1)]
wCO = [wCOin; wCOin*ones(Nz-1)]
wCO2 = [wCO2in; wCO2in*ones(Nz-1)]
wH2 = [wH2in; wH2in*ones(Nz-1)]
wH2O = [wH2Oin; wH2Oin*ones(Nz-1)]
wN2 = [wN2in; wN2in*ones(Nz-1)]
#wCH4  = [wCH4in; zeros(Nz-1)]
#wCO   = [wCOin; zeros(Nz-1)]
#wCO2  = [wCO2in; zeros(Nz-1)]
#wH2   = [wH2in; zeros(Nz-1)]
#wH2O  = [wH2Oin; zeros(Nz-1)]
#wN2   = [wN2in; ones(Nz-1)]

    # Calculate parameters
w   = [wCH4 wCO wCO2 wH2 wH2O wN2]          # Matrix with all the mass fractions
x   = getMolarFractions(w)                 # Matrix with all the molar fractions
M   = getAvgMolarMass(x)                      # Average molar mass [kg mol^{-1}]
const rhoIn   = M[1]*pIn./(R*Tin) #Inlet density [kg/m3]
rho = M.*p./(R*T)
mu  = getViscosity(T,x)                                       # Viscosity [Pa s]
Re  = getReynolds(rho, uz, mu)                                 # Reynolds number
f   = getFrictionFactor(Re)                                    # Friction factor
cp  = getHeatCapacity(T,x)                    # Heat capacity [J K^{-1} kg^{-1}]
dH, reaction = getReaction(T,x,p)      # Enthalpy of reaction and reaction rates
                                       # [J kg^{-1} s^{-1}] [mol kg^{-1} s^{-1}]

# Define A matrices and b vectors
A_p      = [1 zeros(1,Nz-1);
            A[2:end,:]]
b_p      = [pIn;
            -1/dInner*(f.*rho.*uz.^2)[2:end,:]]
A_uz     = [1 zeros(1,Nz-1);
            (rho.*A + (A*rho).*I)[2:end,:]]
b_uz     = [uzIn;
            zeros(Nz-1)]
A_T      = [1 zeros(1,Nz-1);
            (rho.*cp.*uz.*A + 4*U/dInner)[2:end,:]]
b_T      = [Tin;
            4*U/dInner*Ta*ones(Nz-1,1) - ((1-void)*rhoCat*dH)[2:end,:]]
A_w      = [1 zeros(1,Nz-1);
            (rho.*uz.*A + uz.*(A*rho).*I + rho.*(A*uz).*I)[2:end,:]]
b_wCH4   = [wCH4in;
            ((1-void)*rhoCat*molarMass[1]*(reaction*N[:,1]))[2:end,:]]
b_wCO2   = [wCO2in;
            ((1-void)*rhoCat*molarMass[3]*(reaction*N[:,3]))[2:end,:]]
b_wH2    = [wH2in;
            ((1-void)*rhoCat*molarMass[4]*(reaction*N[:,4]))[2:end,:]]
b_wH2O   = [wH2Oin;
            ((1-void)*rhoCat*molarMass[5]*(reaction*N[:,5]))[2:end,:]]
b_wN2    = [wN2in;
            zeros(Nz-1)]
# Under-relaxation factors
gamma_p = 5e-1
gamma_uz = 5e-1
gamma_T = 1e-2
gamma_w0 = 1e-2

converged = false
iter = 1

while (!converged && (iter <= 1000000))
    println("Iteration number: $iter")

    # Solve ergun's equation for pressure
    #println("Solving for pressure")
    p = (1-gamma_p)*p + gamma_p*(A_p\b_p)
    #dp = A_p\b_p - p

    # Solve the ideal gas law for density
    rho = M.*p./(R*T)

    # solve the continuity for velocity
    uz = (1-gamma_uz)*uz + gamma_uz*(A_uz\b_uz)
    uz = max(0.0, uz)
    #println("Solving for velocity")
    #duz = A_uz\b_uz - uz
    #while minimum(uz + duz) <= 0
    #    duz /= 2
    #end
    #uz += duz

    println("Entering w-T loop")
    wResid = 1
    Tresid = 1
    while wResid > 1e-2 && Tresid > 1e-2
        dT = (A_T\b_T - T)*gamma_T
        while maximum(T + dT ) > 1100 || minimum(T + dT) < 500
            dT /= 10
        end
        T += dT
        wMax = 2.0
        wMin = -1.0
        gamma_w = gamma_w0
        while wMax > 1+1e-6 || wMin < -1e-4
            wCH4 = (1-gamma_w)*wCH4 + gamma_w*(A_w\b_wCH4)
            wCO2 = (1-gamma_w)*wCO2 + gamma_w*(A_w\b_wCO2)
            wH2 = (1-gamma_w)*wH2 + gamma_w*(A_w\b_wH2)
            wH2O = (1-gamma_w)*wH2O + gamma_w*(A_w\b_wH2O)
            wN2 = (1-gamma_w)*wN2 + gamma_w*(A_w\b_wN2)
            wCO = 1 - wCH4 - wCO2 - wH2 - wH2O - wN2
            w    = [wCH4 wCO wCO2 wH2 wH2O wN2]
            wMax = maximum(sum(w,2))
            wMin = minimum(w)
            wCH4 = max(0.0, min(wCH4in, wCH4))
            wCO2 = max(0.0, min(1.0-wCH4, wCO2))
            wH2 = max(wH2in, min(1.0-wCH4-wCO2, wH2))
            wH2O = max(0.0, min(wH2Oin, wH2O))
            wN2 = max(0.0, min(1.0-wCH4-wCO2-wH2-wH2O, wN2))
            wCO = max(wCOin, min(1.0, wCO))
            gamma_w /= 10
        end
        w   = [wCH4 wCO wCO2 wH2 wH2O wN2]
        x   = getMolarFractions(w)
        M   = getAvgMolarMass(x)
        rho = M.*p./(R*T)
        mu  = getViscosity(T,x)
        Re  = getReynolds(rho, uz, mu)
        f   = getFrictionFactor(Re)
        cp  = getHeatCapacity(T,x)

        dH, reaction = getReaction(T,x,p)

        A_w      = [1 zeros(1,Nz-1);
                    (rho.*uz.*A + uz.*(A*rho).*I + rho.*(A*uz).*I)[2:end,:]]
        b_wCH4   = [wCH4in;
                    ((1-void)*rhoCat*molarMass[1]*(reaction*N[:,1]))[2:end,:]]
        b_wCO2   = [wCO2in;
                    ((1-void)*rhoCat*molarMass[3]*(reaction*N[:,3]))[2:end,:]]
        b_wH2    = [wH2in;
                    ((1-void)*rhoCat*molarMass[4]*(reaction*N[:,4]))[2:end,:]]
        b_wH2O   = [wH2Oin;
                    ((1-void)*rhoCat*molarMass[5]*(reaction*N[:,5]))[2:end,:]]
        b_wN2    = [wN2in;
                    zeros(Nz-1)]
        A_T      = [1 zeros(1,Nz-1);
                    (rho.*cp.*uz.*A + 4*U/dInner.*I)[2:end,:]]
        b_T      = [Tin;
                    4*U/dInner*Ta*ones(Nz-1,1) - ((1-void)*rhoCat*dH)[2:end,:]]

        wCH4_err    = vec(A_w*wCH4 - b_wCH4)
        wCO2_err    = vec(A_w*wCO2 - b_wCO2)
        wH2_err     = vec(A_w*wH2 - b_wH2)
        wH2O_err    = vec(A_w*wH2O - b_wH2O)
        wN2_err     = vec(A_w*wN2 - b_wN2)
        wResid      = sqrt(dot(wCH4_err,wCH4_err))/mean(wCH4)
        wResid     += sqrt(dot(wCO2_err,wCO2_err))/mean(wCO2)
        wResid     += sqrt(dot(wH2_err,wH2_err))/mean(wH2)
        wResid     += sqrt(dot(wH2O_err,wH2O_err))/mean(wH2O)
        wResid     += sqrt(dot(wN2_err,wN2_err))/mean(wN2)
        
        T_err  = vec(A_T*T - b_T)
        Tresid = sqrt(dot(T_err,T_err))/mean(T)

        println("wResid: $wResid\tTresid: $Tresid")
    end

    # Update matrices and vectors
    b_p      = [pIn;
                -1/dInner*(f.*rho.*uz.^2)[2:end,:]]
    A_uz     = [1 zeros(1,Nz-1);
                (rho.*A + (A*rho).*I)[2:end,:]]
    A_T      = [1 zeros(1,Nz-1);
                (rho.*cp.*uz.*A + 4*U/dInner.*I)[2:end,:]]
    b_T      = [Tin;
                4*U/dInner*Ta*ones(Nz-1,1) - ((1-void)*rhoCat*dH)[2:end,:]]
    A_w      = [1 zeros(1,Nz-1);
                (rho.*uz.*A + uz.*(A*rho).*I + rho.*(A*uz).*I)[2:end,:]]
    b_wCH4   = [wCH4in;
                ((1-void)*rhoCat*molarMass[1]*(reaction*N[:,1]))[2:end,:]]
    b_wCO2   = [wCO2in;
                ((1-void)*rhoCat*molarMass[3]*(reaction*N[:,3]))[2:end,:]]
    b_wH2    = [wH2in;
                ((1-void)*rhoCat*molarMass[4]*(reaction*N[:,4]))[2:end,:]]
    b_wH2O   = [wH2Oin;
                ((1-void)*rhoCat*molarMass[5]*(reaction*N[:,5]))[2:end,:]]
    b_wN2    = [wN2in;
                zeros(Nz-1)]

    # Calculate residuals
    p_err       = vec(A_p*p - b_p)
    uz_err      = vec(A_uz*uz - b_uz)
    T_err       = vec(A_T*T - b_T)
    wCH4_err    = vec(A_w*wCH4 - b_wCH4)
    wCO2_err    = vec(A_w*wCO2 - b_wCO2)
    wH2_err     = vec(A_w*wH2 - b_wH2)
    wH2O_err    = vec(A_w*wH2O - b_wH2O)
    wN2_err     = vec(A_w*wN2 - b_wN2)
    sumResidual = sqrt(dot(p_err, p_err))/mean(p)
    sumResidual+= sqrt(dot(uz_err, uz_err))/mean(uz)
    sumResidual+= sqrt(dot(T_err, T_err))/mean(T)
    sumResidual+= sqrt(dot(wCH4_err,wCH4_err))/mean(wCH4)
    sumResidual+= sqrt(dot(wCO2_err,wCO2_err))/mean(wCO2)
    sumResidual+= sqrt(dot(wH2_err,wH2_err))/mean(wH2)
    sumResidual+= sqrt(dot(wH2O_err,wH2O_err))/mean(wH2O)
    sumResidual+= sqrt(dot(wN2_err,wN2_err))/mean(wN2)
    
    println("Residual: $sumResidual")

    converged = sumResidual < 1e-6

    iter += 1

end

using PyPlot
rc("text", usetex=true)
rc("font", family="serif")
rc("text.latex", preamble="\\usepackage[version=3]{mhchem}\\usepackage{siunitx}")
figure(1)
plot(Z,T)
xlabel(L"$z\quad\SI{}{\meter}$")
ylabel(L"$T\quad\SI{}{\kelvin}$")
savefig("T.pdf",transparent=true,pad_inches=0)
figure(2)
plot(Z,p)
xlabel(L"$z\quad\SI{}{\meter}$")
ylabel(L"$p\quad\SI{}{\pascal}$")
savefig("p.pdf",transparent=true,pad_inches=0)
figure(3)
plot(Z,uz)
xlabel(L"$z\quad\SI{}{\meter}$")
ylabel(L"$u_z\quad\SI{}{\meter\per\second}$")
savefig("uz.pdf",transparent=true,pad_inches=0)
for i = 1:Ncomp
    figure(3+i)
    plot(Z,x[:,i])
    xlabel(L"$z\quad\SI{}{\meter}$")
    ylabel("\$x_\\ce{$(Comp[i])}\$")
    savefig("$(Comp[i]).pdf",transparent=true,pad_inches=0)
end