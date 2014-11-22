function getHeatCoefficients(Re,T,Y,gasViscosity,cpGas,avgMolarMass)
#=
 heatcoef
 The function computes the heat transfer coefficient for radial transport of
 heat from the bed to the surroundings
 Input:
 Re        [=] -               Reynolds number
 T         [=] K               Temperature
 Y         [=] -               Mol fraction
 gasVis 	[=]	  				Gas viscosity
 CPgas     [=] J/kgK           Gas heat capasity

 Output
 Ur        [=] J/m^2sK         Heat coefficient
 LAMBDAer  [=] J/msK           Effective radial conductivity
=#
#  Constants

pconst  = 1.0;
beta    = 1.0;
lambdas = 0.243;
phi     = 0.3;


# Calculates the gas heat conductivity
Tmatrix     = [ones(size(T)) T T.^2 T.^3];
lambdag     = Tmatrix*lambda';
lambdag   .*= Y;
lambdag    *= ones(Ncomp)

# Prandtl number
Pr = gasViscosity.*cpGas./lambdag;

# Radial effective static conduction
alpharv = (0.227/(1+void/(2*(1-void))*(1-pconst)/pconst)*(T/100).^3);
alphars = 0.227*pconst/(2-pconst)*(T/100).^3;
lambdaer0 = lambdag.*(void*(1 + beta*dParticle*alpharv./lambdag) + 
    beta*(1-void)./(1./(1/phi + alphars*dParticle./lambdag) + 2/3*lambdag/lambdas));

# Effective radial conductivity
lambdaer = lambdaer0+0.14*lambdag.*Re.*Pr;

# Heat transfer coefficient near the wall
alphaw0=8.694/(dInner)^(4/3)*lambdaer0[Nr:Nr:end];
alphaw=alphaw0+0.444*Re[Nr:Nr:end].*Pr[Nr:Nr:end].*lambdag[Nr:Nr:end]/dParticle;
Ur=(rInner*log(rOuter/rInner)/lambdaSt+1./alphaw).^(-1);
# Overall heat transfer coefficient
return lambdaer, Ur

end
