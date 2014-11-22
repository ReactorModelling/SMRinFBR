function getReynolds(density, velocity, viscosity)
    #=
    Function that calculates the Reynolds number.
    =# 
    density.*velocity.*dParticle./viscosity;
end