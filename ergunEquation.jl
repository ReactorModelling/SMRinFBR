function ergunEquation!(p, rho, uz, f, b)
    # Get cross sectional averages (skip the inlet)
    rhoAvg = csAverage(rho)[2:end]              # Cross sectional average of rho
    uzAvg = csAverage(uz)[2:end]                 # Cross sectional average of uz
    fAvg = csAverage(f)[2:end]      # Cross sectional average of friction factor
    b[2:end] = [-fAvg.*rhoAvg.*uzAvg.^2/dParticle]         # Insert Ergun into b
end


function ergunEquation!(p, rho, uz, f, A, b)
    A[1,1] = 1.0                                         # Set initial condition
    A[2:end,:] = LagAz[2:end,:]                       # Set the axial derivative
    b[1] = pIn                                           # Set initial condition
    ergunEquation!(p, rho, uz, f, b)                       # Insert ergun into b
end