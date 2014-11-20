function ergunEquation!(p, rho, uz, f, b)
    rhoAvg = csAverage(rho)[2:end]
    uzAvg = csAverage(uz)[2:end]
    fAvg = csAverage(f)[2:end]
    b[:] = [pIn;-fAvg.*rhoAvg.*uzAvg.^2/dParticle]
end


function ergunEquation!(p, rho, uz, f, A, b)
    ergunEquation!(p, rho, uz, f, b)
    A[1,1] = 1.0
    A[2:end,:] = LagAz[2:end,:]
end