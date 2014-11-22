function coupleVelocityDensityPressure!(uz, rho, p, 
                                        M, T, f,
                                        A_uz, b_uz, A_p, b_p,
                                        A, b)
    continuityEquation!(uz, rho, A_uz)
    ergunEquation!(p, rho, uz, f, b_p)
    # Build A matrix for coupled velocity, density and pressure
    A[1:Nglob,1:Nglob]                  = A_uz
    A[Nglob+1:2Nglob,Nglob+1:2Nglob]    = eye(Nglob)          # Velocity
    A[Nglob+1:2Nglob,2Nglob+1:end]      = diagm(-M./(R*T))    # Density
    A[2Nglob+1:end,2Nglob+1:end]        = kron(A_p,eye(Nr))   # Pressure

    b[1:Nglob]      = b_uz               #Initial conditions uz
    b[2Nglob+1:end] = kron(b_p,ones(Nr)) #Initial conditions p
end