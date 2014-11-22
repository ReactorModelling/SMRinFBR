function couple!(w, T, uz, rho, p,
                 f, cp, dH, U, lambdaEff, reaction, D,
                 A_w, b_w, A_T, b_T, A_uz, b_uz, A_p, b_p,
                 A, b)
    # Fill A matrix and B vector for the coupled solution
    # Update matrices and vectors
    continuityEquation!(uz, rho, A_uz)                                # Velocity
    ergunEquation!(p, rho, uz, f, b_p)                                # Pressure
    energyEquation!(T, rho, uz, cp, dH, U, lambdaEff, A_T, b_T)    # Temperature
    speciesMassBalance!(w, rho, uz, reaction, D, A_w, b_w)      # Mass fractions
    # Insert species to be solved with a partial differential equation
    for i = 1:length(CompIndex)
        A[(i-1)*Nglob+1:i*Nglob,(i-1)*Nglob+1:i*Nglob] = A_w[i]
        b[(i-1)*Nglob+1:i*Nglob]                       = b_w[i]
        # Insert identity matrix for solving CO
        A[5Nglob+1:6Nglob,(i-1)*Nglob+1:i*Nglob]       = eye(Nglob)
    end
    # Solve CO using the sum of mass fractions
    A[5Nglob+1:6Nglob,5Nglob+1:6Nglob] = eye(Nglob)
    b[5Nglob+1:6Nglob]                 = ones(Nglob)
    # Temperature
    A[6Nglob+1:7Nglob,6Nglob+1:7Nglob] = A_T
    b[6Nglob+1:7Nglob]                 = b_T
    # Velocity
    A[7Nglob+1:8Nglob,7Nglob+1:8Nglob] = A_uz
    b[7Nglob+1:8Nglob]                 = b_uz
    # Density
    A[8Nglob+1:9Nglob,8Nglob+1:9Nglob] = eye(Nglob)
    A[8Nglob+1:9Nglob,9Nglob+1:end]    = diagm(-M./(R*T))
    # Pressure
    A[9Nglob+1:end,9Nglob+1:end] = kron(A_p,eye(Nr))
    b[9Nglob+1:end]              = kron(b_p,ones(Nr))
end
