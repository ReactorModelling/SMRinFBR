function coupleMassTemperature!(w, T,
                                rho, uz, cp, dH, U, lambdaEff, reaction, D,
                                A_w, b_w, A_T, b_T,
                                A, b)
    
    energyEquation!(T, rho, uz, cp, dH, U, lambdaEff, A_T, b_T)
    speciesMassBalance!(w, rho, uz, reaction, D, A_w, b_w)

    # Species Mass
    for i = 1:length(CompIndex)
        A[(i-1)*Nglob+1:i*Nglob,(i-1)*Nglob+1:i*Nglob] = A_w[i] # Species Mass Balance
        b[(i-1)*Nglob+1:i*Nglob] = b_w[i]                       # Species Mass Balance
        A[5Nglob+1:6Nglob,(i-1)*Nglob+1:i*Nglob] = eye(Nglob)   # Sum fractions to 1
    end
    A[5Nglob+1:6Nglob,5Nglob+1:6Nglob] = eye(Nglob)   # Sum fractions to 1
    b[5Nglob+1:6Nglob] = ones(Nglob)                  # Sum fractions to 1

    # Temperature
    A[6Nglob+1:end,6Nglob+1:end] = A_T  # Energy equation
    b[6Nglob+1:end] = b_T               # Initial conditions
end