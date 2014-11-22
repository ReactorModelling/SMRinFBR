function coupleMassTemperature!(w, T,
                                rho, uz, cp, dH, U, lambdaEff, reaction, D,
                                A_w, b_w, A_T, b_T,
                                A, b)
    
    energyEquation!(T, rho, uz, cp, dH, U, lambdaEff, A_T, b_T)
    speciesMassBalance!(w, rho, uz, reaction, D, A_w, b_w)

    for i = 1:length(CompIndex)
        A[(i-1)*Nglob+1:i*Nglob,(i-1)*Nglob+1:i*Nglob] = A_w[i]
        b[(i-1)*Nglob+1:i*Nglob] = b_w[i]
        A[5Nglob+1:6Nglob,(i-1)*Nglob+1:i*Nglob] = eye(Nglob)
    end
    A[5Nglob+1:6Nglob,5Nglob+1:6Nglob] = eye(Nglob)
    b[5Nglob+1:6Nglob] = ones(Nglob)

    A[6Nglob+1:end,6Nglob+1:end] = A_T
    b[6Nglob+1:end] = b_T
end