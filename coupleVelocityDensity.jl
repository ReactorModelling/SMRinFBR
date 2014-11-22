function coupleVelocityDensity!(uz, rho, M, p, T, A_uz, b_uz, A, b)
    #A = zeros(2Nglob,2Nglob)
    continuityEquation!(uz, rho, A_uz)
    A[1:Nglob,1:Nglob] = A_uz
    A[Nglob+1:end,Nglob+1:end] = eye(Nglob)
    b[1:Nglob] = b_uz
    b[Nglob+1:end] = p.*M./(R*T)
end