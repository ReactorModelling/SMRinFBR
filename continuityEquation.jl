function continuityEquation!(uz, rho, A)
    # Update A for the velocity
    dRhodz = getAxialDerivative(rho)                   # Axial derivative of rho
    A[1:Nr,1:Nr] = eye(Nr)                              # Set initial conditions
    for iZ = 2:Nz
        for iR = 1:Nr
            iGlob = iR + (iZ-1)*Nr                   # Global index i (position)
            for jZ = 1:Nz
                for jR = 1:Nr
                    jGlob = jR + (jZ-1)*Nr         # Global index j (polynomial)
                    # Insert the convective terms of the continuity equation
                    A[iGlob,jGlob] = rho[iGlob]*LagAz[iZ,jZ]*Lagr[iR,jR] +
                                     dRhodz[iGlob]*Lagz[iZ,jZ]*Lagr[iR,jR]
                end
            end
        end
    end
end


function continuityEquation!(uz, rho, A, b)
    # Update A for the velocity
    continuityEquation!(uz, rho, A)
    # Update b for the velocity
    b[1:Nr] = uzIn*ones(Nr)                             # Set initial conditions
end