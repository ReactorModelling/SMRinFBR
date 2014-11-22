function speciesMassBalance!(w, rho, uz, reaction, D, A, b)
    # Update A and b for all the mass fractions
    # The update is done by reference
    for i in 1:length(A) # Loop through all the matrices in the array A
        speciesMassBalance!(w, rho, uz, reaction, D, A, b, i)
    end
end

function speciesMassBalance!(w, rho, uz, reaction, D, A, b, i)
    # Update A and b for all the mass fraction i
    # The update is done by reference
    c = CompIndex[i]                                           # Component index
    dRhodz = getAxialDerivative(rho)                   # Axial derivative of rho
    duzdz  = getAxialDerivative(uz)                     # Axial derivative of uz
    dRhodr = getRadialDerivative(rho)                 # Radial derivative of rho
    dwdr   = getRadialDerivative(w[:,c])              # Radial derivative of w_c
    dw2dr2 = getRadialSecondDerivative(w[:,c]) # Radial second derivative of w_c
    A_w = A[i]                                          # Extract current matrix
    A_w[1:Nr,1:Nr] = eye(Nr)                            # Set initial conditions
    b_w = b[i]                                          # Extract current vector
    b_w[1:Nr] = ones(Nr)*wIn[c]                         # Set initial conditions
    for iZ = 2:Nz
        for iR = 1:Nr
            iGlob = iR + (iZ-1)*Nr                   # Global index i (position)
            if iR == 1 || iR == Nr
                # Insert boundary condition
                # dwdr = 0 at center and wall
                b_w[iGlob] = 0.0
            else
                # Insert species mass balance reaction
                b_w[iGlob] = molarMass[c]*dot(vec(reaction[iGlob,:]),N[:,c])
            end
            for jZ = 1:Nz
                for jR = 1:Nr
                    jGlob = jR + (jZ-1)*Nr         # Global index j (polynomial)
                    if iR == 1 || iR == Nr
                        # Insert boundary condition
                        # dwdr = 0 at center and wall
                        A_w[iGlob,jGlob] = Lagz[iZ,jZ]*LagAr[iR,jR]
                    else
                        # Insert species mass balance convection and diffusion
                        A_w[iGlob,jGlob] = (
                            (
                                rho[iGlob]*uz[iGlob]*LagAz[iZ,jZ]*Lagr[iR,jR]
                              + dRhodz[iGlob]*uz[iGlob]*Lagz[iZ,jZ]*Lagr[iR,jR]
                              + duzdz[iGlob]*rho[iGlob]*Lagz[iZ,jZ]*Lagr[iR,jR]
                            )
                          - D[iGlob]
                          * (
                                rho[iGlob]/r[iR]*Lagz[iZ,jZ]*LagAr[iR,jR]
                              + dRhodr[iGlob]*Lagz[iZ,jZ]*LagAr[iR,jR]
                              + rho[iGlob]*Lagz[iZ,jZ]*LagBr[iR,jR]
                            )
                        )
                    end
                end
            end
        end
    end
end