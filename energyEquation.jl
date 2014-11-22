function energyEquation!(T, rho, uz, cp, dH, U, lambdaEff, A, b)
    # Update A and b for all temperatures
    # The update is done by reference
    A[1:Nr,1:Nr] = eye(Nr)
    for iZ = 2:Nz
        for iR = 1:Nr
            iGlob = iR + (iZ-1)*(Nr)
            for jZ = 1:Nz
                for jR = 1:Nr
                    jGlob = jR + (jZ-1)*(Nr) # Global index j (position)
                    if iR == 1
                        # Insert boundary condition for center symmetry
                        A[iGlob,jGlob] = Lagz[iZ,jZ]*LagAr[iR,jR]

                    elseif iR == Nr
                        # Insert boundary conditon for wall
                        A[iGlob,jGlob] = lambdaEff[iGlob]*Lagz[iZ,jZ]*LagAr[iR,jR]

                    else
                        # Insert energy equation: convection and radial diffusion
                        A[iGlob,jGlob] = 
                        (
                            rho[iGlob]*cp[iGlob]*uz[iGlob]
                          * LagAz[iZ,jZ]*Lagr[iR,jR]
                          - lambdaEff[iGlob]*(1/r[iR]*LagAr[iR,jR]*Lagz[iZ,jZ]
                          + LagBr[iR,jR]*Lagz[iZ,jZ])
                        )

                    end
                end
            end
        end
    end
    b[1:Nr] = ones(Nr)*Tin
    for iZ = 2:Nz
        for iR = 1:Nr
            iGlob = iR + (iZ - 1)*Nr
            if iR == 1
                continue
            elseif iR == Nr
                # Insert heat transfer through wall
                b[iGlob] = -U[iZ]*(T[iGlob] - Ta)

            else
                # Insert heat of reaction
                b[iGlob] = -dH[iGlob]
            end
        end
    end
end