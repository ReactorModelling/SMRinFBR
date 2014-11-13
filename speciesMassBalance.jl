function speciesMassBalance(w, rho, uz, reaction, D, A, b)
    for i in 1:length(A)
        c = CompIndex[i]
        A_w = A[i]
        A_w[1:Nr,1:Nr] = eye(Nr)
        for iZ = 2:Nz
            for iR = 1:Nr
                iGlob = iR + (iZ-1)*(Nr)
                for jZ = 1:Nz
                    for jR = 1:Nr
                        jGlob = jR + (jZ-1)*(Nr)
                        if iR == 1 || iR == Nr
                            A_w[iGlob,jGlob] = Lagz[iZ,jZ]*LagAr[iR,jR]
                        else
                            A_w[iGlob,jGlob] = 
                            (
                                rho[iGlob]*uz[iGlob]*w[jGlob,c]*Lagr[iR,jR]*LagAz[iZ,jZ]
                              + uz[iGlob]*w[iGlob,c]*rho[jGlob]*Lagr[iR,jR]*LagAz[iZ,jZ]
                              + rho[iGlob]*w[iGlob,c]*uz[jGlob]*Lagr[iR,jR]*LagAz[iZ,jZ]
                              - D[iGlob]
                              * (
                                    rho[iGlob]/r[iR]*Lagz[iZ,jZ]*LagAr[iR,jR]
                                  + rho[iGlob]*Lagz[iZ,jZ]*LagBr[iR,jR]
                                  + LagAr[iR,jR]*rho[jGlob]*Lagz[iZ,jZ]*LagAr[iR,jR]
                                )
                            )
                        end
                    end
                end
            end
        end
        b_w = b[i]
        b_w[1:Nr] = ones(Nr)*wIn[c]
        for iZ = 2:Nz
            for iR = 1:Nr
                iGlob = iR + (iZ - 1)*Nr
                if iR == 1 || iR == Nr
                    b_w[iGlob] = 0.0
                else
                    b_w[iGlob] = molarMass[c]*dot(vec(reaction[iGlob,:]),N[:,c])
                end
            end
        end
    end
end
