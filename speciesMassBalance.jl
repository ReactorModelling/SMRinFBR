function speciesMassBalance(w, rho, uz, reaction, D, A, b)
    dRhodz = kron(LagAz, Lagr)*rho
    duzdz  = kron(LagAz, Lagr)*uz
    dRhodr = kron(LagAr, Lagz)*rho
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
                                rho[iGlob]*uz[iGlob]*LagAz[iZ,jZ]*Lagr[iR,jR]
                              + dRhodz[iGlob]*uz[iGlob]*Lagz[iZ,jZ]*Lagr[iR,jR]
                              + duzdz[iGlob]*rho[iGlob]*Lagz[iZ,jZ]*Lagr[iR,jR]
                            )
                            - D*
                            (
                                rho[iGlob]/r[iR]*Lagz[iZ,jZ]*LagAr[iR,jR]
                              + dRhodr[iGlob]*Lagz[iZ,jZ]*LagAr[iR,jR]
                              + rho[iGlob]*Lagz[iZ,jZ]*LagBr[iR,jR]
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
