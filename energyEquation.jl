function energyEquation(T, rho, uz, cp, dH, U, lambda, A, b)
    A[1:Nr,1:Nr] = eye(Nr)
    for iZ = 2:Nz
        for iR = 1:Nr
            iGlob = iR + (iZ-1)*(Nr)
            for jZ = 1:Nz
                for jR = 1:Nr
                    jGlob = jR + (jZ-1)*(Nr)
                    if iR == 1
                        A[iGlob,jGlob] = Lagz[iZ,jZ]*LagAr[iR,jR]
                    elseif iR == Nr
                        A[iGlob,jGlob] = Lagz[iZ,jZ]*LagAr[iR,jR] + 
                                         U[iGlob]/lambda[iGlob]
                    else
                        A[iGlob,jGlob] = rho[iGlob]*cp[iGlob]*uz[iGlob]*LagAz[iZ,jZ]*Lagr[iR,jR] -
                                         lambda[iGlob]*(1/r[iR]*Lagr)
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
                b[iGlob] = 0
            elseif iR == Nr
                b[iGlob] = U[iGlob]/lambda[iGlob]*Ta
            else
                b[iGlob] = dH[iGlob]
            end
        end
    end
end