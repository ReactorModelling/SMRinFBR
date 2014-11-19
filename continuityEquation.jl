function continuityEquation(uz, rho, A)
    dRhodz = kron(LagAz,Lagr)*rho
    A[1:Nr,1:Nr] = eye(Nr)
    for iZ = 2:Nz
        for iR = 1:Nr
            iGlob = iR + (iZ-1)*(Nr)
            for jZ = 1:Nz
                for jR = 1:Nr
                    jGlob = jR + (jZ-1)*(Nr)
                    if iR == 1 || iR == Nr
                        A[iGlob,jGlob] = Lagz[iZ,jZ]*LagAr[iR,jR]
                    else
                        A[iGlob,jGlob] = rho[iGlob]*LagAz[iZ,jZ]*Lagr[iR,jR] +
                                         dRhodz[iGlob]*Lagz[iZ,iR]*Lagr[iR,jR]
                    end
                end
            end
        end
    end
end


function continuityEquation(uz, rho, A, b)
    continuityEquation(uz, rho, A)
    b[1:Nr] = uzIn*ones(Nr)
end