function getAxialDerivative(psi)
    dPsidz = zeros(Nglob)
    for iZ in 1:Nz
        for iR in 1:Nr
            iGlob = iR + (iZ-1)*(Nr)
            for jZ in 1:Nz
                for jR in 1:Nr
                    jGlob = jR + (jZ-1)*Nr
                    dPsidz[iGlob] += psi[jGlob]*Lagr[iR,jR]*LagAz[iZ,jZ]
                end
            end
        end
    end
    return dPsidz
end

function getRadialDerivative(psi)
    dPsidr = zeros(Nglob)
    for iZ in 1:Nz
        for iR in 1:Nr
            iGlob = iR + (iZ-1)*(Nr)
            for jZ in 1:Nz
                for jR in 1:Nr
                    jGlob = jR + (jZ-1)*Nr
                    dPsidr[iGlob] += psi[jGlob]*LagAr[iR,jR]*Lagz[iZ,jZ]
                end
            end
        end
    end
    return dPsidr
end

function getRadialSecondDerivative(psi)
    d2Psidr2 = zeros(Nglob)
    for iZ in 1:Nz
        for iR in 1:Nr
            iGlob = iR + (iZ-1)*(Nr)
            for jZ in 1:Nz
                for jR in 1:Nr
                    jGlob = jR + (jZ-1)*Nr
                    d2Psidr2[iGlob] += psi[jGlob]*LagBr[iR,jR]*Lagz[iZ,jZ]
                end
            end
        end
    end
    return d2Psidr2
end