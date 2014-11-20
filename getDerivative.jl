function getAxialDerivative(psi)
    kron(LagAz,Lagr)*psi
end

function getRadialDerivative(psi)
    kron(Lagz,LagAr)*psi
end

function getRadialSecondDerivative(psi)
    kron(Lagz,LagBr)*psi
end