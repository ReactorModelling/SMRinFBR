function getAxialDerivative(psi)
    # Calculate the axial derivative of psi using the kronecker product
    kron(LagAz,Lagr)*psi
end

function getRadialDerivative(psi)
    # Calculate the radial derivative of psi using the kronecker product
    kron(Lagz,LagAr)*psi
end

function getRadialSecondDerivative(psi)
    # Calculate the radial second derivative of psi using the kronecker product
    kron(Lagz,LagBr)*psi
end