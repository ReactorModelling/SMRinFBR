function getMolarFractions(omega)
    #=
     Converts mass fractions to molar fractions
     Input:
       omega:  mass fractions $\boldsymbol\omega$
               Matrix: Columns correspond to components, rows to points.
     MolarMass:molecular mass $\mathbf{M}_\mathrm{mass}$
               Vector: Each row correspond to each component
     Output:
       y:      mole fractions $\mathbf{y}$
               Matrix: Same as $\boldsymbol\omega$
    =#
    y=(omega./molarMass')./(omega*(1./molarMass))
end