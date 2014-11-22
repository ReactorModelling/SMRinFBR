function getAvgMolarMass(molefraction)
    #=
        Average Molar Mass
         Calculates the average molar mass
         Input: 
           molarMass:      The molecular weights
           molefraction:   Mole fractions of the components
         Output:
           avgmolarmass:   The average molar mass
    =#
    avgmolarmass = molefraction*molarMass;                       # [kg mol^{-1}]
end