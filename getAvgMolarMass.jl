function getAvgMolarMass(molarMass, molefraction)
    #=
    %% Average Molar Mass
    % Calculates the average molar mass
    % Input: 
    %   molarMass:      The molecular weights, $\mathbf{M}_\mathrm{mass} \mathrm{[kg mol^{-1}]}$
    %   molefraction:   Mole fractions of the components, $\mathbf{y}$
    % Output:
    %   avgmolarmass:   The average molar mass, $\bar{M}_\mathrm{mass} \mathrm{[kg mol^{-1}]}$
    
    % $\bar{M}_\mathrm{mass} = \mathbf{M}_\mathrm{mass}^\mathrm{T}\mathbf{y}$
    =#
    avgmolarmass = molarMass'*molefraction; #% $\mathrm{[kg mol^{-1}]}$
end