function getHeatCapacity(temperature, molefraction)
    #=
    %% Heat Capacity
    % Calculates the heat capacities of the components at a given temperature,
    % and returns them as a column vector.
    % Input:
    %   temperature:    The temperature, $T \mathrm{[K]}$
    %   coeffMatrix:    The coefficients $\alpha, \beta, \gamma, \delta$ as a
    %                   matrix
    %   molefraction:   The molefractions (optional)
    % Output: 
    %   cp: The heat capacities
    %       $\mathbf{c}_p = \boldsymbol\alpha + \boldsymbol\beta T + \boldsymbol\gamma T^2 + \boldsymbol\delta T^3$
    %       $\bar{c}_p = \mathbf{y}^\mathrm{T}\mathbf{c}_p$
    =#
    
    cp = [ones(Nz) temperature temperature.^2 temperature.^3]*cpCoeff; 
    # $\mathrm{[J kg^{-1}K^{-1}]}$
    cp.*=molefraction
    cp*=molarMass
end

