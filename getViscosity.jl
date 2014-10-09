function getViscosity(temperature, S, B, molefraction)
    #=
    %% Viscosity
    % Calculates the viscosities of the components at a given temperature
    % and mole fraction
    % Input:
    %   temperature:    The temperature, $T \mathrm{[K]}$
	%	S:              Component viscosity coefficient  $\mathrm{[K]}$
	%	B:              Component viscosity coefficient  $\mathrm{[kg K^{0.5}ms^{-1}]}$
    %   molefraction:   The molefractions (optional)
    % Output: 
    %   mu: The viscosity or the average viscosity depending on the number of
    %       inputs.
    %       $\mu_{i} = \frac{b_{i}T^1.5}{T+S_{i}}$
    %       $\bar{\mu} = \boldsymbol\mu^\mathrm{T}\mathbf{y}$
    =#
    mu = zeros(Nz,Ncomp);
    for i = 1:Ncomp
        mu[:,i] = (B[i]*temperature.^1.5)./(temperature+S[i]);
    end
    mu.*=molefraction
    mu*=ones(Ncomp,1)


    return mu
end