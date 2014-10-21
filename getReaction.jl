function getReaction( T, y, P)

#=
%%  Reaction
%   This function calculates the reaction rates for all the components
%   and the heat of the reaction in all the discretication points.
%   Input:
%       T         [=] K                 Temperature
%       y         [=] -                 Mol fraction vector
%       P         [=] Pa                Total pressure  
%
%   Output
%       rrx     [=] mol/kg(cat)s        Reaction rates for each reaction
%                                       (#r-points x #components-matrix)
%       deltaH  [=] J/m^3s              Reaction heat for each reaction
%                                       (#r-points x #components-matrix)
=#

    #%%  Reaction enthalpies
    #% Three reactions for each point in the r-direction
    deltaH=broadcast(+,ent298',(T-298)/(948-298)*(ent948-ent298)')

    #%%  Partial pressures
    pComp=(p.*y)*1e-5;

    #% Initialization: reaction rate for each point in the r-direction
    rrx = zeros(Nz,3);
    #% Initialization: denominator in the rate expression for each point in 
    #% the r-direction
    denom = zeros(Nz);
    Keq=zeros(3);

    for i=1:Nz
        #% Rate constant
        Krx  = aj.*exp(-actEn./(R*T[i]));
        
        #% Adsorbtion constant
        Kads = ax.*exp(-adEnt./(R*T[i]));
        
        #% Equilibrium constants
        Keq[1] = 10^(-11650/T[i]+13.076);
        Keq[2] = 10^(1910/T[i]-1.784);
        Keq[3] = Keq[1]*Keq[2];

        #% Reaction rates
        denom[i] = 1 + Kads[2]*pComp[i,2] + Kads[3]*pComp[i,4]
            + Kads[1]*pComp[i,1] + Kads[4]*pComp[i,5]/pComp[i,4];

        rrx[i,1]= Krx[1]/(pComp[i,4])^2.5*(pComp[i,1]*pComp[i,5]
            - (pComp[i,4])^3*pComp[i,2]/Keq[1])/(denom[i])^2;
        
        rrx[i,2]= Krx[2]/(pComp[i,4])*(pComp[i,2]*pComp[i,5]
            - (pComp[i,4])*pComp[i,3]/Keq[2])/(denom[i])^2;
        
        rrx[i,3]= Krx[3]/(pComp[i,4])^3.5*(pComp[i,1]*(pComp[i,5])^2
            - (pComp[i,4])^4*pComp[i,3]/Keq[3])/(denom[i])^2;
    end

    return (deltaH, rrx)

end