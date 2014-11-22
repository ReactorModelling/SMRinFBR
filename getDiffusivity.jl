function getDiffusivity(velocity)
#=
     Calculates the diffusivity 
     The Peclet numbers given in the problem formulation
     Input:
       dParticle:  Diameter of the particles $d_\mathrm{p}$
       velocity:   The axial velocity at each radial discretization $u_z$
       dInner:     The inner diameter of the tubes $r_\mathrm{i}$
       Output:
       diffusivity:    The diffusivity $D \mathrm{[m^2 s^{-1}]}$
=#
    
    peNumberrd  = 8*(2-(1-2*dParticle/dInner)^2); 	# Peclet number rd
    peNumbermr  = 1.1*peNumberrd; 				          # Peclet number mr
    diffusivity = velocity*dParticle/peNumbermr;
end
