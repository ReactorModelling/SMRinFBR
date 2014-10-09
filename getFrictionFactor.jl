function getFrictionFactor(Reynolds)
    #=
    Function that calculates the Erguns equation friction factor
    =# 

    f=(1-void)/(void^3)*(1.75+4.2*(Reynolds.^(5/6)).*((1-void)./Reynolds));

end