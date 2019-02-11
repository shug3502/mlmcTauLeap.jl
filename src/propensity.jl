export propensity

function propensity(y, p, batchSize)
# PROPENSITY evaluates the propensity function values for batchSizes
# paths, based on the state vectors give in y, and the parameters in p.
#
# See also EULER, MLMC_L, COUPLED.
terms = adjoint(hcat(p[1]*ones(batchSize),
    p[2]*y[1,:],
    p[3]*y[2,:] .* broadcast(-,y[2,:], 1),
    p[4]*y[1,:],
    p[5]*y[2,:]))
return terms

end
