export tl_tl_coup

function tl_tl_coup(initialConditions, parameters,
    stoichMatrix, terminalTime, targetSpecies, batchSize, K, e)
# TL_TL_COUP provides simultaneous simulations for two consecutive levels 
# in the multi-level algorithm. These are correction terms, and come from
# a batchSize number of pairs.
#
# We require initialConditions, parameters, a stoichMatrix, 
# a terminalTime, and a targetSpecies. Further, K and e determine the 
# time step, which is defined as terminalTime / (K^e). This setup is() 
# convenient for our multi-level estimator.
#
# See also MLMC_FIXED, and PROPENSITY.

# sets initial conditions up
y_f = repeat(initialConditions,outer=[1,batchSize])
y_c = y_f

# set update interval
dt = terminalTime / (K^(e+1))

for i=1:(K^e)
    # update coarse propensity values
    terms_c = propensity(y_c, parameters, batchSize)
    for j = 0:(K-1)
        
        # update fine propensity values
        terms_f = propensity(y_f, parameters, batchSize)
        
        # construction of the virtual reaction channels
        A = zeros(size(terms_f,1),size(terms_f,2),3)
        A[:,:,1] = broadcast(min, terms_f, terms_c)
        A[:,:,2] = broadcast(-, terms_f, A[:,:,1])
        A[:,:,3] = broadcast(-, terms_c, A[:,:,1])
        
        # random number generation
        delta = map(x -> rand(Poisson(x)),A*dt)
        
        # update state vectors
        fineUpdate = sum(delta[:,:,1:2],dims=3)
        coarseUpdate = sum(delta[:,:,1:2:3],dims=3)
        y_f = y_f + stoichMatrix * reshape(fineUpdate,
            size(fineUpdate,1),size(fineUpdate,2))
        y_c = y_c + stoichMatrix * reshape(coarseUpdate,
            size(coarseUpdate,1),size(coarseUpdate,2))
        
        # enforce boundary conditions
        whereNegative = map(x -> Bool(x < 0), y_f)
        if any(any(whereNegative,dims=1),dims=2)[1]
        y_f[whereNegative] = repeat([0],sum(sum(
             whereNegative,dims=1),dims=2)[1])
        end        
    end
    
    # enforce boundary conditions
    whereNegative = map(x -> x < 0, y_c)
    if any(any(whereNegative,dims=1),dims=2)[1]
        y_c[whereNegative] = repeat([0],sum(sum(
             whereNegative,dims=1),dims=2)[1])
    end
end

# output
out = adjoint(hcat(y_f[targetSpecies,:], y_c[targetSpecies,:]))
return out
end
