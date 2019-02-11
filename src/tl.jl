using Distributions

export tl

function tl(initialConditions, stoichMatrix, terminalTime,
    parameters, targetSpecies, batchSize::Int64, K, e)
# TL simultaneously generates sample paths for a biochemical reaction
# network, according to the tau-leaping scheme. It records the terminal
# populations of the targetSpecies species of batchSize paths.
#
# This method requires initialConditions, a stoichMatrix, parameters as
# inputs. Further, K and e determine the time step, which is defined as
# terminalTime / (K^e). This setup is convenient for our multi-level
# estimator.
#
# See also PROPENSITY, and MLMC_FIXED.

y = repeat(initialConditions,outer=[1,batchSize])
dt = terminalTime / (K^e)

for n = 1:(K^e)
    # this steps calculates the propensity functions values
    terms = propensity(y, parameters, batchSize)
    #this step updates the population values
    y = y + stoichMatrix * map(x -> rand(Poisson(x)), terms*dt)

    # This step enforces a boundary condition
    whereNegative = map(x -> x < 0, y)
    if any(any(whereNegative,dims=1),dims=2)[1] 
        y[whereNegative] = repeat([0],sum(sum(
             whereNegative,dims=1),dims=2)[1])
    end
end

out = y[targetSpecies,:]
return out
end

