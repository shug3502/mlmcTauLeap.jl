export dm_tl_coup

function dm_tl_coup(initialConditions, stoichMatrix, 
    terminalTime, parameters, targetSpecies, samplePaths, dt)
# DM_TL_COUP provides a de-biasing estimator for the multi-level
# method. This couples a DM path with a Tau-leaping path.
#
# We require initialConditions, a stoichMatrix, a terminalTime
# parameters, and a targetSpecies. In order to perform the estimate, we
# need a number of samplePaths, as well as the overall tau-leaping time
# step, given as dt.
#
# See also MLMC_FIXED, and PROPENSITY.

out = m = Array{Float64}(undef, 2, 0)

# gets number of reaction channels - will be used below
reactionChannels = size(stoichMatrix,2)

# x stores exact state vectors, and z the biased state vectors
x = repeat(initialConditions, outer=[1, samplePaths])
z = repeat(initialConditions, outer=[1, samplePaths])

# this records the time, t, as well as time of next update, Ttau
t = zeros(samplePaths)
Ttau = repeat([dt], outer=[samplePaths])

# we calculate all biased propensities
zProp = propensity(z,parameters, samplePaths)

while samplePaths > 0
    
    # recalculate the unbiased, DM propensities

    xProp = propensity(x,parameters,samplePaths)
    
    minTerms = broadcast(min, xProp, zProp)
    maxTerms = broadcast(max, xProp, zProp)
    
    cs = cumsum(maxTerms,dims=1)
    
    # time until next reaction
    u = rand(samplePaths)
    delta = [- log(u[x])./cs[end,x] for x in 1:samplePaths]
    
    # detects which paths require their Tau leap propensities
    # updated before continuing

    updateProp = [t[x] + delta[x] >= Ttau[x] for x in 1:samplePaths]
    while any(updateProp)
        updateTodo = sum(updateProp)
        zProp[:,updateProp] = propensity(z[:,updateProp],
            parameters, updateTodo)

        # recalculate virtual channels
        minTerms = broadcast(min, xProp, zProp)
        maxTerms = broadcast(max, xProp, zProp)
        
        cs[:,updateProp] = cumsum(maxTerms[:,updateProp],dims=1)
        
        # recalculate times
        t[updateProp] = broadcast(min,Ttau[updateProp],terminalTime)
        Ttau[updateProp] = broadcast(+,Ttau[updateProp], dt)

        u = rand(length(updateProp))
        delta[updateProp] = [-log(u[x])./cs[end,x] for x in 1:length(updateProp) if updateProp[x]]
        
        updateProp = vec(broadcast(>=, broadcast(+, t, delta), Ttau))
    end
    
    # implements a reaction in each channel

    randReac = reshape(rand(samplePaths).*cs[end,:],1,samplePaths)
    
    # choose which reaction
    reactionIndex = vec(broadcast(+,sum(broadcast(<, cs, repeat(randReac, outer=[reactionChannels,1])),dims=1), 1))
    
    # this implements reaction in both channels
    # and then reverses the reaction if appropriate
    t = broadcast(+,t, delta)
    x = broadcast(+, x, stoichMatrix[:, reactionIndex])
    z = broadcast(+, z, stoichMatrix[:, reactionIndex])
    
    
    q = broadcast(+,reactionIndex,
        (0:reactionChannels:((samplePaths-1)*reactionChannels)))
    
    # test if any reactions need to be reversed
    testRatio = minTerms[q] ./ maxTerms[q]
    sampleRatio = rand(samplePaths)
    
    removeReactions = broadcast(>, sampleRatio, testRatio)
    
    if any(removeReactions)
        zLow = vec((xProp[q] > zProp[q]) .* removeReactions)
        xLow = vec((xProp[q] < zProp[q]) .* removeReactions)
        
        z[:, zLow] = broadcast(-, z[:, zLow], stoichMatrix[:, reactionIndex[zLow]])
        x[:, xLow] = broadcast(-, x[:, xLow], stoichMatrix[:, reactionIndex[xLow]])           
    end
    
    # handle any paths which have reached terminal time
    if any(broadcast(>, t, terminalTime))
        
        completedPaths = broadcast(>,t, terminalTime)
        
        v = adjoint(hcat(x[targetSpecies,completedPaths],
            z[targetSpecies,completedPaths]))
        out = hcat(out, v)
        
        # clears away their internal data, in julia have to reassign to remove from a matrix
        x = x[:, .!completedPaths]
        z = z[:, .!completedPaths]
        
        t = t[.!completedPaths]
        Ttau = Ttau[.!completedPaths]
        zProp = zProp[:,.!completedPaths] 
        
        # number of remaining paths to simulate
        samplePaths = samplePaths - sum(completedPaths)
    end
    
end
return out
end
