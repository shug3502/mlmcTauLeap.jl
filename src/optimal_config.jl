export optimal_config

function optimal_config(sum1, sum2, num, time, err)::Array{Int64,2}

# OPTIMALCONFIGURATION provides an optimal configuration for the
# number of samples to generate for each estimator. This requires
# the sum of the first moments as sum1, the sum of the second moments
# as sum2, the number of simulations as num, the total time taken as
# time, and finally, error as the target estimator variance (distinct
# from a semilength of a confidence interval.

Kl = (sum2./ (num.^2) - (sum1.^2) ./ (num .^ 3)) .* time
targetV = err * Kl.^.5 / sum(Kl.^.5)
targetN = (sum2./num - (sum1.^2)./(num.^2)) ./ targetV
targetN = broadcast(max,broadcast(ceil, broadcast(-, targetN, num)), 0)

return targetN
end
