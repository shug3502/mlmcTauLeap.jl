export mlmc_main

function mlmc_main(K, minLevel, maxLevel,
    epsilon, unbiased)
# MLMC_FIXED is the main method for the multi-level method.
#
# It takes in K, the scaling factor(), a minLevel and a maxLevel, which
# determine the time steps throughout. epsilon determines the statistical
# error target: a 95# confidence interval of semi-length epsilon will be 
# constructed. unbiased is a boolean input as to whether the 
# algorithm is eventually a DM [true] or TL [false] equivalent.
#
# See also tl, tl_tl_coup, dm_tl_coup, ESTIMATE, PROPENSITY.


# ----------------------------------------------------------------------
initialConditions = [0; 0; 0]
parameters = [25; 1000; 0.001; 0.1; 1]
terminalTime = 1

stoichMatrix = [1   0   0  -1   0
                0   1  -2   0  -1
                0   0   1   0   0]

targetSpecies = 3
# ----------------------------------------------------------------------

# This controls the number of paths to use for the trial runs

trials = 1000
trials_final = 100

# ----------------------------------------------------------------------

# output basic information

totalLevels = maxLevel - minLevel + 1 + unbiased

@printf("Discrete state multi-level Monte Carlo estimator\n")
@printf("\nPerforming survey simulations, using %d levels. ", 
    totalLevels)
if unbiased
    @printf("The estimator is unbiased. \n\n")
else
    @printf("The estimator is biased. \n\n")
    
end

# ----------------------------------------------------------------------

# perform trial / initial simulations
# sum1[l] records sum of first moments of samples of the l-th estimator
# similarly, sum2[l] does the same for second moments

# base level

timer = zeros(totalLevels,1)
sum1 = zeros(totalLevels,1)
sum2 = zeros(totalLevels,1)
num = zeros(totalLevels,1)
timer[1] = @elapsed temp = tl(initialConditions,stoichMatrix, terminalTime, 
    parameters, targetSpecies, trials, K, minLevel)

sum1[1] = sum(temp)
sum2[1] = sum(temp.^2)
num[1] = trials

# ----------------------------------------------------------------------

# now we estimate the different, higher levels
for i = 2:(totalLevels - unbiased)
    timer[i] = @elapsed temp = tl_tl_coup(initialConditions, parameters, 
        stoichMatrix, terminalTime, targetSpecies, trials, 
        K, minLevel+i-2)
    
    sum1[i] = sum(temp[1,:]- temp[2,:])
    sum2[i] = sum((temp[1,:]- temp[2,:]).^2)
    num[i] = trials
end

# ----------------------------------------------------------------------

if unbiased
    timer[totalLevels] = @elapsed temp = dm_tl_coup(initialConditions, stoichMatrix, 
        terminalTime, parameters, targetSpecies, trials_final, 
        terminalTime*float(K)^-(maxLevel))
println(size(temp))    
    sum1[totalLevels] = sum(temp[1,:]- temp[2,:])
    sum2[totalLevels] = sum((temp[1,:]- temp[2,:]).^2)
    num[totalLevels] = trials_final
end

# ----------------------------------------------------------------------

# calculate K_l. We have time proportional to [number runs] and also time
# proportional to 1 / (estimator variance). Hence is K_l such that
# Time_l = K_l / V_l, where V_l is variance estimator. Whence

pathsToDo = optimal_config(sum1, sum2, num, timer, (epsilon/1.96)^2)

# ----------------------------------------------------------------------

# outputs plan of action
@printf("We have the following preliminary data: \n")
@printf("Level \t\t Mean \t\t\t Done \t\t\t Remaining \n")
for i=1:totalLevels
   @printf("%2.0f \t\t %6.2f \t\t %6.0f \t\t %6.0f \n", i, 
        sum1[i]./num[i],  num[i], pathsToDo[i])
end

# ----------------------------------------------------------------------

# performs simulations
for i=1:(totalLevels - unbiased)
    if pathsToDo[i] > 0
        # if the base level, use tau leaping, else dm_tl_coup paths
        if i == 1
            
            dmTiming = @elapsed temp = tl(initialConditions,stoichMatrix, terminalTime, 
                parameters, targetSpecies, pathsToDo[i], K,  
                minLevel+i-1)

            timingPart1 = @elapsed sum1[i] = sum1[i] + sum(temp)
            timingPart2 = @elapsed sum2[i] = sum2[i] + sum(temp.^2)
            timingPart3 = @elapsed num[i] = num[i] + length(temp)
        else
            dmTiming = @elapsed temp=tl_tl_coup(initialConditions,parameters,stoichMatrix, 
                terminalTime,targetSpecies,pathsToDo[i],K, 
                minLevel+i-2)

            timingPart1 = @elapsed sum1[i] = sum1[i] + sum(temp[1,:]- temp[2,:])
            timingPart2 = @elapsed sum2[i] = sum2[i] + sum((temp[1,:]- temp[2,:]).^2)
            timingPart3 = @elapsed num[i] = num[i] + length(temp[1,:])
        end
        timer[i] = timer[i] + dmTiming + timingPart1 + timingPart2 + timingPart3
    end
end

if unbiased
    level = totalLevels
    if pathsToDo[level] > 0
    #timing is different to in MATLAB, could rewrite as function for single line
    #or import a different timing method
        dmTiming = @elapsed temp = dm_tl_coup(initialConditions, stoichMatrix, 
            terminalTime, parameters, targetSpecies, 
            pathsToDo[level], terminalTime*float(K)^(-maxLevel) )
println(size(temp))        
        timingPart1 = @elapsed sum1[level] =  sum1[level] + sum(temp[1,:]- temp[2,:])
        timingPart2 = @elapsed sum2[level]  = sum2[level] + sum((temp[1,:]- temp[2,:]).^2)
        timingPart3 = @elapsed num[level] = num[level] + length(temp[1,:])
        
        timer[level] = timer[level] + dmTiming + timingPart1 + 
             timingPart2 + timingPart3
    end
end

# ----------------------------------------------------------------------

# Done with all cases!

@printf("\n\nWe have the following final data: \n")
@printf("Level \t\t Estimate \t\t Samp Var \t\t Performed \n")
for i=1:(totalLevels)
    @printf("%2.0f \t\t %6.2f \t\t  %6.2f \t\t %6.0f \n", i, 
        sum1[i]./num[i], (sum2[i]./num[i] - (sum1[i]./num[i]).^2)  
        , num[i])
end
conf = (sum(broadcast(-, sum2 ./ broadcast(max, (num.^2),1), 
    (sum1.^2)./ broadcast(max,1,(num.^3))), dims=1))[1]^(.5)*1.96

out = sum(sum1 ./ broadcast(max,num,1),dims=1)[1]

# ----------------------------------------------------------------------

@printf("\nWe therefore estimate %0.2f ± %0.2f. \n", out,  
        conf)
@printf("Total time was %0.2f seconds. \nPer level times: ", 
    sum(timer))
map(x -> @printf("%0.2f\t",x), timer)
@printf(".\n")

return out, conf, timer
end
# ----------------------------------------------------------------------

