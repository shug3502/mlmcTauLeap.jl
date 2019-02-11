This is an implementation of Multilevel Monte Carlo for stochastic chemical
kinetics models using tau leap methods. It is based on previous implementations
of MLMC methods from Chris Lester and Mike Giles.  

For further detail on the methods see:

Lester, C., Yates, C. A., & Baker, R. E. (2018).
Robustly simulating biochemical reaction kinetics using multi-level Monte Carlo
approaches.
Journal of Computational Physics, 375, 1401-1423.

Lester, C., Baker, R. E., Giles, M. B., & Yates, C. A. (2016).
Extending the multi-level method for the simulation of stochastic biological
systems.
Bulletin of mathematical biology, 78(8), 1640-1677.

Lester, C., Yates, C. A., Giles, M. B., & Baker, R. E. (2015).
An adaptive multi-level simulation algorithm for stochastic biological systems.
The Journal of chemical physics, 142(2), 01B612_1.

Anderson, D. F., & Higham, D. J. (2012).
Multilevel Monte Carlo for continuous time Markov chains, with applications
in biochemical kinetics.
Multiscale Modeling & Simulation, 10(1), 146-179.

-------------------------------------------------------------------------------

Typical output should look like this:

Discrete state multi-level Monte Carlo estimator

Performing survey simulations, using 5 levels. The estimator is unbiased. 

We have the following preliminary data: 
Level 		 Mean 			 Done 			 Remaining 
 1 		 3212.92 		    100 		 281521 
 2 		 352.87 		    100 		  17725 
 3 		 124.86 		    100 		   4187 
 4 		  38.08 		    100 		   2470 
 5 		  16.32 		    100 		    900 


We have the following final data: 
Level 		 Estimate 		 Samp Var 		 Performed 
 1 		 3189.23 		  1028088.90 		 281621 
 2 		 350.26 		  15805.81 		  17825 
 3 		 118.05 		  2521.04 		   4287 
 4 		  38.55 		  933.27 		   2570 
 5 		  19.67 		  341.16 		   1000 

We therefore estimate 3715.76 ± 4.73. 
Total time was 41.71 seconds. 
Per level times: 15.40, 4.55, 2.45, 3.89, 15.41, .



