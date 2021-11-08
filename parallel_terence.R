
rm(list=ls())
graphics.off()
cat("\014")

# SIMULATION FUNCTION ----------------------------------------------------------

comparison_fun = function(simulation_seed, n)
{
  set.seed(simulation_seed)
  aaa = rnorm(n)
  
  return(aaa)
}


# LOOP -------------------------------------------------------------------------

library(snowfall)

# init cluster parallelization
ncpu = 20
sfInit(par=TRUE,cp=ncpu)

# number of simulations
nsim = 100
seeds = 1:nsim

n = 10
out_mean = sfSapply(x = seeds, 
                    fun = comparison_fun,
                    n = n)

dim(out_mean) # n x nsim

# stop cluster parallelization
sfStop()

