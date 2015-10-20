# Self-adaptive_DE
How parameters in DE evolve along with searching?

In Differential Evolution algorithm, there're two specific parameters __F__ and __CR__(wait a minute! why no one metione __NP__), which control mutation and crossover, finally
the searching direction and convergence performance.

There're several ways to do parameter control

*  __Deterministic__: Adjust the parameters accroding to some deterministic rules, like runing time.
*  __Adaptive__: using the feedback from searching to change parameter values, like "1/5-th rule"
*  __Self-Adaptive__:parameters are encoded in the genome, which gets into evolution as well.

