# interval_censored_regular_vines

Here is some code from my master thesis. 

The file censor_est.R contains some expressions for the copula functions, mainly for the one parameter families. The rest are computed from functions in VineCopula. The copulae are constructed in lists, I didn't bother creating proper classes.

In censor_est.R, you can find the main parts of the interval censored vines, but not the "FullCensored" vines. That can be found in rvine_censor_full_pvalue.R. For simplicity, goodness-of-fit tests are calculated when the vines are constructed. A similar function for simple censoring can be found in rvine_censor_full_pvalue.R. 

censor_gof_parallel.R contains a parallelized bootstrapping scheme.

The simulation experiments can be found in simulations_global_optim_lord_markov.R

