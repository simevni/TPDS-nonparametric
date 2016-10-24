# TPDS-nonparametric
Nonparametric estimation of the welfare change under TPDS policy reform using Bihar and Rajasthan PoC data.
Plot the observed cash take up, estimate and plot the CDF and the demand curve, and compute the welfare change. 

(1) Bihar <bihar_bybpl.m>

Read in price experiment data. 

x: amount of cash offered in lieu of ration entitlement (deviation from the fiscal value of the ration)

y: 0 or 1 cash take up chocie.
 
Estimate the CDF using linear sieve with B-spline basis imposing monotonicity constraint. <bSplineSieve.m>

Find the 95% confidence intervals for the welfare change using the method proposed by Chernozhukov, Newey, and Santos (2015). <CNS2015.m>

Note that bootstrapping cannot be used due to the presence of inequality constraints.


(2) Rajasthan

Read in stated CEV data.

x: CEV


Plot the empirical CDF of CEV.

Estimate the CDF using kernel density estimation method

Bootstrap the 95% confidence intervals for the welfare change. 

