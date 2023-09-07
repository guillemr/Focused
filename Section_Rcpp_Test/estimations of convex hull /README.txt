(1) MAIN_1_EmpirEstimOfCH.R 

This file generates the set of time series (Gaussian and Poisson models without change) with different number of data points in  dimensions p = 1,..5 and estimates the number of facets and vertices of convex hull {P(ts)}. 

For each dimension we write the results in two tables, where the first column is the number of data points, the second is the number of facets (or vertices). These tables we also save in the files. 

(2) MAIN_2_EmpirAndTheorEstimOfCH.R (using the results of MAIN_1)

This file generates Figures with empirical and theretical estimations for Gaussian and Poisson models.
ATTENTION : We add a theoretical estimation from Theorem 6 to the graph. This estimate is constructed using Stirling numbers of the first kind, so the construction of the graphs takes a lot of time (30 mn for each graph).


(3) The results of our testing (MAIN_1) and the figures (MAIN_2) are available in the archive “Results_of_MAIN_1_2.zip”.
