To generate Figure 1 of  Section 2.2 "Bound on The Number of changepoint Candidates"

(and Figure 9 in Appendix D.1):

(1) 3_GET_FIGURE_Stirling_Gaussian_Model_Convex_hull_estimations.R

(2) 4_GET_FIGURE_APPENDIX_Stirling_Poisson_Model_Convex_hull_estimations.R


 we need obtain the results of the following tests:

(1) 1_GET_TEST_1_Gaussian_Model_Convex_hull_estimations.R

(2) 2_GET_APPENDIX_TEST_1_Poisson_Model_Convex_hull_estimations.R

These tests generate the set of time series (Gaussian and Poisson models without change) with different number of data points in  dimensions p = 1,..5 and estimate the number of facets and vertices of convex hull {P(ts)}. 

For each dimension we write the results in two tables, where the first column is the number of data points, the second is the number of facets (or vertices). These tables we also save in the files. 

The results of our testing are available in the archive “5_GET_RES_TEST_1_convex_hull_estimations.zip”.


ATTENTION : We also add a theoretical estimation from Theorem 6 to the graph. This estimate is constructed using Stirling numbers of the first kind, so the construction of the graphs takes a lot of time (30 mn).


You also add plots of these estimates with the theoretical asymptotic estimation from Corollary 1:

(1) not_included_GET_FIGURE_Gaussian_Model_Convex_hull_estimations.R

(2) not_included_GET_APPENDIX_FIGURE_Poisson_Model_Convex_hull_estimations.R

These graphs are not included in the article.

The figures of our testing are available in the archive “6_GET_FIGURES_convex_hull_estimations.zip”.

