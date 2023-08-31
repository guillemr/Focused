To generate Figure 3 of Section 2.3 "A simple online algorithm" (and Figure 10 in Appendix D.2 Empirical run times of MdFOCuS0 and MdFOCuS):

(1) 3_GET_FIGURE_Gaussian_Model_Run_time_alpha_2_beta_1.R

(2) 4_GET_APPENDIX_FIGURE_Poisson_Model_Run_time_step_15_alpha_2_beta_1.R


 we need obtain the results of the following tests:

(1) 1_GET_TEST_2_Gaussian_Model_Run_time_alpha_2_beta_1.R

(2) 2_GET_APPENDIX_TEST_2_Poisson_Model_Run_time_FOCuS_step_14_alpha_2_beta_1.R

These tests generate the set of time series (Gaussian and Poisson models without change) with different number of data points in  dimensions p = 1,..5 and estimate the runtimes of MdFOCuS and MdFOCuS0 (time limit = 3 mn). 

For each dimension we write the results in the table, where the first column is the number of data points, the second is the run time. These tables we also save in the files. 

The results of our testing are available in the archive “5_GET_RES_TEST_2_runtimes_by_Rcpp.zip”.


ATTENTION : We use for Gaussian model first step = p+2, alpha = 2, beta = 1 and for Poisson model first step = p+15, alpha = 2, beta = 1.

The figures of our testing are available in the archive “6_GET_FIGURES_runtimes_byRcpp.zip”.

