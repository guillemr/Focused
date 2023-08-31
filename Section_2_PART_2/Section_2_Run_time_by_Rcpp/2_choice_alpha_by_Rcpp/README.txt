To generate Figure 2 of Section 2.3 "A simple online algorithm" :

(1) 2_GET_FIGURE_Gaussian_Model_Run_time_as_the_function_of_alpha

 we need obtain the results of the following test:

(1) 1_GET_TEST_3_Gaussian_Model_Run_time_as_the_function_of_alpha.R


This test generates the set of time series (Gaussian model without change) with 10^5 data points in  dimensions p = 2 and 3 and estimate the runtimes of MdFOCuS0 for 
alpha = c(0,1,2,3,5,10,20,50,100,200,500). 

For each dimension we write the results in the table, where the first column is the value of alpha, the second is the run time. These tables we also save in the files. 

The results of our testing are available in the archive “3_GET_RES_TEST_3_choice_alpha_by_Rcpp.zip”.


ATTENTION : We use for Gaussian model with first step = p+2 and beta = 1.

The figure of our testing is available in the archive “4_GET_FIGURE_choice_alpha_byRcpp.zip”.

