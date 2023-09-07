(1) MAIN_3_RuntimeAsFunOfAlpha.R

Thise file generates the set of time series (Gaussian model without change with n=10^5 data points) in  dimensions p = 2 and 3 and estimates the runtimes of MdFOCuS0  with different values of alpha. 

For each dimension we write the results in the table, where the first column is the value of alpha, the second is the run time. These tables we also save in the files.



(2) MAIN_4_RuntimeAsFunOfP.R

Thise test generates the set of time series (Gaussian and Poisson models without change) with different number of data points in  dimensions p = 1,..5 and estimate the runtimes of MdFOCuS and MdFOCuS0 (time limit = 3 mn). 

For each dimension we write the results in the table, where the first column is the number of data points, the second is the run time. These tables we also save in the files.

(3) MAIN_5_FiguresOfRuntimeByRcpp.R

This file generates the Figures 3 of this section by the results of MAIN_3 and MAIN_5


(4)The results of our testing and the figures are available in the archive “Results_of_MAIN_3_4_5.zip”.
