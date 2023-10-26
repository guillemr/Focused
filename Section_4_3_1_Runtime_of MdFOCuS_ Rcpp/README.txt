ATTENTION : 

The implementation of the MdFOCuS algorithm in R/C++ uses the Quickhull algorithm implemented in the C++ library "Qhull" (see link http://www.qhull.org), so before proceeding to use this implementation, make sure that this library is installed on your computer (or computing server).

If it is missing, you need to do the following:

(1) to download C++ interface to Qhull from the link: 

https://github.com/qhull/qhull

(2) Install this library on your personal computer (or computing server) using the installation guide on this page ( https://github.com/qhull/qhull ) .

	(a) EXAMPLE Installing Qhull on Unix with gcc

  		To build Qhull, static libraries, shared library, and C++ interface
  		- Download and extract Qhull (either GitHub, .tgz file, or .zip file)
  		- make
  		- export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH
  		- make test
  
  		'make install' installs Qhull at '/usr/local/'.  It installs pkg-config files 
  		at '/usr/local/lib/pkgconfig'.  Change the install directory with DESTDIR and PREFIX.



ATTENTION : 

In case the Qhull library was installed at '/usr/local/' (by default) you can install the implementation of the MdFOCuS algorithm in R/C++ by the following command from package "remotes" (or "devtools") in RStudio:


		install.packages("remotes")
		library(remotes)
		remotes::install_github("lpishchagina/focus", force = TRUE)



otherwise you need:

(1) to download its implementation from the link: 

		https://github.com/lpishchagina/focus

(2) to change in files "../src/Makevars" and "../src/Makevars.win" two addresses:
	
	(a) PKG_CPPFLAGS =  -I"/usr/local/include/libqhulcpp" 
	
	to address of folder "libqhulcpp" (it contains .cpp files)


	(b) PKG_LIBS= -L "/usr/local/lib/" -lqhull_r -lqhullcpp -lqhullstatic_r
	 
	to address of folder "lib" (it contains libraries)


(3) to install this package on your personal computer (or computing server).



(4) R Checking : 

	install.packages("remotes")
	library(remotes)
	remotes::install_github("lpishchagina/focus", force = TRUE)
	
	ts_ <- generate_ts(type = "gauss", p = 2, n = 25, changes = NULL, means = matrix(0, ncol = 1, nrow = 2))
	
	getChangePoints(ts_, method = "FOCuS0", cost = "gauss", common_difference_step = 1, common_ratio_step = 2)
